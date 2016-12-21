/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InterruptedIOException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import org.apache.commons.configuration.CompositeConfiguration;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;

import edu.rit.mp.DoubleBuf;

import ffx.numerics.Potential;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.PDBFilter;

/**
 * An implementation of Transition-Tempered Orthogonal Space Random Walk
 * algorithm.
 *
 * @author Michael J. Schnieders, James Dama, Wei Yang and Pengyu Ren
 */
public class TransitionTemperedOSRW extends AbstractOSRW {

    private static final Logger logger = Logger.getLogger(TransitionTemperedOSRW.class.getName());

    /**
     * The recursion kernel stores the weight of each [lambda][Flambda] bin.
     */
    private double recursionKernel[][];
    /**
     * The recursionWeights stores the [Lambda, FLambda] weight for each
     * process. Therefore the array is of size [number of Processes][2]. Each 2
     * entry array must be wrapped inside a Parallel Java IntegerBuf for the
     * All-Gather communication calls.
     */
    private final double recursionWeights[][];
    private final double myRecursionWeight[];
    /**
     * These DoubleBufs wrap the recusionWeight arrays.
     */
    private final DoubleBuf recursionWeightsBuf[];
    private final DoubleBuf myRecursionWeightBuf;


    /**
     * Total histogram weight.
     */
    private double totalWeight;
    /**
     * A flag to indicate if the transition has been crossed and Dama et al.
     * transition-tempering should begin.
     */
    private boolean tempering = false;
    /**
     * The Dama et al. transition-tempering rate parameter. A reasonable value
     * is about 2 to 4 kT.
     */
    private double deltaT = 4.0 * R * 298.0;
    /**
     * The Dama et al. transition-tempering weight: temperingWeight =
     * exp(-max(G(L,F_L))/deltaT)
     */
    private double temperingWeight = 1.0;
    /**
     * Transition detection flags. The transition is currently defined using the
     * following logic: First, the simulation needs to pass through mid-range
     * [0.45 to 0.55] lambda values. Second, the simulation needs to pass
     * through both low [below 0.05] & high values [above 0.95].
     *
     * The goal is to delay tempering for a landscape that is monotonic. For
     * example, if the simulation starts at the top of the hill (i.e. lambda =
     * 0) and move quickly to the basin (i.e. lambda = 1) this will not count as
     * having found the transition. The transition will not be achieved until
     * the basin is full and lambda = 0 is sampled again.
     */
    private boolean midLambda = false;
    private boolean lowLambda = false;
    private boolean highLambda = false;

    /**
     * The ReceiveThread accumulates OSRW statistics from multiple asynchronous
     * walkers.
     */
    private final ReceiveThread receiveThread;

    /**
     * OSRW Asynchronous MultiWalker Constructor.
     *
     * @param lambdaInterface defines Lambda and dU/dL.
     * @param potential defines the Potential energy.
     * @param lambdaFile contains the current Lambda particle position and
     * velocity.
     * @param histogramFile contains the Lambda and dU/dL histogram.
     * @param properties defines System properties.
     * @param temperature the simulation temperature.
     * @param dt the time step.
     * @param printInterval number of steps between logging updates.
     * @param saveInterval number of steps between restart file updates.
     * @param asynchronous set to true if walkers run asynchronously.
     * @param algorithmListener the AlgorithmListener to be notified of
     * progress.
     */
    public TransitionTemperedOSRW(LambdaInterface lambdaInterface, Potential potential,
            File lambdaFile, File histogramFile, CompositeConfiguration properties,
            double temperature, double dt, double printInterval,
            double saveInterval, boolean asynchronous,
            AlgorithmListener algorithmListener) {
        this(lambdaInterface, potential, lambdaFile, histogramFile, properties,
                temperature, dt, printInterval, saveInterval, asynchronous,
                true, algorithmListener);
    }

    /**
     * OSRW Asynchronous MultiWalker Constructor.
     *
     * @param lambdaInterface defines Lambda and dU/dL.
     * @param potential defines the Potential energy.
     * @param lambdaFile contains the current Lambda particle position and
     * velocity.
     * @param histogramFile contains the Lambda and dU/dL histogram.
     * @param properties defines System properties.
     * @param temperature the simulation temperature.
     * @param dt the time step.
     * @param printInterval number of steps between logging updates.
     * @param saveInterval number of steps between restart file updates.
     * @param asynchronous set to true if walkers run asynchronously.
     * @param resetNumSteps whether to reset energy counts to 0
     * @param algorithmListener the AlgorithmListener to be notified of
     * progress.
     */
    public TransitionTemperedOSRW(LambdaInterface lambdaInterface, Potential potential,
            File lambdaFile, File histogramFile, CompositeConfiguration properties,
            double temperature, double dt, double printInterval,
            double saveInterval, boolean asynchronous, boolean resetNumSteps,
            AlgorithmListener algorithmListener) {
        super(lambdaInterface, potential, lambdaFile, histogramFile, properties, temperature, dt, printInterval, saveInterval, asynchronous, resetNumSteps, algorithmListener);

        /**
         * Allocate space for the recursion kernel that stores weights.
         */
        recursionKernel = new double[lambdaBins][FLambdaBins];

        /**
         * Load the OSRW histogram restart file if it exists.
         */
        boolean readHistogramRestart = false;
        if (histogramFile != null && histogramFile.exists()) {
            try {
                TTOSRWHistogramReader osrwHistogramReader = new TTOSRWHistogramReader(new FileReader(histogramFile));
                osrwHistogramReader.readHistogramFile();
                logger.info(String.format("\n Continuing OSRW histogram from %s.", histogramFile.getName()));
                readHistogramRestart = true;
            } catch (FileNotFoundException ex) {
                logger.info(" Histogram restart file could not be found and will be ignored.");
            }
        }

        /**
         * Load the OSRW lambda restart file if it exists.
         */
        if (lambdaFile != null && lambdaFile.exists()) {
            try {
                TTOSRWLambdaReader osrwLambdaReader = new TTOSRWLambdaReader(new FileReader(lambdaFile));
                osrwLambdaReader.readLambdaFile(resetNumSteps);
                logger.info(String.format("\n Continuing OSRW lambda from %s.", lambdaFile.getName()));
            } catch (FileNotFoundException ex) {
                logger.info(" Lambda restart file could not be found and will be ignored.");
            }
        }

        if (asynchronous) {
            /**
             * Use asynchronous communication.
             */

            myRecursionWeight = new double[3];
            myRecursionWeightBuf = DoubleBuf.buffer(myRecursionWeight);

            receiveThread = new ReceiveThread();
            receiveThread.start();
            recursionWeights = null;
            recursionWeightsBuf = null;
        } else {
            /**
             * Use synchronous communication.
             */
            recursionWeights = new double[numProc][3];
            recursionWeightsBuf = new DoubleBuf[numProc];
            for (int i = 0; i < numProc; i++) {
                recursionWeightsBuf[i] = DoubleBuf.buffer(recursionWeights[i]);
            }
            myRecursionWeight = recursionWeights[rank];
            myRecursionWeightBuf = recursionWeightsBuf[rank];
            receiveThread = null;
        }

        /**
         * Update and print out the recursion slave.
         */
        if (readHistogramRestart) {
            updateFLambda(true);
        }

    }

    @Override
    public double energyAndGradient(double[] x, double[] gradient) {

        double e = potential.energyAndGradient(x, gradient);

        /**
         * OSRW is propagated with the slowly varying terms.
         */
        if (state == STATE.FAST) {
            return e;
        }

        if (osrwOptimization && lambda > osrwOptimizationLambdaCutoff) {
            if (energyCount % osrwOptimizationFrequency == 0) {
                logger.info(String.format(" OSRW Minimization (Step %d)", energyCount));

                // Set Lambda value to 1.0.
                lambdaInterface.setLambda(1.0);

                potential.setEnergyTermState(Potential.STATE.BOTH);

                // Optimize the system.
                Minimize minimize = new Minimize(null, potential, null);
                minimize.minimize(osrwOptimizationEps);

                // Collect the minimum energy.
                double minEnergy = potential.getTotalEnergy();
                // If a new minimum has been found, save its coordinates.
                if (minEnergy < osrwOptimum) {
                    osrwOptimum = minEnergy;
                    logger.info(String.format(" New minimum energy found: %16.8f (Step %d).", osrwOptimum, energyCount));
                    int n = potential.getNumberOfVariables();
                    osrwOptimumCoords = new double[n];
                    osrwOptimumCoords = potential.getCoordinates(osrwOptimumCoords);
                    if (pdbFilter.writeFile(pdbFile, false)) {
                        logger.info(String.format(" Wrote PDB file to " + pdbFile.getName()));
                    }
                }

                // Reset lambda value.
                lambdaInterface.setLambda(lambda);

                // Remove the scaling of coordinates & gradient set by the minimizer.
                potential.setScaling(null);

                // Reset the Potential State
                potential.setEnergyTermState(state);

                // Revert to the coordinates and gradient prior to optimization.
                double eCheck = potential.energyAndGradient(x, gradient);

                if (abs(eCheck - e) > osrwOptimizationTolerance) {
                    logger.warning(String.format(
                            " OSRW optimization could not revert coordinates %16.8f vs. %16.8f.", e, eCheck));
                }
            }
        }

        double biasEnergy = 0.0;
        dEdLambda = lambdaInterface.getdEdL();
        d2EdLambda2 = lambdaInterface.getd2EdL2();
        int lambdaBin = binForLambda(lambda);
        int FLambdaBin = binForFLambda(dEdLambda);
        double dEdU = dEdLambda;

        if (propagateLambda) {
            energyCount++;
            detectTransition();
        }

        /**
         * Calculate recursion kernel G(L, F_L) and its derivatives with respect
         * to L and F_L.
         */
        double dGdLambda = 0.0;
        double dGdFLambda = 0.0;
        double ls2 = (2.0 * dL) * (2.0 * dL);
        double FLs2 = (2.0 * dFL) * (2.0 * dFL);
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int lcenter = lambdaBin + iL;
            double deltaL = lambda - (lcenter * dL);
            double deltaL2 = deltaL * deltaL;
            // Mirror conditions for recursion kernel counts.
            int lcount = lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                mirrorFactor = 2.0;
            } else if (lcount < 0) {
                lcount = -lcount;
            } else if (lcount > lambdaBins - 1) {
                // Number of bins past the last bin
                lcount -= (lambdaBins - 1);
                // Mirror bin
                lcount = lambdaBins - 1 - lcount;
            }
            for (int iFL = -biasCutoff; iFL <= biasCutoff; iFL++) {
                int FLcenter = FLambdaBin + iFL;
                /**
                 * If either of the following FL edge conditions are true, then
                 * there are no counts and we continue.
                 */
                if (FLcenter < 0 || FLcenter >= FLambdaBins) {
                    continue;
                }
                double deltaFL = dEdLambda - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double weight = mirrorFactor * recursionKernel[lcount][FLcenter];
                double bias = weight * biasMag
                        * exp(-deltaL2 / (2.0 * ls2))
                        * exp(-deltaFL2 / (2.0 * FLs2));
                biasEnergy += bias;
                dGdLambda -= deltaL / ls2 * bias;
                dGdFLambda -= deltaFL / FLs2 * bias;
            }
        }

        /**
         * Lambda gradient due to recursion kernel G(L, F_L).
         */
        dEdLambda += dGdLambda + dGdFLambda * d2EdLambda2;

        /**
         * Cartesian coordinate gradient due to recursion kernel G(L, F_L).
         */
        fill(dUdXdL, 0.0);
        lambdaInterface.getdEdXdL(dUdXdL);
        for (int i = 0; i < nVariables; i++) {
            gradient[i] += dGdFLambda * dUdXdL[i];
        }

        if (propagateLambda && energyCount > 0) {
            /**
             * Update free energy F(L) every ~10 steps.
             */
            if (energyCount % 10 == 0) {
                fLambdaUpdates++;
                boolean printFLambda = fLambdaUpdates % fLambdaPrintInterval == 0;
                totalFreeEnergy = updateFLambda(printFLambda);
                /**
                 * Calculating Moving Average & Standard Deviation
                 */
                totalAverage += totalFreeEnergy;
                totalSquare += Math.pow(totalFreeEnergy, 2);
                periodCount++;
                if (periodCount == window - 1) {
                    lastAverage = totalAverage / window;
                    //lastStdDev = Math.sqrt((totalSquare - Math.pow(totalAverage, 2) / window) / window);
                    lastStdDev = Math.sqrt((totalSquare - (totalAverage * totalAverage)) / (window * window));
                    logger.info(String.format(" The running average is %12.4f kcal/mol and the stdev is %8.4f kcal/mol.",
                            lastAverage, lastStdDev));
                    totalAverage = 0;
                    totalSquare = 0;
                    periodCount = 0;
                }
            }
            if (energyCount % saveFrequency == 0) {
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(lambdaOneAssembly);
                }
                /**
                 * Only the rank 0 process writes the histogram restart file.
                 */
                if (rank == 0) {
                    try {
                        TTOSRWHistogramWriter ttOSRWHistogramRestart = new TTOSRWHistogramWriter(
                                new BufferedWriter(new FileWriter(histogramFile)));
                        ttOSRWHistogramRestart.writeHistogramFile();
                        ttOSRWHistogramRestart.flush();
                        ttOSRWHistogramRestart.close();
                        logger.info(String.format(" Wrote TTOSRW histogram restart file to %s.", histogramFile.getName()));
                    } catch (IOException ex) {
                        String message = " Exception writing TTOSRW histogram restart file.";
                        logger.log(Level.INFO, message, ex);
                    }
                }
                /**
                 * All ranks write a lambda restart file.
                 */
                try {
                    TTOSRWLambdaWriter ttOSRWLambdaRestart = new TTOSRWLambdaWriter(new BufferedWriter(new FileWriter(lambdaFile)));
                    ttOSRWLambdaRestart.writeLambdaFile();
                    ttOSRWLambdaRestart.flush();
                    ttOSRWLambdaRestart.close();
                    logger.info(String.format(" Wrote TTOSRW lambda restart file to %s.", lambdaFile.getName()));
                } catch (IOException ex) {
                    String message = " Exception writing TTOSRW lambda restart file.";
                    logger.log(Level.INFO, message, ex);
                }
            }
            /**
             * Write out snapshot upon each full lambda traversal.
             */
            if (writeTraversalSnapshots) {
                double heldTraversalLambda = 0.5;
                if (!traversalInHand.isEmpty()) {
                    heldTraversalLambda = Double.parseDouble(traversalInHand.get(0).split(",")[0]);
                    if ((lambda > 0.2 && traversalSnapshotTarget == 0)
                            || (lambda < 0.8 && traversalSnapshotTarget == 1)) {
                        int snapshotCounts = Integer.parseInt(traversalInHand.get(0).split(",")[1]);
                        traversalInHand.remove(0);
                        File fileToWrite;
                        int numStructures;
                        if (traversalSnapshotTarget == 0) {
                            fileToWrite = lambdaZeroFile;
                            numStructures = ++lambdaZeroStructures;
                        } else {
                            fileToWrite = lambdaOneFile;
                            numStructures = ++lambdaOneStructures;
                        }
                        try {
                            FileWriter fw = new FileWriter(fileToWrite, true);
                            BufferedWriter bw = new BufferedWriter(fw);
                            bw.write(String.format("MODEL        %d          L=%.4f  counts=%d", numStructures, heldTraversalLambda, snapshotCounts));
                            for (int i = 0; i < 50; i++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            for (int i = 0; i < traversalInHand.size(); i++) {
                                bw.write(traversalInHand.get(i));
                                bw.newLine();
                            }
                            bw.write(String.format("ENDMDL"));
                            for (int i = 0; i < 75; i++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            bw.close();
                            logger.info(String.format(" Wrote traversal structure L=%.4f", heldTraversalLambda));
                        } catch (Exception exception) {
                            logger.warning(String.format("Exception writing to file: %s", fileToWrite.getName()));
                        }
                        heldTraversalLambda = 0.5;
                        traversalInHand.clear();
                        traversalSnapshotTarget = 1 - traversalSnapshotTarget;
                    }
                }
                if (((lambda < 0.1 && traversalInHand.isEmpty())
                        || (lambda < heldTraversalLambda - 0.025 && !traversalInHand.isEmpty()))
                        && (traversalSnapshotTarget == 0 || traversalSnapshotTarget == -1)) {
                    if (lambdaZeroFilter == null) {
                        lambdaZeroFilter = new PDBFilter(lambdaZeroFile, lambdaZeroAssembly, null, null);
                        lambdaZeroFilter.setListMode(true);
                    }
                    lambdaZeroFilter.clearListOutput();
                    lambdaZeroFilter.writeFileWithHeader(lambdaFile, new StringBuilder(String.format("%.4f,%d", lambda, totalWeight)));
                    traversalInHand = lambdaZeroFilter.getListOutput();
                    traversalSnapshotTarget = 0;
                } else if (((lambda > 0.9 && traversalInHand.isEmpty()) || (lambda > heldTraversalLambda + 0.025 && !traversalInHand.isEmpty()))
                        && (traversalSnapshotTarget == 1 || traversalSnapshotTarget == -1)) {
                    if (lambdaOneFilter == null) {
                        lambdaOneFilter = new PDBFilter(lambdaOneFile, lambdaOneAssembly, null, null);
                        lambdaOneFilter.setListMode(true);
                    }
                    lambdaOneFilter.clearListOutput();
                    lambdaOneFilter.writeFileWithHeader(lambdaFile, new StringBuilder(String.format("%.4f,%d", lambda, totalWeight)));
                    traversalInHand = lambdaOneFilter.getListOutput();
                    traversalSnapshotTarget = 1;
                }
            }
        }

        /**
         * Compute the energy and gradient for the recursion slave at F(L) using
         * interpolation.
         */
        double freeEnergy = currentFreeEnergy();
        biasEnergy += freeEnergy;

        if (print) {
            logger.info(String.format(" %s %16.8f", "Bias Energy       ", biasEnergy));
            logger.info(String.format(" %s %16.8f  %s",
                    "OSRW Potential    ", e + biasEnergy, "(Kcal/mole)"));
        }

        if (propagateLambda && energyCount > 0) {
            /**
             * Log the current Lambda state.
             */
            if (energyCount % printFrequency == 0) {
                if (lambdaBins < 1000) {
                    logger.info(String.format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f",
                            lambda, lambdaBin, dEdU, dEdLambda - dEdU, dEdLambda));
                } else {
                    logger.info(String.format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f",
                            lambda, lambdaBin, dEdU, dEdLambda - dEdU, dEdLambda));
                }
            }

            /**
             * Metadynamics grid counts (every 'countInterval' steps).
             */
            if (energyCount % countInterval == 0) {
                if (jobBackend != null) {
                    if (world.size() > 1) {
                        jobBackend.setComment(String.format("Overall dG=%10.4f at %7.3e psec, Current: [L=%6.4f, F_L=%10.4f, dG=%10.4f] at %7.3e psec",
                                totalFreeEnergy, totalWeight * dt * countInterval, lambda, dEdU, -freeEnergy, energyCount * dt));
                    } else {
                        jobBackend.setComment(String.format("Overall dG=%10.4f at %7.3e psec, Current: [L=%6.4f, F_L=%10.4f, dG=%10.4f]",
                                totalFreeEnergy, totalWeight * dt * countInterval, lambda, dEdU, -freeEnergy));
                    }
                }
                if (asynchronous) {
                    asynchronousSend(lambda, dEdU);
                } else {
                    synchronousSend(lambda, dEdU);
                }
            }
        }

        /**
         * Propagate the Lambda particle.
         */
        if (propagateLambda) {
            langevin();
        } else {
            equilibrationCounts++;
            if (jobBackend != null) {
                jobBackend.setComment(String.format("Equilibration [L=%6.4f, F_L=%10.4f]", lambda, dEdU));
            }
            if (equilibrationCounts % 10 == 0) {
                logger.info(String.format(" L=%6.4f, F_L=%10.4f", lambda, dEdU));
            }
        }

        totalEnergy = e + biasEnergy;

        return totalEnergy;
    }

    private void detectTransition() {
        if (tempering) {
            return;
        }

        /**
         * Before detecting both lambda extremes, the simulation must pass
         * through intermediate lambda values.
         */
        if (midLambda == false) {
            if (lambda >= 0.45 && lambda <= 0.55) {
                midLambda = true;
            }
        } else {
            /**
             * After passing through intermediate lambda values, detect both
             * extreme values.
             */
            if (lambda >= 0.95) {
                highLambda = true;
            } else if (lambda <= 0.05) {
                lowLambda = true;
            }
            if (lowLambda == true && highLambda == true) {
                tempering = true;
                logger.info(String.format(" Tempering activated at Lambda = %6.4f.", lambda));
            }
        }

    }

    /**
     * Send an OSRW count to all other processes while also receiving an OSRW
     * count from all other processes.
     *
     * @param lambda
     * @param dEdU
     */
    private void synchronousSend(double lambda, double dEdU) {
        /**
         * All-Gather counts from each walker.
         */
        myRecursionWeight[0] = lambda;
        myRecursionWeight[1] = dEdU;
        myRecursionWeight[2] = temperingWeight;
        try {
            world.allGather(myRecursionWeightBuf, recursionWeightsBuf);
        } catch (IOException ex) {
            String message = " Multi-walker OSRW allGather failed.";
            logger.log(Level.SEVERE, message, ex);
        }

        /**
         * Find the minimum and maximum FLambda bin for the gathered counts.
         */
        double minRequired = Double.MAX_VALUE;
        double maxRequired = Double.MIN_VALUE;
        for (int i = 0; i < numProc; i++) {
            minRequired = Math.min(minRequired, recursionWeights[i][1]);
            maxRequired = Math.max(maxRequired, recursionWeights[i][1]);
        }

        /**
         * Check that the FLambda range of the Recursion kernel includes both
         * the minimum and maximum FLambda value.
         */
        checkRecursionKernelSize(minRequired);
        checkRecursionKernelSize(maxRequired);

        /**
         * Increment the Recursion Kernel based on the input of each walker.
         */
        for (int i = 0; i < numProc; i++) {
            int walkerLambda = binForLambda(recursionWeights[i][0]);
            int walkerFLambda = binForFLambda(recursionWeights[i][1]);
            double weight = recursionWeights[i][2];

            /**
             * If the weight is less than 1.0, then a walker has activated
             * tempering.
             */
            if (tempering == false && weight < 1.0) {
                tempering = true;
                logger.info(String.format(" Tempering activated due to recieved weight of (%8.6f)", weight));
            }

            if (resetStatistics && recursionWeights[i][0] > lambdaResetValue) {
                recursionKernel = new double[lambdaBins][FLambdaBins];
                resetStatistics = false;
                logger.info(String.format(" Cleared OSRW histogram (Lambda = %6.4f).", recursionWeights[i][0]));
            }

            recursionKernel[walkerLambda][walkerFLambda] += weight;
        }
    }

    /**
     * Send an OSRW count to all other processes.
     *
     * @param lambda
     * @param dEdU
     */
    private void asynchronousSend(double lambda, double dEdU) {
        myRecursionWeight[0] = lambda;
        myRecursionWeight[1] = dEdU;
        myRecursionWeight[2] = temperingWeight;

        for (int i = 0; i < numProc; i++) {
            try {
                world.send(i, myRecursionWeightBuf);
            } catch (Exception ex) {
                String message = " Asynchronous Multiwalker OSRW send failed.";
                logger.log(Level.SEVERE, message, ex);
            }
        }
    }

    /**
     * If necessary, allocate more space.
     */
    private void checkRecursionKernelSize(double dEdLambda) {
        if (dEdLambda > maxFLambda) {
            logger.info(String.format(" Current F_lambda %8.2f > maximum historgram size %8.2f.",
                    dEdLambda, maxFLambda));

            double origDeltaG = updateFLambda(false);

            int newFLambdaBins = FLambdaBins;
            while (minFLambda + newFLambdaBins * dFL < dEdLambda) {
                newFLambdaBins += 100;
            }
            double newRecursionKernel[][] = new double[lambdaBins][newFLambdaBins];
            /**
             * We have added bins above the indeces of the current counts just
             * copy them into the new array.
             */
            for (int i = 0; i < lambdaBins; i++) {
                System.arraycopy(recursionKernel[i], 0, newRecursionKernel[i], 0, FLambdaBins);
            }
            recursionKernel = newRecursionKernel;
            FLambdaBins = newFLambdaBins;
            maxFLambda = minFLambda + dFL * FLambdaBins;
            logger.info(String.format(" New historgram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            assert (origDeltaG == updateFLambda(false));

        }
        if (dEdLambda < minFLambda) {
            logger.info(String.format(" Current F_lambda %8.2f < minimum historgram size %8.2f.",
                    dEdLambda, minFLambda));

            double origDeltaG = updateFLambda(false);

            int offset = 100;
            while (dEdLambda < minFLambda - offset * dFL) {
                offset += 100;
            }
            int newFLambdaBins = FLambdaBins + offset;
            double newRecursionKernel[][] = new double[lambdaBins][newFLambdaBins];
            /**
             * We have added bins below the current counts, so their indeces
             * must be increased by: offset = newFLBins - FLBins
             */
            for (int i = 0; i < lambdaBins; i++) {
                System.arraycopy(recursionKernel[i], 0, newRecursionKernel[i], offset, FLambdaBins);
            }
            recursionKernel = newRecursionKernel;
            minFLambda = minFLambda - offset * dFL;
            FLambdaBins = newFLambdaBins;
            logger.info(String.format(" New historgram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            assert (origDeltaG == updateFLambda(false));
        }
    }

    /**
     * Eq. 7 from the Xtal Thermodynamics paper.
     *
     * @param print
     * @return the current free energy.
     */
    @Override
    protected double updateFLambda(boolean print) {
        double freeEnergy = 0.0;
        double minFL = Double.MAX_VALUE;
        totalWeight = 0;
        StringBuilder stringBuilder = new StringBuilder();
        if (print) {
            stringBuilder.append(" Weight    Lambda Bins    F_Lambda Bins   <   F_L  >  Max F_L     dG        G\n");
        }
        for (int iL = 0; iL < lambdaBins; iL++) {
            int ulFL = -1;
            int llFL = -1;

            // Find the smallest FL bin.
            for (int jFL = 0; jFL < FLambdaBins; jFL++) {
                double count = recursionKernel[iL][jFL];
                if (count > 0) {
                    llFL = jFL;
                    break;
                }
            }

            // Find the largest FL bin.
            for (int jFL = FLambdaBins - 1; jFL >= 0; jFL--) {
                double count = recursionKernel[iL][jFL];
                if (count > 0) {
                    ulFL = jFL;
                    break;
                }
            }
            
            double lambdaCount = 0;
            // The FL range sampled for lambda bin [iL*dL .. (iL+1)*dL]
            double lla = 0.0;
            double ula = 0.0;
            double maxBias = 0;
            if (ulFL == -1) {
                FLambda[iL] = 0.0;
                minFL = 0.0;
            } else {
                double ensembleAverageFLambda = 0.0;
                double partitionFunction = 0.0;
                for (int jFL = llFL; jFL <= ulFL; jFL++) {
                    double currentFLambda = minFLambda + jFL * dFL + dFL_2;
                    double kernel = evaluateKernel(iL, jFL);
                    if (kernel > maxBias) {
                        maxBias = kernel;
                    }
                    double weight = exp(kernel / (R * temperature));
                    ensembleAverageFLambda += currentFLambda * weight;
                    partitionFunction += weight;
                    lambdaCount += recursionKernel[iL][jFL];
                }
                if (minFL > maxBias) {
                    minFL = maxBias;
                }
                FLambda[iL] = ensembleAverageFLambda / partitionFunction;
                lla = minFLambda + llFL * dFL;
                ula = minFLambda + (ulFL + 1) * dFL;
            }

            // The first and last lambda bins are half size.
            double delta = dL;
            if (iL == 0 || iL == lambdaBins - 1) {
                delta = dL_2;
            }
            double deltaFreeEnergy = FLambda[iL] * delta;
            freeEnergy += deltaFreeEnergy;
            totalWeight += lambdaCount;

            if (print) {
                double llL = iL * dL - dL_2;
                double ulL = llL + dL;
                if (llL < 0.0) {
                    llL = 0.0;
                }
                if (ulL > 1.0) {
                    ulL = 1.0;
                }
                stringBuilder.append(String.format(" %6.2e  %5.3f %5.3f   %7.1f %7.1f   %8.3f  %8.3f  %8.3f %8.3f\n",
                        lambdaCount, llL, ulL, lla, ula,
                        FLambda[iL], maxBias, deltaFreeEnergy, freeEnergy));
            }
        }

        if (tempering) {
            temperingWeight = exp(-minFL / deltaT);
        }

        if (abs(freeEnergy - previousFreeEnergy) > 0.001) {
            if (print) {
                stringBuilder.append(String.format(" Minimum Bias %8.3f", minFL));
                logger.info(stringBuilder.toString());
                previousFreeEnergy = freeEnergy;
            }

        }

        logger.info(String.format(" The free energy is %12.4f kcal/mol (Counts: %6.2e, Weight: %6.4f).",
                freeEnergy, totalWeight, temperingWeight));

        return freeEnergy;
    }


    @Override
    public boolean destroy() {
        if (receiveThread != null) {
            receiveThread.interrupt();
        }
        return true;
    }

    @Override
    protected double evaluateKernel(int cLambda, int cF_Lambda) {
        /**
         * Compute the value of L and FL for the center of the current bin.
         */
        double vL = cLambda * dL;
        double vFL = minFLambda + cF_Lambda * dFL + dFL_2;
        /**
         * Set the variances for the Gaussian bias.
         */
        double Ls2 = 2.0 * dL * 2.0 * dL;
        double FLs2 = 2.0 * dFL * 2.0 * dFL;
        double sum = 0.0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int Lcenter = cLambda + iL;
            double deltaL = vL - Lcenter * dL;
            double deltaL2 = deltaL * deltaL;

            // Mirror condition for Lambda counts.
            int lcount = Lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                /**
                 * The width of the first and last bins is dLambda_2, so the
                 * mirror condition is to double their counts.
                 */
                mirrorFactor = 2.0;
            } else if (lcount < 0) {
                lcount = -lcount;
            } else if (lcount > lambdaBins - 1) {
                // number of bins past the last bin.
                lcount -= (lambdaBins - 1);
                // mirror bin
                lcount = lambdaBins - 1 - lcount;
            }

            for (int jFL = -biasCutoff; jFL <= biasCutoff; jFL++) {
                int FLcenter = cF_Lambda + jFL;
                /**
                 * For FLambda outside the count matrix the weight is 0 so we
                 * continue.
                 */
                if (FLcenter < 0 || FLcenter >= FLambdaBins) {
                    continue;
                }
                double deltaFL = vFL - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double weight = mirrorFactor * recursionKernel[lcount][FLcenter];
                if (weight > 0) {
                    double e = weight * biasMag * exp(-deltaL2 / (2.0 * Ls2))
                            * exp(-deltaFL2 / (2.0 * FLs2));
                    sum += e;
                }
            }
        }
        return sum;
    }

    /**
     * Sets the Dama et al tempering parameter, as a multiple of kbT. T is
     * presently assumed to be 298.0K.
     * @param kbtMult
     */
    public void setDeltaT(double kbtMult) {
        deltaT = R * 298.0 * kbtMult;
    }
    private class TTOSRWHistogramWriter extends PrintWriter {

        public TTOSRWHistogramWriter(Writer writer) {
            super(writer);
        }

        public void writeHistogramFile() {
            printf("Temperature     %15.3f\n", temperature);
            printf("Lambda-Mass     %15.8e\n", thetaMass);
            printf("Lambda-Friction %15.8e\n", thetaFriction);
            printf("Bias-Mag        %15.8e\n", biasMag);
            printf("Bias-Cutoff     %15d\n", biasCutoff);
            printf("Count-Interval  %15d\n", countInterval);
            printf("Lambda-Bins     %15d\n", lambdaBins);
            printf("FLambda-Bins    %15d\n", FLambdaBins);
            printf("Flambda-Min     %15.8e\n", minFLambda);
            printf("Flambda-Width   %15.8e\n", dFL);
            int flag = 0;
            if (tempering) {
                flag = 1;
            }
            printf("Tempering       %15d\n", flag);
            for (int i = 0; i < lambdaBins; i++) {
                printf("%g", recursionKernel[i][0]);
                for (int j = 1; j < FLambdaBins; j++) {
                    printf(" %g", recursionKernel[i][j]);
                }
                println();
            }
        }
    }

    private class TTOSRWLambdaWriter extends PrintWriter {

        public TTOSRWLambdaWriter(Writer writer) {
            super(writer);
        }

        public void writeLambdaFile() {
            printf("Lambda          %15.8f\n", lambda);
            printf("Lambda-Velocity %15.8e\n", halfThetaVelocity);
            printf("Steps-Taken     %15d\n", energyCount);
        }
    }

    private class TTOSRWHistogramReader extends BufferedReader {

        public TTOSRWHistogramReader(Reader reader) {
            super(reader);
        }

        public void readHistogramFile() {
            try {
                temperature = Double.parseDouble(readLine().split(" +")[1]);
                thetaMass = Double.parseDouble(readLine().split(" +")[1]);
                thetaFriction = Double.parseDouble(readLine().split(" +")[1]);
                biasMag = Double.parseDouble(readLine().split(" +")[1]);
                biasCutoff = Integer.parseInt(readLine().split(" +")[1]);
                countInterval = Integer.parseInt(readLine().split(" +")[1]);
                lambdaBins = Integer.parseInt(readLine().split(" +")[1]);
                FLambdaBins = Integer.parseInt(readLine().split(" +")[1]);
                minFLambda = Double.parseDouble(readLine().split(" +")[1]);
                dFL = Double.parseDouble(readLine().split(" +")[1]);
                int flag = Integer.parseInt(readLine().split(" +")[1]);
                if (flag != 0) {
                    tempering = true;
                } else {
                    tempering = false;
                }

                // Allocate memory for the recursion kernel.
                recursionKernel = new double[lambdaBins][FLambdaBins];
                for (int i = 0; i < lambdaBins; i++) {
                    String counts[] = readLine().split(" +");
                    for (int j = 0; j < FLambdaBins; j++) {
                        recursionKernel[i][j] = Double.parseDouble(counts[j]);
                    }
                }
            } catch (Exception e) {
                String message = " Invalid OSRW Histogram file.";
                logger.log(Level.SEVERE, message, e);
            }
        }
    }

    private class TTOSRWLambdaReader extends BufferedReader {

        public TTOSRWLambdaReader(Reader reader) {
            super(reader);
        }

        public void readLambdaFile() {
            readLambdaFile(true);
        }

        public void readLambdaFile(boolean resetEnergyCount) {
            try {
                lambda = Double.parseDouble(readLine().split(" +")[1]);
                halfThetaVelocity = Double.parseDouble(readLine().split(" +")[1]);
                setLambda(lambda);
            } catch (Exception e) {
                String message = " Invalid OSRW Lambda file.";
                logger.log(Level.SEVERE, message, e);
            }
            if (!resetEnergyCount) {
                try {
                    energyCount = Integer.parseUnsignedInt(readLine().split(" +")[1]);
                } catch (Exception e) {
                    logger.log(Level.WARNING, String.format(" Could not find number of steps taken in OSRW Lambda file: %s", e.toString()));
                }
            }
        }
    }

    class ReceiveThread extends Thread {

        final double recursionCount[];
        final DoubleBuf recursionCountBuf;

        public ReceiveThread() {
            recursionCount = new double[3];
            recursionCountBuf = DoubleBuf.buffer(recursionCount);
        }

        @Override
        public void run() {
            while (true) {
                try {
                    world.receive(null, recursionCountBuf);
                } catch (InterruptedIOException ioe) {
                    logger.log(Level.FINE, " ReceiveThread was interrupted at world.receive", ioe);
                    break;
                } catch (Exception e) {
                    String message = e.getMessage();
                    logger.log(Level.WARNING, message, e);
                }
                /**
                 * Check that the FLambda range of the Recursion kernel includes
                 * both the minimum and maximum FLambda value.
                 */
                checkRecursionKernelSize(recursionCount[1]);

                /**
                 * Increment the Recursion Kernel based on the input of current
                 * walker.
                 */
                int walkerLambda = binForLambda(recursionCount[0]);
                int walkerFLambda = binForFLambda(recursionCount[1]);
                double weight = recursionCount[2];

                /**
                 * If the weight is less than 1.0, then a walker has activated
                 * tempering.
                 */
                if (tempering == false && weight < 1.0) {
                    tempering = true;
                    logger.info(String.format(" Tempering activated due to recieved weight of (%8.6f)", weight));
                }

                if (resetStatistics && recursionCount[0] > lambdaResetValue) {
                    recursionKernel = new double[lambdaBins][FLambdaBins];
                    resetStatistics = false;
                    logger.info(String.format(" Cleared OSRW histogram (Lambda = %6.4f).", recursionCount[0]));
                }

                /**
                 * Increase the Recursion Kernel based on the input of current
                 * walker.
                 */
                recursionKernel[walkerLambda][walkerFLambda] += weight;
                if (this.isInterrupted()) {
                    logger.log(Level.FINE, " ReceiveThread was interrupted; ceasing execution");
                    break;
                }
            }
        }
    }
}
