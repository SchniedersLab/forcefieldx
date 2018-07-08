/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.Comm;
import edu.rit.pj.cluster.JobBackend;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;

/**
 * An implementation of the Orthogonal Space Random Walk algorithm.
 *
 * @author Michael J. Schnieders, Wei Yang and Pengyu Ren
 */
public abstract class AbstractOSRW implements CrystalPotential {

    private static final Logger logger = Logger.getLogger(AbstractOSRW.class.getName());
    /**
     * A potential energy that implements the LambdaInterface.
     */
    protected final LambdaInterface lambdaInterface;
    /**
     * The potential energy of the system.
     */
    protected final CrystalPotential potential;
    /**
     * Reference to the Barostat in use; if present this must be turned off
     * during optimization.
     */
    protected final Barostat barostat;
    /**
     * The AlgorithmListener is called each time a count is added.
     */
    protected final AlgorithmListener algorithmListener;
    /**
     * Number of variables.
     */
    protected final int nVariables;
    /**
     * Each walker has a unique lambda restart file.
     */
    protected final File lambdaFile;
    /**
     * Each walker reads the same histogram restart file. Only the walker of
     * rank 0 writes the histogram restart file.
     */
    protected final File histogramFile;
    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    protected double lambda;
    /**
     * Flag to indicate that the Lambda particle should be propagated.
     */
    protected boolean propagateLambda = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "propagateLambda" flag true.
     */
    protected int energyCount;
    protected int biasCount = 0;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005). The final
     * Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     *
     * With this scheme, the maximum of biasing Gaussians is at the edges.
     */
    protected int lambdaBins = 201;
    /**
     * It is useful to have an odd number of bins, so that there is a bin from
     * FL=-dFL/2 to dFL/2 so that as FL approaches zero its contribution to
     * thermodynamic integration goes to zero. Otherwise a contribution of zero
     * from a L bin can only result from equal sampling of the ranges -dFL to 0
     * and 0 to dFL.
     */
    protected int FLambdaBins = 401;
    /**
     * Parallel Java world communicator.
     */
    protected final Comm world;
    /**
     * Number of processes.
     */
    protected final int numProc;
    /**
     * Rank of this process.
     */
    protected final int rank;
    /**
     * A reference to the JobBackend to log messages for the Webserver.
     */
    protected final JobBackend jobBackend;
    /**
     * When evaluating the biasing potential, contributions from Gaussians
     * centered more the "biasCutoff" away will be neglected.
     */
    protected int biasCutoff = 5;
    /**
     * Width of the lambda bin.
     */
    protected double dL = 1.0 / (lambdaBins - 1);
    /**
     * Half the width of a lambda bin.
     */
    protected double dL_2 = dL / 2.0;
    /**
     * The width of the F_lambda bin.
     */
    protected double dFL = 2.0;
    /**
     * Half the width of the F_lambda bin.
     */
    protected double dFL_2 = dFL / 2.0;
    /**
     * The minimum value of the first lambda bin.
     */
    protected double minLambda = -dL_2;
    /**
     * The minimum value of the first F_lambda bin.
     */
    protected double minFLambda = -(dFL * FLambdaBins) / 2.0;
    /**
     * The maximum value of the last F_lambda bin.
     */
    protected double maxFLambda = minFLambda + FLambdaBins * dFL;
    /**
     * Atom gradient for use if "energy" is called.
     */
    private double grad[] = null;
    /**
     * Force Field Potential Energy (i.e. with no bias terms added).
     */
    protected double forceFieldEnergy;
    /**
     * Partial derivative of the force field energy w.r.t. lambda.
     */
    protected double dForceFieldEnergydL;
    /**
     * OSRW Bias energy.
     */
    protected double biasEnergy;
    /**
     * Total partial derivative of the potential (U) being sampled w.r.t.
     * lambda.
     */
    protected double dUdLambda;
    /**
     * 2nd partial derivative of the potential being sampled w.r.t lambda.
     */
    protected double d2UdL2;
    /**
     * Mixed partial derivative with respect to coordinates and lambda.
     */
    protected double dUdXdL[] = null;
    /**
     * Magnitude of each hill (not including tempering) in (kcal/mol).
     */
    protected double biasMag = 0.05;
    /**
     * 1D PMF with respect to lambda F(L).
     */
    protected double FLambda[];
    /**
     * Magnitude of the 2D orthogonal space bias G(L,dE/dL).
     */
    protected double gLdEdL = 0.0;

    /**
     * Gas constant (in Kcal/mol/Kelvin).
     */
    public static final double R = 1.9872066e-3;
    private double theta;
    /**
     * Reasonable thetaFriction is 60 1/ps.
     */
    protected double thetaFriction = 1.0e-19;
    /**
     * Reasonable thetaMass is 100 a.m.u.
     */
    protected double thetaMass = 1.0e-18;
    protected double halfThetaVelocity = 0.0;
    private final Random stochasticRandom;
    /**
     * Random force conversion to kcal/mol/A;
     */
    private static final double randomConvert = sqrt(4.184) / 10e9;
    //private static final double invRandomConvert = 1.0 / randomConvert;
    /**
     * randomConvert squared.
     */
    private static final double randomConvert2 = randomConvert * randomConvert;
    /**
     * Time step in picoseconds.
     */
    protected final double dt;
    /**
     * Temperature in Kelvin.
     */
    protected double temperature = 298.15;
    /**
     * Interval between adding a count to the Recursion kernel in steps.
     */
    protected int countInterval = 10;
    /**
     * Interval between printing information on the lambda particle in steps.
     */
    protected int printFrequency = 100;
    /**
     * Whether to write out lambda-traversal snapshots.
     */
    protected boolean writeTraversalSnapshots = false;
    /**
     * Keeps track of full lambda-traversals for snapshot output.
     */
    protected int traversalSnapshotTarget = -1;
    /**
     * Ensemble files containing full-traversal snapshots, plus an assembly,
     * filter, and counter for each.
     */
    protected File lambdaOneFile;
    protected int lambdaOneStructures = 0;
    protected MolecularAssembly lambdaOneAssembly;
    protected PDBFilter lambdaOneFilter;
    protected File lambdaZeroFile;
    protected int lambdaZeroStructures = 0;
    protected MolecularAssembly lambdaZeroAssembly;
    protected PDBFilter lambdaZeroFilter;
    /**
     * Once the lambda reset value is reached, Transition-Tempered OSRW
     * statistics are reset.
     */
    protected final double lambdaResetValue = 0.99;
    /**
     * Flag set to false once Transition-Tempered OSRW statistics are reset at
     * lambdaResetValue.
     */
    protected boolean resetStatistics = false;
    /**
     * Stores a traversal snapshot that has not yet been written to file.
     */
    protected ArrayList<String> traversalInHand = new ArrayList<>();
    /**
     * Holds the lowest potential-energy parameters from among all visits to
     * lambda < 0.1 and > 0.9
     */
    private double lowEnergyCoordsZero[], lowEnergyCoordsOne[];
    private double lowEnergyZero = Double.MAX_VALUE, lowEnergyOne = Double.MAX_VALUE;

    /**
     * Holds the lowest potential-energy parameters for loopBuilder runs from
     * all visits to lambda > 0.9
     */
    protected double osrwOptimumCoords[];
    protected double osrwOptimum = Double.MAX_VALUE;

    /**
     * Interval between how often the 1D histogram is printed to screen versus
     * silently updated in background.
     */
    protected final int fLambdaPrintInterval = 25;
    protected int fLambdaUpdates = 0;
    /**
     * Interval between writing an OSRW restart file in steps.
     */
    protected int saveFrequency = 1000;
    /**
     * Print detailed energy information.
     */
    protected final boolean print = false;
    /**
     * Total system energy.
     */
    protected double totalEnergy;
    /**
     * Thermodynamic integration from Lambda=0 to Lambda=1.
     */
    protected double totalFreeEnergy;
    /**
     * Save the previous free energy, in order to limit logging to time points
     * where the free energy has changed.
     */
    protected double previousFreeEnergy = 0.0;
    protected double lastAverage = 0.0;
    protected double lastStdDev = 0.0;

    /**
     * Equilibration counts
     */
    protected int equilibrationCounts = 0;
    /**
     * Are FAST varying energy terms being computed, SLOW varying energy terms,
     * or BOTH. OSRW is not active when only FAST varying energy terms are being
     * propagated.
     */
    protected STATE state = STATE.BOTH;
    /**
     * Flag to indicate if OSRW should send and receive counts between processes
     * synchronously or asynchronously. The latter is faster by ~40% because
     * simulation with Lambda > 0.75 must compute two condensed phase
     * self-consistent fields to interpolate polarization.
     */
    protected final boolean asynchronous;

    /**
     * Running average and standard deviation
     */
    protected double totalAverage = 0;
    protected double totalSquare = 0;
    protected int periodCount = 0;
    protected int window = 1000;

    protected boolean osrwOptimization = false;
    protected int osrwOptimizationFrequency = 10000;
    protected double osrwOptimizationLambdaCutoff = 0.5;
    protected double osrwOptimizationEps = 0.1;
    protected double osrwOptimizationTolerance = 1.0e-8;
    protected double osrwOptimizationEnergyWindow = 2.0;

    protected MolecularAssembly molecularAssembly;
    protected SystemFilter systemFilter;
    protected File optFile;

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
     * @param checkpointInterval number of steps between restart file updates.
     * @param asynchronous set to true if walkers run asynchronously.
     * @param algorithmListener the AlgorithmListener to be notified of
     * progress.
     */
    public AbstractOSRW(LambdaInterface lambdaInterface, CrystalPotential potential,
            File lambdaFile, File histogramFile, CompositeConfiguration properties,
            double temperature, double dt, double printInterval,
            double checkpointInterval, boolean asynchronous,
            AlgorithmListener algorithmListener) {
        this(lambdaInterface, potential, lambdaFile, histogramFile, properties,
                temperature, dt, printInterval, checkpointInterval, asynchronous,
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
     * @param checkpointInterval number of steps between restart file updates.
     * @param asynchronous set to true if walkers run asynchronously.
     * @param resetNumSteps whether to reset energy counts to 0
     * @param algorithmListener the AlgorithmListener to be notified of
     * progress.
     */
    public AbstractOSRW(LambdaInterface lambdaInterface, CrystalPotential potential,
            File lambdaFile, File histogramFile, CompositeConfiguration properties,
            double temperature, double dt, double printInterval,
            double checkpointInterval, boolean asynchronous, boolean resetNumSteps,
            AlgorithmListener algorithmListener) {
        this.lambdaInterface = lambdaInterface;
        this.potential = potential;
        this.lambdaFile = lambdaFile;
        this.histogramFile = histogramFile;
        this.temperature = temperature;
        this.asynchronous = asynchronous;
        this.algorithmListener = algorithmListener;

        nVariables = potential.getNumberOfVariables();
        lowEnergyCoordsZero = new double[nVariables];
        lowEnergyCoordsOne = new double[nVariables];

        if (potential instanceof Barostat) {
            barostat = (Barostat) potential;
        } else {
            barostat = null;
        }

        /**
         * Convert the time step to picoseconds.
         */
        this.dt = dt * 0.001;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        /**
         * Convert the save interval to a save frequency.
         */
        saveFrequency = 1000;
        if (checkpointInterval >= this.dt) {
            saveFrequency = (int) (checkpointInterval / this.dt);
        }

        biasCutoff = properties.getInt("lambda-bias-cutoff", 5);
        biasMag = properties.getDouble("bias-gaussian-mag", 0.050);
        dL = properties.getDouble("lambda-bin-width", 0.005);
        dFL = properties.getDouble("flambda-bin-width", 2.0);

        /**
         * Require modest sampling of the lambda path.
         */
        if (dL > 0.1) {
            dL = 0.1;
        }

        /**
         * Many lambda bin widths do not evenly divide into 1.0; here we correct
         * for this by computing an integer number of bins, then re-setting the
         * lambda variable appropriately. Note that we also choose to have an
         * odd number of lambda bins, so that the centers of the first and last
         * bin are at 0 and 1.
         */
        lambdaBins = (int) (1.0 / dL);
        if (lambdaBins % 2 == 0) {
            lambdaBins++;
        }

        /**
         * The initial number of FLambda bins does not really matter, since a
         * larger number is automatically allocated as needed. The center of the
         * central bin is at 0.
         */
        FLambdaBins = 401;
        minFLambda = -(dFL * FLambdaBins) / 2.0;

        energyCount = -1;

        dL = 1.0 / (lambdaBins - 1);
        dL_2 = dL / 2.0;
        minLambda = -dL_2;
        dFL_2 = dFL / 2.0;
        maxFLambda = minFLambda + FLambdaBins * dFL;
        FLambda = new double[lambdaBins];
        dUdXdL = new double[nVariables];
        stochasticRandom = new Random();

        /**
         * Set up the multi-walker communication variables for Parallel Java
         * communication between nodes.
         */
        world = Comm.world();
        numProc = world.size();
        rank = world.rank();
        jobBackend = JobBackend.getJobBackend();

        /**
         * Log OSRW parameters.
         */
        logger.info("\n Orthogonal Space Random Walk Parameters");
        logger.info(String.format("  Gaussian Bias Magnitude:       %6.4f (kcal/mol)", biasMag));
        logger.info(String.format("  Gaussian Bias Cutoff:           %6d bins", biasCutoff));
        logger.info(String.format("  Print Interval:                 %6.3f psec", printInterval));
        logger.info(String.format("  Save Interval:                  %6.3f psec", checkpointInterval));
    }

    public void setPropagateLambda(boolean propagateLambda) {
        this.propagateLambda = propagateLambda;
    }

    public abstract void addBias(double dUdL, double[] x, double[] gradient);

    protected int binForLambda(double lambda) {
        int lambdaBin = (int) floor((lambda - minLambda) / dL);
        if (lambdaBin < 0) {
            lambdaBin = 0;
        }
        if (lambdaBin >= lambdaBins) {
            lambdaBin = lambdaBins - 1;
        }
        return lambdaBin;
    }

    protected int binForFLambda(double dEdLambda) {
        int FLambdaBin = (int) floor((dEdLambda - minFLambda) / dFL);
        if (FLambdaBin == FLambdaBins) {
            FLambdaBin = FLambdaBins - 1;
        }
        assert (FLambdaBin < FLambdaBins);
        assert (FLambdaBin >= 0);
        return FLambdaBin;
    }

    public double getForceFielddEdL() {
        return dForceFieldEnergydL;
    }

    public double getTotaldEdLambda() {
        return dUdLambda;
    }

    public double getForceFieldEnergy() {
        return forceFieldEnergy;
    }

    public double getBiasEnergy() {
        return biasEnergy;
    }

    public Potential getPotentialEnergy() {
        return potential;
    }

    public double getGofLdUdL() {
        return gLdEdL;
    }

    protected double currentFreeEnergy() {
        double biasEnergy = 0.0;
        for (int iL0 = 0; iL0 < lambdaBins - 1; iL0++) {
            int iL1 = iL0 + 1;
            /**
             * Find bin centers and values for interpolation / extrapolation
             * points.
             */
            double L0 = iL0 * dL;
            double L1 = L0 + dL;
            double FL0 = FLambda[iL0];
            double FL1 = FLambda[iL1];
            double deltaFL = FL1 - FL0;
            /**
             * If the lambda is less than or equal to the upper limit, this is
             * the final interval. Set the upper limit to L, compute the partial
             * derivative and break.
             */
            boolean done = false;
            if (lambda <= L1) {
                done = true;
                L1 = lambda;
            }
            /**
             * Upper limit - lower limit of the integral of the extrapolation /
             * interpolation.
             */
            biasEnergy += (FL0 * L1 + deltaFL * L1 * (0.5 * L1 - L0) / dL);
            biasEnergy -= (FL0 * L0 + deltaFL * L0 * (-0.5 * L0) / dL);
            if (done) {
                /**
                 * Compute the gradient d F(L) / dL at L.
                 */
                dUdLambda -= FL0 + (L1 - L0) * deltaFL / dL;
                break;
            }
        }
        return -biasEnergy;
    }

    public double lastFreeEnergy() {
        return totalFreeEnergy;
    }

    public double movingAverageEnergy() {
        return lastAverage;
    }

    public double movingAverageSD() {
        return lastStdDev;
    }

    public void evaluatePMF() {
        StringBuffer sb = new StringBuffer();
        for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
            for (int fLambdaBin = 0; fLambdaBin < FLambdaBins; fLambdaBin++) {
                sb.append(String.format(" %16.8f", evaluateKernel(lambdaBin, fLambdaBin)));
            }
            sb.append("\n");
        }
        logger.info(sb.toString());
    }

    /**
     * Shuts down resources associated with this OSRW, primarily the receive
     * thread.
     *
     * @return Success
     */
    public abstract boolean destroy();

    protected abstract double evaluateKernel(int cLambda, int cF_Lambda);

    /**
     * Evaluates current free energy of the OSRW; intended to be called before
     * any dynamics have been run.
     *
     * @return
     */
    public double evaluateEnergy() {
        return updateFLambda(false);
    }

    protected abstract double updateFLambda(boolean print);

    public void setLambda(double lambda) {
        lambdaInterface.setLambda(lambda);
        this.lambda = lambda;
        theta = asin(sqrt(lambda));
    }

    public double getLambda() {
        return lambda;
    }

    public LambdaInterface getLambdaInterface() {
        return lambdaInterface;
    }

    public void setTraversalOutput(File lambdaOneFile, MolecularAssembly topology1, File lambdaZeroFile, MolecularAssembly topology2) {
        this.writeTraversalSnapshots = true;
        this.lambdaOneFile = lambdaOneFile;
        this.lambdaOneAssembly = topology1;
        this.lambdaZeroFile = lambdaZeroFile;
        this.lambdaZeroAssembly = topology2;
    }

    public void setThetaMass(double thetaMass) {
        this.thetaMass = thetaMass;
    }

    public void setResetStatistics(boolean resetStatistics) {
        this.resetStatistics = resetStatistics;
    }

    public void setThetaFrication(double thetaFriction) {
        this.thetaFriction = thetaFriction;
    }

    /**
     * Set the OSRW Gaussian biasing potential magnitude (kcal/mol).
     *
     * @param biasMag Gaussian biasing potential magnitude (kcal/mol)
     */
    public void setBiasMagnitude(double biasMag) {
        this.biasMag = biasMag;
        logger.info(String.format(" Gaussian Bias Magnitude: %6.4f (kcal/mol)", biasMag));
    }

    /**
     * Set the OSRW count interval. Every 'countInterval' steps the
     * recursionKernel will be incremented based on the current value of the
     * lambda state variable and the derivative of the energy with respect to
     * lambda (dU/dL).
     *
     * @param countInterval Molecular dynamics steps between counts.
     */
    public void setCountInterval(int countInterval) {
        if (countInterval > 0) {
            this.countInterval = countInterval;
        } else {
            logger.info(" OSRW count interval must be greater than 0.");
        }
    }

    /**
     * Propagate Lambda using Langevin dynamics.
     */
    protected void langevin() {
        double rt2 = 2.0 * ffx.algorithms.thermostats.Thermostat.R * temperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / randomConvert;
        //double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() * invRandomConvert;
        double dEdL = -dUdLambda * sin(2.0 * theta);
        halfThetaVelocity = (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
                + randomConvert2 * 2.0 * dt * (dEdL + randomForce))
                / (2.0 * thetaMass + thetaFriction * dt);
        theta = theta + dt * halfThetaVelocity;

        if (theta > PI) {
            theta -= 2.0 * PI;
        } else if (theta <= -PI) {
            theta += 2.0 * PI;
        }

        double sinTheta = sin(theta);
        lambda = sinTheta * sinTheta;
        lambdaInterface.setLambda(lambda);
    }

    @Override
    public void setScaling(double[] scaling) {
        potential.setScaling(scaling);
    }

    @Override
    public double[] getScaling() {
        return potential.getScaling();
    }

    @Override
    public double[] getCoordinates(double[] doubles) {
        return potential.getCoordinates(doubles);
    }

    /**
     * Return a copy of the parameter array containing lowest-energy parameters
     * from amongst visits to the specified end state (either 0 or 1).
     *
     * @param endState
     *
     * @return a double array of parameters
     */
    public double[] getLowEnergyCoordinates(int endState) {
        if (endState == 0) {
            return lowEnergyCoordsZero;
        } else if (endState == 1) {
            return lowEnergyCoordsOne;
        } else {
            logger.severe("Improper function call.");
            return null;
        }
    }

    public void setOSRWOptimum(double prevOSRWOptimum) {
        osrwOptimum = prevOSRWOptimum;
    }

    public double getOSRWOptimum() {
        return osrwOptimum;
    }

    public double[] getLowEnergyLoop() {
        if (osrwOptimum < Double.MAX_VALUE) {
            return osrwOptimumCoords;
        } else {
            logger.info("Lambda > 0.9 was not reached. Try increasing number of timesteps.");
            return null;
        }
    }

    public void setOptimization(boolean osrwOptimization, MolecularAssembly molAss) {
        this.osrwOptimization = osrwOptimization;
        this.molecularAssembly = molAss;
        File file = molecularAssembly.getFile();

        String fileName = FilenameUtils.removeExtension(file.getAbsolutePath());
        String ext = FilenameUtils.getExtension(file.getAbsolutePath());

        if (systemFilter == null) {
            if (ext.toUpperCase().contains("XYZ")) {
                optFile = new File(fileName + "_opt.xyz");
                systemFilter = new XYZFilter(optFile, molecularAssembly, null, null);
            } else {
                optFile = new File(fileName + "_opt.pdb");
                systemFilter = new PDBFilter(optFile, molecularAssembly, null, null);
            }
        }
    }

    /**
     * Returns the number of energy evaluations performed by this ttOSRW,
     * including those picked up in the lambda file.
     *
     * @return Number of energy steps taken by this walker.
     */
    public int getEnergyCount() {
        return energyCount;
    }

    @Override
    public double[] getMass() {
        return potential.getMass();
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return potential.getVariableTypes();
    }

    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    @Override
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
    }

    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        potential.setEnergyTermState(state);
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public double energy(double[] x) {

        boolean origPropagateLambda = propagateLambda;

        propagateLambda = false;

        if (grad == null || grad.length != x.length) {
            grad = new double[x.length];
        }

        double energy = energyAndGradient(x, grad);

        propagateLambda = origPropagateLambda;

        return energy;
    }

    @Override
    public void setVelocity(double[] velocity) {
        potential.setVelocity(velocity);
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        potential.setAcceleration(acceleration);
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        potential.setPreviousAcceleration(previousAcceleration);
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        return potential.getVelocity(velocity);
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        return potential.getAcceleration(acceleration);
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        return potential.getPreviousAcceleration(previousAcceleration);
    }

    @Override
    public void setCrystal(Crystal crystal) {
        potential.setCrystal(crystal);
    }

    @Override
    public Crystal getCrystal() {
        return potential.getCrystal();
    }
}
