/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.numerics.Potential;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;

/**
 * An implementation of the Orthogonal Space Random Walk algorithm.
 * 
 * @author Michael J. Schnieders, Wei Yang and Pengyu Ren
 */
public class OSRW implements Potential {

    private static final Logger logger = Logger.getLogger(OSRW.class.getName());
    /**
     * A potential energy that implements the LambdaInterface. 
     */
    private LambdaInterface lambdaInterface;
    /**
     * The potential energy of the system.
     */
    private Potential potential;
    /**
     * Number of variables.
     */
    private int nVariables;
    /**
     * Restart file.
     */
    private File restartFile;
    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    private double lambda;
    /**
     * Flag to indicate that the Lambda particle should be propogated.
     */
    private boolean propagateLambda = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "propagateLambda" flag true.
     */
    private int energyCount;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005).
     * The final Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     * 
     * With this scheme, the maximum of biasing Gaussians
     * is at the edges. 
     */
    private int lambdaBins = 201;
    /**
     * It is useful to have an odd number of bins, so that there is
     * a bin from FL=-dFL/2 to dFL/2 so that as FL approaches zero its
     * contribution to thermodynamic integration goes to zero. Otherwise 
     * a contribution of zero from a L bin can only result from equal
     * sampling of the ranges -dFL to 0 and 0 to dFL.
     */
    private int FLambdaBins = 401;
    /**
     * The recursion kernel stores the number of visits to 
     * each [lambda][Flambda] bin.
     */
    private int recursionKernel[][];
    /**
     * When evaluating the biasing potential, contributions from Gaussians
     * centered more the "biasCutoff" away will be neglected.
     */
    private int biasCutoff = 5;
    /**
     * Width of the lambda bin.
     */
    private double dL = 1.0 / (lambdaBins - 1);
    /**
     * Half the width of a lambda bin.
     */
    private double dL_2 = dL / 2.0;
    /**
     * The width of the F_lambda bin.
     */
    private double dFL = 2.0;
    /**
     * Half the width of the F_lambda bin.
     */
    private double dFL_2 = dFL / 2.0;
    /**
     * The minimum value of the first lambda bin.  
     */
    private double minLambda = -dL_2;
    /**
     * The minimum value of the first F_lambda bin.
     */
    private double minFLambda = -(dFL * FLambdaBins) / 2.0;
    /**
     * The maximum value of the last F_lambda bin.
     */
    private double maxFLambda = minFLambda + FLambdaBins * dFL;
    /**
     * Total partial derivative of the potential being sampled w.r.t. lambda.
     */
    private double dEdLambda;
    /**
     * 2nd partial derivative of the potential being sampled w.r.t lambda.
     */
    private double d2EdLambda2;
    private double dUdXdL[] = null;
    private double biasMag = 0.005;
    private double FLambda[];
    /**
     * Gas constant (in Kcal/mole/Kelvin).
     */
    public static final double R = 1.9872066e-3;
    private double theta;
    /**
     * Reasonable thetaFriction values are from 20 to 200.
     */
    private double thetaFriction = 1.0e-16;
    private double thetaMass = 1.0e-18;
    private double halfThetaVelocity = 0.0;
    private Random stochasticRandom;
    /**
     * Random force conversion to kcal/mol/A;
     */
    private static final double randomConvert = sqrt(4.184) / 10e9;
    private static final double randomConvert2 = randomConvert * randomConvert;
    /**
     * Time step in picoseconds.
     */
    private double dt;
    /**
     * Temperature in Kelvin.
     */
    private double temperature;
    /**
     * Lambda print frequency (in steps).
     */
    private int lambdaPrintFrequency = 1000;
    /**
     * Free energy print frequency relative to how often the free energy is
     * updated from the count matrix.
     */
    private int fLambdaPrintFrequency = 100;
    private int fLambdaUpdates = 0;
    /**
     * Steps between writing an OSRW restart file.
     */
    private int osrwRestartFrequency = 10000;
    /**
     * Print detailed energy information.
     */
    private boolean print = false;
    /**
     * total system energy
     */
    private double totalEnergy;

    public OSRW(LambdaInterface lambdaInterface, Potential potential,
            File restartFile, CompositeConfiguration properties,
            double temperature, double dt, double restartFrequency) {
        this.lambdaInterface = lambdaInterface;
        this.potential = potential;
        this.restartFile = restartFile;
        this.temperature = temperature;
        nVariables = potential.getNumberOfVariables();

        /**
         * Convert the time step to picoseconds.
         */
        this.dt = dt * 0.001;
        osrwRestartFrequency = (int) (restartFrequency / this.dt);

        biasCutoff = properties.getInt("lambda-bias-cutoff", 5);
        biasMag = properties.getDouble("bias-gaussian-mag", 0.005);
        dL = properties.getDouble("lambda-bin-width", 0.005);
        dFL = properties.getDouble("flambda-bin-width", 2.0);

        /**
         * Require modest sampling of the lambda path. 
         */
        if (dL > 0.1) {
            dL = 0.1;
        }

        /**
         * Many lambda bin widths do not evenly divide into 1.0; here we
         * correct for this by computing an integer number of bins, then
         * re-setting the lambda variable appropriately. Note that we also
         * choose to have an odd number of lambda bins, so that the centers
         * of the first and last bin are at 0 and 1.
         */
        lambdaBins = (int) (1.0 / dL);
        if (lambdaBins % 2 == 0) {
            lambdaBins++;
        }

        /**
         * The initial number of FLambda bins does not really matter, since
         * a larger number is automatically allocated as needed. The center
         * of the central bin is at 0.
         */
        FLambdaBins = 401;
        minFLambda = -(dFL * FLambdaBins) / 2.0;

        /**
         * Allocate space for the recursion kernel that stores counts.
         */
        recursionKernel = new int[lambdaBins][FLambdaBins];

        /**
         * Load the OSRW restart file if it exists. 
         */
        boolean readRestart = false;
        if (restartFile != null && restartFile.exists()) {
            try {
                OSRWRestartReader osrwRestartReader = new OSRWRestartReader(new FileReader(restartFile));
                osrwRestartReader.readRestartFile();
                logger.info(String.format("\n Continuing OSRW from %s.", restartFile.getName()));
                readRestart = true;
            } catch (FileNotFoundException ex) {
                logger.info(" Restart file could not be found and will be ignored.");
            }
        }

        dL = 1.0 / (lambdaBins - 1);
        dL_2 = dL / 2.0;
        minLambda = -dL_2;
        dFL_2 = dFL / 2.0;
        maxFLambda = minFLambda + FLambdaBins * dFL;
        FLambda = new double[lambdaBins];
        dUdXdL = new double[nVariables];
        energyCount = -1;
        stochasticRandom = new Random(0);

        /**
         * Update and print out the recursion slave.
         */
        if (readRestart) {
            updateFLambda(true);
        }

    }

    public void setPropagateLambda(boolean propagateLambda) {
        this.propagateLambda = propagateLambda;
    }

    @Override
    public double energyAndGradient(double[] x, double[] gradient) {
        double e = potential.energyAndGradient(x, gradient);

        if (propagateLambda) {
            energyCount++;
        }

        double biasEnergy = 0.0;
        dEdLambda = lambdaInterface.getdEdL();
        d2EdLambda2 = lambdaInterface.getd2EdL2();

        checkRecursionKernelSize();

        int lambdaBin = (int) floor((lambda - minLambda) / dL);
        if (lambdaBin < 0) {
            lambdaBin = 0;
        }
        if (lambdaBin >= lambdaBins) {
            lambdaBin = lambdaBins - 1;
        }
        int FLambdaBin = (int) floor((dEdLambda - minFLambda) / dFL);
        if (FLambdaBin == FLambdaBins) {
            FLambdaBin = FLambdaBins - 1;
        }
        assert (FLambdaBin < FLambdaBins);
        assert (FLambdaBin >= 0);

        /**
         * Calculate recursion kernel G(L, FL) and its gradient.
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
                 * If either of the following FL edge conditions are true,
                 * then there are no counts and we continue.
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
         * Lambda gradient due to recursion kernel G(L, F(L)).
         */
        double dEdU = dEdLambda;
        dEdLambda += dGdLambda + dGdFLambda * d2EdLambda2;

        /**
         * Gradient due to recursion kernel G(L, F(L)).
         */
        for (int i = 0; i < nVariables; i++) {
            dUdXdL[i] = 0.0;
        }
        lambdaInterface.getdEdXdL(dUdXdL);
        for (int i = 0; i < nVariables; i++) {
            gradient[i] += dGdFLambda * dUdXdL[i];
        }

        /**
         * Update free energy F(L) every ~100 steps.
         */
        if (energyCount > 0 && energyCount % 100 == 0 && propagateLambda) {
            fLambdaUpdates++;
            boolean printFLambda = fLambdaUpdates % fLambdaPrintFrequency == 0;
            updateFLambda(printFLambda);
        }

        if (energyCount > 0 && energyCount % osrwRestartFrequency == 0) {
            try {
                OSRWRestartWriter osrwRestart = new OSRWRestartWriter(new BufferedWriter(new FileWriter(restartFile)));
                osrwRestart.writeRestartFile();
                osrwRestart.flush();
                osrwRestart.close();
                logger.info(String.format(" Wrote OSRW restart file to %s.", restartFile.getName()));
            } catch (IOException ex) {
                String message = " Exception writing OSRW restart file.";
                logger.log(Level.INFO, message, ex);
            }
        }

        /**
         * Compute the energy and gradient for the recursion slave at F(L)
         * using interpolation.
         */
        biasEnergy += computeFreeEnergy();

        if (print) {
            logger.info(String.format(" %s %16.8f", "Bias Energy       ", biasEnergy));
            logger.info(String.format(" %s %16.8f  %s",
                    "OSRW Potential    ", e + biasEnergy, "(Kcal/mole)"));
        }

        /**
         * Log our current state.
         */
        if (energyCount > 0 && energyCount % lambdaPrintFrequency == 0 && propagateLambda) {
            if (lambdaBins < 1000) {
                logger.info(String.format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f",
                        lambda, lambdaBin, dEdU, dEdLambda - dEdU, dEdLambda));
            } else {
                logger.info(String.format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f",
                        lambda, lambdaBin, dEdU, dEdLambda - dEdU, dEdLambda));
            }
        }

        /**
         * Meta-dynamics grid counts (every ~10 steps).
         */
        if (energyCount > 0 && energyCount % 10 == 0 && propagateLambda) {
            recursionKernel[lambdaBin][FLambdaBin]++;
        }

        /**
         * Propagate the Lambda particle.
         */
        if (propagateLambda) {
            langevin();
        }

        totalEnergy = e + biasEnergy;
        return e + biasEnergy;
    }

    public double getTotaldEdLambda() {
        return dEdLambda;
    }

    /**
     * If necessary, allocate more space.
     */
    private void checkRecursionKernelSize() {
        if (dEdLambda > maxFLambda) {
            logger.info(String.format(" Current F_lambda %8.2f > maximum historgram size %8.2f.",
                    dEdLambda, maxFLambda));

            double origDeltaG = updateFLambda(false);

            int newFLambdaBins = FLambdaBins;
            while (minFLambda + newFLambdaBins * dFL < dEdLambda) {
                newFLambdaBins += 100;
            }
            int newRecursionKernel[][] = new int[lambdaBins][newFLambdaBins];
            /**
             * We have added bins above the indeces of the current counts
             * just copy them into the new array.
             */
            for (int i = 0; i < lambdaBins; i++) {
                for (int j = 0; j < FLambdaBins; j++) {
                    newRecursionKernel[i][j] = recursionKernel[i][j];
                }
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
            int newRecursionKernel[][] = new int[lambdaBins][newFLambdaBins];
            /**
             * We have added bins below the current counts,
             * so their indeces must be increased by:
             * offset = newFLBins - FLBins
             */
            for (int i = 0; i < lambdaBins; i++) {
                for (int j = 0; j < FLambdaBins; j++) {
                    newRecursionKernel[i][j + offset] = recursionKernel[i][j];
                }
            }
            recursionKernel = newRecursionKernel;
            minFLambda = minFLambda - offset * dFL;
            FLambdaBins = newFLambdaBins;
            logger.info(String.format(" New historgram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            assert (origDeltaG == updateFLambda(false));
        }
    }

    private double updateFLambda(boolean print) {
        double freeEnergy = 0.0;
        if (print) {
            logger.info(" Count  Lambda Bins    F_Lambda Bins   <   F_L  >       dG        G");
        }
        for (int iL = 0; iL < lambdaBins; iL++) {
            int ulFL = -1;
            int llFL = -1;

            // Find the smallest FL bin.
            for (int jFL = 0; jFL < FLambdaBins; jFL++) {
                int count = recursionKernel[iL][jFL];
                if (count > 0) {
                    llFL = jFL;
                    break;
                }
            }

            // Find the largest FL bin.
            for (int jFL = FLambdaBins - 1; jFL >= 0; jFL--) {
                int count = recursionKernel[iL][jFL];
                if (count > 0) {
                    ulFL = jFL;
                    break;
                }
            }

            int lambdaCount = 0;
            // The FL range sampled for lambda bin [iL*dL .. (iL+1)*dL]
            double lla = 0.0;
            double ula = 0.0;
            if (ulFL == -1) {
                FLambda[iL] = 0.0;
            } else {
                double ensembleAverageFLambda = 0.0;
                double partitionFunction = 0.0;
                for (int jFL = llFL; jFL <= ulFL; jFL++) {
                    double currentFLambda = minFLambda + jFL * dFL + dFL_2;
                    double weight = exp(evaluateKernel(iL, jFL) / (R * temperature));
                    ensembleAverageFLambda += currentFLambda * weight;
                    partitionFunction += weight;
                    lambdaCount += recursionKernel[iL][jFL];
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

            if (print) {
                double llL = iL * dL - dL_2;
                double ulL = llL + dL;
                if (llL < 0.0) {
                    llL = 0.0;
                }
                if (ulL > 1.0) {
                    ulL = 1.0;
                }
                logger.info(String.format(" %5d [%5.3f %5.3f] [%7.1f %7.1f] <%8.3f> %8.3f %8.3f",
                        lambdaCount, llL, ulL, lla, ula,
                        FLambda[iL], deltaFreeEnergy, freeEnergy));
            }
        }

        return freeEnergy;
    }

    private double computeFreeEnergy() {
        double biasEnergy = 0.0;
        for (int iL0 = 0; iL0 < lambdaBins - 1; iL0++) {
            int iL1 = iL0 + 1;
            /**
             * Find bin centers and values for
             * interpolation / extrapolation points.
             */
            double L0 = iL0 * dL;
            double L1 = L0 + dL;
            double FL0 = FLambda[iL0];
            double FL1 = FLambda[iL1];
            double deltaFL = FL1 - FL0;
            /**
             * If the lambda is less than or equal to the upper limit,
             * this is the final interval. Set the upper limit to L,
             * compute the partial derivative and break.
             */
            boolean done = false;
            if (lambda <= L1) {
                done = true;
                L1 = lambda;
            }
            /**
             * Upper limit - lower limit of the integral of the
             * extrapolation / interpolation.
             */
            biasEnergy += (FL0 * L1 + deltaFL * L1 * (0.5 * L1 - L0) / dL);
            biasEnergy -= (FL0 * L0 + deltaFL * L0 * (-0.5 * L0) / dL);
            if (done) {
                /**
                 * Compute the gradient d F(L) / dL at L.
                 */
                dEdLambda -= FL0 + (L1 - L0) * deltaFL / dL;
                break;
            }
        }
        return -biasEnergy;
    }

    public double evaluateKernel(int cLambda, int cF_Lambda) {
        /**
         * Compute the value of L and FL for the
         * center of the current bin.
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
                 * The width of the first and last bins is dLambda_2,
                 * so the mirror condition is to double their counts.
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
                 * For FLambda outside the count matrix the weight is
                 * 0 so we continue.
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

    public void setLambda(double lambda) {
        lambdaInterface.setLambda(lambda);
        this.lambda = lambda;
        theta = Math.asin(Math.sqrt(lambda));
    }

    public void setThetaMass(double thetaMass) {
        this.thetaMass = thetaMass;
    }

    public void setThetaFrication(double thetaFriction) {
        this.thetaFriction = thetaFriction;
    }

    /**
     * Propagate Lambda using Langevin dynamics.
     */
    private void langevin() {
        double rt2 = 2.0 * Thermostat.R * temperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / randomConvert;
        double dEdL = -dEdLambda * sin(2.0 * theta);
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

    @Override
    public double[] getMass() {
        return potential.getMass();
    }
    
    @Override
    public double getTotal() {
        return totalEnergy;
    }

    @Override
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
    }

    private class OSRWRestartWriter extends PrintWriter {

        public OSRWRestartWriter(Writer writer) {
            super(writer);
        }

        public void writeRestartFile() {
            printf("Temperature     %15.3f\n", temperature);
            printf("Lambda          %15.8f\n", lambda);
            printf("Lambda-Mass     %15.8e\n", thetaMass);
            printf("Lambda-Friction %15.8e\n", thetaFriction);
            printf("Lambda-Velocity %15.8e\n", halfThetaVelocity);
            printf("Bias-Mag        %15.8e\n", biasMag);
            printf("Bias-Cutoff     %15d\n", biasCutoff);
            printf("Lambda-Bins     %15d\n", lambdaBins);
            printf("FLambda-Bins    %15d\n", FLambdaBins);
            printf("Flambda-Min     %15.8e\n", minFLambda);
            printf("Flambda-Width   %15.8e\n", dFL);
            for (int i = 0; i < lambdaBins; i++) {
                printf("%d", recursionKernel[i][0]);
                for (int j = 1; j < FLambdaBins; j++) {
                    printf(" %d", recursionKernel[i][j]);
                }
                println();
            }
        }
    }

    private class OSRWRestartReader extends BufferedReader {

        public OSRWRestartReader(Reader reader) {
            super(reader);
        }

        public void readRestartFile() {
            try {
                temperature = Double.parseDouble(readLine().split(" +")[1]);
                lambda = Double.parseDouble(readLine().split(" +")[1]);
                setLambda(lambda);
                thetaMass = Double.parseDouble(readLine().split(" +")[1]);
                thetaFriction = Double.parseDouble(readLine().split(" +")[1]);
                halfThetaVelocity = Double.parseDouble(readLine().split(" +")[1]);
                biasMag = Double.parseDouble(readLine().split(" +")[1]);
                biasCutoff = Integer.parseInt(readLine().split(" +")[1]);
                lambdaBins = Integer.parseInt(readLine().split(" +")[1]);
                FLambdaBins = Integer.parseInt(readLine().split(" +")[1]);
                minFLambda = Double.parseDouble(readLine().split(" +")[1]);
                dFL = Double.parseDouble(readLine().split(" +")[1]);
                // Allocate memory for the recursion kernel.
                recursionKernel = new int[lambdaBins][FLambdaBins];
                for (int i = 0; i < lambdaBins; i++) {
                    String counts[] = readLine().split(" +");
                    for (int j = 0; j < FLambdaBins; j++) {
                        recursionKernel[i][j] = Integer.parseInt(counts[j]);
                    }
                }
            } catch (Exception e) {
                String message = " Invalid OSRW restart file.";
                logger.log(Level.SEVERE, message, e);
            }
        }
    }
}
