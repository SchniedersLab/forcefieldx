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

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import java.util.Random;
import java.util.logging.Logger;

import ffx.numerics.Potential;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import org.apache.commons.configuration.CompositeConfiguration;

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
     * Array of atoms.
     */
    private Atom atoms[];
    /**
     * Number of atoms.
     */
    private int nAtoms;
    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    private double lambda;
    /**
     * Flag to indicate that the number of energy evaluations is currently 
     * being counted.
     */
    private boolean propagateLambda = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "doCounting" flag true.
     */
    private int energyCount;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005).
     * The final Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     * 
     * With this scheme, the maximum of biasing Gaussians
     * is at the edges. 
     */
    private int lambdaBins = 101;
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
    private double dL = 0.01;
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
    private double minLambda = -0.005;
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
    private double dEdLambda = 0.0;
    /**
     * 2nd partial derivative of the potential being sampled w.r.t lambda.
     */
    private double d2EdLambda2 = 0.0;
    private double dUdXdL[] = null;
    private double biasGaussianMag = 0.005;
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
    private boolean print = false;

    public OSRW(LambdaInterface lambdaInterface,
            Potential potential,
            CompositeConfiguration properties,
            Atom atoms[],
            double temperature, double timeStep){
        this.lambdaInterface = lambdaInterface;
        this.potential = potential;
        this.atoms = atoms;
        nAtoms = atoms.length;
        this.temperature = temperature;

        /**
         * Convert the time step to picoseconds.
         */
        this.dt = timeStep * 0.001;

        biasCutoff = properties.getInt("lambda-bias-cutoff", 5);
        biasGaussianMag = properties.getDouble("bias-gaussian-mag", 0.005);
        dL = properties.getDouble("lambda-bin-width", 0.01);

        /**
         * Require modest sampling of the lambda path. 
         */
        if (dL > 0.1) {
            dL = 0.01;
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
        dL = 1.0 / (lambdaBins - 1);
        dL_2 = dL / 2.0;
        minLambda = -dL_2;

        /**
         * The initial number of FLambda bins does not really matter, since
         * a larger number are automatically allocated as needed. The center
         * of the central bin is at 0.
         */
        dFL = properties.getDouble("flambda-bin-width", 2.0);
        dFL_2 = dFL / 2.0;
        FLambdaBins = 401;
        minFLambda = -(dFL * FLambdaBins) / 2.0;
        maxFLambda = minFLambda + FLambdaBins * dFL;

        /**
         * Allocate space for the recursion kernel that stores counts.
         */
        recursionKernel = new int[lambdaBins][FLambdaBins];
        FLambda = new double[lambdaBins];
        dUdXdL = new double[nAtoms * 3];
        energyCount = 0;

        stochasticRandom = new Random(0);

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

        int stateBin = (int) floor((lambda - minLambda) / dL);
        if (stateBin < 0) {
            stateBin = 0;
        }
        if (stateBin >= lambdaBins) {
            stateBin = lambdaBins - 1;
        }
        int FStateBin = (int) floor((dEdLambda - minFLambda) / dFL);
        if (FStateBin == FLambdaBins) {
            FStateBin = FLambdaBins - 1;
        }
        assert (FStateBin < FLambdaBins);
        assert (FStateBin >= 0);

        /**
         * Calculate recursion kernel G(L, dEdL) and gradient.
         */
        double dGdLambda = 0.0;
        double dGdFLambda = 0.0;
        double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
        double FLs2 = dFL * 2.0 * dFL * 2.0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int lcenter = stateBin + iL;
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
                int FLcenter = FStateBin + iFL;
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
                double bias = weight * biasGaussianMag
                        * exp(-deltaL2 / (2.0 * ls2))
                        * exp(-deltaFL2 / (2.0 * FLs2));
                biasEnergy += bias;
                dGdLambda -= deltaL / ls2 * bias;
                dGdFLambda -= deltaFL / FLs2 * bias;
            }
        }

        /**
         * Lambda gradient due to recursion kernel G(L, dEdL).
         */
        dEdLambda += dGdLambda + dGdFLambda * d2EdLambda2;
        
        /**
         * Atomic gradient due to recursion kernel G(L, dEdL).
         */
        for (int i = 0; i < 3 * nAtoms; i++) {
            dUdXdL[i] = 0.0;
        }

        lambdaInterface.getdEdXdL(dUdXdL);

        double grad[] = new double[3];
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            atom.getXYZGradient(grad);
            grad[0] += dGdFLambda * dUdXdL[i * 3];
            grad[1] += dGdFLambda * dUdXdL[i * 3 + 1];
            grad[2] += dGdFLambda * dUdXdL[i * 3 + 2];
            gradient[index++] = grad[0];
            gradient[index++] = grad[1];
            gradient[index++] = grad[2];
            atom.setXYZGradient(grad[0], grad[1], grad[2]);
        }

        /**
         * Update free energy F(L) every ~100 steps.
         */
        if (energyCount % 100 == 0 && propagateLambda) {
            updateFLambda(true);
        }

        /**
         * Compute the energy and gradient for the recursion slave at F(L)
         * using interpolation.
         */
        biasEnergy += computeRecursionSlave();
        
        if (print) {
            logger.info(String.format(" %s %16.8f", "Bias Energy       ", biasEnergy));
            logger.info(String.format(" %s %16.8f  %s",
                    "OSRW Potential    ", e + biasEnergy, "(Kcal/mole)"));
        }
        /**
         * Log our current state.
         */
        logger.info(String.format(" Lambda %8.6f, Bin %d, G %10.4f, dE/dLambda %10.4f",
                lambda, stateBin, biasEnergy, dEdLambda));

        /**
         * Meta-dynamics grid counts (every ~10 steps).
         */
        if (energyCount % 10 == 0 && propagateLambda) {
            recursionKernel[stateBin][FStateBin]++;
        }

        /**
         * Propagate the Lambda particle.
         */
        if (propagateLambda) {
            langevin();
        }

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
            logger.info(String.format("Current F_lambda %8.2f > maximum historgram size %8.2f.",
                    dEdLambda, maxFLambda));

            double origDeltaG = updateFLambda(false);

            int newFStateBins = FLambdaBins;
            while (minFLambda + newFStateBins * dFL < dEdLambda) {
                newFStateBins += 100;
            }
            int newRecursionKernel[][] = new int[lambdaBins][newFStateBins];
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
            FLambdaBins = newFStateBins;
            maxFLambda = minFLambda + dFL * FLambdaBins;
            logger.info(String.format("New historgram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            assert (origDeltaG == updateFLambda(false));

        }
        if (dEdLambda < minFLambda) {
            logger.info(String.format("Current F_lambda %8.2f < minimum historgram size %8.2f.",
                    dEdLambda, minFLambda));
            int offset = 100;
            while (dEdLambda < minFLambda - offset * dFL) {
                offset += 100;
            }
            int newFStateBins = FLambdaBins + offset;
            int newRecursionKernel[][] = new int[lambdaBins][newFStateBins];
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
            FLambdaBins = newFStateBins;
            logger.info(String.format("New historgram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));
        }
    }

    private double updateFLambda(boolean print) {
        double freeEnergy = 0.0;
        if (print) {
            logger.info(" Count  Lambda Bin      F_Lambda Bin   <  F_L  >     dG");
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
            // The FL range that has been sampled for iL*dL to (iL+1)*dL
            double lla = minFLambda + llFL * dFL;
            double ula = minFLambda + ulFL * dFL + dFL;
            if (ulFL == -1) {
                FLambda[iL] = 0.0;
                lla = 0.0;
                ula = 0.0;
            } else {
                double sumFLambda = 0.0;
                double partitionFunction = 0.0;
                for (int jFL = llFL; jFL <= ulFL; jFL++) {
                    double a = minFLambda + jFL * dFL + dFL_2;
                    double e = exp(evaluateKernel(iL, jFL) / (R * 300.0));
                    sumFLambda += a * e;
                    partitionFunction += e;
                    lambdaCount += recursionKernel[iL][jFL];
                }
                FLambda[iL] = sumFLambda / partitionFunction;
            }

            // The first and last bins are half size.
            double delta = dL;
            if (iL == 0 || iL == lambdaBins - 1) {
                delta *= 0.5;
            }
            freeEnergy += FLambda[iL] * delta;

            if (print) {
                double llL = iL * dL - dL_2;
                double ulL = llL + dL;
                if (llL < 0.0) {
                    llL = 0.0;
                }
                if (ulL > 1.0) {
                    ulL = 1.0;
                }

                if (lambdaBins <= 100) {
                    logger.info(String.format(" %5d [%4.2f %4.2f] [%7.1f %7.1f] <%8.3f> %8.3f",
                            lambdaCount, llL, ulL, lla, ula,
                            FLambda[iL], freeEnergy));
                } else {
                    logger.info(String.format(" %5d [%5.3f %5.3f] [%7.1f %7.1f] <%8.3f> %8.3f",
                            lambdaCount, llL, ulL, lla, ula,
                            FLambda[iL], freeEnergy));
                }
            }
        }
        return freeEnergy;
    }

    private double computeRecursionSlave() {
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
                dEdLambda += FL0 + (L1 - L0) * deltaFL / dL;
                break;
            }
        }
        return biasEnergy;
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
                    double e = weight * biasGaussianMag * exp(-deltaL2 / (2.0 * Ls2))
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
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
    }
}
