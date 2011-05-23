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
package ffx.potential;

import ffx.numerics.Potential;
import java.util.logging.Logger;

/**
 * Implements the Orthogonal Space Random Walk method.
 *
 * @author Michael J. Schnieders, Wei Yang & Pengyu Ren
 * @since 1.0
 */
public class OSRW implements Potential {

    private static final Logger logger = Logger.getLogger(OSRW.class.getName());
    protected double lambda;
    protected boolean lambdaGradient = false;
    protected int energyCount;
    protected int lambdaBins = 100;
    protected double λBinWidth = 1.0 / lambdaBins;
    protected double λBinHalfWidth = λBinWidth / 2.0;
    protected double mindEdλ = -500.0, maxdEdλ = 500.0;
    protected double dEdλSpan = maxdEdλ - mindEdλ;
    protected double dEdλWidth = 2.0;
    protected int dEdλBins = (int) Math.floor(dEdλSpan / dEdλWidth);
    protected int recursionKernel[][];
    protected double dEdλ = 0.0;
    protected double d2Edλ2 = 0.0;
    protected double dUdXdL[] = null;
    protected double gaussianMag = 0.002;
    protected double dAdλ[];
    protected double recursionKernelEnergy;
    private final ForceFieldEnergy forceFieldEnergy;

    public OSRW(ForceFieldEnergy forceFieldEnergy) {

        this.forceFieldEnergy = forceFieldEnergy;

        // Zero out the recusion kernel energy.
        recursionKernelEnergy = 0.0;
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {

        energyCount++;

        if (lambdaGradient) {
            dEdλ = 0.0;
            d2Edλ2 = 0.0;
            /* 
            if (vanderWaalsTerm) {
            dEdλ = vanderWaals.getdEdLambda();
            d2Edλ2 = vanderWaals.getd2EdLambda2();
            }
            if (multipoleTerm) {
            dEdλ += particleMeshEwald.getdEdLambda();
            d2Edλ2 += particleMeshEwald.getd2EdLambda2();
            } */
            logger.info(String.format(" Lambda %6.4f dE/dLambda %10.4f", lambda, dEdλ));

            int lambdaBin = (int) Math.floor(lambda * lambdaBins);
            if (lambda == 1.0) {
                lambdaBin = lambdaBins - 1;
            }

            int dUdLBin = (int) Math.floor((dEdλ - mindEdλ) / dEdλWidth);
            if (dEdλ == maxdEdλ) {
                dUdLBin = dUdLBin - 1;
            }

            /**
             * If necessary, allocate more space for dEdλ.
             */
            if (dUdLBin >= dEdλBins) {
                double newMaxdEdλ = maxdEdλ;
                while (newMaxdEdλ < dEdλ) {
                    newMaxdEdλ += 100.0;
                }
                int newdEdλBins = (int) Math.floor((newMaxdEdλ - mindEdλ) / dEdλWidth);
                int newRecursionKernel[][] = new int[lambdaBins][newdEdλBins];
                for (int i = 0; i < lambdaBins; i++) {
                    for (int j = 0; j < dEdλBins; j++) {
                        newRecursionKernel[i][j] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                maxdEdλ = newMaxdEdλ;
                dEdλBins = newdEdλBins;
            }
            if (dUdLBin < 0) {
                double newMindEdλ = mindEdλ;
                while (newMindEdλ > dEdλ) {
                    newMindEdλ -= 100.0;
                }
                int newdEdλBins = (int) Math.floor((maxdEdλ - newMindEdλ) / dEdλWidth);
                int newRecursionKernel[][] = new int[lambdaBins][newdEdλBins];
                for (int i = 0; i < lambdaBins; i++) {
                    for (int j = 0; j < dEdλBins; j++) {
                        newRecursionKernel[i][j] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                mindEdλ = newMindEdλ;
                dEdλBins = newdEdλBins;
            }

            /**
             * Calculate recursion kernel G(λ, dEdλ) and gradient.
             */
            double dgdL = 0.0;
            double dgdEdL = 0.0;
            double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
            double dEdLs2 = dEdλWidth * 2.0 * dEdλWidth * 2.0;
            for (int l = -5; l <= 5; l++) {
                int lcenter = lambdaBin + 1;
                double lv = lcenter / lambdaBins + 0.5 / lambdaBins;
                double lv2 = (lambda - lv) * (lambda - lv);
                // Mirror conditions for recursion kernel counts.
                int lcount = lcenter;
                if (lcount < 0) {
                    lcount = -lcount;
                }
                if (lcount >= lambdaBins) {
                    lcount = lambdaBins - (lcount - lambdaBins) - 1;
                }
                for (int dl = -5; dl <= 5; dl++) {
                    int dlcenter = dUdLBin + dl;
                    double dlv = dlcenter * dEdλWidth + dEdλWidth / 2.0;
                    double dlv2 = (dEdλ - dlv) * (dEdλ - dlv);
                    int weight = recursionKernel[lcount][dlcenter];
                    double e = weight * gaussianMag * Math.exp(-lv2 / (2.0 * ls2)) * Math.exp(-dlv2 / (2.0 * dEdLs2));
                    recursionKernelEnergy += e;
                    dgdL += -lv / ls2 * e;
                    dgdEdL += -dlv / dEdLs2 * e;
                }
            }

            /**
             * λ gradient due to recursion kernel G(λ, dEdλ).
             */
            dEdλ += dgdL + dgdEdL * d2Edλ2;

            /**
             * Atomic gradient due to recursion kernel G(λ, dEdλ).
             */
            /*
            for (int i=0; i<3*nAtoms; i++) {
            dUdXdL[i] = 0.0;
            }
            if (vanderWaalsTerm) {
            vanderWaals.getdEdLambdaGradient(dUdXdL);
            }
            if (multipoleTerm) {
            particleMeshEwald.getdEdLambdaGradient(dUdXdL);
            } 
             * 
            double grad[] = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.getXYZGradient(grad);
                grad[0] += dgdEdL * dUdXdL[i * 3];
                grad[1] += dgdEdL * dUdXdL[i * 3 + 1];
                grad[2] += dgdEdL * dUdXdL[i * 3 + 2];
                atom.setXYZGradient(grad[0], grad[1], grad[2]);
            } */

            // Update free energy F(λ) every ~100 steps.
            if (energyCount % 100 == 0) {
                double freeEnergy = 0.0;
                double binHalf = 0.5 / lambdaBins;
                logger.info("λ              dAdλ           Cumulative");
                for (int i = 0; i < lambdaBins; i++) {
                    double lc = binHalf + i * λBinWidth;

                    int ul = -1;
                    int ll = -1;
                    // Find the smallest dUdL bin.
                    for (int j = 0; j < dEdλBins; j++) {
                        int count = recursionKernel[i][j];
                        if (count != 0 && ll == -1) {
                            ll = j;
                            break;
                        }
                    }
                    // Find the largest dUdL bin.
                    for (int j = dEdλBins - 1; j >= 0; j--) {
                        int count = recursionKernel[i][j];
                        if (count != 0 && ul == -1) {
                            ul = j;
                            break;
                        }
                    }

                    if (ul == -1) {
                        dAdλ[i] = 0.0;
                    } else {
                        double wdUdL = 0.0;
                        double part = 0.0;
                        for (int j = ll; j <= ul; j++) {
                            double dUdLc = mindEdλ + (j + 0.5) * dEdλWidth;
                            double e = Math.exp(gKernel(lc, dUdLc) / (R * 300.0));
                            wdUdL += dUdLc * e;
                            part += e;
                        }
                        dAdλ[i] = wdUdL / part;
                    }
                    freeEnergy += dAdλ[i] * λBinWidth;
                    logger.info(String.format("%15.8f %15.8f %15.8f",
                                              i * λBinWidth + λBinHalfWidth, dAdλ[i], freeEnergy));
                }
            }

            /**
             * Force interpolation due to recursion slave F(λ).
             */
            if (lambda > λBinHalfWidth && lambda < 1.0 - λBinHalfWidth) {
                double binCenter = lambdaBin * λBinWidth + λBinHalfWidth;
                int lb;
                int ub;
                if (lambda > binCenter) {
                    lb = lambdaBin;
                } else {
                    lb = lambdaBin - 1;
                }
                ub = lb + 1;
                double dAdLm = dAdλ[lb];
                double dAdLp = dAdλ[ub];
                double m1c = lb * λBinWidth + λBinHalfWidth;
                double p1c = ub * λBinWidth + λBinHalfWidth;
                dEdλ -= ((lambda - m1c) * dAdLp + (p1c - lambda) * dAdLm) / λBinWidth;
            } else if (lambda <= λBinHalfWidth) {
                double mlc = λBinHalfWidth;
                double plc = λBinWidth + λBinHalfWidth;
                dEdλ -= ((lambda - mlc) * dAdλ[1] + (plc - lambda) * dAdλ[0]) / λBinWidth;
            } else {
                double mlc = 1.0 - 1.5 * λBinWidth;
                double plc = 1.0 - λBinHalfWidth;
                dEdλ -= ((lambda - mlc) * dAdλ[lambdaBins - 1] + (plc - lambda) * dAdλ[lambdaBins - 2]) / λBinWidth;
            }

            /**
             * Meta-dynamic grid counts (every ~10 steps).
             */
            if (energyCount % 10 == 0) {
                recursionKernel[lambdaBin][dUdLBin]++;
            }
        }

        return recursionKernelEnergy;
    }
    public static final double R = 1.9872066e-3;

    public double gKernel(double lambda, double dUdL) {

        int lambdaBin = (int) Math.ceil(lambda * lambdaBins) - 1;
        int dUdLBin = (int) Math.ceil((dUdL - mindEdλ) / dEdλSpan * dEdλBins) - 1;

        double sum = 0.0;
        double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
        double dUdLs2 = dEdλWidth * 2.0 * dEdλWidth * 2.0;

        for (int l = -5; l <= 5; l++) {
            int lcenter = lambdaBin + 1;
            double lv = lcenter / lambdaBins + 0.5 / lambdaBins;
            double lv2 = (lambda - lv) * (lambda - lv);

            // Mirror condition for Lambda counts.
            int lcount = lcenter;
            if (lcount < 0) {
                lcount = -lcount;
            }
            if (lcount >= lambdaBins) {
                lcount = lambdaBins - (lcount - lambdaBins) - 1;
            }

            for (int dl = -5; dl <= 5; dl++) {
                int dlcenter = dUdLBin + dl;
                double dlv = dlcenter * dEdλWidth + dEdλWidth / 2.0;
                double dlv2 = (dUdL - dlv) * (dUdL - dlv);
                int weight = recursionKernel[lcount][dlcenter];
                double e = weight * gaussianMag * Math.exp(-lv2 / (2.0 * ls2)) * Math.exp(-dlv2 / (2.0 * dUdLs2));
                sum += e;
            }
        }
        return sum;
    }

    @Override
    public void setScaling(double[] scaling) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getScaling() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getMass() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getNumberOfVariables() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
