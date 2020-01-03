//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.xray;

import java.util.logging.Logger;
import static java.lang.Double.isNaN;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.Potential;
import ffx.numerics.math.ComplexNumber;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import static ffx.crystal.Crystal.invressq;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.mat3Mat3;
import static ffx.numerics.math.VectorMath.mat3SymVec6;
import static ffx.numerics.math.VectorMath.transpose3;
import static ffx.numerics.math.VectorMath.vec3Mat3;

/**
 * Fit bulk solvent and aniso B scaling terms to correct calculated structure
 * factors against data
 *
 * @author Timothy D. Fenn
 * @see <a href="http://dx.doi.org/10.1107/S0907444905007894" target="_blank">
 * P. V. Afonine, R. W. Grosse-Kunstleve and P. D. Adams, Acta Cryst. (2005).
 * D61, 850-855</a>
 * @see <a href="http://dx.doi.org/10.1107/S0021889802008580" target="_blank">
 * R. W. Grosse-Kunstleve and P. D. Adams, J. Appl. Cryst. (2002). 35,
 * 477-480.</a>
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 * @since 1.0
 */
public class ScaleBulkEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(ScaleBulkEnergy.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private static final double[][] u11 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double[][] u22 = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double[][] u33 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    private static final double[][] u12 = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double[][] u13 = {{0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    private static final double[][] u23 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
    private final double[][] recipt;
    private final double[][] j11;
    private final double[][] j22;
    private final double[][] j33;
    private final double[][] j12;
    private final double[][] j13;
    private final double[][] j23;

    private final ReflectionList reflectionList;
    private final Crystal crystal;
    private final DiffractionRefinementData refinementData;
    private final double[][] fc;
    private final double[][] fcTot;
    private final double[][] fSigF;
    private final int n;
    private final int solventN;
    private double[] optimizationScaling = null;
    private final ParallelTeam parallelTeam;
    private final ScaleBulkEnergyRegion scaleBulkEnergyRegion;
    private double totalEnergy;
    private STATE state = STATE.BOTH;

    /**
     * <p>
     * Constructor for ScaleBulkEnergy.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
     * @param n              a int.
     * @param parallelTeam   the ParallelTeam to execute the ScaleBulkEnergy.
     */
    ScaleBulkEnergy(ReflectionList reflectionList, DiffractionRefinementData refinementData,
                    int n, ParallelTeam parallelTeam) {
        this.reflectionList = reflectionList;
        this.crystal = reflectionList.crystal;
        this.refinementData = refinementData;
        this.fc = refinementData.fc;
        this.fcTot = refinementData.fcTot;
        this.fSigF = refinementData.fSigF;
        this.n = n;
        this.solventN = n - refinementData.nScale;

        recipt = transpose3(crystal.A);
        j11 = mat3Mat3(mat3Mat3(crystal.A, u11), recipt);
        j22 = mat3Mat3(mat3Mat3(crystal.A, u22), recipt);
        j33 = mat3Mat3(mat3Mat3(crystal.A, u33), recipt);
        j12 = mat3Mat3(mat3Mat3(crystal.A, u12), recipt);
        j13 = mat3Mat3(mat3Mat3(crystal.A, u13), recipt);
        j23 = mat3Mat3(mat3Mat3(crystal.A, u23), recipt);

        int threadCount = parallelTeam.getThreadCount();
        this.parallelTeam = parallelTeam;
        scaleBulkEnergyRegion = new ScaleBulkEnergyRegion(threadCount);
    }


    private class ScaleBulkEnergyRegion extends ParallelRegion {

        boolean gradient = true;
        double[] x;
        double[] g;
        double solventK;
        double modelK;
        double solventUEq;
        private final double[] modelB = new double[6];
        private final double[][] uStar = new double[3][3];
        private final double[][] resM = new double[3][3];
        SharedDouble r;
        SharedDouble rf;
        SharedDouble rFree;
        SharedDouble rFreeF;
        SharedDouble sum;
        SharedDouble sumFo;
        SharedDoubleArray grad;
        ScaleBulkEnergyLoop[] scaleBulkEnergyLoop;

        ScaleBulkEnergyRegion(int nThreads) {
            scaleBulkEnergyLoop = new ScaleBulkEnergyLoop[nThreads];
            r = new SharedDouble();
            rf = new SharedDouble();
            rFree = new SharedDouble();
            rFreeF = new SharedDouble();
            sum = new SharedDouble();
            sumFo = new SharedDouble();
        }

        public void init(double x[], double g[], boolean gradient) {
            this.x = x;
            this.g = g;
            this.gradient = gradient;
        }

        @Override
        public void start() {
            r.set(0.0);
            rf.set(0.0);
            rFree.set(0.0);
            rFreeF.set(0.0);
            sum.set(0.0);
            sumFo.set(0.0);

            for (int i = 0; i < 6; i++) {
                if (crystal.scaleB[i] >= 0) {
                    modelB[i] = x[solventN + crystal.scaleB[i]];
                }
            }

            modelK = x[0];
            solventK = refinementData.bulkSolventK;
            solventUEq = refinementData.bulkSolventUeq;
            if (solventN > 1) {
                solventK = x[1];
                solventUEq = x[2];
            }

            // Generate Ustar.
            mat3SymVec6(crystal.A, modelB, resM);
            mat3Mat3(resM, recipt, uStar);

            if (gradient) {
                if (grad == null) {
                    grad = new SharedDoubleArray(g.length);
                }
                for (int i = 0; i < g.length; i++) {
                    grad.set(i, 0.0);
                }
            }
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();
            if (scaleBulkEnergyLoop[ti] == null) {
                scaleBulkEnergyLoop[ti] = new ScaleBulkEnergyLoop();
            }

            try {
                execute(0, reflectionList.hkllist.size() - 1, scaleBulkEnergyLoop[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        @Override
        public void finish() {
            if (gradient) {
                for (int i = 0; i < g.length; i++) {
                    g[i] = grad.get(i);
                }
            }
        }

        private class ScaleBulkEnergyLoop extends IntegerForLoop {

            private final double[] resv = new double[3];
            private final double[] ihc = new double[3];
            private final ComplexNumber resc = new ComplexNumber();
            private final ComplexNumber fcc = new ComplexNumber();
            private final ComplexNumber fsc = new ComplexNumber();
            private final ComplexNumber fct = new ComplexNumber();
            private final ComplexNumber kfct = new ComplexNumber();
            private final double[] lgrad;
            private double lr;
            private double lrf;
            private double lrfree;
            private double lrfreef;
            private double lsum;
            private double lsumfo;

            ScaleBulkEnergyLoop() {
                lgrad = new double[g.length];
            }

            @Override
            public void start() {
                lr = 0.0;
                lrf = 0.0;
                lrfree = 0.0;
                lrfreef = 0.0;
                lsum = 0.0;
                lsumfo = 0.0;
                fill(lgrad, 0.0);
            }

            @Override
            public void run(int lb, int ub) throws Exception {

                for (int j = lb; j <= ub; j++) {
                    HKL ih = reflectionList.hkllist.get(j);
                    int i = ih.index();
                    if (isNaN(fc[i][0]) || isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                        continue;
                    }

                    // Constants
                    double s = invressq(crystal, ih);
                    ihc[0] = ih.h();
                    ihc[1] = ih.k();
                    ihc[2] = ih.l();
                    vec3Mat3(ihc, uStar, resv);
                    double u = modelK - dot(resv, ihc);
                    double expBS = exp(-twopi2 * solventUEq * s);
                    double ksExpBS = solventK * expBS;
                    double expU = exp(0.25 * u);

                    // Structure Factors
                    refinementData.getFcIP(i, fcc);
                    refinementData.getFsIP(i, fsc);
                    fct.copy(fcc);
                    if (refinementData.crystalReciprocalSpaceFs.solventModel != SolventModel.NONE) {
                        resc.copy(fsc);
                        resc.timesIP(ksExpBS);
                        fct.plusIP(resc);
                    }
                    kfct.copy(fct);
                    kfct.timesIP(expU);

                    // Total structure factor (for refinement)
                    fcTot[i][0] = kfct.re();
                    fcTot[i][1] = kfct.im();

                    // Target
                    double f1 = refinementData.getF(i);
                    double akfct = kfct.abs();
                    double af1 = abs(f1);
                    double d = f1 - akfct;
                    double d2 = d * d;
                    double dr = -2.0 * d;

                    lsum += d2;
                    lsumfo += f1 * f1;

                    if (refinementData.isFreeR(i)) {
                        lrfree += abs(af1 - abs(akfct));
                        lrfreef += af1;
                    } else {
                        lr += abs(af1 - abs(akfct));
                        lrf += af1;
                    }

                    if (gradient) {
                        // modelK/modelB - common derivative element
                        double dfm = 0.25 * akfct * dr;
                        // bulk solvent - common derivative element
                        double afsc = fsc.abs();
                        double dfb = expBS * (fcc.re() * fsc.re() + fcc.im() * fsc.im() + ksExpBS * afsc * afsc);

                        // modelK derivative
                        lgrad[0] += dfm;
                        if (solventN > 1) {
                            double iafct = 1.0 / fct.abs();
                            // solventK derivative
                            lgrad[1] += expU * dfb * dr * iafct;
                            // solventUeq derivative
                            lgrad[2] += expU * -twopi2 * s * solventK * dfb * dr * iafct;
                        }

                        for (int jj = 0; jj < 6; jj++) {
                            if (crystal.scaleB[jj] >= 0) {
                                switch (jj) {
                                    case (0):
                                        // B11
                                        vec3Mat3(ihc, j11, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                    case (1):
                                        // B22
                                        vec3Mat3(ihc, j22, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                    case (2):
                                        // B33
                                        vec3Mat3(ihc, j33, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                    case (3):
                                        // B12
                                        vec3Mat3(ihc, j12, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                    case (4):
                                        // B13
                                        vec3Mat3(ihc, j13, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                    case (5):
                                        // B23
                                        vec3Mat3(ihc, j23, resv);
                                        lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                                        break;
                                }
                            }
                        }
                    }
                }
            }

            @Override
            public void finish() {
                r.addAndGet(lr);
                rf.addAndGet(lrf);
                rFree.addAndGet(lrfree);
                rFreeF.addAndGet(lrfreef);
                sum.addAndGet(lsum);
                sumFo.addAndGet(lsumfo);
                for (int i = 0; i < lgrad.length; i++) {
                    grad.getAndAdd(i, lgrad[i]);
                }
            }
        }
    }

    /**
     * <p>
     * target</p>
     *
     * @param x        an array of double.
     * @param g        an array of double.
     * @param gradient a boolean.
     * @param print    a boolean.
     * @return a double.
     */
    public double target(double[] x, double[] g, boolean gradient, boolean print) {

        try {
            scaleBulkEnergyRegion.init(x, g, gradient);
            parallelTeam.execute(scaleBulkEnergyRegion);
        } catch (Exception e) {
            logger.info(e.toString());
        }

        double sum = scaleBulkEnergyRegion.sum.get();
        double sumfo = scaleBulkEnergyRegion.sumFo.get();
        double r = scaleBulkEnergyRegion.r.get();
        double rf = scaleBulkEnergyRegion.rf.get();
        double rfree = scaleBulkEnergyRegion.rFree.get();
        double rfreef = scaleBulkEnergyRegion.rFreeF.get();

        if (gradient) {
            double isumfo = 1.0 / sumfo;
            for (int i = 0; i < g.length; i++) {
                g[i] *= isumfo;
            }
        }

        if (print) {
            StringBuilder sb = new StringBuilder("\n");
            sb.append("Bulk solvent and scale fit\n");
            sb.append(String.format("   residual:  %8.3f\n", sum / sumfo));
            sb.append(String.format("   R:  %8.3f  Rfree:  %8.3f\n",
                    (r / rf) * 100.0, (rfree / rfreef) * 100.0));
            sb.append("x: ");
            for (double x1 : x) {
                sb.append(String.format("%8g ", x1));
            }
            sb.append("\ng: ");
            for (double v : g) {
                sb.append(String.format("%8g ", v));
            }
            sb.append("\n");
            logger.info(sb.toString());
        }
        totalEnergy = sum / sumfo;
        return sum / sumfo;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {
        unscaleCoordinates(x);
        double sum = target(x, null, false, false);
        scaleCoordinates(x);
        return sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        unscaleCoordinates(x);
        double sum = target(x, g, true, false);
        scaleCoordinatesAndGradient(x, g);
        return sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling) {
        if (scaling != null && scaling.length == n) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     * <p>
     * Return a reference to each variables type.
     */
    @Override
    public Potential.VARIABLE_TYPE[] getVariableTypes() {
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setEnergyTermState(Potential.STATE state) {
        this.state = state;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        // The parallelTeam should have been passed in by DiffractionData, which handles destroying it.
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
