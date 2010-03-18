/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.xray;

import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.PI;
import static java.lang.Math.pow;

import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.mat3mat3;
import static ffx.numerics.VectorMath.mat3symvec6;
import static ffx.numerics.VectorMath.transpose3;
import static ffx.numerics.VectorMath.vec3mat3;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.ComplexNumber;
import ffx.numerics.Optimizable;

/**
 *
 * Fit bulk solvent and aniso B scaling terms to correct calculated structure
 * factors against data
 *
 * @author Tim Fenn<br>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444905007894" target="_blank">
 * P. V. Afonine, R. W. Grosse-Kunstleve and P. D. Adams,
 * Acta Cryst. (2005). D61, 850-855
 *
 * @see <a href="http://dx.doi.org/10.1107/S0021889802008580" target="_blank">
 * R. W. Grosse-Kunstleve and P. D. Adams, J. Appl. Cryst. (2002). 35, 477-480.
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank">
 * J. A. Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank">
 * J. S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.
 */
public class ScaleBulkEnergy implements Optimizable {

    private static final Logger logger = Logger.getLogger(ScaleBulkEnergy.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private static final double eightpi2 = 8.0 * PI * PI;
    private static final double u11[][] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u22[][] = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u33[][] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    private static final double u12[][] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u13[][] = {{0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    private static final double u23[][] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
    private final ReflectionList reflectionlist;
    private final Crystal crystal;
    private final RefinementData refinementdata;
    private final double fc[][];
    private final double fs[][];
    private final double fctot[][];
    private final double fo[][];
    private final int freer[];
    private final int n;
    protected double[] optimizationScaling = null;
    private final double recipt[][];
    private final double j11[][];
    private final double j22[][];
    private final double j33[][];
    private final double j12[][];
    private final double j13[][];
    private final double j23[][];

    public ScaleBulkEnergy(ReflectionList reflectionlist, RefinementData refinementdata, int n) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.refinementdata = refinementdata;
        this.fc = refinementdata.fc;
        this.fs = refinementdata.fs;
        this.fctot = refinementdata.fctot;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;
        this.n = n;

        recipt = transpose3(crystal.recip);
        j11 = mat3mat3(mat3mat3(crystal.recip, u11), recipt);
        j22 = mat3mat3(mat3mat3(crystal.recip, u22), recipt);
        j33 = mat3mat3(mat3mat3(crystal.recip, u33), recipt);
        j12 = mat3mat3(mat3mat3(crystal.recip, u12), recipt);
        j13 = mat3mat3(mat3mat3(crystal.recip, u13), recipt);
        j23 = mat3mat3(mat3mat3(crystal.recip, u23), recipt);
    }

    public double target(double x[], double g[],
            boolean gradient, boolean print) {
        double r, rf, rfree, rfreef, sum, sumfo;

        // zero out the gradient
        if (gradient) {
            for (int i = 0; i < g.length; i++) {
                g[i] = 0.0;
            }
        }

        int scale_n = crystal.scale_n;
        double model_k = x[0];
        double solvent_k = 0.0;
        double solvent_ueq = 0.0;
        if (refinementdata.solvent_n > 1) {
            solvent_k = x[1];
            solvent_ueq = x[2];
        }
        double model_b[] = new double[6];
        for (int i = 0; i < 6; i++) {
            if (crystal.scale_b[i] >= 0) {
                model_b[i] = x[refinementdata.solvent_n + crystal.scale_b[i]];
            }
        }
        double ustar[][] = mat3mat3(mat3symvec6(crystal.recip, model_b), recipt);

        r = rf = rfree = rfreef = sum = sumfo = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            //  constants
            double s = Crystal.invressq(crystal, ih);
            double ihc[] = {ih.h(), ih.k(), ih.l()};
            double u = model_k - dot(vec3mat3(ihc, ustar), ihc);
            double ebs = exp(-twopi2 * solvent_ueq * s);
            double ksebs = solvent_k * ebs;
            double kmebm = exp(0.25 * u);

            // structure factors
            ComplexNumber fcc = refinementdata.fc(i);
            ComplexNumber fsc = refinementdata.fs(i);
            ComplexNumber fct = (refinementdata.solvent_n > 1)
                    ? fcc.plus(fsc.times(ksebs)) : fcc;
            ComplexNumber kfct = fct.times(kmebm);

            // total structure factor (for refinement)
            fctot[i][0] = kfct.re();
            fctot[i][1] = kfct.im();

            // target
            double f1 = refinementdata.f(i);
            double f2 = kfct.abs();
            double d = f1 - f2;
            double d2 = d * d;
            double dr = -2.0 * d;

            sum += d2;
            sumfo += f1 * f1;

            if (refinementdata.isfreer(i)) {
                rfree += abs(abs(f1) - abs(kfct.abs()));
                rfreef += abs(f1);
            } else {
                r += abs(abs(f1) - abs(kfct.abs()));
                rf += abs(f1);
            }

            if (gradient) {
                // model_k/model_b - common derivative element
                double dfm = 0.25 * kfct.abs() * dr;
                // bulk solvent - common derivative element
                double dfb = ebs * (fcc.re() * fsc.re()
                        + fcc.im() * fsc.im()
                        + ksebs * pow(fsc.abs(), 2.0));

                // model_k derivative
                g[0] += dfm;
                if (refinementdata.solvent_n > 1) {
                    // solvent_k derivative
                    g[1] += kmebm * dfb * dr / fct.abs();
                    // solvent_ueq derivative
                    g[2] += kmebm * -twopi2 * s * solvent_k * dfb * dr / fct.abs();
                }

                int sn = refinementdata.solvent_n;
                for (int j = 0; j < 6; j++) {
                    if (crystal.scale_b[j] >= 0) {
                        switch (j) {
                            case (0):
                                // B11
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j11), ihc);
                                break;
                            case (1):
                                // B22
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j22), ihc);
                                break;
                            case (2):
                                // B33
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j33), ihc);
                                break;
                            case (3):
                                // B12
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j12), ihc);
                                break;
                            case (4):
                                // B13
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j13), ihc);
                                break;
                            case (5):
                                // B23
                                // g[sn + crystal.scale_b[j]] += 0.25 * kfct.abs() * -2.0 * ihf[1] * ihf[2] * dr;
                                g[sn + crystal.scale_b[j]] += -dfm
                                        * dot(vec3mat3(ihc, j23), ihc);
                                break;
                        }
                    }
                }
            }
        }

        if (gradient) {
            for (int i = 0; i < g.length; i++) {
                g[i] /= sumfo;
            }
        }

        if (print) {
            StringBuffer sb = new StringBuffer("\n");
            sb.append(" Bulk solvent and scale fit\n");
            sb.append(String.format("   residual:  %8.3f\n", sum / sumfo));
            sb.append(String.format("   R:  %8.3f  Rfree:  %8.3f\n",
                    (r / rf) * 100.0, (rfree / rfreef) * 100.0));
            sb.append("x: ");
            for (int i = 0; i < x.length; i++) {
                sb.append(String.format("%8g ", x[i]));
            }
            sb.append("\ng: ");
            for (int i = 0; i < g.length; i++) {
                sb.append(String.format("%8g ", g[i]));
            }
            sb.append("\n");
            logger.info(sb.toString());
        }

        return sum / sumfo;
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, g, true, true);

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }

        return sum;
    }

    @Override
    public void setOptimizationScaling(double[] scaling) {
        if (scaling != null && scaling.length == n) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getOptimizationScaling() {
        return optimizationScaling;
    }
}
