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

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
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
 */
public class ScaleBulkOptimizer implements Optimizable {

    private static final Logger logger = Logger.getLogger(ScaleBulkOptimizer.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private static final double eightpi2 = 8.0 * PI * PI;
    private final ReflectionList reflectionlist;
    private final ReflectionSpline spline;
    private final Crystal crystal;
    private final RefinementData refinementdata;
    private final double fc[][];
    private final double fs[][];
    private final double fctot[][];
    private final double fo[][];
    private final int freer[];
    protected double[] optimizationScaling = null;

    public ScaleBulkOptimizer(ReflectionList reflectionlist, RefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.refinementdata = refinementdata;
        this.fc = refinementdata.fc;
        this.fs = refinementdata.fs;
        this.fctot = refinementdata.fctot;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;

        // initialize params
        this.spline = new ReflectionSpline(reflectionlist,
                refinementdata.nparams);
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

        int scale_flag = crystal.scale_flag;
        double model_k = x[0];
        double model_b[] = {x[1], x[2], x[3], x[4], x[5], x[6]};
        double solvent_k = x[7];
        double solvent_ueq = x[8];
        r = rf = rfree = rfreef = sum = sumfo = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (ih.allowed() == 0.0
                    || Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])) {
                continue;
            }

            //  constants
            double ihc[] = {ih.h(), ih.k(), ih.l()};
            double ihf[] = new double[3];
            crystal.toFractionalCoordinates(ihc, ihf);
            double u = exp(-0.25
                    * (pow(ihf[0], 2.0) * model_b[0]
                    + pow(ihf[1], 2.0) * model_b[1]
                    + pow(ihf[2], 2.0) * model_b[2]
                    + 2.0 * ihf[0] * ihf[1] * model_b[3]
                    + 2.0 * ihf[0] * ihf[2] * model_b[4]
                    + 2.0 * ihf[1] * ihf[2] * model_b[5]));
            double s = Crystal.invressq(crystal, ih);
            double ebs = exp(-twopi2 * solvent_ueq * s);
            double ksebs = solvent_k * ebs;
            double kmems = model_k * u;

            // spline setup
            double fh = spline.f(s, refinementdata.spline);

            // structure factors
            ComplexNumber fcc = new ComplexNumber(fc[i][0], fc[i][1]);
            ComplexNumber fsc = new ComplexNumber(fs[i][0], fs[i][1]);
            ComplexNumber fct = fcc.plus(fsc.times(ksebs));
            ComplexNumber kfct = fct.times(kmems);

            // total structure factor (for refinement)
            fctot[i][0] = kfct.re();
            fctot[i][1] = kfct.im();

            // target
            double f1 = fo[i][0];
            double f2 = kfct.abs();
            double d = f1 - fh * f2;
            double d2 = d * d;
            double dr = -2.0 * fh * d;

            sum += d2;
            sumfo += f1 * f1;

            if (freer[i] == 0) {
                rfree += abs(abs(fo[i][0]) - abs(fh * kfct.abs()));
                rfreef += abs(fo[i][0]);
            } else {
                r += abs(abs(fo[i][0]) - abs(fh * kfct.abs()));
                rf += abs(fo[i][0]);
            }

            if (gradient) {
                // common derivative element
                double dfp = fcc.re() * fsc.re()
                        + fcc.im() * fsc.im()
                        + ksebs * pow(fsc.abs(), 2.0);

                g[0] += fct.abs() * u * dr;

                if ((scale_flag & Crystal.SCALE_B11) == Crystal.SCALE_B11) {
                    g[1] += kfct.abs() * -0.25 * pow(ihf[0], 2.0) * dr;
                }
                if ((scale_flag & Crystal.SCALE_B22) == Crystal.SCALE_B22) {
                    g[2] += kfct.abs() * -0.25 * pow(ihf[1], 2.0) * dr;
                }
                if ((scale_flag & Crystal.SCALE_B33) == Crystal.SCALE_B33) {
                    g[3] += kfct.abs() * -0.25 * pow(ihf[2], 2.0) * dr;
                }
                if ((scale_flag & Crystal.SCALE_B12) == Crystal.SCALE_B12) {
                    g[4] += kfct.abs() * -0.5 * ihf[0] * ihf[1] * dr;
                }
                if ((scale_flag & Crystal.SCALE_B13) == Crystal.SCALE_B13) {
                    g[5] += kfct.abs() * -0.5 * ihf[0] * ihf[2] * dr;
                }
                if ((scale_flag & Crystal.SCALE_B23) == Crystal.SCALE_B23) {
                    g[6] += kfct.abs() * -0.5 * ihf[1] * ihf[2] * dr;
                }

                g[7] += model_k * u * (ebs * dfp) * dr / fct.abs();
                g[8] += model_k * u * (-twopi2 * s * ksebs * dfp) * dr / fct.abs();
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
        if (scaling != null && scaling.length == 9) {
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
