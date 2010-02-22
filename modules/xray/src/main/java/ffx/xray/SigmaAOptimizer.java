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
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.numerics.Optimizable;

/**
 *
 * Fit structure factors using spline coefficients
 *
 * @author Tim Fenn<br>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0021889804031474" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2005). 38, 193-198
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444992007352" target="_blank">
 * A. T. Brunger, Acta Cryst. (1993). D49, 24-36
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444996012255" target="_blank">
 * G. N. Murshudov, A. A. Vagin and E. J. Dodson,
 * Acta Cryst. (1997). D53, 240-255
 * 
 */
public class SigmaAOptimizer implements Optimizable {

    private static final Logger logger = Logger.getLogger(SigmaAOptimizer.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private static final double sim_a = 1.639294;
    private static final double sim_b = 3.553967;
    private static final double sim_c = 2.228716;
    private static final double sim_d = 3.524142;
    private static final double sim_e = 7.107935;
    private static final double sim_A = -1.28173889903;
    private static final double sim_B = 0.69231689903;
    private static final double sim_C = -1.33099462667;
    private static final double sim_g = 2.13643992379;
    private static final double sim_p = 0.04613803811;
    private static final double sim_q = 1.82167089029;
    private static final double sim_r = -0.74817947490;
    private final ReflectionList reflectionlist;
    private final ReflectionSpline spline;
    private final int n;
    private final Crystal crystal;
    private final RefinementData refinementdata;
    private final double fo[][];
    private final int freer[];
    private final double fc[][];
    private final double fs[][];
    private final double fctot[][];
    private final double fofc1[][];
    private final double fofc2[][];
    private final double fomphi[][];
    private final double dfc[][];
    private final double dfs[][];
    protected double[] optimizationScaling = null;

    public SigmaAOptimizer(ReflectionList reflectionlist,
            RefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.refinementdata = refinementdata;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;
        this.fc = refinementdata.fc;
        this.fs = refinementdata.fs;
        this.fctot = refinementdata.fctot;
        this.fomphi = refinementdata.fomphi;
        this.fofc1 = refinementdata.fofc1;
        this.fofc2 = refinementdata.fofc2;
        this.dfc = refinementdata.dfc;
        this.dfs = refinementdata.dfs;
        this.n = refinementdata.nparams;

        // initialize params
        this.spline = new ReflectionSpline(reflectionlist, n);
    }

    public double target(double x[], double g[],
            boolean gradient, boolean print) {
        int nsum, nsumr;
        double sum, sumr;

        // zero out the gradient
        if (gradient) {
            for (int i = 0; i < g.length; i++) {
                g[i] = 0.0;
            }
        }

        double model_k = refinementdata.model_k;
        double solvent_k = refinementdata.solvent_k;
        double solvent_ueq = refinementdata.solvent_ueq;
        double model_b[] = new double[6];
        for (int i = 0; i < 6; i++) {
            model_b[i] = refinementdata.aniso_b[i];
        }
        double sa[] = new double[n];
        double wa[] = new double[n];
        for (int i = 0; i < n; i++) {
            sa[i] = 1.0 + x[i];
            wa[i] = x[n + i];
        }
        nsum = nsumr = 0;
        sum = sumr = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (ih.allowed() == 0.0) {
                continue;
            }

            // constants
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
            double u2 = u * u;
            double k2 = model_k * model_k;
            double s = Crystal.invressq(crystal, ih);
            double ebs = exp(-twopi2 * solvent_ueq * s);
            double ksebs = solvent_k * ebs;
            double kmems = model_k * u;
            double km2 = kmems * kmems;
            double epsc = ih.epsilonc();

            // spline setup
            double ecscale = spline.f(s, refinementdata.fcesq);
            double eoscale = spline.f(s, refinementdata.foesq);
            double sai = spline.f(s, sa);
            double wai = spline.f(s, wa);
            double sa2 = pow(sai, 2.0);

            // structure factors
            ComplexNumber fcc = new ComplexNumber(fc[i][0], fc[i][1]);
            ComplexNumber fsc = new ComplexNumber(fs[i][0], fs[i][1]);
            ComplexNumber fct = fcc.plus(fsc.times(ksebs));
            ComplexNumber kfct = fct.times(kmems);

            ComplexNumber ecc = fcc.times(sqrt(ecscale));
            ComplexNumber esc = fsc.times(sqrt(ecscale));
            ComplexNumber ect = fct.times(sqrt(ecscale));
            ComplexNumber kect = kfct.times(sqrt(ecscale));
            double eo = fo[i][0] * sqrt(eoscale);
            double sigeo = fo[i][1] * sqrt(eoscale);
            double eo2 = pow(eo, 2.0);
            double ect2 = pow(ect.abs(), 2.0);
            double kect2 = pow(kect.abs(), 2.0);

            double d = 2.0 * sigeo * sigeo + epsc * wai;
            double d2 = d * d;
            double fomx = 2.0 * eo * sai * kect.abs() / d;
            double dfbulk = 2.0 * ebs * (fcc.re() * fsc.re()
                    + fcc.im() * fsc.im()
                    + ksebs * pow(fsc.abs(), 2.0));

            double inot, dinot, cf;
            if (ih.centric()) {
                inot = (abs(fomx) < 10.0) ? log(cosh(fomx)) : abs(fomx) + log(0.5);
                dinot = tanh(fomx);
                cf = 0.5;
            } else {
                inot = sim_integ(fomx);
                dinot = sim(fomx);
                cf = 1.0;
            }
            double llk = cf * log(d) + (eo2 + sa2 * kect2) / d - inot;

            // derivatives (scary)
            double dfcr = (2.0 * sa2 * km2 * ect.re()) / d - ((2.0 * eo * sai * kmems * ect.re()) / (d * ect.abs())) * dinot;
            double dfci = (2.0 * sa2 * km2 * ect.im()) / d - ((2.0 * eo * sai * kmems * ect.im()) / (d * ect.abs())) * dinot;
            double dfsr = (2.0 * sa2 * km2 * ksebs * ect.re()) / d - ((2.0 * eo * sai * kmems * ksebs * ect.re()) / (d * ect.abs())) * dinot;
            double dfsi = (2.0 * sa2 * km2 * ksebs * ect.im()) / d - ((2.0 * eo * sai * kmems * ksebs * ect.im()) / (d * ect.abs())) * dinot;
            double dfsa = 2.0 * sai * kect2 / d - (2.0 * eo * kect.abs() / d) * dinot;
            double dfwa = epsc * (cf / d - (eo2 + sa2 * kect2) / d2 + (2.0 * eo * sai * kect.abs() / d2) * dinot);
            /*
            double dfkm = (2.0 * u2 * model_k * sa2 * ect2) / d - ((2.0 * u * eo * sai * ect.abs()) / d) * dinot;
            double dfks = (km2 * sa2 * ecscale * dfbulk) / d - ((kmems * eo * sai * ecscale * dfbulk) / (d * ect.abs())) * dinot;
            double dfbs = (km2 * sa2 * -twopi2 * s * solvent_k * ecscale * dfbulk) / d - ((kmems * eo * sai * -twopi2 * s * solvent_k * ecscale * dfbulk) / (d * ect.abs())) * dinot;
            double dfb11 = (2.0 * kmems * model_k * sa2 * ect2 * -0.25 * pow(ihf[0], 2.0)) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.25 * pow(ihf[0], 2.0)) / d) * dinot;
            double dfb22 = (2.0 * kmems * model_k * sa2 * ect2 * -0.25 * pow(ihf[1], 2.0)) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.25 * pow(ihf[1], 2.0)) / d) * dinot;
            double dfb33 = (2.0 * kmems * model_k * sa2 * ect2 * -0.25 * pow(ihf[2], 2.0)) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.25 * pow(ihf[2], 2.0)) / d) * dinot;
            double dfb12 = (2.0 * kmems * model_k * sa2 * ect2 * -0.5 * ihf[0] * ihf[1]) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.5 * ihf[0] * ihf[1]) / d) * dinot;
            double dfb13 = (2.0 * kmems * model_k * sa2 * ect2 * -0.5 * ihf[0] * ihf[2]) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.5 * ihf[0] * ihf[2]) / d) * dinot;
            double dfb23 = (2.0 * kmems * model_k * sa2 * ect2 * -0.5 * ihf[1] * ihf[2]) / d - ((2.0 * kmems * eo * sai * ect.abs() * -0.5 * ihf[1] * ihf[2]) / d) * dinot;
             */

            double f = dinot * eo;
            double phi = kect.phase();
            fomphi[i][0] = dinot;
            fomphi[i][1] = phi;
            ComplexNumber mfo = new ComplexNumber(f * cos(phi), f * sin(phi));
            ComplexNumber mfo2 = new ComplexNumber(2.0 * f * cos(phi), 2.0 * f * sin(phi));
            ComplexNumber dfcc = new ComplexNumber(sai * kect.abs() * cos(phi), sai * kect.abs() * sin(phi));
            // map and derivative coefficients
            if (Double.isNaN(fctot[i][0])) {
                fofc1[i][0] = 0.0;
                fofc1[i][1] = 0.0;
                fofc2[i][0] = mfo.re();
                fofc2[i][1] = mfo.im();
                dfc[i][0] = 0.0;
                dfc[i][1] = 0.0;
                dfs[i][0] = 0.0;
                dfs[i][1] = 0.0;
                continue;
            }
            if (Double.isNaN(fo[i][0])) {
                fofc1[i][0] = 0.0;
                fofc1[i][1] = 0.0;
                if (Double.isNaN(fctot[i][0])) {
                    fofc2[i][0] = 0.0;
                    fofc2[i][1] = 0.0;
                } else {
                    fofc2[i][0] = dfcc.re();
                    fofc2[i][1] = dfcc.im();
                }
                dfc[i][0] = 0.0;
                dfc[i][1] = 0.0;
                dfs[i][0] = 0.0;
                dfs[i][1] = 0.0;
                continue;
            }
            fofc1[i][0] = mfo.minus(dfcc).re();
            fofc1[i][1] = mfo.minus(dfcc).im();
            fofc2[i][0] = mfo2.minus(dfcc).re();
            fofc2[i][1] = mfo2.minus(dfcc).im();
            dfc[i][0] = dfcr;
            dfc[i][1] = dfci;
            dfs[i][0] = dfsr;
            dfs[i][1] = dfsi;

            // only use freeR flagged reflections in overall sum
            if (freer[i] == refinementdata.rfreeflag) {
                sum += llk;
                nsum++;
            } else {
                sumr += llk;
                dfsa = dfwa = 0.0;
                nsumr++;
            }

            if (gradient) {
                int i0 = spline.i0();
                int i1 = spline.i1();
                int i2 = spline.i2();
                double g0 = spline.dfi0();
                double g1 = spline.dfi1();
                double g2 = spline.dfi2();

                // s derivative
                g[i0] += dfsa * g0;
                g[i1] += dfsa * g1;
                g[i2] += dfsa * g2;
                // w derivative
                g[n + i0] += dfwa * g0;
                g[n + i1] += dfwa * g1;
                g[n + i2] += dfwa * g2;
            }
        }
        if (print) {
            StringBuffer sb = new StringBuffer("\n");
            sb.append(" sigmaA[s and w] fit using ONLY Rfree reflections\n");
            sb.append(String.format("      # HKL: %10d (free set) %10d (working set) %10d (total)\n",
                    nsum, nsumr, nsum + nsumr));
            sb.append(String.format("   residual: %10g (free set) %10g (working set) %10g (total)\n",
                    sum, sumr, sum + sumr));
            sb.append(String.format("    x: "));
            for (int i = 0; i < x.length; i++) {
                sb.append(String.format("%g ", x[i]));
            }
            sb.append("\n    g: ");
            for (int i = 0; i < g.length; i++) {
                sb.append(String.format("%g ", g[i]));
            }
            sb.append("\n");
            logger.info(sb.toString());
        }
        return sum;
    }

    @Override
    public double energyAndGradient(double x[], double g[]) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, g, true, true);

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i
                    < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }

        return sum;
    }

    @Override
    public void setOptimizationScaling(double[] scaling) {
        if (scaling != null && scaling.length == n * 2) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getOptimizationScaling() {
        return optimizationScaling;
    }

    public static double sim(double x) {
        if (x >= 0.0) {
            return (((x + sim_a) * x + sim_b) * x)
                    / (((x + sim_c) * x + sim_d) * x + sim_e);
        } else {
            return -(-(-(-x + sim_a) * x + sim_b) * x)
                    / (-(-(-x + sim_c) * x + sim_d) * x + sim_e);
        }
    }

    public static double sim_integ(double x0) {
        double x = abs(x0);
        double z = (x + sim_p) / sim_q;

        return sim_A * log(x + sim_g)
                + 0.5 * sim_B * log(z * z + 1.0)
                + sim_r * atan(z)
                + x + 1.0;
    }
}
