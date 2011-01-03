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
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.numerics.Potential;

/**
 *
 * Fit structure factors using spline coefficients
 *
 * @author Tim Fenn<br>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2002). 35, 655-663</a>
 *
 */
public class SplineEnergy implements Potential {

    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static interface Type {

        public static final int FOFC = 1;
        public static final int F1F2 = 2;
        public static final int FCTOESQ = 3;
        public static final int FOTOESQ = 4;
    }
    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private final ReflectionList reflectionlist;
    private final ReflectionSpline spline;
    private final int nparams;
    private final int type;
    private final Crystal crystal;
    private final RefinementData refinementdata;
    private final double fc[][];
    private final double fs[][];
    private final double fctot[][];
    private final double fo[][];
    private final int freer[];
    protected double[] optimizationScaling = null;
    private ComplexNumber fct = new ComplexNumber();

    public SplineEnergy(ReflectionList reflectionlist,
            RefinementData refinementdata, int nparams, int type) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.refinementdata = refinementdata;
        this.type = type;
        this.fc = refinementdata.fc;
        this.fs = refinementdata.fs;
        this.fctot = refinementdata.fctot;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;

        // initialize params
        this.spline = new ReflectionSpline(reflectionlist, nparams);
        this.nparams = nparams;
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

        r = rf = rfree = rfreef = sum = sumfo = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            if (type == Type.FOTOESQ
                    && fo[i][0] <= 0.0) {
                continue;
            }

            double eps = ih.epsilon();
            double s = Crystal.invressq(crystal, ih);

            // spline setup
            double fh = spline.f(s, x);

            refinementdata.get_fctot_ip(i, fct);

            double f1, f2, d, d2, dr, w;
            f1 = f2 = d = d2 = dr = w = 0.0;
            switch (type) {
                case Type.FOFC:
                    w = 1.0;
                    f1 = refinementdata.get_f(i);
                    f2 = fct.abs();
                    d = f1 - fh * f2;
                    d2 = d * d;
                    dr = -2.0 * f2 * d;
                    sumfo += f1 * f1;
                    break;
                case Type.F1F2:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(fct.abs(), 2.0) / eps;
                    f2 = pow(refinementdata.get_f(i), 2.0) / eps;
                    d = fh * f1 - f2;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
                case Type.FCTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(fct.abs() / sqrt(eps), 2.0);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
                case Type.FOTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(refinementdata.get_f(i) / sqrt(eps), 2.0);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
            }

            sum += w * d2;

            if (refinementdata.isfreer(i)) {
                rfree += abs(abs(fo[i][0]) - abs(fh * fct.abs()));
                rfreef += abs(fo[i][0]);
            } else {
                r += abs(abs(fo[i][0]) - abs(fh * fct.abs()));
                rf += abs(fo[i][0]);
            }

            if (gradient) {
                int i0 = spline.i0();
                int i1 = spline.i1();
                int i2 = spline.i2();
                double g0 = spline.dfi0();
                double g1 = spline.dfi1();
                double g2 = spline.dfi2();

                g[i0] += w * dr * g0;
                g[i1] += w * dr * g1;
                g[i2] += w * dr * g2;
            }
        }

        if (gradient) {
            for (int i = 0; i < g.length; i++) {
                g[i] /= sumfo;
            }
        }

        if (print) {
            StringBuilder sb = new StringBuilder("\n");
            sb.append(" Computed Potential Energy\n");
            sb.append(String.format("   residual:  %8.3f\n",
                    sum / sumfo));
            if (type == Type.FOFC || type == Type.F1F2) {
                sb.append(String.format("   R:  %8.3f  Rfree:  %8.3f\n",
                        (r / rf) * 100.0, (rfree / rfreef) * 100.0));
            }
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
    public double energyAndGradient(double x[], double g[]) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, g, true, false);

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
    public void setScaling(double[] scaling) {
        if (scaling != null && scaling.length == nparams) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getScaling() {
        return optimizationScaling;
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
