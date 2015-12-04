/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.xray;

import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
 * @author Timothy D. Fenn<br>
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2002). 35, 655-663</a>
 *
 */
public class SplineEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private final ReflectionList reflectionlist;
    private final ReflectionSpline spline;
    private final int nparams;
    private final int type;
    private final Crystal crystal;
    private final DiffractionRefinementData refinementdata;
    private final double fc[][];
    private final double fo[][];
    protected double[] optimizationScaling = null;
    private final ComplexNumber fct = new ComplexNumber();
    private double totalEnergy;
    private STATE state = STATE.BOTH;

    /**
     * <p>
     * Constructor for SplineEnergy.</p>
     *
     * @param reflectionlist a {@link ffx.crystal.ReflectionList} object.
     * @param refinementdata a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param nparams a int.
     * @param type a int.
     */
    public SplineEnergy(ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, int nparams, int type) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.refinementdata = refinementdata;
        this.type = type;
        this.fc = refinementdata.fc;
        this.fo = refinementdata.fsigf;

        // initialize params
        this.spline = new ReflectionSpline(reflectionlist, nparams);
        this.nparams = nparams;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * <p>
     * target</p>
     *
     * @param x an array of double.
     * @param g an array of double.
     * @param gradient a boolean.
     * @param print a boolean.
     * @return a double.
     */
    public double target(double x[], double g[], boolean gradient, boolean print) {

        double r = 0.0;
        double rf = 0.0;
        double rfree = 0.0;
        double rfreef = 0.0;
        double sum = 0.0;
        double sumfo = 0.0;

        // Zero out the gradient.
        if (gradient) {
            fill(g, 0.0);
        }

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
            refinementdata.getFcTotIP(i, fct);

            double d2 = 0.0;
            double dr = 0.0;
            double w = 0.0;
            switch (type) {
                case Type.FOFC:
                    w = 1.0;
                    double f1 = refinementdata.getF(i);
                    double f2 = fct.abs();
                    double d = f1 - fh * f2;
                    d2 = d * d;
                    dr = -2.0 * f2 * d;
                    sumfo += f1 * f1;
                    break;
                case Type.F1F2:
                    w = 2.0 / ih.epsilonc();
                    double ieps = 1.0 / eps;
                    f1 = pow(fct.abs(), 2.0) * ieps;
                    f2 = pow(refinementdata.getF(i), 2) * ieps;
                    d = fh * f1 - f2;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
                case Type.FCTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(fct.abs() / sqrt(eps), 2);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
                case Type.FOTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(refinementdata.getF(i) / sqrt(eps), 2);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    sumfo = 1.0;
                    break;
            }

            sum += w * d2;

            double afo = abs(fo[i][0]);
            double afh = abs(fh * fct.abs());
            if (refinementdata.isFreeR(i)) {
                rfree += abs(afo - afh);
                rfreef += afo;
            } else {
                r += abs(afo - afh);
                rf += afo;
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

        /**
         * Tim - should this only be done for Type.FOFC??
         */
        if (gradient) {
            double isumfo = 1.0 / sumfo;
            for (int i = 0; i < g.length; i++) {
                g[i] *= isumfo;
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

        totalEnergy = sum / sumfo;
        return sum / sumfo;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double x[]) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, null, false, false);

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }

        return sum;
    }

    /**
     * {@inheritDoc}
     */
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

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling) {
        if (scaling != null && scaling.length == nparams) {
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
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return null;
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
    }

    @Override
    public void setVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static interface Type {

        public static final int FOFC = 1;
        public static final int F1F2 = 2;
        public static final int FCTOESQ = 3;
        public static final int FOTOESQ = 4;
    }
}
