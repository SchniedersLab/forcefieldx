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

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.Potential;
import ffx.numerics.math.ComplexNumber;
import static ffx.crystal.Crystal.invressq;

/**
 * Fit structure factors using spline coefficients
 *
 * @author Timothy D. Fenn<br>
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2002). 35, 655-663</a>
 * @since 1.0
 */
public class SplineEnergy implements Potential {

    public interface Type {

        int FOFC = 1;
        int F1F2 = 2;
        int FCTOESQ = 3;
        int FOTOESQ = 4;
    }

    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private final ReflectionList reflectionList;
    private final ReflectionSpline spline;
    private final int nParams;
    private final int type;
    private final Crystal crystal;
    private final DiffractionRefinementData refinementData;
    private final double[][] fc;
    private final double[][] fo;
    private double[] optimizationScaling = null;
    private final ComplexNumber fct = new ComplexNumber();
    private double totalEnergy;
    private STATE state = STATE.BOTH;

    /**
     * <p>
     * Constructor for SplineEnergy.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
     * @param nParams        a int.
     * @param type           a int.
     */
    SplineEnergy(ReflectionList reflectionList,
                        DiffractionRefinementData refinementData, int nParams, int type) {
        this.reflectionList = reflectionList;
        this.crystal = reflectionList.crystal;
        this.refinementData = refinementData;
        this.type = type;
        this.fc = refinementData.fc;
        this.fo = refinementData.fSigF;

        // initialize params
        this.spline = new ReflectionSpline(reflectionList, nParams);
        this.nParams = nParams;
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

        double r = 0.0;
        double rf = 0.0;
        double rfree = 0.0;
        double rfreef = 0.0;
        double sum = 0.0;
        double sumfo = 1.0;

        // Zero out the gradient.
        if (gradient) {
            fill(g, 0.0);
        }

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            if (isNaN(fc[i][0]) || isNaN(fo[i][0]) || fo[i][1] <= 0.0) {
                continue;
            }

            if (type == Type.FOTOESQ
                    && fo[i][0] <= 0.0) {
                continue;
            }

            double eps = ih.epsilon();
            double s = invressq(crystal, ih);
            // Spline setup
            double fh = spline.f(s, x);
            refinementData.getFcTotIP(i, fct);

            double d2 = 0.0;
            double dr = 0.0;
            double w = 0.0;
            switch (type) {
                case Type.FOFC:
                    w = 1.0;
                    double f1 = refinementData.getF(i);
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
                    f2 = pow(refinementData.getF(i), 2) * ieps;
                    d = fh * f1 - f2;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    break;
                case Type.FCTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(fct.abs() / sqrt(eps), 2);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    break;
                case Type.FOTOESQ:
                    w = 2.0 / ih.epsilonc();
                    f1 = pow(refinementData.getF(i) / sqrt(eps), 2);
                    d = f1 * fh - 1.0;
                    d2 = d * d / f1;
                    dr = 2.0 * d;
                    break;
            }

            sum += w * d2;

            double afo = abs(fo[i][0]);
            double afh = abs(fh * fct.abs());
            if (refinementData.isFreeR(i)) {
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

        if (gradient && type == Type.FOFC) {
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
        if (scaling != null && scaling.length == nParams) {
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
     * <p>
     * Return a reference to each variables type.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
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
    public void setEnergyTermState(STATE state) {
        this.state = state;
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

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        // Should be handled upstream.
        return true;
    }

}
