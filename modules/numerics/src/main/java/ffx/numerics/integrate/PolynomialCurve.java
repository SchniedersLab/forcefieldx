/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.numerics.integrate;

/**
 * A PolynomialCurve describes points along a polynomial function.
 *
 * @author Jacob M. Litman
 */
public class PolynomialCurve extends FunctionDataCurve {
    
    private final double[] coeff;

    /**
     * Default constructor, assumes constant-width bins. Functional form will
     * be a0 + a1x + a2x^2 + a3x^3 + ... + anx^n.
     *
     * @param x an array of {@link double} objects.
     * @param coefficients Lowest-order coefficients first
     */
    public PolynomialCurve(double[] x, double[] coefficients) {
        this(x, false, coefficients);
    }
    

    /**
     * Default constructor, assumes constant-width bins. Functional form will
     * be a0 + a1x + a2x^2 + a3x^3 + ... + anx^n.
     *
     * @param x an array of {@link double} objects.
     * @param halfWidthEnds Specifies that first and last bins are half-width.l
     * @param coefficients Lowest-order coefficients first
     */
    public PolynomialCurve(double[] x, boolean halfWidthEnds, double[] coefficients) {
        int npoints = x.length;
        points = new double[npoints];
        coeff = new double[coefficients.length];
        System.arraycopy(coefficients, 0, coeff, 0, coefficients.length);
        this.halfWidthEnd = halfWidthEnds;
        
        for (int i = 0; i < points.length; i++) {
            points[i] = polynomialAt(x[i]);
        }
        lb = x[0];
        ub = x[npoints-1];
        assertXIntegrity(x);
        this.x = new double[x.length];
        System.arraycopy(x, 0, this.x, 0, x.length);
    }
    
    /** {@inheritDoc} */
    @Override
    public double integralAt(double x) {
        //double total = x * coeff[0];
        double total = 0;
        for (int i = 0; i < coeff.length; i++) {
            double val = 1.0 / ((double) i+1);
            val *= coeff[i];
            for (int j = 0; j <= i; j++) {
                val *= x;
            }
            total += val;
        }
        return total;
    }
    
    /** {@inheritDoc} */
    @Override
    public double fX(double x) {
        return polynomialAt(x);
    }
    
    // Private, non-overrideable method for use in the constructor.
    private double polynomialAt(double x) {
        double total = 0.0;
        for (int i = 0; i < coeff.length; i++) {
            double val = coeff[i];
            for (int j = 0; j < i; j++) {
                val *= x;
            }
            total += val;
        }
        return total;
    }
    
    /** {@inheritDoc} */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Polynomial curve of degree ");
        sb.append(coeff.length).append(String.format(" with %d points from lower bound %9.3g and upper bound %9.3g", points.length, lb, ub));
        if (halfWidthEnd) {
            sb.append(" and half-width start/end bins");
        }
        sb.append(".\nCoefficients: ");
        if (coeff.length > 0) {
            sb.append(coeff[0]);
        }
        for (int i = 1; i < coeff.length; i++) {
            sb.append(",").append(coeff[i]);
        }
        
        return sb.toString();
    }
}
