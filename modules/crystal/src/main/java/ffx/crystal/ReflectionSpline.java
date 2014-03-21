/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.crystal;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * The ReflectionSpline class represents a reflection spline basis.
 *
 * @author Timothy D. Fenn<br>
 *
 * @see
 * <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2002). 35, 655-663
 * </a>
 *
 * @since 1.0
 */
public class ReflectionSpline {

    private final ReflectionList reflectionList;
    private final Crystal crystal;
    private final int nParams;
    private double f;
    private int i0, i1, i2;
    private double dfi0, dfi1, dfi2;

    /**
     * <p>
     * Constructor for ReflectionSpline.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param nParams a int.
     */
    public ReflectionSpline(ReflectionList reflectionList, int nParams) {
        this.reflectionList = reflectionList;
        this.crystal = reflectionList.crystal;
        this.nParams = nParams;
    }

    /**
     * <p>
     * f</p>
     *
     * @return a double.
     */
    public double f() {
        return f;
    }

    /**
     * <p>
     * i0</p>
     *
     * @return a int.
     */
    public int i0() {
        return i0;
    }

    /**
     * <p>
     * i1</p>
     *
     * @return a int.
     */
    public int i1() {
        return i1;
    }

    /**
     * <p>
     * i2</p>
     *
     * @return a int.
     */
    public int i2() {
        return i2;
    }

    /**
     * <p>
     * dfi0</p>
     *
     * @return a double.
     */
    public double dfi0() {
        return dfi0;
    }

    /**
     * <p>
     * dfi1</p>
     *
     * @return a double.
     */
    public double dfi1() {
        return dfi1;
    }

    /**
     * <p>
     * dfi2</p>
     *
     * @return a double.
     */
    public double dfi2() {
        return dfi2;
    }

    /**
     * Evaluate basis function and derivative at a given resolution (Equations
     * 24 and 25 in Cowtan et al.).
     *
     * @param invressq resolution of desired spline interpolation
     * @param params current spline parameters
     * @return value at invressq
     */
    public double f(double invressq, double params[]) {
        double s = nParams * reflectionList.ordinal(invressq);
        int i = (int) floor(s);
        double ds = s - i - 0.5;
        i0 = min(max(0, i - 1), nParams - 1);
        i1 = min(max(0, i), nParams - 1);
        i2 = min(max(0, i + 1), nParams - 1);

        f = params[i0] * 0.5 * (ds - 0.5) * (ds - 0.5)
                + params[i1] * (0.75 - ds * ds)
                + params[i2] * 0.5 * (ds + 0.5) * (ds + 0.5);

        dfi0 = 0.5 * (ds - 0.5) * (ds - 0.5);
        dfi1 = 0.75 - ds * ds;
        dfi2 = 0.5 * (ds + 0.5) * (ds + 0.5);

        return f;
    }
}
