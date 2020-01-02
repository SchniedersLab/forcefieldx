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
package ffx.potential.parameters;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.sort;

import static org.apache.commons.math3.util.FastMath.abs;

import static ffx.potential.parameters.ForceField.ForceFieldType.TORTORS;

/**
 * The TorsionTorsionType class defines a Torsion-Torsion spline.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class TorsionTorsionType extends BaseType implements Comparator<String> {

    private static final Logger logger = Logger.getLogger(TorsionTorsionType.class.getName());

    /**
     * Atom classes that form this Torsion-Torsion type.
     */
    public final int[] atomClasses;
    /**
     * Energy values.
     */
    public final double[] energy;
    /**
     * Number of points along x.
     */
    public final int nx;
    /**
     * Number of point along y.
     */
    public final int ny;
    /**
     * Torsion values along x.
     */
    public final double[] tx;
    /**
     * Torsion values along y.
     */
    public final double[] ty;
    /**
     * First derivative along x.
     */
    public final double[] dx;
    /**
     * First derivative along y.
     */
    public final double[] dy;
    /**
     * Second derivatives.
     */
    public final double[] dxy;
    /**
     * Grid points.
     */
    private final int[] gridPoints;
    /**
     * Convert Torsion-Torsion energy to kcal/mole.
     */
    public static final double units = 1.0;

    /**
     * <p>
     * Constructor for TorsionTorsionType.</p>
     *
     * @param atomClasses an array of int.
     * @param gridPoints  an array of int.
     * @param torsion1    an array of double.
     * @param torsion2    an array of double.
     * @param energy      an array of double.
     */
    public TorsionTorsionType(int[] atomClasses, int[] gridPoints,
                              double[] torsion1, double[] torsion2, double[] energy) {
        super(TORTORS, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        nx = gridPoints[0];
        ny = gridPoints[1];
        if (nx != ny) {
            logger.severe("Untested TORTOR parameters: nx != ny: " + nx + ", " + ny);
        }

        this.energy = energy;
        this.gridPoints = gridPoints;
        tx = new double[nx];
        ty = new double[ny];
        dx = new double[nx * ny];
        dy = new double[nx * ny];
        dxy = new double[nx * ny];
        sort(torsion1);
        sort(torsion2);
        tx[0] = torsion1[0];
        ty[0] = torsion2[0];
        int j1 = 1;
        int j2 = 1;
        for (int i = 1; i < nx; i++) {
            while (torsion1[j1] == tx[i - 1]) {
                j1++;
            }
            while (torsion2[j2] == ty[i - 1]) {
                j2++;
            }
            tx[i] = torsion1[j1];
            ty[i] = torsion2[j2];
        }

        // Check for cyclic energy.
        boolean isCyclic = true;
        double eps = 0.0001;
        if (abs(tx[0] - tx[nx - 1]) - 360.0 > eps) {
            isCyclic = false;
            if (logger.isLoggable(Level.FINEST)) {
                logger.finest(" tortor is aperiodic: " + tx[0] + ", " + tx[nx - 1]);
            }
        }
        if (isCyclic) {
            for (int i = 0; i < ny; i++) {
                int k = i * nx;
                if (abs(energy[k] - energy[k + nx - 1]) > eps) {
                    isCyclic = false;
                    if (logger.isLoggable(Level.FINEST)) {
                        logger.finest(" tortor is apreriodic: " + k + ", " + (k + nx - 1) + ": " + abs(energy[k] - energy[k + nx - 1]));
                    }
                    break;
                }
            }
        }
        if (isCyclic) {
            int k = (ny - 1) * nx;
            for (int i = 0; i < nx; i++) {
                if (abs(energy[i] - energy[i + k]) > eps) {
                    if (logger.isLoggable(Level.FINEST)) {
                        logger.fine(" tortor is aperiodic: " + i + ", " + i + k + ": " + abs(energy[i] - energy[i + k]));
                    }
                    isCyclic = false;
                    break;
                }
            }
        }

        boolean cyclic = isCyclic;
        double[] tmp1 = new double[nx];
        double[] tmp2 = new double[nx];
        double[] tmp3 = new double[nx];
        double[] tmp4 = new double[nx];
        double[] tmp5 = new double[nx];
        double[] tmp6 = new double[nx];
        double[] tmp7 = new double[nx];
        double[] bs = new double[nx];
        double[] cs = new double[nx];
        double[] ds = new double[nx];

        // Spline fit the derivatives about the first torsion.
        arraycopy(tx, 0, tmp1, 0, nx);

        int m = 0;
        for (int j = 0; j < ny; j++) {
            arraycopy(energy, m, tmp2, 0, nx);
            if (cyclic) {
                cspline(nx - 1, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, tmp5, tmp6, tmp7);
            } else {
                nspline(nx - 1, tmp1, tmp2, 0.0, 0.0, bs, cs, tmp3, tmp4, tmp5, tmp6, tmp7);
            }
            arraycopy(bs, 0, dx, m, nx);
            m = m + nx;
        }

        // Spline fit the derivatives about the second torsion.
        arraycopy(ty, 0, tmp1, 0, ny);

        m = 0;
        for (int j = 0; j < nx; j++) {
            for (int k = 0; k < ny; k++) {
                tmp2[k] = energy[m + k * nx];
            }
            if (cyclic) {
                cspline(ny - 1, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, tmp5, tmp6, tmp7);
            } else {
                nspline(ny - 1, tmp1, tmp2, 0.0, 0.0, bs, cs, tmp3, tmp4, tmp5, tmp6, tmp7);
            }
            for (int k = 0; k < ny; k++) {
                dy[m + k * nx] = bs[k];
            }
            m = m + 1;
        }

        // Spline fit the cross derivatives about both torsions.
        m = 0;
        for (int j = 0; j < nx; j++) {
            for (int k = 0; k < ny; k++) {
                tmp2[k] = dx[m + k * nx];
            }
            if (cyclic) {
                cspline(ny - 1, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, tmp5, tmp6, tmp7);
            } else {
                nspline(ny - 1, tmp1, tmp2, 0.0, 0.0, bs, cs, tmp3, tmp4, tmp5, tmp6, tmp7);
            }
            for (int k = 0; k < ny; k++) {
                dxy[m + k * nx] = bs[k];
            }
            m = m + 1;
        }
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    public void incrementClasses(int increment) {
        for (int i = 0; i < atomClasses.length; i++) {
            atomClasses[i] += increment;
        }
        setKey(sortKey(atomClasses));
    }

    /**
     * Remap new atom classes to known internal ones.
     *
     * @param typeMap a lookup between new atom types and known atom types.
     */
    public void patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        for (AtomType newType : typeMap.keySet()) {

            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
                    count++;
                }
            }
        }
        if (count > 0 && count < atomClasses.length) {
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < atomClasses.length; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        atomClasses[i] = knownType.atomClass;
                    }
                }

            }
            setKey(sortKey(atomClasses));
        }
    }

    /**
     * Computes the coefficients for an aperiodic interpolating cubic spline.
     *
     * @param n
     * @param x0
     * @param y0
     * @param y21
     * @param y2n
     * @param s1
     * @param s2
     * @param h
     * @param g
     * @param dy
     * @param dla
     * @param dmu
     */
    private void nspline(int n, double[] x0, double[] y0, double y21, double y2n, double[] s1,
                         double[] s2, double[] h, double[] g, double[] dy, double[] dla, double[] dmu) {

        // Calculate the intervals.
        for (int i = 0; i < n; i++) {
            h[i] = x0[i + 1] - x0[i];
            dy[i] = (y0[i + 1] - y0[i]) / h[i];
        }

        // Get the coeffcients.
        for (int i = 1; i < n; i++) {
            dla[i] = h[i] / (h[i] + h[i - 1]);
            dmu[i] = 1.0 - dla[i];
            g[i] = 3.0 * (dla[i] * dy[i - 1] + dmu[i] * dy[i]);
        }

        // Set the initial value via natural boundary condition.
        dla[n] = 1.0;
        dla[0] = 0.0;
        dmu[n] = 0.0;
        dmu[0] = 1.0;
        g[0] = 3.0 * dy[0] - 0.5 * h[0] * y21;
        g[n] = 3.0 * dy[n - 1] + 0.5 * h[n - 1] * y2n;

        // Solve the triagonal system of linear equations.
        dmu[0] = 0.5 * dmu[0];
        g[0] = 0.5 * g[0];
        for (int i = 1; i <= n; i++) {
            double t = 2.0 - dmu[i - 1] * dla[i];
            dmu[i] = dmu[i] / t;
            g[i] = (g[i] - g[i - 1] * dla[i]) / t;
        }

        for (int i = n - 1; i >= 0; i--) {
            g[i] = g[i] - dmu[i] * g[i + 1];
        }

        // Get the first derivative at each grid point.
        arraycopy(g, 0, s1, 0, n + 1);

        // Get the second derivative at each grid point.
        s2[0] = y21;
        s2[n] = y2n;
        for (int i = 1; i < n; i++) {
            s2[i] = 6.0 * (y0[i + 1] - y0[i]) / (h[i] * h[i])
                    - 4.0 * s1[i] / h[i] - 2.0 * s1[i + 1] / h[i];
        }

    }

    /**
     * Computes the coefficients for a periodic interpolating cubic spline.
     *
     * @param n
     * @param xn
     * @param fn
     * @param b
     * @param d
     * @param h
     * @param du
     * @param dm
     * @param rc
     * @param rs
     */
    private void cspline(int n, double[] xn, double[] fn, double[] b,
                         double[] c, double[] d, double[] h, double[] du, double[] dm,
                         double[] rc, double[] rs) {
        double eps = 0.000001;
        if (abs(fn[n] - fn[0]) > eps) {
            logger.severe("TORTOR values are not periodic.");
        }

        // Make the spline exactly periodic making the ends equal to their mean.
        double mean = 0.5 * (fn[0] + fn[n]);
        fn[0] = mean;
        fn[n] = mean;

        // Set up auxiliary variables and matrix elements on first call.
        for (int i = 0; i < n; i++) {
            h[i] = xn[i + 1] - xn[i];
        }
        h[n] = h[0];

        if (n - 1 >= 0) arraycopy(h, 1, du, 1, n - 1);

        du[n] = h[0];
        for (int i = 1; i <= n; i++) {
            dm[i] = 2.0 * (h[i - 1] + h[i]);
        }

        // Compute the RHS.
        double temp1 = (fn[1] - fn[0]) / h[0];
        for (int i = 1; i < n; i++) {
            double temp2 = (fn[i + 1] - fn[i]) / h[i];
            rs[i] = 3.0 * (temp2 - temp1);
            temp1 = temp2;
        }
        rs[n] = 3.0 * ((fn[1] - fn[0]) / h[0] - temp1);

        // Solve the linear system with factorization.
        if (cytsy(n, dm, du, rc, rs, c) != 1) {
            return;
        }

        // Compute remaining spline coefficients.
        c[0] = c[n];
        for (int i = 0; i < n; i++) {
            b[i] = (fn[i + 1] - fn[i]) / h[i] - h[i] / 3.0 * (c[i + 1] + 2.0 * c[i]);
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
        b[n] = (fn[1] - fn[n]) / h[n] - h[n] / 3.0 * (c[1] + 2.0 * c[n]);
    }

    /**
     * <p>average.</p>
     *
     * @param torsionTorsionType1 a {@link ffx.potential.parameters.TorsionTorsionType} object.
     * @param torsionTorsionType2 a {@link ffx.potential.parameters.TorsionTorsionType} object.
     * @param atomClasses         an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.TorsionTorsionType} object.
     */
    public static TorsionTorsionType average(TorsionTorsionType torsionTorsionType1,
                                             TorsionTorsionType torsionTorsionType2, int[] atomClasses) {
        if (torsionTorsionType1 == null || torsionTorsionType2 == null || atomClasses == null) {
            return null;
        }
        return null;
    }

    /**
     * Solves a system of linear equations for a cyclically tridiagonal,
     * symmetric, positive definite matrix.
     *
     * @param n
     * @param dm
     * @param du
     * @param rs
     * @param c
     * @return positive or negative 1.
     */
    private static int cytsy(int n, double[] dm, double[] du, double[] cr, double[] rs, double[] c) {
        if (n < 3) {
            return -2;
        }

        // Factorization of the input matrix.
        if (cytsyp(n, dm, du, cr) == 1) {
            // Update and back-substitute as necessary.
            cytsys(n, dm, du, cr, rs, c);
            return 1;
        }
        return -1;
    }

    /**
     * Finds the Cholesky factors of a cyclically tridiagonal symmetric,
     * positive definite matrix given by two vectors.
     *
     * @param n
     * @param dm
     * @param du
     * @param cr
     * @return an integer.
     */
    private static int cytsyp(int n, double[] dm, double[] du, double[] cr) {
        double eps = 0.00000001;
        // Check for n < 3.
        if (n < 3) {
            return -2;
        }

        // Check if the matrix is positive definite.
        double row = abs(dm[1]) + abs(du[1]) + abs(du[n]);
        if (row == 0.0) {
            return 0;
        }
        double d = 1.0 / row;
        if (dm[1] < 0.0) {
            return -1;
        } else if (abs(dm[1]) * d < eps) {
            return 0;
        }

        // Factoring while checking for a positive definite and strong non-singular matrix.
        double temp1 = du[1];
        double temp2 = 0.0;
        du[1] = du[1] / dm[1];
        cr[1] = du[n] / dm[1];
        for (int i = 2; i < n; i++) {
            row = abs(dm[i]) + abs(du[i]) + abs(temp1);
            if (row == 0.0) {
                return 0;
            }
            d = 1.0 / row;
            dm[i] = dm[i] - temp1 * du[i - 1];
            if (dm[i] < 0.0) {
                return -1;
            } else if (abs(dm[i]) * d < eps) {
                return 0;
            }
            if (i < (n - 1)) {
                cr[i] = -temp1 * cr[i - 1] / dm[i];
                temp1 = du[i];
                du[i] = du[i] / dm[i];
            } else {
                temp2 = du[i];
                du[i] = (du[i] - temp1 * cr[i - 1]) / dm[i];
            }
        }
        row = abs(du[n]) + abs(dm[n]) + abs(temp2);
        if (row == 0.0) {
            return 0;
        }
        d = 1.0 / row;
        dm[n] = dm[n] - dm[n - 1] * du[n - 1] * du[n - 1];
        temp1 = 0.0;
        for (int i = 1; i < n - 1; i++) {
            temp1 = temp1 + dm[i] * cr[i] * cr[i];
        }
        dm[n] = dm[n] - temp1;
        if (dm[n] < 0.0) {
            return -1;
        } else if (abs(dm[n]) * d < eps) {
            return 0;
        }
        return 1;
    }

    /**
     * Solves a cyclically tridiagonal linear system given the Cholesky factors.
     *
     * @param n
     * @param dm
     * @param du
     * @param cr
     * @param rs
     * @param c
     */
    private static void cytsys(int n, double[] dm, double[] du, double[] cr, double[] rs, double[] c) {
        // Updating phase.
        double temp = rs[1];
        rs[1] = temp / dm[1];
        double sum = cr[1] * temp;
        for (int i = 2; i < n; i++) {
            temp = rs[i] - du[i - 1] * temp;
            rs[i] = temp / dm[i];
            if (i != (n - 1)) {
                sum = sum + cr[i] * temp;
            }
        }
        temp = rs[n] - du[n - 1] * temp;
        temp = temp - sum;
        rs[n] = temp / dm[n];

        // Back-substitution phase.
        c[n] = rs[n];
        c[n - 1] = rs[n - 1] - du[n - 1] * c[n];
        for (int i = n - 2; i >= 1; i--) {
            c[i] = rs[i] - du[i] * c[i + 1] - cr[i] * c[n];
        }
    }

    /**
     * No sorting is done for the Torsion-Torsion lookup.
     *
     * @param c atomClasses
     * @return lookup key
     */
    public static String sortKey(int[] c) {
        return c[0] + " " + c[1] + " " + c[2] + " " + c[3] + " " + c[4];
    }

    /**
     * Reversed key for the Torsion-Torsion lookup.
     *
     * @param c atomClasses
     * @return lookup key
     */
    public static String reverseKey(int[] c) {
        return c[4] + " " + c[3] + " " + c[2] + " " + c[1] + " " + c[0];
    }

    /**
     * Construct a TorsionTorsionType from multiple input lines.
     *
     * @param input  The overall input String.
     * @param tokens The input String tokenized.
     * @param br     a BufferedReader instance.
     * @return a TorsionTorsionType instance.
     */
    public static TorsionTorsionType parse(String input, String[] tokens, BufferedReader br) {
        if (tokens.length < 8) {
            logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[5];
                for (int i = 0; i < 5; i++) {
                    atomClasses[i] = parseInt(tokens[i + 1]);
                }
                int[] gridPoints = new int[2];
                gridPoints[0] = parseInt(tokens[6]);
                gridPoints[1] = parseInt(tokens[7]);
                int points = gridPoints[0] * gridPoints[1];
                double[] torsion1 = new double[points];
                double[] torsion2 = new double[points];
                double[] energy = new double[points];
                for (int i = 0; i < points; i++) {
                    input = br.readLine();
                    tokens = input.trim().split(" +");
                    if (tokens.length != 3) {
                        logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
                        return null;
                    }
                    torsion1[i] = parseDouble(tokens[0]);
                    torsion2[i] = parseDouble(tokens[1]);
                    energy[i] = parseDouble(tokens[2]);
                }
                return new TorsionTorsionType(atomClasses, gridPoints, torsion1, torsion2, energy);
            } catch (NumberFormatException | IOException e) {
                String message = "Exception parsing TORTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * Construct a TorsionTorsionType from a single input line.
     *
     * @param input  The overall input String.
     * @param tokens The input String tokenized.
     * @return a TorsionTorsionType instance.
     */
    public static TorsionTorsionType parse(String input, String[] tokens) {
        if (tokens.length < 8) {
            logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[5];
                for (int i = 0; i < 5; i++) {
                    atomClasses[i] = parseInt(tokens[i + 1]);
                }
                int[] gridPoints = new int[2];
                gridPoints[0] = parseInt(tokens[6]);
                gridPoints[1] = parseInt(tokens[7]);
                int points = gridPoints[0] * gridPoints[1];
                int numTokens = points * 3 + 8;
                if (tokens.length < numTokens) {
                    logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
                    return null;
                }
                double[] torsion1 = new double[points];
                double[] torsion2 = new double[points];
                double[] energy = new double[points];
                int index = 8;
                for (int i = 0; i < points; i++) {
                    torsion1[i] = parseDouble(tokens[index++]);
                    torsion2[i] = parseDouble(tokens[index++]);
                    energy[i] = parseDouble(tokens[index++]);
                }
                return new TorsionTorsionType(atomClasses, gridPoints, torsion1, torsion2, energy);
            } catch (NumberFormatException e) {
                String message = "Exception parsing TORTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted torsion-torsion type.
     */
    @Override
    public String toString() {
        StringBuilder tortorBuffer = new StringBuilder("tortors");
        for (int i : atomClasses) {
            tortorBuffer.append(format("  %5d", i));
        }
        tortorBuffer.append(format("  %2d  %2d", gridPoints[0],
                gridPoints[1]));
        for (int i = 0; i < energy.length; i++) {
            int nxi = i % nx;
            int nyi = i / ny;
            tortorBuffer.append(format(" \\\n  % 6.1f  % 6.1f  % 8.5f",
                    tx[nxi], ty[nyi], energy[i]));
        }
        return tortorBuffer.toString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String key1, String key2) {
        String[] keys1 = key1.split(" ");
        String[] keys2 = key2.split(" ");
        int c1 = parseInt(keys1[2]);
        int c2 = parseInt(keys2[2]);
        if (c1 < c2) {
            return -1;
        } else if (c1 > c2) {
            return 1;
        }
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        TorsionTorsionType torsionTorsionType = (TorsionTorsionType) o;
        return Arrays.equals(atomClasses, torsionTorsionType.atomClasses);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(atomClasses);
    }
}
