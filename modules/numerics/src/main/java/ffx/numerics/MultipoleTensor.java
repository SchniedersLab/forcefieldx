/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.numerics;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.PI;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.numerics.Erf.erfc;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.doubleFactorial;
import static ffx.numerics.VectorMath.norm;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;

/**
 * The MultipoleTensor class computes derivatives of 1/|<b>r</b>| via recursion
 * to arbitrary order for Cartesian multipoles in either a global frame or a
 * quasi-internal frame.
 *
 * @author Michael J. Schnieders
 *
 * @see
 * <a href="http://doi.org/10.1142/9789812830364_0002"
 * target="_blank">
 * Matt Challacombe, Eric Schwegler and Jan Almlof, Modern developments in
 * Hartree-Fock theory: Fast methods for computing the Coulomb matrix.
 * Computational Chemistry: Review of Current Trends. pp. 53-107, Ed. J.
 * Leczszynski, World Scientifc, 1996.
 * </a>
 *
 * @since 1.0
 */
public class MultipoleTensor {

    private static final Logger logger = Logger.getLogger(MultipoleTensor.class.getName());

    /**
     * Operators that are supported.
     */
    public enum OPERATOR {
        COULOMB, SCREENED_COULOMB, THOLE_FIELD
    };

    private OPERATOR operator;
    private final int order;
    /**
     * These are the "source" terms for the recursion. Source terms exist for
     * for:<br>
     * 1) the Coulomb operator (1/R).<br>
     * 2) the screened Coulomb operator erfc(R)/R.
     *
     */
    private final double T000j[];
    private double beta;
    private double damp;
    private double aiak;
    private final int o1;
    private final int il;
    private final int im;
    private final int in;
    private final int size;
    private final double sqrtPI = sqrt(PI);
    /**
     * Store the auxillary tensor memory to avoid memory consumption.
     */
    private final double T000[];
    /**
     * Store the work array to avoid memory consumption. Note that rather than
     * use an array for intermediate values, a 4D matrix was tried. It was
     * approximately 50% slower than the linear work array.
     */
    private final double work[];

    /**
     * <p>
     * Constructor for MultipoleTensor.</p>
     *
     * @param operator The tensor operator.
     * @param order The order of the tensor.
     * @param beta The screening parameter.
     */
    public MultipoleTensor(OPERATOR operator, int order, double beta) {
        assert (order > 0);
        o1 = order + 1;
        il = o1;
        im = il * o1;
        in = im * o1;
        size = (order + 1) * (order + 2) * (order + 3) / 6;
        work = new double[in * o1];
        T000j = new double[o1];
        this.order = order;
        this.operator = operator;
        this.beta = beta;

        setOperator(operator);

        T000 = new double[order + 1];
        // l + m + n = 0 (1)
        t000 = ti(0, 0, 0, order);
        // l + m + n = 1 (3)   4
        t100 = ti(1, 0, 0, order);
        t010 = ti(0, 1, 0, order);
        t001 = ti(0, 0, 1, order);
        // l + m + n = 2 (6)  10
        t200 = ti(2, 0, 0, order);
        t020 = ti(0, 2, 0, order);
        t002 = ti(0, 0, 2, order);
        t110 = ti(1, 1, 0, order);
        t101 = ti(1, 0, 1, order);
        t011 = ti(0, 1, 1, order);
        // l + m + n = 3 (10) 20
        t300 = ti(3, 0, 0, order);
        t030 = ti(0, 3, 0, order);
        t003 = ti(0, 0, 3, order);
        t210 = ti(2, 1, 0, order);
        t201 = ti(2, 0, 1, order);
        t120 = ti(1, 2, 0, order);
        t021 = ti(0, 2, 1, order);
        t102 = ti(1, 0, 2, order);
        t012 = ti(0, 1, 2, order);
        t111 = ti(1, 1, 1, order);
        // l + m + n = 4 (15) 35
        t400 = ti(4, 0, 0, order);
        t040 = ti(0, 4, 0, order);
        t004 = ti(0, 0, 4, order);
        t310 = ti(3, 1, 0, order);
        t301 = ti(3, 0, 1, order);
        t130 = ti(1, 3, 0, order);
        t031 = ti(0, 3, 1, order);
        t103 = ti(1, 0, 3, order);
        t013 = ti(0, 1, 3, order);
        t220 = ti(2, 2, 0, order);
        t202 = ti(2, 0, 2, order);
        t022 = ti(0, 2, 2, order);
        t211 = ti(2, 1, 1, order);
        t121 = ti(1, 2, 1, order);
        t112 = ti(1, 1, 2, order);
        // l + m + n = 5 (21) 56
        t500 = ti(5, 0, 0, order);
        t050 = ti(0, 5, 0, order);
        t005 = ti(0, 0, 5, order);
        t410 = ti(4, 1, 0, order);
        t401 = ti(4, 0, 1, order);
        t140 = ti(1, 4, 0, order);
        t041 = ti(0, 4, 1, order);
        t104 = ti(1, 0, 4, order);
        t014 = ti(0, 1, 4, order);
        t320 = ti(3, 2, 0, order);
        t302 = ti(3, 0, 2, order);
        t230 = ti(2, 3, 0, order);
        t032 = ti(0, 3, 2, order);
        t203 = ti(2, 0, 3, order);
        t023 = ti(0, 2, 3, order);
        t311 = ti(3, 1, 1, order);
        t131 = ti(1, 3, 1, order);
        t113 = ti(1, 1, 3, order);
        t221 = ti(2, 2, 1, order);
        t212 = ti(2, 1, 2, order);
        t122 = ti(1, 2, 2, order);
    }

    /**
     * Set the Operator.
     *
     * @param operator
     */
    public final void setOperator(OPERATOR operator) {
        this.operator = operator;
        switch (operator) {
            case SCREENED_COULOMB:
                double prefactor = 2.0 * beta / sqrtPI;
                double twoBeta2 = -2.0 * beta * beta;
                for (int n = 0; n <= order; n++) {
                    T000j[n] = prefactor * pow(twoBeta2, n);
                }
                break;
            default:
            case THOLE_FIELD:
            case COULOMB:
                for (int n = 0; n <= order; n++) {
                    /**
                     * Math.pow(-1.0, j) returns positive for all j, with -1.0
                     * as the // argument rather than -1. This is a bug?
                     */
                    T000j[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
                }
        }
    }

    /**
     * Set the Ewald beta parameters.
     *
     * @param beta
     */
    public void setEwaldBeta(double beta) {
        this.beta = beta;
    }

    /**
     * Set the Thole damping parameters.
     *
     * @param damp
     * @param aiak
     */
    public void setTholeDamping(double damp, double aiak) {
        this.damp = damp;
        this.aiak = aiak;
    }

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param r Cartesian vector.
     * @param T000 Location to store the source terms.
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    private void source(double r[], double T000[]) {
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        final double r2 = (x * x + y * y + z * z);

        switch (operator) {
            case SCREENED_COULOMB:
                // Sagui et al. Eq. 2.28.
                double R = sqrt(r2);
                double betaR = beta * R;
                double betaR2 = betaR * betaR;
                double iBetaR2 = 1.0 / (2.0 * betaR2);
                double expBR2 = exp(-betaR2);
                // Fnc(x^2) = Sqrt(PI) * erfc(x) / (2*x)
                // where x = Beta*R
                double Fnc = sqrtPI * erfc(betaR) / (2.0 * betaR);
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * Fnc;
                    // Generate F(n+1)c from Fnc (Eq. 2.24 in Sagui et al.)
                    // F(n+1)c = [(2*n+1) Fnc(x) + exp(-x)] / 2x
                    // where x = (Beta*R)^2
                    Fnc = ((2.0 * n + 1.0) * Fnc + expBR2) * iBetaR2;
                }
                break;
            case THOLE_FIELD:
                assert (order <= 4);
                double ir2 = 1.0 / r2;
                double ir = sqrt(ir2);
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * ir;
                    ir *= ir2;
                }
                /**
                 * Add the Thole damping terms: edamp = exp(-damp*u^3).
                 */
                double u = sqrt(r2) / aiak;
                double u3 = damp * pow(u, 3);
                double u6 = u3 * u3;
                double u9 = u6 * u3;
                double expU3 = exp(-u3);
                T000[0] = 0.0; // The zeroth order term is not calculated for Thole damping.
                T000[1] *= expU3;
                T000[2] *= (1.0 + u3) * expU3;
                T000[3] *= (1.0 + u3 + (3.0 / 5.0) * u6) * expU3;
                T000[4] *= (1.0 + u3 + (18.0 * u6 + 9.0 * u9) / 35.0) * expU3;
                break;
            case COULOMB:
            default:
                // Challacombe et al. Equation 40.
                ir2 = 1.0 / r2;
                ir = sqrt(ir2);
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * ir;
                    ir *= ir2;
                }
        }
    }

    /**
     * Log the tensors.
     *
     * @param tensor
     */
    public void log(double tensor[]) {

        StringBuilder sb = new StringBuilder();

        sb.append(String.format("\n %s Operator to order %d:", operator, order));
        sb.append(String.format("\n%5s %4s %4s %4s %12s\n", "Index", "d/dx", "d/dy", "d/dz", "Tensor"));
        sb.append(String.format("%5d %4d %4d %4d %12.8f\n", 0, 0, 0, 0, tensor[0]));
        int count = 1;
        // Print (d/dx)^l for l = 1..order (m = 0, n = 0)
        for (int l = 1; l <= order; l++) {
            double value = tensor[ti(l, 0, 0)];
            if (value != 0.0) {
                sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, 0, 0), l, 0, 0, value));
                count++;
            }
        }
        // Print (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
        for (int l = 0; l <= o1; l++) {
            for (int m = 1; m <= order - l; m++) {
                double value = tensor[ti(l, m, 0)];
                if (value != 0.0) {
                    sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, m, 0), l, m, 0, value));
                    count++;
                }
            }
        }
        // Print (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l <= o1; l++) {
            for (int m = 0; m <= o1 - l; m++) {
                for (int n = 1; n <= order - l - m; n++) {
                    double value = tensor[ti(l, m, n)];
                    if (value != 0.0) {
                        sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, m, n), l, m, n, value));
                        count++;
                    }
                }
            }
        }
        sb.append(String.format("\n Total number of active tensors: %d\n", count));
        logger.log(Level.INFO, sb.toString());
    }

    /**
     * Returns the number of tensors for derivatives to the given order.
     *
     * @param order maximum number of derivatives.
     * @return the number of tensors.
     * @since 1.0
     */
    public static int tensorCount(int order) {
        long ret = VectorMath.binomial(order + 3, 3);
        assert (ret < Integer.MAX_VALUE);
        return (int) ret;
    }

    /**
     * The index is based on the idea of filling tetrahedron.
     * <p>
     * 1/r has an index of 0
     * <br>
     * derivatives of x are first; indeces from 1..o for d/dx..(d/dx)^o
     * <br>
     * derivatives of x and y are second; base triangle of size (o+1)(o+2)/2
     * <br>
     * derivatives of x, y and z are last; total size (o+1)*(o+2)*(o+3)/6
     * <br>
     * <p>
     * This function is useful to set up masking constants:
     * <br>
     * static int Tlmn = ti(l,m,n,order)
     * <br>
     * For example the (d/dy)^2 (1/R) storage location: <br> static int T020 =
     * ti(0,2,0,order)
     * <p>
     *
     * @param dx int The number of d/dx operations.
     * @param dy int The number of d/dy operations.
     * @param dz int The number of d/dz operations.
     * @param order int The maximum tensor order (0 .LE. dx + dy + dz .LE.
     * order).
     *
     * @return int in the range (0..binomial(order + 3, 3) - 1)
     */
    public static int ti(int dx, int dy, int dz, int order) {
        int size = (order + 1) * (order + 2) * (order + 3) / 6;
        /**
         * We only get to the top of the tetrahedron if dz = order, otherwise
         * subtract off the top, including the level of the requested tensor
         * index.
         */
        int top = order + 1 - dz;
        top = top * (top + 1) * (top + 2) / 6;
        int zindex = size - top;
        /**
         * Given the "dz level", dy can range from 0..order - dz) To get to the
         * row for a specific value of dy, dy*(order + 1) - dy*(dy-1)/2 indeces
         * are skipped. This is an operation that looks like the area of
         * rectangle, minus the area of an empty triangle.
         */
        int yindex = dy * (order - dz) - (dy - 1) * (dy - 2) / 2 + 1;
        /**
         * Given the dz level and dy row, dx can range from (0..order - dz - dy)
         * The dx index is just walking down the dy row for "dx" steps.
         */
        int ret = dx + yindex + zindex;
        return ret;
    }

    /**
     * The index is based on the idea of filling tetrahedron.
     * <p>
     * 1/r has an index of 0
     * <br>
     * derivatives of x are first; indeces from 1..o for d/dx..(d/dx)^o
     * <br>
     * derivatives of x and y are second; base triangle of size (o+1)(o+2)/2
     * <br>
     * derivatives of x, y and z are last; total size (o+1)*(o+2)*(o+3)/6
     * <br>
     * <p>
     * This function is useful to set up masking constants:
     * <br>
     * static int Tlmn = ti(l,m,n,order)
     * <br>
     * For example the (d/dy)^2 (1/R) storage location: <br> static int T020 =
     * ti(0,2,0,order)
     * <p>
     *
     * @param dx int The number of d/dx operations.
     * @param dy int The number of d/dy operations.
     * @param dz int The number of d/dz operations.
     *
     * @return int in the range (0..binomial(order + 3, 3) - 1)
     */
    public int ti(int dx, int dy, int dz) {
        return ti(dx, dy, dz, order);
    }

    /**
     * This method is a driver to collect elements of the Cartesian multipole
     * tensor given the recursion relationships implemented by the method
     * "Tlmnj", which can be called directly to get a single tensor element. It
     * does not store intermediate values of the recursion, causing it to scale
     * O(order^8). For order = 5, this approach is a factor of 10 slower than
     * recursion.
     *
     * @param r double[] vector between two sites.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void noStorageRecursion(double r[], double tensor[]) {
        source(r, T000);
        // 1/r
        tensor[0] = T000[0];
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        for (int l = 1; l <= order; l++) {
            tensor[ti(l, 0, 0)] = Tlmnj(l, 0, 0, 0, r, T000);
        }
        // Find (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
        for (int l = 0; l <= o1; l++) {
            for (int m = 1; m <= order - l; m++) {
                tensor[ti(l, m, 0)] = Tlmnj(l, m, 0, 0, r, T000);
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l <= o1; l++) {
            for (int m = 0; m <= o1 - l; m++) {
                for (int n = 1; n <= order - l - m; n++) {
                    tensor[ti(l, m, n)] = Tlmnj(l, m, n, 0, r, T000);
                }
            }
        }
    }

    /**
     * This method is a driver to collect elements of the Cartesian multipole
     * tensor given the recursion relationships implemented by the method
     * "Tlmnj", which can be called directly to get a single tensor element. It
     * does not store intermediate values of the recursion, causing it to scale
     * O(order^8). For order = 5, this approach is a factor of 10 slower than
     * recursion.
     *
     * @param r double[] vector between two sites. r[0] and r[1] must equal 0.0.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void noStorageRecursionQI(double r[], double tensor[]) {
        assert (r[0] == 0.0 && r[1] == 0.0);
        source(r, T000);
        // 1/r
        tensor[0] = T000[0];
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        for (int l = 1; l <= order; l++) {
            tensor[ti(l, 0, 0)] = TlmnjQI(l, 0, 0, 0, r, T000);
        }
        // Find (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
        for (int l = 0; l <= o1; l++) {
            for (int m = 1; m <= order - l; m++) {
                tensor[ti(l, m, 0)] = TlmnjQI(l, m, 0, 0, r, T000);
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l <= o1; l++) {
            for (int m = 0; m <= o1 - l; m++) {
                for (int n = 1; n <= order - l - m; n++) {
                    tensor[ti(l, m, n)] = TlmnjQI(l, m, n, 0, r, T000);
                }
            }
        }
    }

    /**
     * This routine implements the recurrence relations for computation of any
     * Cartesian multipole tensor in ~O(L^8) time, where L is the total order l
     * + m + n, given the auxiliary elements T0000.
     * <br>
     * It implements the recursion relationships in brute force fashion, without
     * saving intermediate values. This is useful for finding a single tensor,
     * rather than all binomial(L + 3, 3).
     * <br>
     * The specific recursion equations (41-43) and set of auxiliary tensor
     * elements from equation (40) can be found in Challacombe et al.
     *
     * @param l int The number of (d/dx) operations.
     * @param m int The number of (d/dy) operations.
     * @param n int The number of (d/dz) operations.
     * @param j int j = 0 is the Tlmn tensor, j .GT. 0 is an intermediate.
     * @param r double[] The {x,y,z} coordinates.
     * @param T000 double[] Initial auxiliary tensor elements from Eq. (40).
     *
     * @return double The requested Tensor element (intermediate if j .GT. 0).
     *
     * @since 1.0
     */
    public static double Tlmnj(final int l, final int m, final int n,
            final int j, final double[] r, final double[] T000) {
        if (m == 0 && n == 0) {
            if (l > 1) {
                return r[0] * Tlmnj(l - 1, 0, 0, j + 1, r, T000) + (l - 1) * Tlmnj(l - 2, 0, 0, j + 1, r, T000);
            } else if (l == 1) { // l == 1, d/dx is done.
                return r[0] * Tlmnj(0, 0, 0, j + 1, r, T000);
            } else { // l = m = n = 0. Recursion is done.
                return T000[j];
            }
        } else if (n == 0) { // m >= 1
            if (m > 1) {
                return r[1] * Tlmnj(l, m - 1, 0, j + 1, r, T000) + (m - 1) * Tlmnj(l, m - 2, 0, j + 1, r, T000);
            }
            return r[1] * Tlmnj(l, 0, 0, j + 1, r, T000);
        } else { // n >= 1
            if (n > 1) {
                return r[2] * Tlmnj(l, m, n - 1, j + 1, r, T000) + (n - 1) * Tlmnj(l, m, n - 2, j + 1, r, T000);
            }
            return r[2] * Tlmnj(l, m, 0, j + 1, r, T000);
        }
    }

    /**
     * This routine implements the recurrence relations for computation of any
     * Cartesian multipole tensor in ~O(L^8) time, where L is the total order l
     * + m + n, given the auxiliary elements T0000.
     * <br>
     * It implements the recursion relationships in brute force fashion, without
     * saving intermediate values. This is useful for finding a single tensor,
     * rather than all binomial(L + 3, 3).
     * <br>
     * The specific recursion equations (41-43) and set of auxiliary tensor
     * elements from equation (40) can be found in Challacombe et al.
     *
     * @param l int The number of (d/dx) operations.
     * @param m int The number of (d/dy) operations.
     * @param n int The number of (d/dz) operations.
     * @param j int j = 0 is the Tlmn tensor, j .GT. 0 is an intermediate.
     * @param r double[] The {x,y,z} coordinates.
     * @param T000 double[] Initial auxiliary tensor elements from Eq. (40).
     *
     * @return double The requested Tensor element (intermediate if j .GT. 0).
     *
     * @since 1.0
     */
    public static double TlmnjQI(final int l, final int m, final int n,
            final int j, final double[] r, final double[] T000) {

        double z = r[2];
        assert (r[0] == 0.0 && r[1] == 0.0);

        if (m == 0 && n == 0) {
            if (l > 1) {
                return (l - 1) * Tlmnj(l - 2, 0, 0, j + 1, r, T000);
            } else if (l == 1) { // l == 1, d/dx is done.
                return 0.0;
            } else { // l = m = n = 0. Recursion is done.
                return T000[j];
            }
        } else if (n == 0) { // m >= 1
            if (m > 1) {
                return (m - 1) * Tlmnj(l, m - 2, 0, j + 1, r, T000);
            }
            return 0.0;
        } else { // n >= 1
            if (n > 1) {
                return z * Tlmnj(l, m, n - 1, j + 1, r, T000) + (n - 1) * Tlmnj(l, m, n - 2, j + 1, r, T000);
            }
            return z * Tlmnj(l, m, 0, j + 1, r, T000);
        }
    }

    /**
     * This function is a driver to collect elements of the Cartesian multipole
     * tensor. Collecting all tensors scales slightly better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example, up
     * to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1
     * <br>
     * tensor[0] = 1/|r| <br>
     * tensor[1] = -x/|r|^3 <br>
     * tensor[2] = -y/|r|^3 <br>
     * tensor[3] = -z/|r|^3 <br>
     * <p>
     * @param r double[] vector between two sites.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     * @since 1.0
     */
    public void recursion(final double r[], final double tensor[]) {
        source(r, work);
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        tensor[0] = work[0];
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        // Any (d/dx) term can be formed as
        // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
        // All intermediate terms are indexed as l*il + m*im + n*in + j;
        double current;
        double previous = work[1];
        // Store the l=1 tensor T100 (d/dx)
        tensor[ti(1, 0, 0)] = x * previous;
        // Starting the loop at l=2 avoids an if statement.
        for (int l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=2 avoid an if statement.
            // T100(l-1) = x * T000(l)
            current = x * work[l];
            int iw = il + l - 1;
            work[iw] = current;
            for (int a = 1; a < l - 1; a++) {
                // T200(l-2) = x * T100(l-1) + (2 - 1) * T000(l-1)
                // T300(l-3) = x * T200(l-2) + (3 - 1) * T100(l-2)
                // ...
                // T(l-1)001 = x * T(l-2)002 + (l - 2) * T(l-3)002
                current = x * current + a * work[iw - il];
                iw += il - 1;
                work[iw] = current;
            }
            // Store the Tl00 tensor (d/dx)^l
            // Tl00 = x * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
            tensor[ti(l, 0, 0)] = x * current + (l - 1) * previous;
            previous = current;
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (int l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[ti(l, 1, 0)] = y * previous;
            for (int m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = y * work[iw];
                iw += im - 1;
                work[iw] = current;
                for (int a = 1; a < m - 1; a++) {
                    // Tl20(m-2) = Y * Tl10(m-1) + (2 - 1) * T100(m-1)
                    // Tl30(m-3) = Y * Tl20(m-2) + (3 - 1) * Tl10(m-2)
                    // ...
                    // Tl(m-1)01 = Y * Tl(m-2)02 + (m - 2) * T(m-3)02
                    current = y * current + a * work[iw - im];
                    iw += im - 1;
                    work[iw] = current;
                }
                // Store the tensor (d/dx)^l * (d/dy)^m
                // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
                tensor[ti(l, m, 0)] = y * current + (m - 1) * previous;
                previous = current;
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
        // Any (d/dz) term can be formed as:
        // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m + l < order; m++) {
                // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
                // Tlmn = z * Tlm01
                final int lm = m + l;
                final int lilmim = l * il + m * im;
                previous = work[lilmim + 1];
                tensor[ti(l, m, 1)] = z * previous;
                for (int n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    final int n1 = n - 1;
                    for (int a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
                    tensor[ti(l, m, n)] = z * current + n1 * previous;
                    previous = current;
                }
            }
        }
    }

    /**
     * This function is a driver to collect elements of the Cartesian multipole
     * tensor for Quasi-Internal coordinate evaluation.
     *
     * Thus, we assume r[0] and r[1] (i.e. dx and dy) are zero.
     *
     * Collecting all tensors scales better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example, up
     * to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1<br> tensor[0] = 1/|r|
     * <br>
     * tensor[1] = -x/|r|^3 <br>
     * tensor[2] = -y/|r|^3 <br>
     * tensor[3] = -z/|r|^3 <br>
     * <p>
     *
     * @param r double[] vector between two sites (assumes r[0] and r[1] = 0.0).
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     * @since 1.0
     */
    public void recursionQI(final double r[], final double tensor[]) {
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        assert (x == 0.0 && y == 0.0);
        source(r, work);
        tensor[0] = work[0];
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        // Any (d/dx) term can be formed as
        // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
        // All intermediate terms are indexed as l*il + m*im + n*in + j;
        // Store the l=1 tensor T100 (d/dx)
        // Starting the loop at l=2 avoids an if statement.
        double current;
        double previous = work[1];
        tensor[ti(1, 0, 0)] = 0.0;
        for (int l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=2 avoid an if statement.
            // T100(l-1) = 0.0 * T000(l)
            current = 0.0;
            int iw = il + l - 1;
            work[iw] = current;
            for (int a = 1; a < l - 1; a++) {
                // T200(l-2) = 0.0 * T100(l-1) + (2 - 1) * T000(l-1)
                // T300(l-3) = 0.0 * T200(l-2) + (3 - 1) * T100(l-2)
                // ...
                // T(l-1)001 = 0.0 * T(l-2)002 + (l - 2) * T(l-3)002
                current = a * work[iw - il];
                iw += il - 1;
                work[iw] = current;
            }
            // Store the Tl00 tensor (d/dx)^l
            // Tl00 = 0.0 * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
            tensor[ti(l, 0, 0)] = (l - 1) * previous;
            previous = current;
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (int l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[ti(l, 1, 0)] = 0.0;
            for (int m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = 0.0;
                iw += im - 1;
                work[iw] = current;
                for (int a = 1; a < m - 1; a++) {
                    // Tl20(m-2) = 0.0 * Tl10(m-1) + (2 - 1) * T100(m-1)
                    // Tl30(m-3) = 0.0 * Tl20(m-2) + (3 - 1) * Tl10(m-2)
                    // ...
                    // Tl(m-1)01 = 0.0 * Tl(m-2)02 + (m - 2) * Tl(m-3)02
                    current = a * work[iw - im];
                    iw += im - 1;
                    work[iw] = current;
                }
                // Store the tensor (d/dx)^l * (d/dy)^m
                // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
                tensor[ti(l, m, 0)] = (m - 1) * previous;
                previous = current;
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
        // Any (d/dz) term can be formed as:
        // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m + l < order; m++) {
                // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
                // Tlmn = z * Tlm01
                final int lm = m + l;
                final int lilmim = l * il + m * im;
                previous = work[lilmim + 1];
                tensor[ti(l, m, 1)] = z * previous;
                for (int n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    final int n1 = n - 1;
                    for (int a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
                    tensor[ti(l, m, n)] = z * current + n1 * previous;
                    previous = current;
                }
            }
        }
    }

    /**
     * This function is a driver to collect elements of the Cartesian multipole
     * tensor. Collecting all tensors scales slightly better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example, up
     * to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1<br> tensor[0] = 1/|r|
     * <br>
     * tensor[1] = -x/|r|^3 <br> tensor[2] = -y/|r|^3 <br> tensor[3] = -z/|r|^3
     * <br>
     * <p>
     *
     * @param r double[] vector between two sites.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     * @return Java code for the tensor recursion.
     *
     * @since 1.0
     */
    public String codeTensorRecursion(final double r[], final double tensor[]) {
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        source(r, work);
        StringBuilder sb = new StringBuilder();
        tensor[0] = work[0];
        if (work[0] > 0) {
            sb.append(format("tensor[%s] = %s;\n", tlmn(0, 0, 0), term(0, 0, 0, 0)));
        }
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        // Any (d/dx) term can be formed as
        // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
        // All intermediate terms are indexed as l*il + m*im + n*in + j;
        double current;
        double previous = work[1];
        // Store the l=1 tensor T100 (d/dx)
        tensor[ti(1, 0, 0)] = x * previous;
        sb.append(format("tensor[%s] = x * %s;\n", tlmn(1, 0, 0), term(0, 0, 0, 1)));
        // Starting the loop at l=2 avoids an if statement.
        for (int l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=2 avoid an if statement.
            // T100(l-1) = x * T000(l)
            current = x * work[l];
            int iw = il + l - 1;
            work[iw] = current;
            sb.append(format("double %s = x * %s;\n", term(1, 0, 0, l - 1), term(0, 0, 0, l)));
            for (int a = 1; a < l - 1; a++) {
                // T200(l-2) = x * T100(l-1) + (2 - 1) * T000(l-1)
                // T300(l-3) = x * T200(l-2) + (3 - 1) * T100(l-2)
                // ...
                // T(l-1)001 = x * T(l-2)002 + (l - 2) * T(l-3)002
                current = x * current + a * work[iw - il];
                iw += il - 1;
                work[iw] = current;
                if (a > 1) {
                    sb.append(format("double %s = x * %s + %d * %s;\n",
                            term(a + 1, 0, 0, l - a - 1), term(a, 0, 0, l - a), a, term(a - 1, 0, 0, l - a)));
                } else {
                    sb.append(format("double %s = x * %s + %s;\n",
                            term(a + 1, 0, 0, l - a - 1), term(a, 0, 0, l - a), term(a - 1, 0, 0, l - a)));
                }
            }
            // Store the Tl00 tensor (d/dx)^l
            // Tl00 = x * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
            tensor[ti(l, 0, 0)] = x * current + (l - 1) * previous;
            previous = current;
            if (l > 2) {
                sb.append(format("tensor[%s] = x * %s + %d * %s;\n",
                        tlmn(l, 0, 0), term(l - 1, 0, 0, 1), (l - 1), term(l - 2, 0, 0, 1)));
            } else {
                sb.append(format("tensor[%s] = x * %s + %s;\n",
                        tlmn(l, 0, 0), term(l - 1, 0, 0, 1), term(l - 2, 0, 0, 1), l, 0, 0));
            }
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (int l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[ti(l, 1, 0)] = y * previous;
            sb.append(format("tensor[%s] = y * %s;\n", tlmn(l, 1, 0), term(l, 0, 0, 1)));
            for (int m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = y * work[iw];
                iw += im - 1;
                work[iw] = current;
                sb.append(format("double %s = y * %s;\n",
                        term(l, 1, 0, m - 1), term(l, 0, 0, m)));
                for (int a = 1; a < m - 1; a++) {
                    // Tl20(m-2) = Y * Tl10(m-1) + (2 - 1) * T100(m-1)
                    // Tl30(m-3) = Y * Tl20(m-2) + (3 - 1) * Tl10(m-2)
                    // ...
                    // Tl(m-1)01 = Y * Tl(m-2)02 + (m - 2) * T(m-3)02
                    current = y * current + a * work[iw - im];
                    iw += im - 1;
                    work[iw] = current;
                    if (a > 1) {
                        sb.append(format("double %s = y * %s + %d * %s;\n",
                                term(l, a + 1, 0, m - a - 1), term(l, a, 0, m - a), a, term(l, a - 1, 0, m - a)));
                    } else {
                        sb.append(format("double %s = y * %s + %s;\n",
                                term(l, a + 1, 0, m - a - 1), term(l, a, 0, m - a), term(l, a - 1, 0, m - a)));
                    }
                }
                // Store the tensor (d/dx)^l * (d/dy)^m
                // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
                tensor[ti(l, m, 0)] = y * current + (m - 1) * previous;
                previous = current;
                if (m > 2) {
                    sb.append(format("tensor[%s] = y * %s + %d * %s;\n",
                            tlmn(l, m, 0), term(l, m - 1, 0, 1), (m - 1), term(l, m - 2, 0, 1)));
                } else {
                    sb.append(format("tensor[%s] = y * %s + %s;\n",
                            tlmn(l, m, 0), term(l, m - 1, 0, 1), term(l, m - 2, 0, 1)));
                }

            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
        // Any (d/dz) term can be formed as:
        // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m + l < order; m++) {
                // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
                // Tlmn = z * Tlm01
                final int lm = m + l;
                final int lilmim = l * il + m * im;
                previous = work[lilmim + 1];
                tensor[ti(l, m, 1)] = z * previous;
                sb.append(format("tensor[%s] = z * %s;\n", tlmn(l, m, 1), term(l, m, 0, 1), l, m, 1));
                for (int n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    sb.append(format("double %s = z * %s;\n", term(l, m, 1, n - 1), term(l, m, 0, n)));
                    final int n1 = n - 1;
                    for (int a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                        if (a > 1) {
                            sb.append(format("double %s = z * %s + %d * %s;\n",
                                    term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), a, term(l, m, a - 1, n - a)));
                        } else {
                            sb.append(format("double %s = z * %s + %s;\n",
                                    term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), term(l, m, a - 1, n - a)));
                        }
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
                    tensor[ti(l, m, n)] = z * current + n1 * previous;
                    previous = current;
                    if (n > 2) {
                        sb.append(format("tensor[%s] = z * %s + %d * %s;\n",
                                tlmn(l, m, n), term(l, m, n - 1, 1), (n - 1), term(l, m, n - 2, 1), l, m, n));
                    } else {
                        sb.append(format("tensor[%s] = z * %s + %s;\n",
                                tlmn(l, m, n), term(l, m, n - 1, 1), term(l, m, n - 2, 1), l, m, n));
                    }
                }
            }
        }
        return sb.toString();
    }

    /**
     * This function write Java code to collect elements of the Cartesian
     * multipole tensor for Quasi-Internal coordinate evaluation.
     *
     * The r[0] and r[1] (i.e. dx and dy) must be zero.
     *
     * Collecting all tensors scales better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example, up
     * to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1<br> tensor[0] = 1/|r|
     * <br>
     * tensor[1] = -x/|r|^3 <br>
     * tensor[2] = -y/|r|^3 <br>
     * tensor[3] = -z/|r|^3 <br>
     * <p>
     *
     * @param r double[] vector between two sites (assumes r[0] and r[1] = 0.0).
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @return Java code to build the tensors.
     *
     * @since 1.0
     */
    public String codeTensorRecursionQI(final double r[], final double tensor[]) {
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        assert (x == 0.0 && y == 0.0);
        source(r, work);
        StringBuilder sb = new StringBuilder();
        tensor[0] = work[0];
        if (work[0] > 0) {
            sb.append(format("tensor[%s] = %s;\n", tlmn(0, 0, 0), term(0, 0, 0, 0)));
        }
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        // Any (d/dx) term can be formed as
        // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
        // All intermediate terms are indexed as l*il + m*im + n*in + j;
        // Store the l=1 tensor T100 (d/dx)
        // Starting the loop at l=2 avoids an if statement.
        double current;
        double previous = work[1];
        tensor[ti(1, 0, 0)] = 0.0;
        for (int l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=1 avoids an if statement.
            // T100(l-1) = 0.0 * T000(l)
            current = 0.0;
            int iw = il + (l - 1);
            work[iw] = current;
            // sb.append(format("double %s = 0.0;\n", term(1, 0, 0, l - 1)));
            for (int a = 1; a < l - 1; a++) {
                // T200(l-2) = 0.0 * T100(l-1) + (2 - 1) * T000(l-1)
                // T300(l-3) = 0.0 * T200(l-2) + (3 - 1) * T100(l-2)
                // ...
                // T(l-1)001 = 0.0 * T(l-2)002 + (l - 2) * T(l-3)002
                // iw = (a - 1) * il + (l - a)
                current = a * work[iw - il];
                iw += il - 1;
                // iw = (a + 1) * il + (l - a - 1)
                work[iw] = current;
                if (current != 0) {
                    if (a > 2) {
                        sb.append(format("double %s = %d * %s;\n",
                                term(a + 1, 0, 0, l - a - 1), a, term(a - 1, 0, 0, l - a)));
                    } else {
                        sb.append(format("double %s = %s;\n",
                                term(a + 1, 0, 0, l - a - 1), term(a - 1, 0, 0, l - a)));
                    }
                }
            }
            // Store the Tl00 tensor (d/dx)^l
            // Tl00 = 0.0 * T(l-1)001 + (l - 1) * T(l-2)001
            int index = ti(l, 0, 0);
            tensor[index] = (l - 1) * previous;
            previous = current;
            if (tensor[index] != 0) {
                if (l > 2) {
                    sb.append(format("tensor[%s] = %d * %s;\n", tlmn(l, 0, 0), (l - 1), term(l - 2, 0, 0, 1)));
                } else {
                    sb.append(format("tensor[%s] = %s;\n", tlmn(l, 0, 0), term(l - 2, 0, 0, 1)));
                }
            }
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (int l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[ti(l, 1, 0)] = 0.0;
            for (int m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = 0.0;
                iw += im - 1;
                work[iw] = current;
                for (int a = 1; a < m - 1; a++) {
                    // Tl20(m-2) = 0.0 * Tl10(m-1) + (2 - 1) * Tl00(m-1)
                    // Tl30(m-3) = 0.0 * Tl20(m-2) + (3 - 1) * Tl10(m-2)
                    // ...
                    // Tl(m-1)01 = 0.0 * Tl(m-2)02 + (m - 2) * Tl(m-3)02
                    current = a * work[iw - im];
                    iw += im - 1;
                    work[iw] = current;
                    if (current != 0) {
                        if (a > 1) {
                            sb.append(format("double %s = %d * %s;\n",
                                    term(l, a + 1, 0, m - a - 1), a, term(l, a - 1, 0, m - a)));
                        } else {
                            sb.append(format("double %s = %s;\n",
                                    term(l, a + 1, 0, m - a - 1), term(l, a - 1, 0, m - a)));
                        }

                    }
                }
                // Store the tensor (d/dx)^l * (d/dy)^m
                // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
                int index = ti(l, m, 0);
                tensor[index] = (m - 1) * previous;
                previous = current;
                if (tensor[index] != 0) {
                    if (m > 2) {
                        sb.append(format("tensor[%s] = %d * %s;\n", tlmn(l, m, 0), (m - 1), term(l, m - 2, 0, 1)));
                    } else {
                        sb.append(format("tensor[%s] = %s;\n", tlmn(l, m, 0), term(l, m - 2, 0, 1)));
                    }
                }

            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
        // Any (d/dz) term can be formed as:
        // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m + l < order; m++) {
                // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
                // Tlmn = z * Tlm01
                final int lm = m + l;
                final int lilmim = l * il + m * im;
                previous = work[lilmim + 1];
                int index = ti(l, m, 1);
                tensor[index] = z * previous;
                if (tensor[index] != 0) {
                    sb.append(format("tensor[%s] = z * %s;\n", tlmn(l, m, 1), term(l, m, 0, 1)));
                }
                for (int n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    if (current != 0) {
                        sb.append(format("double %s = z * %s;\n",
                                term(l, m, 1, n - 1), term(l, m, 0, n)));
                    }
                    final int n1 = n - 1;
                    for (int a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                        if (current != 0) {
                            if (a > 1) {
                                sb.append(format("double %s = z * %s + %d * %s;\n",
                                        term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), a, term(l, m, a - 1, n - a)));
                            } else {
                                sb.append(format("double %s = z * %s + %s;\n",
                                        term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), term(l, m, a - 1, n - a)));
                            }
                        }
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
                    index = ti(l, m, n);
                    tensor[index] = z * current + n1 * previous;
                    previous = current;
                    if (tensor[index] != 0) {
                        if (n > 2) {
                            sb.append(format("tensor[%s] = z * %s + %d * %s;\n",
                                    tlmn(l, m, n), term(l, m, n - 1, 1), (n - 1), term(l, m, n - 2, 1)));
                        } else {
                            sb.append(format("tensor[%s] = z * %s + %s;\n",
                                    tlmn(l, m, n), term(l, m, n - 1, 1), term(l, m, n - 2, 1)));
                        }
                    }
                }
            }
        }
        return sb.toString();
    }

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 4th
     * order, in the global frame, which is sufficient for quadrupole-induced
     * dipole forces.
     *
     * @param r Separation between particles i and k
     * @param tensor the Cartesian multipole tensors in the global frame.
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void order4(final double r[]) {
        source(r, work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        R000 = term0000;
        R100 = x * term0001;
        double term1001 = x * term0002;
        R200 = x * term1001 + term0001;
        double term1002 = x * term0003;
        double term2001 = x * term1002 + term0002;
        R300 = x * term2001 + 2 * term1001;
        double term1003 = x * term0004;
        double term2002 = x * term1003 + term0003;
        double term3001 = x * term2002 + 2 * term1002;
        R400 = x * term3001 + 3 * term2001;
        R010 = y * term0001;
        double term0101 = y * term0002;
        R020 = y * term0101 + term0001;
        double term0102 = y * term0003;
        double term0201 = y * term0102 + term0002;
        R030 = y * term0201 + 2 * term0101;
        double term0103 = y * term0004;
        double term0202 = y * term0103 + term0003;
        double term0301 = y * term0202 + 2 * term0102;
        R040 = y * term0301 + 3 * term0201;
        R110 = y * term1001;
        double term1101 = y * term1002;
        R120 = y * term1101 + term1001;
        double term1102 = y * term1003;
        double term1201 = y * term1102 + term1002;
        R130 = y * term1201 + 2 * term1101;
        R210 = y * term2001;
        double term2101 = y * term2002;
        R220 = y * term2101 + term2001;
        R310 = y * term3001;
        R001 = z * term0001;
        double term0011 = z * term0002;
        R002 = z * term0011 + term0001;
        double term0012 = z * term0003;
        double term0021 = z * term0012 + term0002;
        R003 = z * term0021 + 2 * term0011;
        double term0013 = z * term0004;
        double term0022 = z * term0013 + term0003;
        double term0031 = z * term0022 + 2 * term0012;
        R004 = z * term0031 + 3 * term0021;
        R011 = z * term0101;
        double term0111 = z * term0102;
        R012 = z * term0111 + term0101;
        double term0112 = z * term0103;
        double term0121 = z * term0112 + term0102;
        R013 = z * term0121 + 2 * term0111;
        R021 = z * term0201;
        double term0211 = z * term0202;
        R022 = z * term0211 + term0201;
        R031 = z * term0301;
        R101 = z * term1001;
        double term1011 = z * term1002;
        R102 = z * term1011 + term1001;
        double term1012 = z * term1003;
        double term1021 = z * term1012 + term1002;
        R103 = z * term1021 + 2 * term1011;
        R111 = z * term1101;
        double term1111 = z * term1102;
        R112 = z * term1111 + term1101;
        R121 = z * term1201;
        R201 = z * term2001;
        double term2011 = z * term2002;
        R202 = z * term2011 + term2001;
        R211 = z * term2101;
        R301 = z * term3001;
    }

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 4th
     * order, based on a quasi-internal frame, which is sufficient for
     * quadrupole-induced dipole forces.
     *
     * @param r Separation between particles i and k.
     * @param tensor the Cartesian multipole tensors in the quasi-internal
     * frame.
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void order4QI(final double r[]) {
        source(r, work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        R000 = term0000;
        R200 = term0001;
        double term2001 = term0002;
        double term2002 = term0003;
        R400 = 3 * term2001;
        R020 = term0001;
        double term0201 = term0002;
        double term0202 = term0003;
        R040 = 3 * term0201;
        R220 = term2001;
        R001 = z * term0001;
        double term0011 = z * term0002;
        R002 = z * term0011 + term0001;
        double term0012 = z * term0003;
        double term0021 = z * term0012 + term0002;
        R003 = z * term0021 + 2 * term0011;
        double term0013 = z * term0004;
        double term0022 = z * term0013 + term0003;
        double term0031 = z * term0022 + 2 * term0012;
        R004 = z * term0031 + 3 * term0021;
        R021 = z * term0201;
        double term0211 = z * term0202;
        R022 = z * term0211 + term0201;
        R201 = z * term2001;
        double term2011 = z * term2002;
        R202 = z * term2011 + term2001;
    }

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 5th
     * order, in the global frame, which is sufficient for quadrupole-quadrupole
     * forces.
     *
     * @param r Separation between particles i and k
     * @param tensor the Cartesian multipole tensors in the global frame.
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void order5(final double r[]) {
        source(r, work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
        double term0005 = work[5];
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        R000 = term0000;
        R100 = x * term0001;
        double term1001 = x * term0002;
        R200 = x * term1001 + term0001;
        double term1002 = x * term0003;
        double term2001 = x * term1002 + term0002;
        R300 = x * term2001 + 2 * term1001;
        double term1003 = x * term0004;
        double term2002 = x * term1003 + term0003;
        double term3001 = x * term2002 + 2 * term1002;
        R400 = x * term3001 + 3 * term2001;
        double term1004 = x * term0005;
        double term2003 = x * term1004 + term0004;
        double term3002 = x * term2003 + 2 * term1003;
        double term4001 = x * term3002 + 3 * term2002;
        R500 = x * term4001 + 4 * term3001;
        R010 = y * term0001;
        double term0101 = y * term0002;
        R020 = y * term0101 + term0001;
        double term0102 = y * term0003;
        double term0201 = y * term0102 + term0002;
        R030 = y * term0201 + 2 * term0101;
        double term0103 = y * term0004;
        double term0202 = y * term0103 + term0003;
        double term0301 = y * term0202 + 2 * term0102;
        R040 = y * term0301 + 3 * term0201;
        double term0104 = y * term0005;
        double term0203 = y * term0104 + term0004;
        double term0302 = y * term0203 + 2 * term0103;
        double term0401 = y * term0302 + 3 * term0202;
        R050 = y * term0401 + 4 * term0301;
        R110 = y * term1001;
        double term1101 = y * term1002;
        R120 = y * term1101 + term1001;
        double term1102 = y * term1003;
        double term1201 = y * term1102 + term1002;
        R130 = y * term1201 + 2 * term1101;
        double term1103 = y * term1004;
        double term1202 = y * term1103 + term1003;
        double term1301 = y * term1202 + 2 * term1102;
        R140 = y * term1301 + 3 * term1201;
        R210 = y * term2001;
        double term2101 = y * term2002;
        R220 = y * term2101 + term2001;
        double term2102 = y * term2003;
        double term2201 = y * term2102 + term2002;
        R230 = y * term2201 + 2 * term2101;
        R310 = y * term3001;
        double term3101 = y * term3002;
        R320 = y * term3101 + term3001;
        R410 = y * term4001;
        R001 = z * term0001;
        double term0011 = z * term0002;
        R002 = z * term0011 + term0001;
        double term0012 = z * term0003;
        double term0021 = z * term0012 + term0002;
        R003 = z * term0021 + 2 * term0011;
        double term0013 = z * term0004;
        double term0022 = z * term0013 + term0003;
        double term0031 = z * term0022 + 2 * term0012;
        R004 = z * term0031 + 3 * term0021;
        double term0014 = z * term0005;
        double term0023 = z * term0014 + term0004;
        double term0032 = z * term0023 + 2 * term0013;
        double term0041 = z * term0032 + 3 * term0022;
        R005 = z * term0041 + 4 * term0031;
        R011 = z * term0101;
        double term0111 = z * term0102;
        R012 = z * term0111 + term0101;
        double term0112 = z * term0103;
        double term0121 = z * term0112 + term0102;
        R013 = z * term0121 + 2 * term0111;
        double term0113 = z * term0104;
        double term0122 = z * term0113 + term0103;
        double term0131 = z * term0122 + 2 * term0112;
        R014 = z * term0131 + 3 * term0121;
        R021 = z * term0201;
        double term0211 = z * term0202;
        R022 = z * term0211 + term0201;
        double term0212 = z * term0203;
        double term0221 = z * term0212 + term0202;
        R023 = z * term0221 + 2 * term0211;
        R031 = z * term0301;
        double term0311 = z * term0302;
        R032 = z * term0311 + term0301;
        R041 = z * term0401;
        R101 = z * term1001;
        double term1011 = z * term1002;
        R102 = z * term1011 + term1001;
        double term1012 = z * term1003;
        double term1021 = z * term1012 + term1002;
        R103 = z * term1021 + 2 * term1011;
        double term1013 = z * term1004;
        double term1022 = z * term1013 + term1003;
        double term1031 = z * term1022 + 2 * term1012;
        R104 = z * term1031 + 3 * term1021;
        R111 = z * term1101;
        double term1111 = z * term1102;
        R112 = z * term1111 + term1101;
        double term1112 = z * term1103;
        double term1121 = z * term1112 + term1102;
        R113 = z * term1121 + 2 * term1111;
        R121 = z * term1201;
        double term1211 = z * term1202;
        R122 = z * term1211 + term1201;
        R131 = z * term1301;
        R201 = z * term2001;
        double term2011 = z * term2002;
        R202 = z * term2011 + term2001;
        double term2012 = z * term2003;
        double term2021 = z * term2012 + term2002;
        R203 = z * term2021 + 2 * term2011;
        R211 = z * term2101;
        double term2111 = z * term2102;
        R212 = z * term2111 + term2101;
        R221 = z * term2201;
        R301 = z * term3001;
        double term3011 = z * term3002;
        R302 = z * term3011 + term3001;
        R311 = z * term3101;
        R401 = z * term4001;
    }

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 5th
     * order, based on a quasi-internal frame, which is sufficient for
     * quadrupole-quadrupole forces.
     *
     * @param r Separation between particles i and k.
     * @param tensor the Cartesian multipole tensors in the quasi-internal
     * frame.
     * @param damp Thole damping for this interaction.
     * @param aiak Polarizability of (site i * site k)^6.
     */
    public void order5QI(final double r[]) {
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        assert (x == 0.0 && y == 0.0);
        source(r, work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
        double term0005 = work[5];
        R000 = term0000;
        R200 = term0001;
        double term2001 = term0002;
        double term2002 = term0003;
        R400 = 3 * term2001;
        double term2003 = term0004;
        double term4001 = 3 * term2002;
        R020 = term0001;
        double term0201 = term0002;
        double term0202 = term0003;
        R040 = 3 * term0201;
        double term0203 = term0004;
        double term0401 = 3 * term0202;
        R220 = term2001;
        double term2201 = term2002;
        R001 = z * term0001;
        double term0011 = z * term0002;
        R002 = z * term0011 + term0001;
        double term0012 = z * term0003;
        double term0021 = z * term0012 + term0002;
        R003 = z * term0021 + 2 * term0011;
        double term0013 = z * term0004;
        double term0022 = z * term0013 + term0003;
        double term0031 = z * term0022 + 2 * term0012;
        R004 = z * term0031 + 3 * term0021;
        double term0014 = z * term0005;
        double term0023 = z * term0014 + term0004;
        double term0032 = z * term0023 + 2 * term0013;
        double term0041 = z * term0032 + 3 * term0022;
        R005 = z * term0041 + 4 * term0031;
        R021 = z * term0201;
        double term0211 = z * term0202;
        R022 = z * term0211 + term0201;
        double term0212 = z * term0203;
        double term0221 = z * term0212 + term0202;
        R023 = z * term0221 + 2 * term0211;
        R041 = z * term0401;
        R201 = z * term2001;
        double term2011 = z * term2002;
        R202 = z * term2011 + term2001;
        double term2012 = z * term2003;
        double term2021 = z * term2012 + term2002;
        R203 = z * term2021 + 2 * term2011;
        R221 = z * term2201;
        R401 = z * term4001;
    }

    public void Ek5() {
        double term000 = 0.0;
        term000 += qi * R000;
        term000 += dxi * R100;
        term000 += dyi * R010;
        term000 += dzi * R001;
        term000 += qxxi * R200;
        term000 += qyyi * R020;
        term000 += qzzi * R002;
        double term0002 = 0.0;
        term0002 += qxyi * R110;
        term0002 += qxzi * R101;
        term0002 += qyzi * R011;
        term000 += 2.0 * term0002;
        E000 = term000;
        double term100 = 0.0;
        term100 += qi * R100;
        term100 += dxi * R200;
        term100 += dyi * R110;
        term100 += dzi * R101;
        term100 += qxxi * R300;
        term100 += qyyi * R120;
        term100 += qzzi * R102;
        double term1002 = 0.0;
        term1002 += qxyi * R210;
        term1002 += qxzi * R201;
        term1002 += qyzi * R111;
        term100 += 2.0 * term1002;
        E100 = term100;
        double term010 = 0.0;
        term010 += qi * R010;
        term010 += dxi * R110;
        term010 += dyi * R020;
        term010 += dzi * R011;
        term010 += qxxi * R210;
        term010 += qyyi * R030;
        term010 += qzzi * R012;
        double term0102 = 0.0;
        term0102 += qxyi * R120;
        term0102 += qxzi * R111;
        term0102 += qyzi * R021;
        term010 += 2.0 * term0102;
        E010 = term010;
        double term001 = 0.0;
        term001 += qi * R001;
        term001 += dxi * R101;
        term001 += dyi * R011;
        term001 += dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        double term0012 = 0.0;
        term0012 += qxyi * R111;
        term0012 += qxzi * R102;
        term0012 += qyzi * R012;
        term001 += 2.0 * term0012;
        E001 = term001;
        double term200 = 0.0;
        term200 += qi * R200;
        term200 += dxi * R300;
        term200 += dyi * R210;
        term200 += dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        double term2002 = 0.0;
        term2002 += qxyi * R310;
        term2002 += qxzi * R301;
        term2002 += qyzi * R211;
        term200 += 2.0 * term2002;
        E200 = term200;
        double term020 = 0.0;
        term020 += qi * R020;
        term020 += dxi * R120;
        term020 += dyi * R030;
        term020 += dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        double term0202 = 0.0;
        term0202 += qxyi * R130;
        term0202 += qxzi * R121;
        term0202 += qyzi * R031;
        term020 += 2.0 * term0202;
        E020 = term020;
        double term002 = 0.0;
        term002 += qi * R002;
        term002 += dxi * R102;
        term002 += dyi * R012;
        term002 += dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        double term0022 = 0.0;
        term0022 += qxyi * R112;
        term0022 += qxzi * R103;
        term0022 += qyzi * R013;
        term002 += 2.0 * term0022;
        E002 = term002;
        double term110 = 0.0;
        term110 += qi * R110;
        term110 += dxi * R210;
        term110 += dyi * R120;
        term110 += dzi * R111;
        term110 += qxxi * R310;
        term110 += qyyi * R130;
        term110 += qzzi * R112;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term1102 += qxzi * R211;
        term1102 += qyzi * R121;
        term110 += 2.0 * term1102;
        E110 = term110;
        double term101 = 0.0;
        term101 += qi * R101;
        term101 += dxi * R201;
        term101 += dyi * R111;
        term101 += dzi * R102;
        term101 += qxxi * R301;
        term101 += qyyi * R121;
        term101 += qzzi * R103;
        double term1012 = 0.0;
        term1012 += qxyi * R211;
        term1012 += qxzi * R202;
        term1012 += qyzi * R112;
        term101 += 2.0 * term1012;
        E101 = term101;
        double term011 = 0.0;
        term011 += qi * R011;
        term011 += dxi * R111;
        term011 += dyi * R021;
        term011 += dzi * R012;
        term011 += qxxi * R211;
        term011 += qyyi * R031;
        term011 += qzzi * R013;
        double term0112 = 0.0;
        term0112 += qxyi * R121;
        term0112 += qxzi * R112;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E011 = term011;
    }

    public void Ei5() {
        double term000 = 0.0;
        term000 += qk * R000;
        term000 += dxk * R100;
        term000 += dyk * R010;
        term000 += dzk * R001;
        term000 += qxxk * R200;
        term000 += qyyk * R020;
        term000 += qzzk * R002;
        double term0002 = 0.0;
        term0002 += qxyk * R110;
        term0002 += qxzk * R101;
        term0002 += qyzk * R011;
        term000 += 2.0 * term0002;
        E000 = term000;
        double term100 = 0.0;
        term100 += qk * R100;
        term100 += dxk * R200;
        term100 += dyk * R110;
        term100 += dzk * R101;
        term100 += qxxk * R300;
        term100 += qyyk * R120;
        term100 += qzzk * R102;
        double term1002 = 0.0;
        term1002 += qxyk * R210;
        term1002 += qxzk * R201;
        term1002 += qyzk * R111;
        term100 += 2.0 * term1002;
        E100 = term100;
        double term010 = 0.0;
        term010 += qk * R010;
        term010 += dxk * R110;
        term010 += dyk * R020;
        term010 += dzk * R011;
        term010 += qxxk * R210;
        term010 += qyyk * R030;
        term010 += qzzk * R012;
        double term0102 = 0.0;
        term0102 += qxyk * R120;
        term0102 += qxzk * R111;
        term0102 += qyzk * R021;
        term010 += 2.0 * term0102;
        E010 = term010;
        double term001 = 0.0;
        term001 += qk * R001;
        term001 += dxk * R101;
        term001 += dyk * R011;
        term001 += dzk * R002;
        term001 += qxxk * R201;
        term001 += qyyk * R021;
        term001 += qzzk * R003;
        double term0012 = 0.0;
        term0012 += qxyk * R111;
        term0012 += qxzk * R102;
        term0012 += qyzk * R012;
        term001 += 2.0 * term0012;
        E001 = term001;
        double term200 = 0.0;
        term200 += qk * R200;
        term200 += dxk * R300;
        term200 += dyk * R210;
        term200 += dzk * R201;
        term200 += qxxk * R400;
        term200 += qyyk * R220;
        term200 += qzzk * R202;
        double term2002 = 0.0;
        term2002 += qxyk * R310;
        term2002 += qxzk * R301;
        term2002 += qyzk * R211;
        term200 += 2.0 * term2002;
        E200 = term200;
        double term020 = 0.0;
        term020 += qk * R020;
        term020 += dxk * R120;
        term020 += dyk * R030;
        term020 += dzk * R021;
        term020 += qxxk * R220;
        term020 += qyyk * R040;
        term020 += qzzk * R022;
        double term0202 = 0.0;
        term0202 += qxyk * R130;
        term0202 += qxzk * R121;
        term0202 += qyzk * R031;
        term020 += 2.0 * term0202;
        E020 = term020;
        double term002 = 0.0;
        term002 += qk * R002;
        term002 += dxk * R102;
        term002 += dyk * R012;
        term002 += dzk * R003;
        term002 += qxxk * R202;
        term002 += qyyk * R022;
        term002 += qzzk * R004;
        double term0022 = 0.0;
        term0022 += qxyk * R112;
        term0022 += qxzk * R103;
        term0022 += qyzk * R013;
        term002 += 2.0 * term0022;
        E002 = term002;
        double term110 = 0.0;
        term110 += qk * R110;
        term110 += dxk * R210;
        term110 += dyk * R120;
        term110 += dzk * R111;
        term110 += qxxk * R310;
        term110 += qyyk * R130;
        term110 += qzzk * R112;
        double term1102 = 0.0;
        term1102 += qxyk * R220;
        term1102 += qxzk * R211;
        term1102 += qyzk * R121;
        term110 += 2.0 * term1102;
        E110 = term110;
        double term101 = 0.0;
        term101 += qk * R101;
        term101 += dxk * R201;
        term101 += dyk * R111;
        term101 += dzk * R102;
        term101 += qxxk * R301;
        term101 += qyyk * R121;
        term101 += qzzk * R103;
        double term1012 = 0.0;
        term1012 += qxyk * R211;
        term1012 += qxzk * R202;
        term1012 += qyzk * R112;
        term101 += 2.0 * term1012;
        E101 = term101;
        double term011 = 0.0;
        term011 += qk * R011;
        term011 += dxk * R111;
        term011 += dyk * R021;
        term011 += dzk * R012;
        term011 += qxxk * R211;
        term011 += qyyk * R031;
        term011 += qzzk * R013;
        double term0112 = 0.0;
        term0112 += qxyk * R121;
        term0112 += qxzk * R112;
        term0112 += qyzk * R022;
        term011 += 2.0 * term0112;
        E011 = term011;
    }

    public void Ex5() {
        double term100 = 0.0;
        term100 += qi * R100;
        term100 += dxi * R200;
        term100 += dyi * R110;
        term100 += dzi * R101;
        term100 += qxxi * R300;
        term100 += qyyi * R120;
        term100 += qzzi * R102;
        double term1002 = 0.0;
        term1002 += qxyi * R210;
        term1002 += qxzi * R201;
        term1002 += qyzi * R111;
        term100 += 2.0 * term1002;
        E000 = term100;
        double term200 = 0.0;
        term200 += qi * R200;
        term200 += dxi * R300;
        term200 += dyi * R210;
        term200 += dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        double term2002 = 0.0;
        term2002 += qxyi * R310;
        term2002 += qxzi * R301;
        term2002 += qyzi * R211;
        term200 += 2.0 * term2002;
        E100 = term200;
        double term110 = 0.0;
        term110 += qi * R110;
        term110 += dxi * R210;
        term110 += dyi * R120;
        term110 += dzi * R111;
        term110 += qxxi * R310;
        term110 += qyyi * R130;
        term110 += qzzi * R112;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term1102 += qxzi * R211;
        term1102 += qyzi * R121;
        term110 += 2.0 * term1102;
        E010 = term110;
        double term101 = 0.0;
        term101 += qi * R101;
        term101 += dxi * R201;
        term101 += dyi * R111;
        term101 += dzi * R102;
        term101 += qxxi * R301;
        term101 += qyyi * R121;
        term101 += qzzi * R103;
        double term1012 = 0.0;
        term1012 += qxyi * R211;
        term1012 += qxzi * R202;
        term1012 += qyzi * R112;
        term101 += 2.0 * term1012;
        E001 = term101;
        double term300 = 0.0;
        term300 += qi * R300;
        term300 += dxi * R400;
        term300 += dyi * R310;
        term300 += dzi * R301;
        term300 += qxxi * R500;
        term300 += qyyi * R320;
        term300 += qzzi * R302;
        double term3002 = 0.0;
        term3002 += qxyi * R410;
        term3002 += qxzi * R401;
        term3002 += qyzi * R311;
        term300 += 2.0 * term3002;
        E200 = term300;
        double term120 = 0.0;
        term120 += qi * R120;
        term120 += dxi * R220;
        term120 += dyi * R130;
        term120 += dzi * R121;
        term120 += qxxi * R320;
        term120 += qyyi * R140;
        term120 += qzzi * R122;
        double term1202 = 0.0;
        term1202 += qxyi * R230;
        term1202 += qxzi * R221;
        term1202 += qyzi * R131;
        term120 += 2.0 * term1202;
        E020 = term120;
        double term102 = 0.0;
        term102 += qi * R102;
        term102 += dxi * R202;
        term102 += dyi * R112;
        term102 += dzi * R103;
        term102 += qxxi * R302;
        term102 += qyyi * R122;
        term102 += qzzi * R104;
        double term1022 = 0.0;
        term1022 += qxyi * R212;
        term1022 += qxzi * R203;
        term1022 += qyzi * R113;
        term102 += 2.0 * term1022;
        E002 = term102;
        double term210 = 0.0;
        term210 += qi * R210;
        term210 += dxi * R310;
        term210 += dyi * R220;
        term210 += dzi * R211;
        term210 += qxxi * R410;
        term210 += qyyi * R230;
        term210 += qzzi * R212;
        double term2102 = 0.0;
        term2102 += qxyi * R320;
        term2102 += qxzi * R311;
        term2102 += qyzi * R221;
        term210 += 2.0 * term2102;
        E110 = term210;
        double term201 = 0.0;
        term201 += qi * R201;
        term201 += dxi * R301;
        term201 += dyi * R211;
        term201 += dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        double term2012 = 0.0;
        term2012 += qxyi * R311;
        term2012 += qxzi * R302;
        term2012 += qyzi * R212;
        term201 += 2.0 * term2012;
        E101 = term201;
        double term111 = 0.0;
        term111 += qi * R111;
        term111 += dxi * R211;
        term111 += dyi * R121;
        term111 += dzi * R112;
        term111 += qxxi * R311;
        term111 += qyyi * R131;
        term111 += qzzi * R113;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term1112 += qxzi * R212;
        term1112 += qyzi * R122;
        term111 += 2.0 * term1112;
        E011 = term111;
    }

    public void Ey5() {
        double term010 = 0.0;
        term010 += qi * R010;
        term010 += dxi * R110;
        term010 += dyi * R020;
        term010 += dzi * R011;
        term010 += qxxi * R210;
        term010 += qyyi * R030;
        term010 += qzzi * R012;
        double term0102 = 0.0;
        term0102 += qxyi * R120;
        term0102 += qxzi * R111;
        term0102 += qyzi * R021;
        term010 += 2.0 * term0102;
        E000 = term010;
        double term110 = 0.0;
        term110 += qi * R110;
        term110 += dxi * R210;
        term110 += dyi * R120;
        term110 += dzi * R111;
        term110 += qxxi * R310;
        term110 += qyyi * R130;
        term110 += qzzi * R112;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term1102 += qxzi * R211;
        term1102 += qyzi * R121;
        term110 += 2.0 * term1102;
        E100 = term110;
        double term020 = 0.0;
        term020 += qi * R020;
        term020 += dxi * R120;
        term020 += dyi * R030;
        term020 += dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        double term0202 = 0.0;
        term0202 += qxyi * R130;
        term0202 += qxzi * R121;
        term0202 += qyzi * R031;
        term020 += 2.0 * term0202;
        E010 = term020;
        double term011 = 0.0;
        term011 += qi * R011;
        term011 += dxi * R111;
        term011 += dyi * R021;
        term011 += dzi * R012;
        term011 += qxxi * R211;
        term011 += qyyi * R031;
        term011 += qzzi * R013;
        double term0112 = 0.0;
        term0112 += qxyi * R121;
        term0112 += qxzi * R112;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E001 = term011;
        double term210 = 0.0;
        term210 += qi * R210;
        term210 += dxi * R310;
        term210 += dyi * R220;
        term210 += dzi * R211;
        term210 += qxxi * R410;
        term210 += qyyi * R230;
        term210 += qzzi * R212;
        double term2102 = 0.0;
        term2102 += qxyi * R320;
        term2102 += qxzi * R311;
        term2102 += qyzi * R221;
        term210 += 2.0 * term2102;
        E200 = term210;
        double term030 = 0.0;
        term030 += qi * R030;
        term030 += dxi * R130;
        term030 += dyi * R040;
        term030 += dzi * R031;
        term030 += qxxi * R230;
        term030 += qyyi * R050;
        term030 += qzzi * R032;
        double term0302 = 0.0;
        term0302 += qxyi * R140;
        term0302 += qxzi * R131;
        term0302 += qyzi * R041;
        term030 += 2.0 * term0302;
        E020 = term030;
        double term012 = 0.0;
        term012 += qi * R012;
        term012 += dxi * R112;
        term012 += dyi * R022;
        term012 += dzi * R013;
        term012 += qxxi * R212;
        term012 += qyyi * R032;
        term012 += qzzi * R014;
        double term0122 = 0.0;
        term0122 += qxyi * R122;
        term0122 += qxzi * R113;
        term0122 += qyzi * R023;
        term012 += 2.0 * term0122;
        E002 = term012;
        double term120 = 0.0;
        term120 += qi * R120;
        term120 += dxi * R220;
        term120 += dyi * R130;
        term120 += dzi * R121;
        term120 += qxxi * R320;
        term120 += qyyi * R140;
        term120 += qzzi * R122;
        double term1202 = 0.0;
        term1202 += qxyi * R230;
        term1202 += qxzi * R221;
        term1202 += qyzi * R131;
        term120 += 2.0 * term1202;
        E110 = term120;
        double term111 = 0.0;
        term111 += qi * R111;
        term111 += dxi * R211;
        term111 += dyi * R121;
        term111 += dzi * R112;
        term111 += qxxi * R311;
        term111 += qyyi * R131;
        term111 += qzzi * R113;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term1112 += qxzi * R212;
        term1112 += qyzi * R122;
        term111 += 2.0 * term1112;
        E101 = term111;
        double term021 = 0.0;
        term021 += qi * R021;
        term021 += dxi * R121;
        term021 += dyi * R031;
        term021 += dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        double term0212 = 0.0;
        term0212 += qxyi * R131;
        term0212 += qxzi * R122;
        term0212 += qyzi * R032;
        term021 += 2.0 * term0212;
        E011 = term021;
    }

    public void Ez5() {
        double term001 = 0.0;
        term001 += qi * R001;
        term001 += dxi * R101;
        term001 += dyi * R011;
        term001 += dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        double term0012 = 0.0;
        term0012 += qxyi * R111;
        term0012 += qxzi * R102;
        term0012 += qyzi * R012;
        term001 += 2.0 * term0012;
        E000 = term001;
        double term101 = 0.0;
        term101 += qi * R101;
        term101 += dxi * R201;
        term101 += dyi * R111;
        term101 += dzi * R102;
        term101 += qxxi * R301;
        term101 += qyyi * R121;
        term101 += qzzi * R103;
        double term1012 = 0.0;
        term1012 += qxyi * R211;
        term1012 += qxzi * R202;
        term1012 += qyzi * R112;
        term101 += 2.0 * term1012;
        E100 = term101;
        double term011 = 0.0;
        term011 += qi * R011;
        term011 += dxi * R111;
        term011 += dyi * R021;
        term011 += dzi * R012;
        term011 += qxxi * R211;
        term011 += qyyi * R031;
        term011 += qzzi * R013;
        double term0112 = 0.0;
        term0112 += qxyi * R121;
        term0112 += qxzi * R112;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E010 = term011;
        double term002 = 0.0;
        term002 += qi * R002;
        term002 += dxi * R102;
        term002 += dyi * R012;
        term002 += dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        double term0022 = 0.0;
        term0022 += qxyi * R112;
        term0022 += qxzi * R103;
        term0022 += qyzi * R013;
        term002 += 2.0 * term0022;
        E001 = term002;
        double term201 = 0.0;
        term201 += qi * R201;
        term201 += dxi * R301;
        term201 += dyi * R211;
        term201 += dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        double term2012 = 0.0;
        term2012 += qxyi * R311;
        term2012 += qxzi * R302;
        term2012 += qyzi * R212;
        term201 += 2.0 * term2012;
        E200 = term201;
        double term021 = 0.0;
        term021 += qi * R021;
        term021 += dxi * R121;
        term021 += dyi * R031;
        term021 += dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        double term0212 = 0.0;
        term0212 += qxyi * R131;
        term0212 += qxzi * R122;
        term0212 += qyzi * R032;
        term021 += 2.0 * term0212;
        E020 = term021;
        double term003 = 0.0;
        term003 += qi * R003;
        term003 += dxi * R103;
        term003 += dyi * R013;
        term003 += dzi * R004;
        term003 += qxxi * R203;
        term003 += qyyi * R023;
        term003 += qzzi * R005;
        double term0032 = 0.0;
        term0032 += qxyi * R113;
        term0032 += qxzi * R104;
        term0032 += qyzi * R014;
        term003 += 2.0 * term0032;
        E002 = term003;
        double term111 = 0.0;
        term111 += qi * R111;
        term111 += dxi * R211;
        term111 += dyi * R121;
        term111 += dzi * R112;
        term111 += qxxi * R311;
        term111 += qyyi * R131;
        term111 += qzzi * R113;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term1112 += qxzi * R212;
        term1112 += qyzi * R122;
        term111 += 2.0 * term1112;
        E110 = term111;
        double term102 = 0.0;
        term102 += qi * R102;
        term102 += dxi * R202;
        term102 += dyi * R112;
        term102 += dzi * R103;
        term102 += qxxi * R302;
        term102 += qyyi * R122;
        term102 += qzzi * R104;
        double term1022 = 0.0;
        term1022 += qxyi * R212;
        term1022 += qxzi * R203;
        term1022 += qyzi * R113;
        term102 += 2.0 * term1022;
        E101 = term102;
        double term012 = 0.0;
        term012 += qi * R012;
        term012 += dxi * R112;
        term012 += dyi * R022;
        term012 += dzi * R013;
        term012 += qxxi * R212;
        term012 += qyyi * R032;
        term012 += qzzi * R014;
        double term0122 = 0.0;
        term0122 += qxyi * R122;
        term0122 += qxzi * R113;
        term0122 += qyzi * R023;
        term012 += 2.0 * term0122;
        E011 = term012;
    }

    public void EkQI5() {
        double term000 = 0.0;
        term000 += qi * R000;
        term000 += dzi * R001;
        term000 += qxxi * R200;
        term000 += qyyi * R020;
        term000 += qzzi * R002;
        E000 = term000;
        double term100 = 0.0;
        term100 += dxi * R200;
        double term1002 = 0.0;
        term1002 += qxzi * R201;
        term100 += 2.0 * term1002;
        E100 = term100;
        double term010 = 0.0;
        term010 += dyi * R020;
        double term0102 = 0.0;
        term0102 += qyzi * R021;
        term010 += 2.0 * term0102;
        E010 = term010;
        double term001 = 0.0;
        term001 += qi * R001;
        term001 += dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        E001 = term001;
        double term200 = 0.0;
        term200 += qi * R200;
        term200 += dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        E200 = term200;
        double term020 = 0.0;
        term020 += qi * R020;
        term020 += dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        E020 = term020;
        double term002 = 0.0;
        term002 += qi * R002;
        term002 += dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        E002 = term002;
        double term110 = 0.0;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term110 += 2.0 * term1102;
        E110 = term110;
        double term101 = 0.0;
        term101 += dxi * R201;
        double term1012 = 0.0;
        term1012 += qxzi * R202;
        term101 += 2.0 * term1012;
        E101 = term101;
        double term011 = 0.0;
        term011 += dyi * R021;
        double term0112 = 0.0;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E011 = term011;
    }

    public void EiQI5() {
        double term000 = 0.0;
        term000 += qk * R000;
        term000 += dzk * R001;
        term000 += qxxk * R200;
        term000 += qyyk * R020;
        term000 += qzzk * R002;
        E000 = term000;
        double term100 = 0.0;
        term100 += dxk * R200;
        double term1002 = 0.0;
        term1002 += qxzk * R201;
        term100 += 2.0 * term1002;
        E100 = term100;
        double term010 = 0.0;
        term010 += dyk * R020;
        double term0102 = 0.0;
        term0102 += qyzk * R021;
        term010 += 2.0 * term0102;
        E010 = term010;
        double term001 = 0.0;
        term001 += qk * R001;
        term001 += dzk * R002;
        term001 += qxxk * R201;
        term001 += qyyk * R021;
        term001 += qzzk * R003;
        E001 = term001;
        double term200 = 0.0;
        term200 += qk * R200;
        term200 += dzk * R201;
        term200 += qxxk * R400;
        term200 += qyyk * R220;
        term200 += qzzk * R202;
        E200 = term200;
        double term020 = 0.0;
        term020 += qk * R020;
        term020 += dzk * R021;
        term020 += qxxk * R220;
        term020 += qyyk * R040;
        term020 += qzzk * R022;
        E020 = term020;
        double term002 = 0.0;
        term002 += qk * R002;
        term002 += dzk * R003;
        term002 += qxxk * R202;
        term002 += qyyk * R022;
        term002 += qzzk * R004;
        E002 = term002;
        double term110 = 0.0;
        double term1102 = 0.0;
        term1102 += qxyk * R220;
        term110 += 2.0 * term1102;
        E110 = term110;
        double term101 = 0.0;
        term101 += dxk * R201;
        double term1012 = 0.0;
        term1012 += qxzk * R202;
        term101 += 2.0 * term1012;
        E101 = term101;
        double term011 = 0.0;
        term011 += dyk * R021;
        double term0112 = 0.0;
        term0112 += qyzk * R022;
        term011 += 2.0 * term0112;
        E011 = term011;
    }

    public void ExQI5() {
        double term100 = 0.0;
        term100 += dxi * R200;
        double term1002 = 0.0;
        term1002 += qxzi * R201;
        term100 += 2.0 * term1002;
        E000 = term100;
        double term200 = 0.0;
        term200 += qi * R200;
        term200 += dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        E100 = term200;
        double term110 = 0.0;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term110 += 2.0 * term1102;
        E010 = term110;
        double term101 = 0.0;
        term101 += dxi * R201;
        double term1012 = 0.0;
        term1012 += qxzi * R202;
        term101 += 2.0 * term1012;
        E001 = term101;
        double term300 = 0.0;
        term300 += dxi * R400;
        double term3002 = 0.0;
        term3002 += qxzi * R401;
        term300 += 2.0 * term3002;
        E200 = term300;
        double term120 = 0.0;
        term120 += dxi * R220;
        double term1202 = 0.0;
        term1202 += qxzi * R221;
        term120 += 2.0 * term1202;
        E020 = term120;
        double term102 = 0.0;
        term102 += dxi * R202;
        double term1022 = 0.0;
        term1022 += qxzi * R203;
        term102 += 2.0 * term1022;
        E002 = term102;
        double term210 = 0.0;
        term210 += dyi * R220;
        double term2102 = 0.0;
        term2102 += qyzi * R221;
        term210 += 2.0 * term2102;
        E110 = term210;
        double term201 = 0.0;
        term201 += qi * R201;
        term201 += dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        E101 = term201;
        double term111 = 0.0;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term111 += 2.0 * term1112;
        E011 = term111;
    }

    public void EyQI5() {
        double term010 = 0.0;
        term010 += dyi * R020;
        double term0102 = 0.0;
        term0102 += qyzi * R021;
        term010 += 2.0 * term0102;
        E000 = term010;
        double term110 = 0.0;
        double term1102 = 0.0;
        term1102 += qxyi * R220;
        term110 += 2.0 * term1102;
        E100 = term110;
        double term020 = 0.0;
        term020 += qi * R020;
        term020 += dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        E010 = term020;
        double term011 = 0.0;
        term011 += dyi * R021;
        double term0112 = 0.0;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E001 = term011;
        double term210 = 0.0;
        term210 += dyi * R220;
        double term2102 = 0.0;
        term2102 += qyzi * R221;
        term210 += 2.0 * term2102;
        E200 = term210;
        double term030 = 0.0;
        term030 += dyi * R040;
        double term0302 = 0.0;
        term0302 += qyzi * R041;
        term030 += 2.0 * term0302;
        E020 = term030;
        double term012 = 0.0;
        term012 += dyi * R022;
        double term0122 = 0.0;
        term0122 += qyzi * R023;
        term012 += 2.0 * term0122;
        E002 = term012;
        double term120 = 0.0;
        term120 += dxi * R220;
        double term1202 = 0.0;
        term1202 += qxzi * R221;
        term120 += 2.0 * term1202;
        E110 = term120;
        double term111 = 0.0;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term111 += 2.0 * term1112;
        E101 = term111;
        double term021 = 0.0;
        term021 += qi * R021;
        term021 += dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        E011 = term021;
    }

    public void EzQI5() {
        double term001 = 0.0;
        term001 += qi * R001;
        term001 += dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        E000 = term001;
        double term101 = 0.0;
        term101 += dxi * R201;
        double term1012 = 0.0;
        term1012 += qxzi * R202;
        term101 += 2.0 * term1012;
        E100 = term101;
        double term011 = 0.0;
        term011 += dyi * R021;
        double term0112 = 0.0;
        term0112 += qyzi * R022;
        term011 += 2.0 * term0112;
        E010 = term011;
        double term002 = 0.0;
        term002 += qi * R002;
        term002 += dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        E001 = term002;
        double term201 = 0.0;
        term201 += qi * R201;
        term201 += dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        E200 = term201;
        double term021 = 0.0;
        term021 += qi * R021;
        term021 += dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        E020 = term021;
        double term003 = 0.0;
        term003 += qi * R003;
        term003 += dzi * R004;
        term003 += qxxi * R203;
        term003 += qyyi * R023;
        term003 += qzzi * R005;
        E002 = term003;
        double term111 = 0.0;
        double term1112 = 0.0;
        term1112 += qxyi * R221;
        term111 += 2.0 * term1112;
        E110 = term111;
        double term102 = 0.0;
        term102 += dxi * R202;
        double term1022 = 0.0;
        term1022 += qxzi * R203;
        term102 += 2.0 * term1022;
        E101 = term102;
        double term012 = 0.0;
        term012 += dyi * R022;
        double term0122 = 0.0;
        term0122 += qyzi * R023;
        term012 += 2.0 * term0122;
        E011 = term012;
    }

    public void torqueI(double torque[]) {
        // Torque on dipole moments due to the field.
        double dx = dyi * E001 - dzi * E010;
        double dy = dzi * E100 - dxi * E001;
        double dz = dxi * E010 - dyi * E100;
        // Torque on quadrupole moments due to the gradient of the field.
        double qx = qxyi * E101 + qyyi * E011 + qyzi * E002
                - (qxzi * E110 + qyzi * E020 + qzzi * E011);
        double qy = qxzi * E200 + qyzi * E110 + qzzi * E101
                - (qxxi * E101 + qxyi * E011 + qxzi * E002);
        double qz = qxxi * E110 + qxyi * E020 + qxzi * E011
                - (qxyi * E200 + qyyi * E110 + qyzi * E101);
        torque[0] = dx + 2.0 * qx;
        torque[1] = dy + 2.0 * qy;
        torque[2] = dz + 2.0 * qz;
    }

    public void torqueK(double torque[]) {
        // Torque on dipole moments due to the field.
        double dx = dyk * E001 - dzk * E010;
        double dy = dzk * E100 - dxk * E001;
        double dz = dxk * E010 - dyk * E100;
        // Torque on quadrupole moments due to the gradkent of the fkeld.
        double qx = qxyk * E101 + qyyk * E011 + qyzk * E002
                - (qxzk * E110 + qyzk * E020 + qzzk * E011);
        double qy = qxzk * E200 + qyzk * E110 + qzzk * E101
                - (qxxk * E101 + qxyk * E011 + qxzk * E002);
        double qz = qxxk * E110 + qxyk * E020 + qxzk * E011
                - (qxyk * E200 + qyyk * E110 + qyzk * E101);
        torque[0] = dx + 2.0 * qx;
        torque[1] = dy + 2.0 * qy;
        torque[2] = dz + 2.0 * qz;
    }

    /**
     * Contract multipole moments with their respective electrostatic potential
     * derivatives.
     *
     * @param Q array of Cartesian multipole moments
     * @param T array of electrostatic potential and partial derivatives
     * @param l apply (d/dx)^l to the potential
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     * @return the contracted interaction.
     */
    public double contract(double T[], int l, int m, int n) {
        double total = 0.0;
        double total2 = 0.0;
        total += qi * T[ti(l, m, n)];
        total += dxi * T[ti(l + 1, m, n)];
        total += dyi * T[ti(l, m + 1, n)];
        total += dzi * T[ti(l, m, n + 1)];
        total += qxxi * T[ti(l + 2, m, n)];
        total += qyyi * T[ti(l, m + 2, n)];
        total += qzzi * T[ti(l, m, n + 2)];
        total2 += qxyi * T[ti(l + 1, m + 1, n)];
        total2 += qxzi * T[ti(l + 1, m, n + 1)];
        total2 += qyzi * T[ti(l, m + 1, n + 1)];
        return total + 2.0 * total2;
    }

    /**
     * Contract multipole moments with their respective electrostatic potential
     * derivatives.
     *
     * @param Q array of Cartesian multipole moments
     * @param T array of electrostatic potential and partial derivatives
     * @param l apply (d/dx)^l to the potential
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     * @param sb the code will be appended to the StringBuilfer.
     * @return the contracted interaction.
     */
    public double codeContract(double T[], int l, int m, int n, StringBuilder sb) {
        double total = 0.0;
        String name = term(l, m, n);
        sb.append(format("double %s = 0.0;\n", name));
        StringBuilder sb1 = new StringBuilder();
        double term = qi * T[ti(l, m, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q000i * T[%s];\n", name, tlmn(l, m, n)));
        }
        term = dxi * T[ti(l + 1, m, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q100i * T[%s];\n", name, tlmn(l + 1, m, n)));
        }
        term = dyi * T[ti(l, m + 1, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q010i * T[%s];\n", name, tlmn(l, m + 1, n)));
        }
        term = dzi * T[ti(l, m, n + 1)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q001i * T[%s];\n", name, tlmn(l, m, n + 1)));
        }
        StringBuilder traceSB = new StringBuilder();
        double trace = 0.0;
        term = qxxi * T[ti(l + 2, m, n)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qxx: %16.15f T: %16.15f Term: %16.15f", q200, T[ti(l + 2, m, n)], term));
            traceSB.append(format("%s += q200i * T[%s];\n", name, tlmn(l + 2, m, n)));
        }
        term = qyyi * T[ti(l, m + 2, n)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qyy: %16.15f T: %16.15f Term: %16.15f", q020, T[ti(l, m + 2, n)], term));
            traceSB.append(format("%s += q020i * T[%s];\n", name, tlmn(l, m + 2, n)));
        }
        term = qzzi * T[ti(l, m, n + 2)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qzz: %16.15f T: %16.15f Term: %16.15f", q002, T[ti(l, m, n + 2)], term));
            traceSB.append(format("%s += q002i * T[%s];\n", name, tlmn(l, m, n + 2)));
        }
        total += trace;
        if (total != 0) {
            sb.append(sb1.toString());
            if (trace != 0) {
                //logger.info(format(" Trace: %16.15f", trace));
                sb.append(traceSB);
            }
        }
        StringBuilder sb2 = new StringBuilder();
        double total2 = 0.0;
        term = qxyi * T[ti(l + 1, m + 1, n)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q110i * T[%s];\n", name, tlmn(l + 1, m + 1, n)));
        }
        term = qxzi * T[ti(l + 1, m, n + 1)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q101i * T[%s];\n", name, tlmn(l + 1, m, n + 1)));
        }
        term = qyzi * T[ti(l, m + 1, n + 1)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q011i * T[%s];\n", name, tlmn(l, m + 1, n + 1)));
        }
        if (total2 != 0.0) {
            sb.append(format("double %s2 = 0.0;\n", name));
            sb.append(sb2);
            total += 2.0 * total2;
            sb.append(format("%s += 2.0 * %s2;\n", name, name));
        }
        return total;
    }

    /**
     * Collect the field at R due to Q multipole moments at the origin.
     *
     * @param Q Cartesian multipole moments at the origin.
     * @param T Electrostatic potential and partial derivatives
     * @param E The components of the field at R.
     * @param l apply (d/dx)^l to the potential
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     */
    public void field(double T[], int l, int m, int n) {
        E000 = contract(T, l, m, n);
        E100 = contract(T, l + 1, m, n);
        E010 = contract(T, l, m + 1, n);
        E001 = contract(T, l, m, n + 1);
        E200 = contract(T, l + 2, m, n);
        E020 = contract(T, l, m + 2, n);
        E002 = contract(T, l, m, n + 2);
        E110 = contract(T, l + 1, n + 1, m);
        E101 = contract(T, l + 1, m, n + 1);
        E011 = contract(T, l, m + 1, n + 1);
        return;
    }

    /**
     * Collect the field at R due to Q multipole moments at the origin.
     *
     * @param Q Cartesian multipole moments at the origin.
     * @param T Electrostatic potential and partial derivatives
     * @param E The components of the field at R.
     * @param l apply (d/dx)^l to the potential
     * @param sb
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     */
    public void codeField(double T[], int l, int m, int n, StringBuilder sb) {
        E000 = codeContract(T, l, m, n, sb);
        if (E000 != 0) {
            sb.append(format("e000 = %s;\n", term(l, m, n)));
        }
        E100 = codeContract(T, l + 1, m, n, sb);
        if (E100 != 0) {
            sb.append(format("e100 = %s;\n", term(l + 1, m, n)));
        }
        E010 = codeContract(T, l, m + 1, n, sb);
        if (E100 != 0) {
            sb.append(format("e010 = %s;\n", term(l, m + 1, n)));
        }
        E001 = codeContract(T, l, m, n + 1, sb);
        if (E001 != 0) {
            sb.append(format("e001 = %s;\n", term(l, m, n + 1)));
        }
        E200 = codeContract(T, l + 2, m, n, sb);
        if (E200 != 0) {
            sb.append(format("e200 = %s;\n", term(l + 2, m, n)));
        }
        E020 = codeContract(T, l, m + 2, n, sb);
        if (E020 != 0) {
            sb.append(format("e020 = %s;\n", term(l, m + 2, n)));
        }
        E002 = codeContract(T, l, m, n + 2, sb);
        if (E002 != 0) {
            sb.append(format("e002 = %s;\n", term(l, m, n + 2)));
        }
        E110 = codeContract(T, l + 1, m + 1, n, sb);
        if (E110 != 0) {
            sb.append(format("e110 = %s;\n", term(l + 1, m + 1, n)));
        }
        E101 = codeContract(T, l + 1, m, n + 1, sb);
        if (E101 != 0) {
            sb.append(format("e101 = %s;\n", term(l + 1, m, n + 1)));
        }
        E011 = codeContract(T, l, m + 1, n + 1, sb);
        if (E011 != 0) {
            sb.append(format("e011 = %s;\n", term(l, m + 1, n + 1)));
        }
    }

    private double dotK() {
        double total = qk * E000;
        total += dxk * E100;
        total += dyk * E010;
        total += dzk * E001;
        total += qxxk * E200;
        total += qyyk * E020;
        total += qzzk * E002;
        double total2 = qxyi * E110;
        total2 += qxzk * E101;
        total2 += qyzk * E011;
        return total + 2.0 * total2;
    }

    public double codeInteract5(double r[], double Qi[], double Qk[],
            double Fi[], double Fk[], double Ti[], double Tk[]) {
        double T[] = new double[tensorCount(5)];
        recursion(r, T);

        setMultipoleI(Qi);
        setMultipoleK(Qk);

        StringBuilder sb = new StringBuilder("\n\npublic void E5(double T[]) {\n");
        codeField(T, 0, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ex5(double T[]) {\n");
        codeField(T, 1, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ey5(double T[]) {\n");
        codeField(T, 0, 1, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ez5(double T[]) {\n");
        codeField(T, 0, 0, 1, sb);
        sb.append("}\n");
        logger.log(Level.INFO, sb.toString());

        return 0.0;
    }

    private void setMultipoleI(double Qi[]) {
        qi = Qi[0];
        dxi = Qi[1];
        dyi = Qi[2];
        dzi = Qi[3];
        qxxi = Qi[4];
        qyyi = Qi[5];
        qzzi = Qi[6];
        qxyi = Qi[7];
        qxzi = Qi[8];
        qyzi = Qi[9];
    }

    private void setMultipoleK(double Qk[]) {
        qk = Qk[0];
        dxk = Qk[1];
        dyk = Qk[2];
        dzk = Qk[3];
        qxxk = Qk[4];
        qyyk = Qk[5];
        qzzk = Qk[6];
        qxyk = Qk[7];
        qxzk = Qk[8];
        qyzk = Qk[9];
    }

    public double codeInteractQI5(double r[], double Qi[], double Qk[],
            double Fi[], double Fk[], double Ti[], double Tk[]) {
        double T[] = new double[tensorCount(5)];
        assert (r[0] == 0.0 && r[1] == 0.0);

        recursionQI(r, T);
        setMultipoleI(Qi);
        setMultipoleK(Qk);

        StringBuilder sb = new StringBuilder("\n\npublic void EQI5(double T[]) {\n");
        codeField(T, 0, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void ExQI5(double T[]) {\n");
        codeField(T, 1, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void EyQI5(double T[]) {\n");
        codeField(T, 0, 1, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void EzQI5(double T[]) {\n");
        codeField(T, 0, 0, 1, sb);
        sb.append("}\n");
        logger.log(Level.INFO, sb.toString());

        return 0.0;
    }

    public double interact5(double r[], double Qi[], double Qk[],
            double Fi[], double Fk[], double Ti[], double Tk[]) {
        order5(r);
        setMultipoleI(Qi);
        setMultipoleK(Qk);

        Ek5();
        double energy = dotK();

        // Torques
        torqueK(Tk);
        Ei5();
        torqueI(Tk);

        // Forces
        Ex5();
        Fi[0] = -dotK();
        Ey5();
        Fi[1] = -dotK();
        Ez5();
        Fi[2] = -dotK();

        Fk[0] = -Fi[0];
        Fk[1] = -Fi[1];
        Fk[2] = -Fi[2];

        return energy;
    }

    public double interact5QI(double r[], double Qi[], double Qk[],
            double Fi[], double Fk[], double Ti[], double Tk[]) {

        // Find the QI Tensor.
        double z[] = {0.0, 0.0, r(r)};
        order5QI(z);

        // Find the rotation matrix from the Global frame to the QI frame.
        setQIRotationMatrix(r);

        // Rotate multipole I and K.
        multipoleItoQI(Qi);
        multipoleKtoQI(Qk);

        // Compute the potential due to site I at site K.
        EiQI5();

        // Dot the potential, field, field gradient with multipole K.
        double energy = dotK();

        // Compute the torque on site K due to the field from site I.
        torqueK(Tk);

        // Compute the field at site I due to site K.
        EkQI5();
        // Compute the torque on site I due to the field from site K.
        torqueI(Ti);

        // Compute the force on site I F = {-dE/dx, -dE/dy, -dE/dz}.
        ExQI5();
        Fi[0] = -dotK();
        EyQI5();
        Fi[1] = -dotK();
        EzQI5();
        Fi[2] = -dotK();

        // Rotate the force and torques from the QI frame into the Global frame.
        vectorsToGlobal(Fi, Ti, Tk);

        // Set the force on site K as -Fi.
        Fk[0] = -Fi[0];
        Fk[1] = -Fi[1];
        Fk[2] = -Fi[2];

        return energy;
    }

    public static void main(String args[]) {
        if (args == null || args.length < 4) {
            logger.info(" Usage: java ffx.numerics.TensorRecursion order dx dy dz");
        }
        int order = Integer.parseInt(args[0]);
        double dx = Double.parseDouble(args[1]);
        double dy = Double.parseDouble(args[2]);
        double dz = Double.parseDouble(args[3]);
        double r[] = {dx, dy, dz};

        double n2 = 710643;
        double cycles = 10;

        MultipoleTensor multipoleTensor;
        multipoleTensor = new MultipoleTensor(OPERATOR.COULOMB, order, 0.0);

        double Fi[] = new double[3];
        double Fk[] = new double[3];
        double Ti[] = new double[3];
        double Tk[] = new double[3];
        double Qi[] = {0.11,
            0.21, 0.31, 0.41,
            -0.51, -0.61, 1.12, 0.71, 0.81, 0.91};
        double Qk[] = {0.11,
            0.21, 0.31, 0.41,
            -0.51, -0.61, 1.12, 0.7, 0.81, 0.91};

        for (int j = 0; j < cycles; j++) {
            long timeGlobalT = -System.nanoTime();
            for (int i = 0; i < n2; i++) {
                r[0] = Math.random();
                r[1] = Math.random();
                r[2] = Math.random();
                multipoleTensor.order5(r);
            }
            timeGlobalT += System.nanoTime();

            long timeGlobal = -System.nanoTime();
            for (int i = 0; i < n2; i++) {
                r[0] = Math.random();
                r[1] = Math.random();
                r[2] = Math.random();
                multipoleTensor.interact5(r, Qi, Qk, Fi, Fk, Ti, Tk);
            }
            timeGlobal += System.nanoTime();

            long timeQIT = -System.nanoTime();
            for (int i = 0; i < n2; i++) {
                r[0] = Math.random();
                r[1] = Math.random();
                r[2] = Math.random();
                multipoleTensor.order5QI(r);
            }
            timeQIT += System.nanoTime();

            long timeQI = -System.nanoTime();
            for (int i = 0; i < n2; i++) {
                r[0] = 0.0;
                r[1] = 0.0;
                r[2] = Math.random();
                multipoleTensor.interact5QI(r, Qi, Qk, Fi, Fk, Ti, Tk);
            }
            timeQI += System.nanoTime();

            logger.info(String.format("\n Global Frame: %6.4f %6.4f\n QI:           %6.4f %6.4f\n",
                    timeGlobalT * 1.0e-9, timeGlobal * 1.0e-9, timeQIT * 1.0e-9, timeQI * 1.0e-9));
        }

        // String string = tensorRecursion.codeTensorRecursion(r, tensors, damp, aiak);
        // logger.info(" Java Code:\n" + string);
        // string = tensorRecursion.codeTensorRecursionQI(r, tensors, damp, aiak);
        // logger.info(" Java Code:\n" + string);
    }

    /**
     * Convenience method for writing out intermediate terms in the recursion.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @param j the jth intermediate term.
     * @return a String of the form <code>termlmnj</code>.
     */
    private static String term(int l, int m, int n, int j) {
        return String.format("term%d%d%d%d", l, m, n, j);
    }

    /**
     * Convenience method for writing out intermediate terms in the recursion.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @param j the jth intermediate term.
     * @return a String of the form <code>termlmnj</code>.
     */
    private static String term(int l, int m, int n) {
        return String.format("term%d%d%d", l, m, n);
    }

    /**
     * Convenience method for writing out tensor indeces.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @return a String of the form <code>tlmn</code>.
     */
    private static String tlmn(int l, int m, int n) {
        return String.format("t%d%d%d", l, m, n);
    }

    public void getTensor(double T[]) {
        if (T == null || order < 0) {
            return;
        }
        // l + m + n = 0 (1)
        T[t000] = R000;
        if (order < 1) {
            return;
        }
        // l + m + n = 1 (3)   4
        T[t100] = R100;
        T[t010] = R010;
        T[t001] = R001;
        if (order < 2) {
            return;
        }
        // l + m + n = 2 (6)  10
        T[t200] = R200;
        T[t020] = R020;
        T[t002] = R002;
        T[t110] = R110;
        T[t101] = R101;
        T[t011] = R011;
        if (order < 3) {
            return;
        }
        // l + m + n = 3 (10) 20
        T[t300] = R300;
        T[t030] = R030;
        T[t003] = R003;
        T[t210] = R210;
        T[t201] = R201;
        T[t120] = R120;
        T[t021] = R021;
        T[t102] = R102;
        T[t012] = R012;
        T[t111] = R111;
        if (order < 4) {
            return;
        }
        // l + m + n = 4 (15) 35
        T[t400] = R400;
        T[t040] = R040;
        T[t004] = R004;
        T[t310] = R310;
        T[t301] = R301;
        T[t130] = R130;
        T[t031] = R031;
        T[t103] = R103;
        T[t013] = R013;
        T[t220] = R220;
        T[t202] = R202;
        T[t022] = R022;
        T[t211] = R211;
        T[t121] = R121;
        T[t112] = R112;
        if (order < 5) {
            return;
        }
        // l + m + n = 5 (21) 56
        T[t500] = R500;
        T[t050] = R050;
        T[t005] = R005;
        T[t410] = R410;
        T[t401] = R401;
        T[t140] = R140;
        T[t041] = R041;
        T[t104] = R104;
        T[t014] = R014;
        T[t320] = R320;
        T[t302] = R302;
        T[t230] = R230;
        T[t032] = R032;
        T[t203] = R203;
        T[t023] = R023;
        T[t311] = R311;
        T[t131] = R131;
        T[t113] = R113;
        T[t221] = R221;
        T[t212] = R212;
        T[t122] = R122;
    }

    public void setTensor(double T[]) {
        if (T == null || order < 0) {
            return;
        }
        // l + m + n = 0 (1)
        R000 = T[t000];
        if (order < 1) {
            return;
        }
        // l + m + n = 1 (3)   4
        R100 = T[t100];
        R010 = T[t010];
        R001 = T[t001];
        if (order < 2) {
            return;
        }
        // l + m + n = 2 (6)  10
        R200 = T[t200];
        R020 = T[t020];
        R002 = T[t002];
        R110 = T[t110];
        R101 = T[t101];
        R011 = T[t011];
        if (order < 3) {
            return;
        }
        // l + m + n = 3 (10) 20
        R300 = T[t300];
        R030 = T[t030];
        R003 = T[t003];
        R210 = T[t210];
        R201 = T[t201];
        R120 = T[t120];
        R021 = T[t021];
        R102 = T[t102];
        R012 = T[t012];
        R111 = T[t111];
        if (order < 4) {
            return;
        }
        // l + m + n = 4 (15) 35
        R400 = T[t400];
        R040 = T[t040];
        R004 = T[t004];
        R310 = T[t310];
        R301 = T[t301];
        R130 = T[t130];
        R031 = T[t031];
        R103 = T[t103];
        R013 = T[t013];
        R220 = T[t220];
        R202 = T[t202];
        R022 = T[t022];
        R211 = T[t211];
        R121 = T[t121];
        R112 = T[t112];
        if (order < 5) {
            return;
        }
        // l + m + n = 5 (21) 56
        R500 = T[t500];
        R050 = T[t050];
        R005 = T[t005];
        R410 = T[t410];
        R401 = T[t401];
        R140 = T[t140];
        R041 = T[t041];
        R104 = T[t104];
        R014 = T[t014];
        R320 = T[t320];
        R302 = T[t302];
        R230 = T[t230];
        R032 = T[t032];
        R203 = T[t203];
        R023 = T[t023];
        R311 = T[t311];
        R131 = T[t131];
        R113 = T[t113];
        R221 = T[t221];
        R212 = T[t212];
        R122 = T[t122];
    }

    private void setQIRotationMatrix(double r[]) {
        double zAxis[] = {r[0], r[1], r[2]};
        double xAxis[] = new double[3];
        xAxis[0] = r[0] + 1.0;
        xAxis[1] = r[1];
        xAxis[2] = r[2];

        norm(zAxis, zAxis);
        ir02 = zAxis[0];
        ir12 = zAxis[1];
        ir22 = zAxis[2];

        double dot = dot(xAxis, zAxis);
        scalar(zAxis, dot, zAxis);
        diff(xAxis, zAxis, xAxis);
        norm(xAxis, xAxis);

        ir00 = xAxis[0];
        ir10 = xAxis[1];
        ir20 = xAxis[2];
        ir01 = ir20 * ir12 - ir10 * ir22;
        ir11 = ir00 * ir22 - ir20 * ir02;
        ir21 = ir10 * ir02 - ir00 * ir12;

        // Set the forward elements as the transpose of the inverse matrix.
        r00 = ir00;
        r11 = ir11;
        r22 = ir22;
        r01 = ir10;
        r02 = ir20;
        r10 = ir01;
        r12 = ir21;
        r20 = ir02;
        r21 = ir12;

    }

    private void multipoleItoQI(double Qi[]) {

        qi = Qi[0];

        double dx = Qi[1];
        double dy = Qi[2];
        double dz = Qi[3];

        dxi = r00 * dx + r01 * dy + r02 * dz;
        dyi = r10 * dx + r11 * dy + r12 * dz;
        dzi = r20 * dx + r21 * dy + r22 * dz;

        double qxx = Qi[4];
        double qyy = Qi[5];
        double qzz = Qi[6];
        double qxy = Qi[7];
        double qxz = Qi[8];
        double qyz = Qi[9];

        // i=0, j=0
        // qij   r0k *  r00 * qkx + r01 * qky + r02 * qkz
        qxxi = r00 * (r00 * qxx + r01 * qxy + r02 * qxz)
                + r01 * (r00 * qxy + r01 * qyy + r02 * qyz)
                + r02 * (r00 * qxz + r01 * qyz + r02 * qzz);

        // i=0, j=1
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxyi = r00 * (r10 * qxx + r11 * qxy + r12 * qxz)
                + r01 * (r10 * qxy + r11 * qyy + r12 * qyz)
                + r02 * (r10 * qxz + r11 * qyz + r12 * qzz);

        // i=0, j=2
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxzi = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);

        // i=1, j=1
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qyyi = r10 * (r10 * qxx + r11 * qxy + r12 * qxz)
                + r11 * (r10 * qxy + r11 * qyy + r12 * qyz)
                + r12 * (r10 * qxz + r11 * qyz + r12 * qzz);

        // i=1, j=2
        // qij   r1k *  r20 * qkx + r21 * qky + r22 * qkz
        qyzi = r10 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r11 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r12 * (r20 * qxz + r21 * qyz + r22 * qzz);

        // i=2, j=2
        // qij   r2k *  r20 * qkx + r21 * qky + r22 * qkz
        qzzi = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);
    }

    private void multipoleKtoQI(double Qk[]) {

        qk = Qk[0];

        double dx = Qk[1];
        double dy = Qk[2];
        double dz = Qk[3];

        dxk = r00 * dx + r01 * dy + r02 * dz;
        dyk = r10 * dx + r11 * dy + r12 * dz;
        dzk = r20 * dx + r21 * dy + r22 * dz;

        double qxx = Qk[4];
        double qyy = Qk[5];
        double qzz = Qk[6];
        double qxy = Qk[7];
        double qxz = Qk[8];
        double qyz = Qk[9];

        // i=0, j=0
        // qij   r0k *  r00 * qkx + r01 * qky + r02 * qkz
        qxxk = r00 * (r00 * qxx + r01 * qxy + r02 * qxz)
                + r01 * (r00 * qxy + r01 * qyy + r02 * qyz)
                + r02 * (r00 * qxz + r01 * qyz + r02 * qzz);

        // i=0, j=1
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxyk = r00 * (r10 * qxx + r11 * qxy + r12 * qxz)
                + r01 * (r10 * qxy + r11 * qyy + r12 * qyz)
                + r02 * (r10 * qxz + r11 * qyz + r12 * qzz);

        // i=0, j=2
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxzk = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);

        // i=1, j=1
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qyyk = r10 * (r10 * qxx + r11 * qxy + r12 * qxz)
                + r11 * (r10 * qxy + r11 * qyy + r12 * qyz)
                + r12 * (r10 * qxz + r11 * qyz + r12 * qzz);

        // i=1, j=2
        // qij   r1k *  r20 * qkx + r21 * qky + r22 * qkz
        qyzk = r10 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r11 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r12 * (r20 * qxz + r21 * qyz + r22 * qzz);

        // i=2, j=2
        // qij   r2k *  r20 * qkx + r21 * qky + r22 * qkz
        qzzk = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);
    }

    private void vectorsToGlobal(double v1[], double v2[],
            double v3[]) {
        double vx = v1[0];
        double vy = v1[1];
        double vz = v1[2];
        v1[0] = ir00 * vx + ir01 * vy + ir02 * vz;
        v1[1] = ir10 * vx + ir11 * vy + ir12 * vz;
        v1[2] = ir20 * vx + ir21 * vy + ir22 * vz;

        vx = v2[0];
        vy = v2[1];
        vz = v2[2];
        v2[0] = ir00 * vx + ir01 * vy + ir02 * vz;
        v2[1] = ir10 * vx + ir11 * vy + ir12 * vz;
        v2[2] = ir20 * vx + ir21 * vy + ir22 * vz;

        vx = v3[0];
        vy = v3[1];
        vz = v3[2];
        v3[0] = ir00 * vx + ir01 * vy + ir02 * vz;
        v3[1] = ir10 * vx + ir11 * vy + ir12 * vz;
        v3[2] = ir20 * vx + ir21 * vy + ir22 * vz;
    }

    // Multipole components for atom i.
    private double qi;
    private double dxi;
    private double dyi;
    private double dzi;
    private double qxxi;
    private double qyyi;
    private double qzzi;
    private double qxyi;
    private double qxzi;
    private double qyzi;

    // Multipole components for atom k.
    private double qk;
    private double dxk;
    private double dyk;
    private double dzk;
    private double qxxk;
    private double qyyk;
    private double qzzk;
    private double qxyk;
    private double qxzk;
    private double qyzk;

    // Rotation Matrix from Global to QI.
    private double r00, r01, r02;
    private double r10, r11, r12;
    private double r20, r21, r22;

    // Rotation Matrix from QI to Global.
    private double ir00, ir01, ir02;
    private double ir10, ir11, ir12;
    private double ir20, ir21, ir22;

    // Cartesian tensor elements (for 1/R, erfc(Beta*R)/R or Thole damping.
    // l + m + n = 0 (1)
    private double R000;
    // l + m + n = 1 (3)   4
    private double R100;
    private double R010;
    private double R001;
    // l + m + n = 2 (6)  10
    private double R200;
    private double R020;
    private double R002;
    private double R110;
    private double R101;
    private double R011;
    // l + m + n = 3 (10) 20
    private double R300;
    private double R030;
    private double R003;
    private double R210;
    private double R201;
    private double R120;
    private double R021;
    private double R102;
    private double R012;
    private double R111;
    // l + m + n = 4 (15) 35
    private double R400;
    private double R040;
    private double R004;
    private double R310;
    private double R301;
    private double R130;
    private double R031;
    private double R103;
    private double R013;
    private double R220;
    private double R202;
    private double R022;
    private double R211;
    private double R121;
    private double R112;
    // l + m + n = 5 (21) 56
    private double R500;
    private double R050;
    private double R005;
    private double R410;
    private double R401;
    private double R140;
    private double R041;
    private double R104;
    private double R014;
    private double R320;
    private double R302;
    private double R230;
    private double R032;
    private double R203;
    private double R023;
    private double R311;
    private double R131;
    private double R113;
    private double R221;
    private double R212;
    private double R122;

    // Components of the potential, field and field gradient.
    private double E000; // Potential
    private double E100; // X Component of the Field
    private double E010; // Y Component of the Field
    private double E001; // Z Component of the Field
    private double E200; // XX Component of the Field Gradient
    private double E020; // YY Component of the Field Gradient
    private double E002; // ZZ Component of the Field Gradient
    private double E110; // XY Component of the Field Gradient
    private double E101; // XZ Component of the Field Gradient
    private double E011; // YZ Component of the Field Gradient

    // l + m + n = 0 (1)
    public final int t000;
    // l + m + n = 1 (3)   4
    public final int t100;
    public final int t010;
    public final int t001;
    // l + m + n = 2 (6)  10
    public final int t200;
    public final int t020;
    public final int t002;
    public final int t110;
    public final int t101;
    public final int t011;
    // l + m + n = 3 (10) 20
    public final int t300;
    public final int t030;
    public final int t003;
    public final int t210;
    public final int t201;
    public final int t120;
    public final int t021;
    public final int t102;
    public final int t012;
    public final int t111;
    // l + m + n = 4 (15) 35
    public final int t400;
    public final int t040;
    public final int t004;
    public final int t310;
    public final int t301;
    public final int t130;
    public final int t031;
    public final int t103;
    public final int t013;
    public final int t220;
    public final int t202;
    public final int t022;
    public final int t211;
    public final int t121;
    public final int t112;
    // l + m + n = 5 (21) 56
    public final int t500;
    public final int t050;
    public final int t005;
    public final int t410;
    public final int t401;
    public final int t140;
    public final int t041;
    public final int t104;
    public final int t014;
    public final int t320;
    public final int t302;
    public final int t230;
    public final int t032;
    public final int t203;
    public final int t023;
    public final int t311;
    public final int t131;
    public final int t113;
    public final int t221;
    public final int t212;
    public final int t122;

}
