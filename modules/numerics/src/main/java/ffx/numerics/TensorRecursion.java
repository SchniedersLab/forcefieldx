/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.numerics;

import static java.lang.Math.sqrt;

/**
 * The TensorRecusion class compute derivates of 1/|<b>r</b>| via recursion to
 * arbitrary order.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TensorRecursion {

    /**
     * The index is based on the idea of filling tetrahedron.
     * <p>
     * 1/r has an index of 0 <br>
     * derivatives of x are first -> indeces from 1..o for d/dx..do/dox) <br>
     * derivatives of x & y are second -> base triangle of size (o+1)(o+2)/2 <br>
     * derivatives of x & y & z are last -> total size (o+1)*(o+2)*(o+3)/6 <br>
     * <p>
     * This function is useful to set up masking constants:<br>
     * static int Tlmn = tensorIndex(l,m,n,order) <br>
     * For example the (d/dy)^2 (1/R) storage location: <br>
     * static int T020 = tensorIndex(0,2,0,order)
     * <p>
     *
     * @param dx int The number of d/dx operations.
     * @param dy int The number of d/dy operations.
     * @param dz int The number of d/dz operations.
     * @param order int The maximum tensor order (0 <= dx + dy + dz <= order).
     *
     * @return int in the range (0..binomial(order + 3, 3) - 1)
     * @since 1.0
     */
    public static int tensorIndex(int dx, int dy, int dz, int order) {
        int size = (order + 1) * (order + 2) * (order + 3) / 6;
        // We only get to the top of the tetrahedron if dz = order,
        // otherwise subtract off the top, including the level of the requested
        // tensor index.
        int top = order + 1 - dz;
        top = top * (top + 1) * (top + 2) / 6;
        int zindex = size - top;
        // Given the "dz level", dy can range from 0..order - dz)
        // To get to the row for a specific value of dy,
        // dy*(order + 1) - dy*(dy-1)/2 indeces are skipped.
        // This is an operation that looks like the area of rectangle, minus
        // the area of an empty triangle.
        int yindex = dy * (order - dz) - (dy - 1) * (dy - 2) / 2 + 1;
        // Given the dz level and dy row, dx can range from (0..order - dz - dy)
        // The dx index is just walking down the dy row for "dx" steps.
        int ret = dx + yindex + zindex;
        return ret;
    }

    /**
     * Returns the number of tensors for derivatives to the given order.
     * @param order maximum number of derivatives.
     * @return the number of tensors.
     *
     * @since 1.0
     */
    public static int tensorCount(int order) {
        long ret = VectorMath.binomial(order + 3, 3);
        assert (ret < Integer.MAX_VALUE);
        return (int) ret;
    }

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
    private final int order;
    private final double t000j_Constants[];
    private final int o1;
    private final int il;
    private final int im;
    private final int in;

    public TensorRecursion(int order) {
        assert(order > 0);
        o1 = order + 1;
        il = o1;
        im = il * o1;
        in = im * o1;
        work = new double[in * o1];
        t000j_Constants = new double[o1];
        for (int j = 0; j <= order; j++) {
            // Math.pow(-1.0, j) returns positive for all j, with -1.0 as the
            // arguement rather than -1. This is a bug.
            t000j_Constants[j] = Math.pow(-1, j) * VectorMath.doublefactorial(2 * j - 1);
        }
        this.order = order;
        T000 = new double[order + 1];
    }

    /**
     * This method is a driver to collect elements of the Cartesion multipole
     * tensor given the recursion relationships implemented by the method
     * "Tlmnj", which can be called directly to get a single tensor element. It
     * does not store intermediate values of the recursion, causing it to scale
     * O(order^8). For order = 5, this approach is a factor of 10 slower than
     * tensorRecursion.
     *
     * @param r
     *            double[] vector between two sites.
     * @param tensor
     *            double[] length must be at least binomial(order + 3, 3).
     */
    public void noStorageTensorRecursion(double r[], double tensor[]) {
        // 1/r
        double rr = 1.0 / VectorMath.r(r);
        // 1/r^2
        double rr2 = rr * rr;
        // Create the auxillary tensors elements (equation 40).
        for (int j = 0; j <= order; j++) {
            T000[j] = t000j_Constants[j] * rr;
            rr = rr * rr2;
        }
        // 1/r
        tensor[0] = T000[0];
        // Start the tensor index rolling.
        int index = 1;
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        for (int l = 1; l <= order; l++) {
            tensor[index++] = Tlmnj(l, 0, 0, 0, r, T000);
        }
        // Find (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
        for (int l = 0; l <= o1; l++) {
            for (int m = 1; m <= order - l; m++) {
                tensor[index++] = Tlmnj(l, m, 0, 0, r, T000);
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l <= o1; l++) {
            for (int m = 0; m <= o1 - l; m++) {
                for (int n = 1; n <= order - l - m; n++) {
                    tensor[index++] = Tlmnj(l, m, n, 0, r, T000);
                }
            }
        }
    }

    /**
     * This function is a driver to collect elements of the Cartesion multipole
     * tensor. Collecting all tensors scales slightly better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example,
     * up to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1<br>
     * tensor[0] = 1/|r| <br>
     * tensor[1] = -x/|r|^3 <br>
     * tensor[2] = -y/|r|^3 <br>
     * tensor[3] = -z/|r|^3 <br>
     * <p>
     * @param r
     *            double[] vector between two sites.
     * @param tensor
     *            double[] length must be at least binomial(order + 3, 3).
     * @since 1.0
     */
    public void tensorRecursion(final double r[], final double tensor[]) {
        int a;
        final double x = r[0];
        final double y = r[1];
        final double z = r[2];
        int n;
        int m;
        int l;
        final double rr2 = 1.0 / (x * x + y * y + z * z);
        double rr = sqrt(rr2);
        // Create the auxillary tensors elements (equation 40).
        for (int j = 0; j < o1; j++) {
            work[j] = t000j_Constants[j] * rr;
            rr *= rr2;
        }
        tensor[0] = work[0];
        // Start the tensor index rolling.
        int index = 1;
        // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
        // Any (d/dx) term can be formed as
        // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
        // All intermediate terms are indexed as l*il + m*im + n*in + j;
        // Store the l=1 tensor T100 (d/dx)
        // Starting the loop at l=2 avoids an if statement.
        double current;
        double previous = work[1];
        tensor[index++] = x * previous;
        for (l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=2 avoid an if statement.
            // T100(l-1) = x * T000(l)
            current = x * work[l];
            int iw = il + l - 1;
            work[iw] = current;
            for (a = 1; a < l - 1; a++) {
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
            tensor[index++] = x * current + (l - 1) * previous;
            previous = current;
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[index++] = y * previous;
            for (m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = y * work[iw];
                iw += im - 1;
                work[iw] = current;
                for (a = 1; a < m - 1; a++) {
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
                tensor[index++] = y * current + (m - 1) * previous;
                previous = current;
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
        // Any (d/dz) term can be formed as:
        // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
        for (l = 0; l < order; l++) {
            for (m = 0; m + l < order; m++) {
                // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
                // Tlmn = z * Tlm01
                final int lm = m + l;
                final int lilmim = l * il + m * im;
                previous = work[lilmim + 1];
                tensor[index++] = z * previous;
                for (n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    final int n1 = n - 1;
                    for (a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 n - 1) * Tlm(n-2)1
                    tensor[index++] = z * current + n1 * previous;
                    previous = current;
                }
            }
        }
    }
    
    /**
     * This routine implements the recurrence relations for computation of
     * any Cartesion multipole tensor in ~O(L^8) time, where L is the total
     * order l + m + n, given the auxillary elements T0000.
     * <p>
     * It implements the recursion relationships from the reference below
     * in brute force fashion, without saving intermediate values. This is
     * useful for finding a single tensor, rather than all binomial(L + 3, 3).
     * <p>
     * The specific recursion equations (41-43) and set of auxillary tensor
     * elements from equation (40) can be found in:
     * <p>
     * Matt Challacombe, Eric Schwegler and Jan Almlof, Modern developments in
     * Hartree-Fock theory: Fast methods for computing the
     * Coulomb matrix. Computational Chemistry: Review of Current Trends pp.
     * 53-107, Ed. J. Leczszynski, World Scientifc, 1996.
     *
     * @param l int The number of (d/dx) operations.
     * @param m int The number of (d/dy) operations.
     * @param n int The number of (d/dz) operations.
     * @param j int j = 0 is the Tlmn tensor, j > 0 is an intermediate.
     * @param r double[] The {x,y,z} coordinates.
     * @param T000 double[] Initial auxillary tensor elements from Eq. (40).
     * @return double The requested Tensor element (intermediate if j > 0).
     * @since 1.0
     */
    private static double Tlmnj(final int l, final int m, final int n,
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

}
