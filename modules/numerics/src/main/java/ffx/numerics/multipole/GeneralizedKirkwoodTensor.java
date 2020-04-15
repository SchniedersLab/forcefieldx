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
package ffx.numerics.multipole;

import static java.lang.System.arraycopy;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.numerics.math.ScalarMath.doubleFactorial;

/**
 * The GeneralizedKirkwoodTensor class contains utilities for generated
 * Generalized Kirkwood interaction tensors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GeneralizedKirkwoodTensor {

    /**
     * Compute the potential auxiliary function for a multipole of order n
     *
     * @param n  Multipole order.
     * @param r2 Separation distance squared.
     * @param ai Born radius on atom i.
     * @param aj Born radius on atom j.
     * @param gc Generalized Kirkwood constant.
     * @param Eh Homogeneous dielectric constant.
     * @param Es Solvent dielectric constant.
     * @return The potential auxiliary function for a multipole of order n.
     */
    public static double an0(int n, double r2, double ai, double aj, double gc, double Eh, double Es) {
        var expTerm = exp(-r2 / (gc * ai * aj));
        var f = sqrt(r2 + ai * aj * expTerm);
        return cn(n, Eh, Es) * pow(-1, n) * doubleFactorial(2 * n - 1) / pow(f, 2 * n + 1);
    }

    /**
     * Compute the mth potential gradient auxiliary function for a multipole of order n.
     *
     * @param n  Multipole order.
     * @param m  Mth potential gradient auxiliary function.
     * @param r2 Separation distance squared.
     * @param ai Born radius on atom i.
     * @param aj Born radius on atom j.
     * @param gc Generalized Kirkwood constant.
     * @param Eh Homogeneous dielectric constant.
     * @param Es Solvent dielectric constant.
     * @return Returns the mth potential gradient auxiliary function for a multipole of order n.
     */
    public static double anm(int n, int m, double r2, double ai, double aj, double gc, double Eh, double Es) {
        if (m == 0) {
            return an0(n, r2, ai, aj, gc, Eh, Es);
        }
        var ret = 0.0;
        var coef = anmc(n);
        for (int i = 1; i <= m; i++) {
            ret += coef[i - 1] * fn(i, r2, ai, aj, gc) * anm(n + 1, m - i, r2, ai, aj, gc, Eh, Es);
        }
        return ret;
    }

    /**
     * Return coefficients needed when taking derivatives of auxiliary functions.
     *
     * @param n Multipole order.
     * @return Returns coefficients needed when taking derivatives of auxiliary functions.
     */
    public static double[] anmc(int n) {
        double[] ret = new double[n];
        double[] prev = new double[n];
        ret[0] = 1.0;
        if (n == 1) {
            return ret;
        }
        ret[1] = 1.0;
        if (n == 2) {
            return ret;
        }
        prev[0] = 1.0;
        prev[1] = 1.0;
        for (int i = 3; i <= n; i++) {
            for (int j = 2; j <= i - 1; j++) {
                ret[j - 1] = prev[j - 2] + prev[j - 1];
            }
            ret[i - 1] = 0.1e1;
            arraycopy(ret, 0, prev, 0, i);
        }
        return ret;
    }

    /**
     * Returns nth value of the function b, which are chain rule terms from
     * differentiating zeroth order auxiliary functions (an0) with respect to Ai or Aj.
     *
     * @param n  Multipole order.
     * @param r2 Separation distance squared.
     * @param Ai Born radius on atom i.
     * @param Aj Born radius on atom j.
     * @param gc Generalized Kirkwood constant.
     * @return Returns the nth value of the function f.
     */
    public static double bn(int n, double r2, double Ai, double Aj, double gc) {
        var gcAiAj = gc * Ai * Aj;
        var ratio = -r2 / gcAiAj;
        var expTerm = exp(ratio);
        if (n == 0) {
            return 0.5 * expTerm * (1.0 - ratio);
        }
        if (n == 1) {
            return -r2 * expTerm / (gcAiAj * gcAiAj);
        }
        var b2 = 2.0 * expTerm / (gcAiAj * gcAiAj) * (-ratio - 1.0);
        var br = 2.0 / (gcAiAj * Ai * Aj);
        var f2 = 2.0 / (gc * gcAiAj) * expTerm;
        var fr = -2.0 / (gcAiAj);
        return (n - 2) * pow(fr, n - 3) * br * f2 + pow(fr, n - 2) * b2;
    }

    /**
     * Compute the derivative with respect to a Born radius of the mth potential gradient auxiliary function for a multipole of order n.
     *
     * @param n  Multipole order.
     * @param m  Mth potential gradient auxiliary function.
     * @param r2 Separation distance squared.
     * @param ai Born radius on atom i.
     * @param aj Born radius on atom j.
     * @param gc Generalized Kirkwood constant.
     * @param Eh Homogeneous dielectric constant.
     * @param Es Solvent dielectric constant.
     * @return Returns the derivative with respect to a Born radius of the mth potential gradient auxiliary function for a multipole of order n.
     */
    public static double bnm(int n, int m, double r2, double ai, double aj, double gc, double Eh, double Es) {
        if (m == 0) {
            return bn(0, r2, ai, aj, gc) * an0(n + 1, r2, ai, aj, gc, Eh, Es);
        }
        var ret = 0;
        var coef = anmc(n);
        for (int i = 1; i <= m; i++) {
            ret += coef[i - 1] * bn(i, r2, ai, aj, gc) * anm(n + 1, m - i, r2, ai, aj, gc, Eh, Es);
            ret += coef[i - 1] * fn(i, r2, ai, aj, gc) * bnm(n + 1, m - i, r2, ai, aj, gc, Eh, Es);
        }
        return ret;
    }

    /**
     * Compute the Kirkwood dieletric function for a multipole of order n.
     *
     * @param n  Multipole order.
     * @param Eh Homogeneous dieletric.
     * @param Es Solvent dieletric.
     * @return Returns (n+1)*(Eh-Es)/((n+1)*Es + n*Eh))
     */
    public static double cn(int n, double Eh, double Es) {
        return (n + 1) * (Eh - Es) / ((n + 1) * Es + n * Eh);
    }

    /**
     * Returns nth value of the function f, which are chain rule terms from
     * differentiating zeroth order auxiliary functions (an0) with respect to x, y or z.
     *
     * @param n  Multipole order.
     * @param r2 Separation distance squared.
     * @param Ai Born radius on atom i.
     * @param Aj Born radius on atom j.
     * @param gc Generalized Kirkwood constant.
     * @return Returns the nth value of the function f.
     */
    public static double fn(int n, double r2, double Ai, double Aj, double gc) {
        var gcAiAj = gc * Ai * Aj;
        var ratio = -r2 / gcAiAj;
        var expTerm = exp(ratio);
        if (n == 0) {
            return sqrt(r2 + Ai * Aj * expTerm);
        }
        if (n == 1) {
            return 1.0 - expTerm / gc;
        }
        var f2 = 2.0 * expTerm / (gc * gcAiAj);
        var fr = -2.0 / gcAiAj;
        return pow(fr, n - 2) * f2;
    }
}
