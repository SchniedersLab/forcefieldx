//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.numerics.special;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Implementation of the modified Bessel function of the first kind using
 * Chebyshev polynomials.
 * <p>
 * Adapted from CERN's cern.jet.math.tdouble package as included with
 * ParallelColt: http://sourceforge.net/projects/parallelcolt
 *
 * @author Stephen LuCore
 */
public class ModifiedBessel {

    /**
     * Compute log(i0(x)).
     *
     * @param x input parameter
     * @return the log(i0(x))
     */
    public static double lnI0(double x) {
        return log(i0(x));
    }

    /**
     * Compute the ration of i1(x) to i0(x).
     *
     * @param x input parameter
     * @return i1(x) / i0(x)
     */
    public static double i1OverI0(double x) {
        return (i1(x) / i0(x));
    }

    /**
     * Modified zero-order Bessel function. The function is defined as i0(x) =
     * J0( ix ). The range is partitioned into the two intervals [0,8] and (8,
     * infinity). Chebyshev polynomial expansions are employed in each interval.
     *
     * @param x input parameter
     * @return i0(x)
     */
    public static double i0(double x) {
        double y;
        if (x < 0) {
            x = -x;
        }
        if (x <= 8.0) {
            y = (x * 0.5) - 2.0;
            return (eToThe(x) * chbevl(y, A_i0, 30));
        }
        double ix = 1.0 / x;
        return (eToThe(x) * chbevl(32.0 * ix - 2.0, B_i0, 25) * sqrt(ix));
    }

    /**
     * Modified 1st-order ModifiedBessel function. The function is defined as
     * i1(x) = -i j1( ix ). The range is partitioned into the two intervals
     * [0,8] and (8, infinity). Chebyshev polynomial expansions are employed in
     * each interval.
     *
     * @param x input parameter
     * @return i1(x)
     */
    public static double i1(double x) {
        double y, z;

        z = abs(x);
        if (z <= 8.0) {
            y = (z * 0.5) - 2.0;
            z = chbevl(y, A_i1, 29) * z * eToThe(z);
        } else {
            double iz = 1.0 / z;
            z = eToThe(z) * chbevl(32.0 * iz - 2.0, B_i1, 25) * sqrt(iz);
        }
        if (x < 0.0) {
            z = -z;
        }
        return (z);
    }

    /**
     * Returns Double.MAX_VALUE in place of Double.POSITIVE_INFINITY; Returns
     * Double.MIN_VALUE in place of Double.NEGATIVE_INFINITY
     *
     * @param x input parameter
     * @return exp(x)
     */
    private static double eToThe(double x) {
        double res = exp(x);
        if (res == Double.POSITIVE_INFINITY) {
            return Double.MAX_VALUE;
        } else if (res == Double.NEGATIVE_INFINITY) {
            return Double.MIN_VALUE;
        }
        return res;
    }

    /**
     * Chebyshev coefficients for exp(-x) i0(x) in the interval [0,8].
     * lim(x->0){ exp(-x) i0(x) } = 1.
     */
    private static final double[] A_i0 = {-4.41534164647933937950E-18, 3.33079451882223809783E-17,
            -2.43127984654795469359E-16, 1.71539128555513303061E-15, -1.16853328779934516808E-14,
            7.67618549860493561688E-14, -4.85644678311192946090E-13, 2.95505266312963983461E-12,
            -1.72682629144155570723E-11, 9.67580903537323691224E-11, -5.18979560163526290666E-10,
            2.65982372468238665035E-9, -1.30002500998624804212E-8, 6.04699502254191894932E-8,
            -2.67079385394061173391E-7, 1.11738753912010371815E-6, -4.41673835845875056359E-6,
            1.64484480707288970893E-5, -5.75419501008210370398E-5, 1.88502885095841655729E-4,
            -5.76375574538582365885E-4, 1.63947561694133579842E-3, -4.32430999505057594430E-3,
            1.05464603945949983183E-2, -2.37374148058994688156E-2, 4.93052842396707084878E-2,
            -9.49010970480476444210E-2, 1.71620901522208775349E-1, -3.04682672343198398683E-1,
            6.76795274409476084995E-1};

    /**
     * Chebyshev coefficients for exp(-x) sqrt(x) i0(x) in the inverted interval
     * [8,infinity]. lim(x->inf){ exp(-x) sqrt(x) i0(x) } = 1/sqrt(2pi).
     */
    private static final double[] B_i0 = {-7.23318048787475395456E-18, -4.83050448594418207126E-18,
            4.46562142029675999901E-17, 3.46122286769746109310E-17, -2.82762398051658348494E-16,
            -3.42548561967721913462E-16, 1.77256013305652638360E-15, 3.81168066935262242075E-15,
            -9.55484669882830764870E-15, -4.15056934728722208663E-14, 1.54008621752140982691E-14,
            3.85277838274214270114E-13, 7.18012445138366623367E-13, -1.79417853150680611778E-12,
            -1.32158118404477131188E-11, -3.14991652796324136454E-11, 1.18891471078464383424E-11,
            4.94060238822496958910E-10, 3.39623202570838634515E-9, 2.26666899049817806459E-8,
            2.04891858946906374183E-7, 2.89137052083475648297E-6, 6.88975834691682398426E-5, 3.36911647825569408990E-3,
            8.04490411014108831608E-1};

    /**
     * Chebyshev coefficients for exp(-x) i1(x) / x in the interval [0,8].
     * lim(x->0){ exp(-x) i1(x) / x } = 1/2.
     */
    private static final double[] A_i1 = {2.77791411276104639959E-18, -2.11142121435816608115E-17,
            1.55363195773620046921E-16, -1.10559694773538630805E-15, 7.60068429473540693410E-15,
            -5.04218550472791168711E-14, 3.22379336594557470981E-13, -1.98397439776494371520E-12,
            1.17361862988909016308E-11, -6.66348972350202774223E-11, 3.62559028155211703701E-10,
            -1.88724975172282928790E-9, 9.38153738649577178388E-9, -4.44505912879632808065E-8,
            2.00329475355213526229E-7, -8.56872026469545474066E-7, 3.47025130813767847674E-6,
            -1.32731636560394358279E-5, 4.78156510755005422638E-5, -1.61760815825896745588E-4,
            5.12285956168575772895E-4, -1.51357245063125314899E-3, 4.15642294431288815669E-3,
            -1.05640848946261981558E-2, 2.47264490306265168283E-2, -5.29459812080949914269E-2,
            1.02643658689847095384E-1, -1.76416518357834055153E-1, 2.52587186443633654823E-1};

    /*
     * Chebyshev coefficients for exp(-x) sqrt(x) i1(x) in the inverted interval [8,infinity].
     * lim(x->inf){ exp(-x) sqrt(x) i1(x) } = 1/sqrt(2pi).
     */
    private static final double[] B_i1 = {7.51729631084210481353E-18, 4.41434832307170791151E-18,
            -4.65030536848935832153E-17, -3.20952592199342395980E-17, 2.96262899764595013876E-16,
            3.30820231092092828324E-16, -1.88035477551078244854E-15, -3.81440307243700780478E-15,
            1.04202769841288027642E-14, 4.27244001671195135429E-14, -2.10154184277266431302E-14,
            -4.08355111109219731823E-13, -7.19855177624590851209E-13, 2.03562854414708950722E-12,
            1.41258074366137813316E-11, 3.25260358301548823856E-11, -1.89749581235054123450E-11,
            -5.58974346219658380687E-10, -3.83538038596423702205E-9, -2.63146884688951950684E-8,
            -2.51223623787020892529E-7, -3.88256480887769039346E-6, -1.10588938762623716291E-4,
            -9.76109749136146840777E-3, 7.78576235018280120474E-1};

    /**
     * Evaluates Chebyshev polynomials (first kind) at x/2.
     * <p>
     * NOTE: Argument x must first be transformed to the interval (-1,1).
     * <p>
     * NOTE: Coefficients are in reverse; zero-order term is last.
     *
     * @param x    argument to the polynomial.
     * @param coef the coefficients of the polynomial.
     * @param N    the number of coefficients.
     * @return the result
     */
    private static double chbevl(double x, double[] coef, int N) {
        double b0, b1, b2;

        int p = 0;
        int i;

        b0 = coef[p++];
        b1 = 0.0;
        i = N - 1;

        do {
            b2 = b1;
            b1 = b0;
            b0 = x * b1 - b2 + coef[p++];
        } while (--i > 0);

        return (0.5 * (b0 - b2));
    }
}
