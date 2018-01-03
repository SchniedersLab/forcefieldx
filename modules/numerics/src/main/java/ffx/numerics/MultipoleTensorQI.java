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
package ffx.numerics;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.norm;
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
public class MultipoleTensorQI extends MultipoleTensor {

    private static final Logger logger = Logger.getLogger(MultipoleTensorQI.class.getName());

    private double dEdZ = 0.0;
    private double d2EdZ2 = 0.0;

    public MultipoleTensorQI(OPERATOR operator, int order, double aewald) {
        super(operator, COORDINATES.QI, order, aewald);
    }

    @Override
    public double multipoleEnergy(double Fi[], double Ti[], double Tk[]) {
        // Compute the potential due to site I at site K.
        multipoleIField();

        // Dot the potential, field, field gradient with multipole K.
        double energy = dotMultipoleK();

        // Compute the torque on site K due to the field from site I.
        multipoleKTorque(Tk);

        // Compute the field at site I due to site K.
        multipoleKField();
        // Compute the torque on site I due to the field from site K.
        multipoleITorque(Ti);

        // Compute the force on site I F = {-dE/dx, -dE/dy, -dE/dz}.
        // Force on site K, Fk == -Fi.
        multipoleIdX();
        Fi[0] = -dotMultipoleK();
        multipoleIdY();
        Fi[1] = -dotMultipoleK();
        multipoleIdZ();
        Fi[2] = -dotMultipoleK();

        dEdZ = -Fi[2];
        if (order >= 6) {
            multipoleIdZ2();
            d2EdZ2 = dotMultipoleK();
        }

        /* Rotate the force and torques from the QI frame into the Global frame.
         * Has no effect if already in global frame. */
        qiToGlobal(Fi, Ti, Tk);
        return energy;
    }

    @Override
    public double polarizationEnergy(double scaleField, double scaleEnergy, double scaleMutual,
            double Fi[], double Ti[], double Tk[]) {
        // Find the potential, field, etc at k due to the induced dipole i.
        inducedIField();
        // Energy of multipole k in the field of induced dipole i.
        double energy = scaleEnergy * dotMultipoleK();

        /**
         * Get the induced-induced portion of the force.
         */
        Fi[0] = -0.5 * scaleMutual * (pxk * E200 + pyk * E110 + pzk * E101);
        Fi[1] = -0.5 * scaleMutual * (pxk * E110 + pyk * E020 + pzk * E011);
        Fi[2] = -0.5 * scaleMutual * (pxk * E101 + pyk * E011 + pzk * E002);

        // Find the potential, field, etc at i due to the induced dipole k.
        inducedKField();
        // Energy of multipole i in the field of induced dipole k.
        energy += scaleEnergy * dotMultipoleI();

        /**
         * Get the induced-induced portion of the force.
         */
        Fi[0] += 0.5 * scaleMutual * (pxi * E200 + pyi * E110 + pzi * E101);
        Fi[1] += 0.5 * scaleMutual * (pxi * E110 + pyi * E020 + pzi * E011);
        Fi[2] += 0.5 * scaleMutual * (pxi * E101 + pyi * E011 + pzi * E002);

        /**
         * Apply scale factors directly to induced dipole components for
         * efficiency and convenience in computing remaining force terms and
         * torques.
         */
        scaleInduced(scaleField, scaleEnergy);

        // Find the potential, field, etc at k due to (ind + indCR) at i.
        inducedIFieldForTorque();
        // inducedIFieldCR();
        // Torque on multipole k.
        multipoleKTorque(Tk);

        // Find the potential, field, etc at i due to (ind + indCR) at k.
        inducedKFieldForTorque();
        // inducedKFieldCR();
        // Torque on multipole i.
        multipoleITorque(Ti);

        // Forces
        inducedIdX();
        Fi[0] -= dotMultipoleK();
        inducedIdY();
        Fi[1] -= dotMultipoleK();
        inducedIdZ();
        Fi[2] -= dotMultipoleK();

        inducedKdX();
        Fi[0] -= dotMultipoleI();
        inducedKdY();
        Fi[1] -= dotMultipoleI();
        inducedKdZ();
        Fi[2] -= dotMultipoleI();

        qiToGlobal(Fi, Ti, Tk);
        return energy;
    }

    @Override
    public double getdEdZ() {
        return dEdZ;
    }

    @Override
    public double getd2EdZ2() {
        if (order < 6) {
//            logger.warning("Use sixth-order tensor to get second lambda derivative from QI.");
        }
        return d2EdZ2;
    }

    /**
     * @return Whether anything changed as a result.
     */
    @Override
    protected boolean setR(double r[], double lambdaFunction) {
        if (rprev != null && r[0] == rprev[0] && r[1] == rprev[1]
                && r[2] == rprev[2] && lambdaFunction == rprev[3]) {
            return true;
        }
        rprev[0] = r[0];
        rprev[1] = r[1];
        rprev[2] = r[2];
        rprev[3] = lambdaFunction;
        x = 0.0;
        y = 0.0;
        double zl = r[2] + lambdaFunction;
        r2 = r[0] * r[0] + r[1] * r[1] + zl * zl;
        z = sqrt(r2);
        R = z;
        setQIRotationMatrix(r[0], r[1], r[2] + lambdaFunction);
        return false;
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
     */
    @Override
    protected void noStorageRecursion(double r[], double tensor[]) {
        setR(r);
        r = new double[]{x, y, z};
        assert (r[0] == 0.0 && r[1] == 0.0);
        source(T000);
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
     * @see {@code MultipoleTensor#recursion(double[],double[])
     */
    @Override
    protected void recursion(final double r[], final double tensor[]) {
        setR(r);
        assert (x == 0.0 && y == 0.0);
        source(work);
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

    @Override
    protected double Tlmnj(final int l, final int m, final int n,
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
     * Hard coded computation of all Cartesian multipole tensors up to 4th
     * order, based on a quasi-internal frame, which is sufficient for
     * quadrupole-induced dipole forces.
     */
    @Override
    protected void order4() {
        source(work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
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
     * order, based on a quasi-internal frame, which is sufficient for
     * quadrupole-quadrupole forces.
     */
    @Override
    protected void order5() {
        source(work);
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

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 6th
     * order, based on a quasi-internal frame, which is sufficient for
     * quadrupole-quadrupole forces and orthogonal space sampling.
     */
    @Override
    protected void order6() {
        source(work);
        double term0000 = work[0];
        double term0001 = work[1];
        double term0002 = work[2];
        double term0003 = work[3];
        double term0004 = work[4];
        double term0005 = work[5];
        double term0006 = work[6];
        R000 = term0000;
        R200 = term0001;
        double term2001 = term0002;
        double term2002 = term0003;
        R400 = 3 * term2001;
        double term2003 = term0004;
        double term4001 = 3 * term2002;
        double term2004 = term0005;
        double term4002 = 3 * term2003;
        R600 = 5 * term4001;
        R020 = term0001;
        double term0201 = term0002;
        double term0202 = term0003;
        R040 = 3 * term0201;
        double term0203 = term0004;
        double term0401 = 3 * term0202;
        double term0204 = term0005;
        double term0402 = 3 * term0203;
        R060 = 5 * term0401;
        R220 = term2001;
        double term2201 = term2002;
        double term2202 = term2003;
        R240 = 3 * term2201;
        R420 = term4001;
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
        double term0015 = z * term0006;
        double term0024 = z * term0015 + term0005;
        double term0033 = z * term0024 + 2 * term0014;
        double term0042 = z * term0033 + 3 * term0023;
        double term0051 = z * term0042 + 4 * term0032;
        R006 = z * term0051 + 5 * term0041;
        R021 = z * term0201;
        double term0211 = z * term0202;
        R022 = z * term0211 + term0201;
        double term0212 = z * term0203;
        double term0221 = z * term0212 + term0202;
        R023 = z * term0221 + 2 * term0211;
        double term0213 = z * term0204;
        double term0222 = z * term0213 + term0203;
        double term0231 = z * term0222 + 2 * term0212;
        R024 = z * term0231 + 3 * term0221;
        R041 = z * term0401;
        double term0411 = z * term0402;
        R042 = z * term0411 + term0401;
        R201 = z * term2001;
        double term2011 = z * term2002;
        R202 = z * term2011 + term2001;
        double term2012 = z * term2003;
        double term2021 = z * term2012 + term2002;
        R203 = z * term2021 + 2 * term2011;
        double term2013 = z * term2004;
        double term2022 = z * term2013 + term2003;
        double term2031 = z * term2022 + 2 * term2012;
        R204 = z * term2031 + 3 * term2021;
        R221 = z * term2201;
        double term2211 = z * term2202;
        R222 = z * term2211 + term2201;
        R401 = z * term4001;
        double term4011 = z * term4002;
        R402 = z * term4011 + term4001;
    }

    @Override
    protected void multipoleIField() {
        double term000 = qi * R000;
        term000 -= dzi * R001;
        term000 += qxxi * R200;
        term000 += qyyi * R020;
        term000 += qzzi * R002;
        E000 = term000;
        double term100 = -dxi * R200;
        term100 += qxzi * R201;
        E100 = term100;
        double term010 = -dyi * R020;
        term010 += qyzi * R021;
        E010 = term010;
        double term001 = qi * R001;
        term001 -= dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        E001 = term001;
        double term200 = qi * R200;
        term200 -= dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        E200 = term200;
        double term020 = qi * R020;
        term020 -= dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        E020 = term020;
        double term002 = qi * R002;
        term002 -= dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        E002 = term002;
        double term110 = qxyi * R220;
        E110 = term110;
        double term101 = -dxi * R201;
        term101 += qxzi * R202;
        E101 = term101;
        double term011 = -dyi * R021;
        term011 += qyzi * R022;
        E011 = term011;
    }

    @Override
    protected void multipoleKField() {
        double term000 = qk * R000;
        term000 += dzk * R001;
        term000 += qxxk * R200;
        term000 += qyyk * R020;
        term000 += qzzk * R002;
        E000 = term000;
        double term100 = dxk * R200;
        term100 += qxzk * R201;
        E100 = term100;
        double term010 = dyk * R020;
        term010 += qyzk * R021;
        E010 = term010;
        double term001 = qk * R001;
        term001 += dzk * R002;
        term001 += qxxk * R201;
        term001 += qyyk * R021;
        term001 += qzzk * R003;
        E001 = term001;
        double term200 = qk * R200;
        term200 += dzk * R201;
        term200 += qxxk * R400;
        term200 += qyyk * R220;
        term200 += qzzk * R202;
        E200 = term200;
        double term020 = qk * R020;
        term020 += dzk * R021;
        term020 += qxxk * R220;
        term020 += qyyk * R040;
        term020 += qzzk * R022;
        E020 = term020;
        double term002 = qk * R002;
        term002 += dzk * R003;
        term002 += qxxk * R202;
        term002 += qyyk * R022;
        term002 += qzzk * R004;
        E002 = term002;
        double term110 = qxyk * R220;
        E110 = term110;
        double term101 = dxk * R201;
        term101 += qxzk * R202;
        E101 = term101;
        double term011 = dyk * R021;
        term011 += qyzk * R022;
        E011 = term011;
    }

    @Override
    protected void multipoleIdX() {
        double term100 = -dxi * R200;
        term100 += qxzi * R201;
        E000 = term100;
        double term200 = qi * R200;
        term200 -= dzi * R201;
        term200 += qxxi * R400;
        term200 += qyyi * R220;
        term200 += qzzi * R202;
        E100 = term200;
        double term110 = qxyi * R220;
        E010 = term110;
        double term101 = -dxi * R201;
        term101 += qxzi * R202;
        E001 = term101;
        double term300 = -dxi * R400;
        term300 += qxzi * R401;
        E200 = term300;
        double term120 = -dxi * R220;
        term120 += qxzi * R221;
        E020 = term120;
        double term102 = -dxi * R202;
        term102 += qxzi * R203;
        E002 = term102;
        double term210 = -dyi * R220;
        term210 += qyzi * R221;
        E110 = term210;
        double term201 = qi * R201;
        term201 -= dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        E101 = term201;
        double term111 = qxyi * R221;
        E011 = term111;
    }

    @Override
    protected void multipoleIdY() {
        double term010 = -dyi * R020;
        term010 += qyzi * R021;
        E000 = term010;
        double term110 = qxyi * R220;
        E100 = term110;
        double term020 = qi * R020;
        term020 -= dzi * R021;
        term020 += qxxi * R220;
        term020 += qyyi * R040;
        term020 += qzzi * R022;
        E010 = term020;
        double term011 = -dyi * R021;
        term011 += qyzi * R022;
        E001 = term011;
        double term210 = -dyi * R220;
        term210 += qyzi * R221;
        E200 = term210;
        double term030 = -dyi * R040;
        term030 += qyzi * R041;
        E020 = term030;
        double term012 = -dyi * R022;
        term012 += qyzi * R023;
        E002 = term012;
        double term120 = -dxi * R220;
        term120 += qxzi * R221;
        E110 = term120;
        double term111 = qxyi * R221;
        E101 = term111;
        double term021 = qi * R021;
        term021 -= dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        E011 = term021;
    }

    @Override
    protected void multipoleIdZ() {
        double term001 = qi * R001;
        term001 -= dzi * R002;
        term001 += qxxi * R201;
        term001 += qyyi * R021;
        term001 += qzzi * R003;
        E000 = term001;
        double term101 = -dxi * R201;
        term101 += qxzi * R202;
        E100 = term101;
        double term011 = -dyi * R021;
        term011 += qyzi * R022;
        E010 = term011;
        double term002 = qi * R002;
        term002 -= dzi * R003;
        term002 += qxxi * R202;
        term002 += qyyi * R022;
        term002 += qzzi * R004;
        E001 = term002;
        double term201 = qi * R201;
        term201 -= dzi * R202;
        term201 += qxxi * R401;
        term201 += qyyi * R221;
        term201 += qzzi * R203;
        E200 = term201;
        double term021 = qi * R021;
        term021 -= dzi * R022;
        term021 += qxxi * R221;
        term021 += qyyi * R041;
        term021 += qzzi * R023;
        E020 = term021;
        double term003 = qi * R003;
        term003 -= dzi * R004;
        term003 += qxxi * R203;
        term003 += qyyi * R023;
        term003 += qzzi * R005;
        E002 = term003;
        double term111 = qxyi * R221;
        E110 = term111;
        double term102 = -dxi * R202;
        term102 += qxzi * R203;
        E101 = term102;
        double term012 = -dyi * R022;
        term012 += qyzi * R023;
        E011 = term012;
    }

    /**
     * WARNING: Requires 6th-order tensors!
     */
    @Override
    protected void multipoleIdZ2() {
        double term001 = qi * R002;
        term001 -= dzi * R003;
        term001 += qxxi * R202;
        term001 += qyyi * R022;
        term001 += qzzi * R004;
        E000 = term001;
        double term101 = -dxi * R202;
        term101 += qxzi * R203;
        E100 = term101;
        double term011 = -dyi * R022;
        term011 += qyzi * R023;
        E010 = term011;
        double term002 = qi * R003;
        term002 -= dzi * R004;
        term002 += qxxi * R203;
        term002 += qyyi * R023;
        term002 += qzzi * R005;
        E001 = term002;
        double term201 = qi * R202;
        term201 -= dzi * R203;
        term201 += qxxi * R402;
        term201 += qyyi * R222;
        term201 += qzzi * R204;
        E200 = term201;
        double term021 = qi * R022;
        term021 -= dzi * R023;
        term021 += qxxi * R222;
        term021 += qyyi * R042;
        term021 += qzzi * R024;
        E020 = term021;
        double term003 = qi * R004;
        term003 -= dzi * R005;
        term003 += qxxi * R204;
        term003 += qyyi * R024;
        term003 += qzzi * R006;
        E002 = term003;
        double term111 = qxyi * R222;
        E110 = term111;
        double term102 = -dxi * R203;
        term102 += qxzi * R204;
        E101 = term102;
        double term012 = -dyi * R023;
        term012 += qyzi * R024;
        E011 = term012;
    }

    @Override
    protected void inducedIField() {
        E000 = -uzi * R001;
        E100 = -uxi * R200;
        E010 = -uyi * R020;
        E001 = -uzi * R002;
        E200 = -uzi * R201;
        E020 = -uzi * R021;
        E002 = -uzi * R003;
        E110 = 0.0;
        E101 = -uxi * R201;
        E011 = -uyi * R021;
    }

    @Override
    protected void inducedKField() {
        E000 = uzk * R001;
        E100 = uxk * R200;
        E010 = uyk * R020;
        E001 = uzk * R002;
        E200 = uzk * R201;
        E020 = uzk * R021;
        E002 = uzk * R003;
        E110 = 0.0;
        E101 = uxk * R201;
        E011 = uyk * R021;
    }

    @Override
    protected void inducedIFieldCR() {
        E000 = -pzi * R001;
        E100 = -pxi * R200;
        E010 = -pyi * R020;
        E001 = -pzi * R002;
        E200 = -pzi * R201;
        E020 = -pzi * R021;
        E002 = -pzi * R003;
        E110 = 0.0;
        E101 = -pxi * R201;
        E011 = -pyi * R021;
    }

    @Override
    protected void inducedKFieldCR() {
        E000 = pzk * R001;
        E100 = pxk * R200;
        E010 = pyk * R020;
        E001 = pzk * R002;
        E200 = pzk * R201;
        E020 = pzk * R021;
        E002 = pzk * R003;
        E110 = 0.0;
        E101 = pxk * R201;
        E011 = pyk * R021;
    }

    @Override
    protected void inducedIFieldForTorque() {
        E000 = -szi * R001;
        E100 = -sxi * R200;
        E010 = -syi * R020;
        E001 = -szi * R002;
        E200 = -szi * R201;
        E020 = -szi * R021;
        E002 = -szi * R003;
        E110 = 0.0;
        E101 = -sxi * R201;
        E011 = -syi * R021;
    }

    @Override
    protected void inducedKFieldForTorque() {
        E000 = szk * R001;
        E100 = sxk * R200;
        E010 = syk * R020;
        E001 = szk * R002;
        E200 = szk * R201;
        E020 = szk * R021;
        E002 = szk * R003;
        E110 = 0.0;
        E101 = sxk * R201;
        E011 = syk * R021;
    }

    @Override
    protected void inducedIdX() {
        E000 = -sxi * R200;
        E100 = -szi * R201;
        E010 = 0.0;
        E001 = -sxi * R201;
        E200 = -sxi * R400;
        E020 = -sxi * R220;
        E002 = -sxi * R202;
        E110 = -syi * R220;
        E101 = -szi * R202;
        E011 = 0.0;
    }

    @Override
    protected void inducedIdY() {
        E000 = -syi * R020;
        E100 = 0.0;
        E010 = -szi * R021;
        E001 = -syi * R021;
        E200 = -syi * R220;
        E020 = -syi * R040;
        E002 = -syi * R022;
        E110 = -sxi * R220;
        E101 = 0.0;
        E011 = -szi * R022;
    }

    @Override
    protected void inducedIdZ() {
        E000 = -szi * R002;
        E100 = -sxi * R201;
        E010 = -syi * R021;
        E001 = -szi * R003;
        E200 = -szi * R202;
        E020 = -szi * R022;
        E002 = -szi * R004;
        E110 = 0.0;
        E101 = -sxi * R202;
        E011 = -syi * R022;
    }

    @Override
    protected void inducedKdX() {
        E000 = sxk * R200;
        E100 = szk * R201;
        E010 = 0.0;
        E001 = sxk * R201;
        E200 = sxk * R400;
        E020 = sxk * R220;
        E002 = sxk * R202;
        E110 = syk * R220;
        E101 = szk * R202;
        E011 = 0.0;
    }

    @Override
    protected void inducedKdY() {
        E000 = syk * R020;
        E100 = 0.0;
        E010 = szk * R021;
        E001 = syk * R021;
        E200 = syk * R220;
        E020 = syk * R040;
        E002 = syk * R022;
        E110 = sxk * R220;
        E101 = 0.0;
        E011 = szk * R022;
    }

    @Override
    protected void inducedKdZ() {
        E000 = szk * R002;
        E100 = sxk * R201;
        E010 = syk * R021;
        E001 = szk * R003;
        E200 = szk * R202;
        E020 = szk * R022;
        E002 = szk * R004;
        E110 = 0.0;
        E101 = sxk * R202;
        E011 = syk * R022;
    }

    /**
     * Specific to QI; sets transform to rotate multipoles to (and from)
     * quasi-internal frame.
     */
    private void setQIRotationMatrix(double dx, double dy, double dz) {

        double zAxis[] = {dx, dy, dz};
        double xAxis[] = {dx, dy, dz};
        if (dy != 0.0 || dz != 0.0) {
            xAxis[0] += 1.0;
        } else {
            xAxis[1] += 1.0;
        }

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

    @Override
    protected void setDipoleI(double Ui[], double UiCR[]) {
        double dx = Ui[0];
        double dy = Ui[1];
        double dz = Ui[2];
        uxi = r00 * dx + r01 * dy + r02 * dz;
        uyi = r10 * dx + r11 * dy + r12 * dz;
        uzi = r20 * dx + r21 * dy + r22 * dz;
        dx = UiCR[0];
        dy = UiCR[1];
        dz = UiCR[2];
        pxi = r00 * dx + r01 * dy + r02 * dz;
        pyi = r10 * dx + r11 * dy + r12 * dz;
        pzi = r20 * dx + r21 * dy + r22 * dz;
    }

    @Override
    protected void setDipoleK(double Uk[], double UkCR[]) {
        double dx = Uk[0];
        double dy = Uk[1];
        double dz = Uk[2];
        uxk = r00 * dx + r01 * dy + r02 * dz;
        uyk = r10 * dx + r11 * dy + r12 * dz;
        uzk = r20 * dx + r21 * dy + r22 * dz;
        dx = UkCR[0];
        dy = UkCR[1];
        dz = UkCR[2];
        pxk = r00 * dx + r01 * dy + r02 * dz;
        pyk = r10 * dx + r11 * dy + r12 * dz;
        pzk = r20 * dx + r21 * dy + r22 * dz;
    }

    @Override
    protected void setMultipoleI(double Qi[]) {

        qi = Qi[0];

        double dx = Qi[1];
        double dy = Qi[2];
        double dz = Qi[3];

        dxi = r00 * dx + r01 * dy + r02 * dz;
        dyi = r10 * dx + r11 * dy + r12 * dz;
        dzi = r20 * dx + r21 * dy + r22 * dz;

        double qxx = Qi[4] * oneThird;
        double qyy = Qi[5] * oneThird;
        double qzz = Qi[6] * oneThird;
        double qxy = Qi[7] * oneThird;
        double qxz = Qi[8] * oneThird;
        double qyz = Qi[9] * oneThird;

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
        qxyi *= 2.0;

        // i=0, j=2
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxzi = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);
        qxzi *= 2.0;

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
        qyzi *= 2.0;

        // i=2, j=2
        // qij   r2k *  r20 * qkx + r21 * qky + r22 * qkz
        qzzi = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);

    }

    @Override
    protected void setMultipoleK(double Qk[]) {

        qk = Qk[0];

        double dx = Qk[1];
        double dy = Qk[2];
        double dz = Qk[3];

        dxk = r00 * dx + r01 * dy + r02 * dz;
        dyk = r10 * dx + r11 * dy + r12 * dz;
        dzk = r20 * dx + r21 * dy + r22 * dz;

        double qxx = Qk[4] * oneThird;
        double qyy = Qk[5] * oneThird;
        double qzz = Qk[6] * oneThird;
        double qxy = Qk[7] * oneThird;
        double qxz = Qk[8] * oneThird;
        double qyz = Qk[9] * oneThird;

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
        qxyk *= 2.0;

        // i=0, j=2
        // qij   rik *  rj0 * qkx + rj1 * qky + rj2 * qkz
        qxzk = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);
        qxzk *= 2.0;

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
        qyzk *= 2.0;

        // i=2, j=2
        // qij   r2k *  r20 * qkx + r21 * qky + r22 * qkz
        qzzk = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
                + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
                + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);
    }

    @Override
    protected void qiToGlobal(double v1[], double v2[],
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

    // Rotation Matrix from Global to QI.
    protected double r00, r01, r02;
    protected double r10, r11, r12;
    protected double r20, r21, r22;

    // Rotation Matrix from QI to Global.
    protected double ir00, ir01, ir02;
    protected double ir10, ir11, ir12;
    protected double ir20, ir21, ir22;
}
