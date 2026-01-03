// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
// ******************************************************************************
package ffx.xray.scatter;

import ffx.crystal.HKL;
import ffx.potential.bonded.Atom;
import ffx.xray.refine.RefinementMode;

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.DoubleMath.length2;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.numerics.math.MatrixMath.mat3Determinant;
import static ffx.numerics.math.MatrixMath.mat3Inverse;
import static ffx.numerics.math.MatrixMath.mat3Mat3Multiply;
import static ffx.numerics.math.MatrixMath.vec3Mat3;
import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.numerics.math.ScalarMath.quadForm;
import static ffx.numerics.math.ScalarMath.u2b;
import static ffx.xray.scatter.XRayScatteringParameters.getFormFactor;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * This implementation uses either coefficients from Su and Coppens or 3 coefficient parameters derived
 * from CCTBX.
 *
 * @author Timothy D. Fenn<br>
 * @see <a href="http://www.iucr.org/resources/commissions/crystallographic-computing/newsletters/3"
 * target="_blank"> R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams. Newsletter of the IUCr
 * Commission on Crystallographic Computing. (2004). 3, 22-31.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0907444909022707" target="_blank"> M. J. Schnieders, T.
 * D. Fenn, V. S. Pande and A. T. Brunger, Acta Cryst. (2009). D65 952-965.</a>
 * @since 1.0
 */
public final class XRayFormFactor implements FormFactor {

    private static final Logger logger = Logger.getLogger(XRayFormFactor.class.getName());

    private final Atom atom;
    private final double[] xyz = new double[3];
    private final double[] dxyz = new double[3];
    private final double[] resv = new double[3];
    private final double[] a = new double[6];
    private final double[] b = new double[6];
    private final double[] ainv = new double[6];
    private final double[] binv = new double[6];
    private final double[] gradu = new double[6];
    private final double[][][] u = new double[6][3][3];
    private final double[][][] uInv = new double[6][3][3];
    private final double[][] jmat = new double[3][3];
    private final int nGaussians;
    private boolean hasAnisou;
    private double[] anisou = null;
    private double uAdd;
    private double occupancy;

    /**
     * Constructor for XRayFormFactor.
     *
     * @param atom  a {@link ffx.potential.bonded.Atom} object.
     * @param use3G a boolean.
     * @param badd  a double.
     */
    public XRayFormFactor(Atom atom, boolean use3G, double badd) {
        this(atom, use3G, badd, atom.getXYZ(null));
    }

    /**
     * Constructor for XRayFormFactor.
     *
     * @param atom  a {@link ffx.potential.bonded.Atom} object.
     * @param use3G a boolean.
     * @param badd  a double.
     * @param xyz   an array of double.
     */
    public XRayFormFactor(Atom atom, boolean use3G, double badd, double[] xyz) {
        this.atom = atom;
        this.uAdd = b2u(badd);

        // Assign scattering factors.
        int charge = 0;
        if (atom.getMultipoleType() != null) {
            charge = (int) atom.getMultipoleType().getCharge();
        }

        XRayScatteringParameters parameters = getFormFactor(atom.getAtomicNumber(), charge, use3G);
        double[][] formFactor = parameters.formFactor();

        int i;
        for (i = 0; i < formFactor[1].length; i++) {
            if (formFactor[1][i] < 0.01) {
                break;
            }
            a[i] = formFactor[1][i];
            b[i] = formFactor[2][i];
        }
        nGaussians = i;
        assert (nGaussians > 0);
        occupancy = atom.getOccupancy();

        if (occupancy <= 0.0 && logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, " Zero occupancy for atom: {0}", atom.toString());
        }

        update(xyz, uAdd);
    }

    /**
     * f
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a double.
     */
    public double f(HKL hkl) {
        return fN(hkl, nGaussians);
    }

    /**
     * f_n
     *
     * @param hkl        a {@link ffx.crystal.HKL} object.
     * @param nGaussians a int.
     * @return a double.
     */
    public double fN(HKL hkl, int nGaussians) {
        double sum = 0.0;

        for (int i = 0; i < nGaussians; i++) {
            sum += a[i] * exp(-twoPI2 * hkl.quadForm(u[i]));
        }
        return occupancy * sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double[] xyz) {
        return rhoN(f, lambda, xyz, nGaussians);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementMode) {
        rhoGradN(xyz, nGaussians, dfc, refinementMode);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double[] xyz) {
        update(xyz, u2b(uAdd));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double[] xyz, double bAdd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
        uAdd = b2u(bAdd);
        occupancy = atom.getOccupancy();
        double bIso = atom.getTempFactor();

        // check occ is valid
        if (occupancy < 0.0) {
            logger.warning(format(" %s negative occupancy reset to 0.0.", atom));
            occupancy = 0.0;
            atom.setOccupancy(0.0);
        }

        // check if anisou changed
        if (atom.getAnisou(null) == null) {
            if (anisou == null) {
                anisou = new double[6];
            }
            hasAnisou = false;
        } else {
            hasAnisou = true;
        }

        if (hasAnisou) {
            // Check that the ANISOU is valid.
            anisou = atom.getAnisou(null);
            double det = mat3Determinant(anisou);

            if (det <= 1e-14) {
                String message = format(" %s non-positive definite ANISOU\n  Reset to isotropic B: %6.2f", atom, bIso);
                logger.warning(message);

                double uIso = b2u(bIso);
                anisou[0] = uIso;
                anisou[1] = uIso;
                anisou[2] = uIso;
                anisou[3] = 0.0;
                anisou[4] = 0.0;
                anisou[5] = 0.0;
                atom.setAnisou(anisou);
            }
        } else {
            if (bIso < 0.0) {
                StringBuilder sb = new StringBuilder();
                sb.append(" Negative B factor for atom: ").append(atom);
                sb.append("\n Resetting B to 5.0\n");
                logger.warning(sb.toString());
                bIso = 5.0;
                atom.setTempFactor(5.0);
            }
            double uIso = b2u(bIso);
            anisou[0] = uIso;
            anisou[1] = uIso;
            anisou[2] = uIso;
            anisou[3] = 0.0;
            anisou[4] = 0.0;
            anisou[5] = 0.0;
        }

        for (int i = 0; i < nGaussians; i++) {
            double uIso = b2u(b[i]);
            u[i][0][0] = anisou[0] + uIso + uAdd;
            u[i][1][1] = anisou[1] + uIso + uAdd;
            u[i][2][2] = anisou[2] + uIso + uAdd;
            u[i][0][1] = anisou[3];
            u[i][1][0] = anisou[3];
            u[i][0][2] = anisou[4];
            u[i][2][0] = anisou[4];
            u[i][1][2] = anisou[5];
            u[i][2][1] = anisou[5];
            mat3Inverse(u[i], uInv[i]);
            double det = mat3Determinant(u[i]);
            ainv[i] = a[i] / sqrt(det);
            det = mat3Determinant(uInv[i]);
            binv[i] = pow(det, oneThird);
        }
    }

    /**
     * rho_n
     *
     * @param f          a double.
     * @param lambda     a double.
     * @param xyz        an array of double.
     * @param nGaussians a int.
     * @return a double.
     */
    private double rhoN(double f, double lambda, double[] xyz, int nGaussians) {
        assert (nGaussians > 0 && nGaussians <= this.nGaussians);
        sub(this.xyz, xyz, xyz);

        // Compare r^2 to form factor width^2 to avoid expensive sqrt.
        if (length2(xyz) > atom.getFormFactorWidth2()) {
            return f;
        }

        double sum = 0.0;
        for (int i = 0; i < nGaussians; i++) {
            sum += ainv[i] * exp(-0.5 * quadForm(xyz, uInv[i]));
        }
        return f + (lambda * occupancy * inverseTwoPI32 * sum);
    }

    /**
     * rho_grad_n
     *
     * @param xyz            an array of double.
     * @param nGaussians     a int.
     * @param dfc            a double.
     * @param refinementMode a {@link RefinementMode} object.
     */
    private void rhoGradN(double[] xyz, int nGaussians, double dfc, RefinementMode refinementMode) {
        assert (nGaussians > 0 && nGaussians <= this.nGaussians);
        sub(this.xyz, xyz, dxyz);
        double r2 = length2(dxyz);

        // Compare r^2 to form factor width^2 to avoid expensive sqrt.
        if (r2 > atom.getFormFactorWidth2()) {
            return;
        }

        boolean refineXYZ = refinementMode.includesCoordinates();
        boolean refineBFactor = refinementMode.includesBFactors();
        boolean refineOccupancy = refinementMode.includesOccupancies();

        for (int i = 0; i < nGaussians; i++) {
            double aex = ainv[i] * exp(-0.5 * quadForm(dxyz, uInv[i]));

            if (refineXYZ) {
                vec3Mat3(dxyz, uInv[i], resv);
                double scale = aex * dfc * occupancy * inverseTwoPI32;
                double dex = -scale * resv[0];
                double dey = -scale * resv[1];
                double dez = -scale * resv[2];
                atom.addToXYZGradient(dex, dey, dez);
            }

            if (refineOccupancy) {
                atom.addToOccupancyGradient(dfc * inverseTwoPI32 * aex);
            }

            if (refineBFactor) {
                double scale = 0.5 * aex * dfc * occupancy * inverseTwoPI32;
                double deb = scale * (r2 * binv[i] * binv[i] - 3.0 * binv[i]);
                atom.addToTempFactorGradient(b2u(deb));
                if (hasAnisou) {
                    // U11
                    mat3Mat3Multiply(uInv[i], dUdU11, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[0] = scale * (quadForm(dxyz, jmat) - uInv[i][0][0]);
                    // U22
                    mat3Mat3Multiply(uInv[i], dUdU22, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[1] = scale * (quadForm(dxyz, jmat) - uInv[i][1][1]);
                    // U33
                    mat3Mat3Multiply(uInv[i], dUdU33, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[2] = scale * (quadForm(dxyz, jmat) - uInv[i][2][2]);
                    // U12
                    mat3Mat3Multiply(uInv[i], dUdU12, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[3] = scale * (quadForm(dxyz, jmat) - uInv[i][0][1] * 2.0);
                    // U13
                    mat3Mat3Multiply(uInv[i], dUdU13, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[4] = scale * (quadForm(dxyz, jmat) - uInv[i][0][2] * 2.0);
                    // U23
                    mat3Mat3Multiply(uInv[i], dUdU23, jmat);
                    mat3Mat3Multiply(jmat, uInv[i], jmat);
                    gradu[5] = scale * (quadForm(dxyz, jmat) - uInv[i][1][2] * 2.0);
                    // Store the anisotropic B-factor gradient.
                    atom.addToAnisouGradient(gradu);
                }
            }
        }
    }
}
