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
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Neutron scattering form factor.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public final class NeutronFormFactor implements FormFactor {

    private static final Logger logger = Logger.getLogger(NeutronFormFactor.class.getName());

    /**
     * The scattering atom.
     */
    private final Atom atom;
    /**
     * Current atomic coordinates of this atom.
     */
    private final double[] xyz = new double[3];
    /**
     * The location, relative to the atomic center, where the neutron form factor is being evaluated.
     */
    private final double[] dxyz = new double[3];
    /**
     * Neutron scattering length.
     */
    private final double[] a = new double[2];
    /**
     * Real portion of the neutron scattering length normalized by the determinant of the U matrix.
     */
    private double ainv;
    /**
     * Anisotropic U matrix with B_add included.
     */
    private final double[][] u = new double[3][3];
    /**
     * Inverse of U.
     */
    private final double[][] uInv = new double[3][3];
    /**
     * Determinant of the inverse U matrix to the 1/3 power.
     */
    private double binv;
    /**
     * Work array for atomic coordinate gradient.
     */
    private final double[] resv = new double[3];
    /**
     * If true, anisotropic B-factor is present.
     */
    private boolean hasAnisou;
    /**
     * Anisotropic B-factor.
     */
    private double[] uaniso = null;
    /**
     * Work array for anisotropic B-factor gradient.
     */
    private final double[] dU = new double[6];
    /**
     * Chain-rule term for anisotropic B-factor gradient.
     */
    private final double[][] jmat = new double[3][3];
    /**
     * U version of B_add.
     */
    private double uadd;
    /**
     * Occupancy of this atom.
     */
    private double occ;

    /**
     * Constructor for NeutronFormFactor.
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param badd a double.
     */
    public NeutronFormFactor(Atom atom, double badd) {
        this(atom, badd, atom.getXYZ(null));
    }

    /**
     * Constructor for NeutronFormFactor.
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param badd a double.
     * @param xyz  an array of double.
     */
    public NeutronFormFactor(Atom atom, double badd, double[] xyz) {
        this.atom = atom;
        this.uadd = b2u(badd);

        String key;
        if (atom.getAtomicNumber() == 1) {
            if (atom.isDeuterium()) {
                key = atom.getAtomicNumber() + "_2";
            } else {
                key = atom.getAtomicNumber() + "_1";
            }
        } else {
            key = "" + atom.getAtomicNumber();
        }

        NeutronScatteringParameters parameters = NeutronScatteringParameters.getFormFactor(key);
        double[] ffactor = parameters.formFactor();
        arraycopy(ffactor, 0, a, 0, ffactor.length);
        if (a[1] != 0.0) {
            logger.severe(" Complex neutron form factor method not supported");
        }

        occ = atom.getOccupancy();

        if (occ <= 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Zero occupancy: ").append(atom);
            logger.fine(sb.toString());
        }

        update(xyz, uadd);
    }

    /**
     * Returns the string representation of the NeutronFormFactor object.
     * The string includes details about the associated atom and its neutron
     * scattering magnitude.
     *
     * @return A formatted string containing the atom's information and its
     * neutron scattering magnitude.
     */
    @Override
    public String toString() {
        return atom.toString() + format(" Neutron mag: %9.6f", a[0]);
    }

    /**
     * f
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a double.
     */
    public double f(HKL hkl) {
        double sum = a[0] * exp(-twoPI2 * hkl.quadForm(u));
        return occ * sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double[] xyz) {
        // Compute the vector to the location of interest.
        sub(this.xyz, xyz, xyz);

        // Check if this is beyond the form factor cutoff.
        if (length2(xyz) > atom.getFormFactorWidth2()) {
            return f;
        }

        double sum = ainv * exp(-0.5 * quadForm(xyz, uInv));
        return f + (lambda * occ * inverseTwoPI32 * sum);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementMode) {
        // Compute the vector to the location of interest.
        sub(this.xyz, xyz, dxyz);

        // Check if this is beyond the form factor cutoff.
        double r2 = length2(dxyz);
        if (r2 > atom.getFormFactorWidth2()) {
            return;
        }

        double aex = ainv * exp(-0.5 * quadForm(dxyz, uInv));

        if (refinementMode.includesCoordinates()) {
            vec3Mat3(dxyz, uInv, resv);
            double scale = aex * dfc * occ * inverseTwoPI32;
            double dex = -scale * resv[0];
            double dey = -scale * resv[1];
            double dez = -scale * resv[2];
            atom.addToXYZGradient(dex, dey, dez);
        }

        if (refinementMode.includesOccupancies()) {
            atom.addToOccupancyGradient(dfc * inverseTwoPI32 * aex);
        }

        if (refinementMode.includesBFactors()) {
            double scale = 0.5 * aex * dfc * occ * inverseTwoPI32;
            double deb = scale * (r2 * binv * binv - 3.0 * binv);
            atom.addToTempFactorGradient(b2u(deb));
            if (hasAnisou) {
                // U11
                mat3Mat3Multiply(uInv, dUdU11, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[0] = scale * (quadForm(dxyz, jmat) - uInv[0][0]);
                // U22
                mat3Mat3Multiply(uInv, dUdU22, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[1] = scale * (quadForm(dxyz, jmat) - uInv[1][1]);
                // U33
                mat3Mat3Multiply(uInv, dUdU33, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[2] = scale * (quadForm(dxyz, jmat) - uInv[2][2]);
                // U12
                mat3Mat3Multiply(uInv, dUdU12, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[3] = scale * (quadForm(dxyz, jmat) - uInv[0][1] * 2.0);
                // U13
                mat3Mat3Multiply(uInv, dUdU13, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[4] = scale * (quadForm(dxyz, jmat) - uInv[0][2] * 2.0);
                // U23
                mat3Mat3Multiply(uInv, dUdU23, jmat);
                mat3Mat3Multiply(jmat, uInv, jmat);
                dU[5] = scale * (quadForm(dxyz, jmat) - uInv[1][2] * 2.0);
                // Store the anisotropic B-factor gradient.
                atom.addToAnisouGradient(dU);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double[] xyz) {
        update(xyz, u2b(uadd));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double[] xyz, double badd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
        double biso = atom.getTempFactor();
        uadd = b2u(badd);
        occ = atom.getOccupancy();

        // Check occ is valid.
        if (occ < 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Negative occupancy for atom: ").append(atom).append("\n");
            sb.append(" Resetting to 0.0\n");
            logger.warning(sb.toString());
            occ = 0.0;
            atom.setOccupancy(0.0);
        }

        // Check if anisou changed.
        if (atom.getAnisou(null) == null) {
            if (uaniso == null) {
                uaniso = new double[6];
            }
            hasAnisou = false;
        } else {
            hasAnisou = true;
        }

        if (hasAnisou) {
            // first check the ANISOU is valid
            uaniso = atom.getAnisou(null);
            double det = mat3Determinant(uaniso);
            if (det <= 1e-14) {
                StringBuilder sb = new StringBuilder();
                sb.append(" Non-positive definite ANISOU for atom: ").append(atom).append("\n");
                sb.append(" Resetting ANISOU based on isotropic B: (").append(biso).append(")\n");
                logger.warning(sb.toString());
                double val = b2u(biso);
                uaniso[0] = val;
                uaniso[1] = val;
                uaniso[2] = val;
                uaniso[3] = 0.0;
                uaniso[4] = 0.0;
                uaniso[5] = 0.0;
                atom.setAnisou(uaniso);
            }
        } else {
            if (biso < 0.0) {
                StringBuilder sb = new StringBuilder();
                sb.append(" Negative B factor for atom: ").append(atom).append("\n");
                sb.append(" Resetting B to 1.0\n");
                logger.warning(sb.toString());
                atom.setTempFactor(1.0);
                double val = b2u(1.0);
                uaniso[0] = val;
                uaniso[1] = val;
                uaniso[2] = val;
            } else {
                double val = b2u(biso);
                uaniso[0] = val;
                uaniso[1] = val;
                uaniso[2] = val;
            }
            uaniso[3] = 0.0;
            uaniso[4] = 0.0;
            uaniso[5] = 0.0;
        }

        u[0][0] = uaniso[0] + uadd;
        u[1][1] = uaniso[1] + uadd;
        u[2][2] = uaniso[2] + uadd;
        u[0][1] = uaniso[3];
        u[1][0] = uaniso[3];
        u[0][2] = uaniso[4];
        u[2][0] = uaniso[4];
        u[1][2] = uaniso[5];
        u[2][1] = uaniso[5];
        mat3Inverse(u, uInv);

        double det = mat3Determinant(u);
        ainv = a[0] / sqrt(det);
        det = mat3Determinant(uInv);
        binv = pow(det, oneThird);
    }
}
