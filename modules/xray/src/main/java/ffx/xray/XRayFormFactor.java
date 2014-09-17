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
package ffx.xray;

import java.util.HashMap;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.numerics.VectorMath.b2u;
import static ffx.numerics.VectorMath.determinant3;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.mat3Inverse;
import static ffx.numerics.VectorMath.mat3Mat3;
import static ffx.numerics.VectorMath.rsq;
import static ffx.numerics.VectorMath.scalarMat3Mat3;
import static ffx.numerics.VectorMath.u2b;
import static ffx.numerics.VectorMath.vec3Mat3;

/**
 * This implementation uses the coefficients from Su and Coppens and 3
 * coefficient parameters derived from CCTBX.
 *
 * @author Timothy D. Fenn<br>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767397004558" target="_blank">
 * Z. Su and P. Coppens, Acta Cryst. (1997). A53, 749-762</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S010876739800124X" target="_blank">
 * Z. Su and P. Coppens, Acta Cryst. (1998). A54, 357</a>
 *
 * @see <a href="http://harker.chem.buffalo.edu/group/groupindex.html"
 * target="_blank"> The Coppens lab website (Source data)</a>
 *
 * @see <a
 * href="http://www.iucr.org/resources/commissions/crystallographic-computing/newsletters/3"
 * target="_blank"> R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams.
 * Newsletter of the IUCr Commission on Crystallographic Computing. (2004). 3,
 * 22-31.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444909022707" target="_blank">
 * M. J. Schnieders, T. D. Fenn, V. S. Pande and A. T. Brunger, Acta Cryst.
 * (2009). D65 952-965.</a>
 *
 */
public final class XRayFormFactor implements FormFactor {

    private static final Logger logger = Logger.getLogger(ffx.xray.XRayFormFactor.class.getName());

    private final Atom atom;
    private final double xyz[] = new double[3];
    private final double dxyz[] = new double[3];
    private final double resv[] = new double[3];
    private final double a[] = new double[6];
    private final double b[] = new double[6];
    private final double ainv[] = new double[6];
    private final double binv[] = new double[6];
    private final double gradp[] = new double[6];
    private final double gradu[] = new double[6];
    private final double resm[][] = new double[3][3];
    private final double u[][][] = new double[6][3][3];
    private final double uinv[][][] = new double[6][3][3];
    private final double jmat[][][] = new double[6][3][3];
    private double anisou[] = null;
    protected final int ffIndex;
    private double bIso;
    private double uAdd;
    private double occupancy;
    private boolean hasAnisou;
    private final int nGaussians;

    /**
     * <p>
     * Constructor for XRayFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     */
    public XRayFormFactor(Atom atom) {
        this(atom, true, 0.0, atom.getXYZ());
    }

    /**
     * <p>
     * Constructor for XRayFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param use3G a boolean.
     */
    public XRayFormFactor(Atom atom, boolean use3G) {
        this(atom, use3G, 0.0, atom.getXYZ());
    }

    /**
     * <p>
     * Constructor for XRayFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param use3G a boolean.
     * @param badd a double.
     */
    public XRayFormFactor(Atom atom, boolean use3G, double badd) {
        this(atom, use3G, badd, atom.getXYZ());
    }

    /**
     * <p>
     * Constructor for XRayFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param use3G a boolean.
     * @param badd a double.
     * @param xyz an array of double.
     */
    public XRayFormFactor(Atom atom, boolean use3G, double badd, double xyz[]) {
        this.atom = atom;
        this.uAdd = b2u(badd);
        double ffactor[][];
        String key = "" + atom.getAtomicNumber();
        int charge = 0;
        if (atom.getMultipoleType() != null) {
            charge = (int) atom.getMultipoleType().charge;
        }

        int atomindex = atom.getFormFactorIndex();
        if (atomindex < 0) {
            // if it has a charge, first try to find Su&Coppens 6G params
            if (formfactors.containsKey(key + "_" + charge)) {
                ffactor = getFormFactor(key + "_" + charge);
            } else {
                // if not, use 3G params if requested
                if (use3G) {
                    // first look for charged form
                    if (formfactors.containsKey(key + "_" + charge + "_3g")) {
                        ffactor = getFormFactor(key + "_" + charge + "_3g");
                    } else {
                        // if this fails, we don't have the SFs
                        ffactor = getFormFactor(key + "_3g");
                    }
                } else {
                    ffactor = getFormFactor(key);
                }
            }
            ffIndex = (int) ffactor[0][0];
            atom.setFormFactorIndex(ffIndex);
        } else {
            ffIndex = atomindex;
            ffactor = ffactors[atomindex];
        }

        int i;
        for (i = 0; i < ffactor[1].length; i++) {
            if (ffactor[1][i] < 0.01) {
                break;
            }
            a[i] = ffactor[1][i];
            b[i] = ffactor[2][i];
        }
        nGaussians = i;
        assert (nGaussians > 0);
        occupancy = atom.getOccupancy();

        if (occupancy <= 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Zero occupancy for atom: " + atom.toString());
            logger.info(sb.toString());
        }

        update(xyz, uAdd);
    }

    /**
     * <p>
     * getFormFactorIndex</p>
     *
     * @param atom a {@link java.lang.String} object.
     * @return a int.
     */
    public static int getFormFactorIndex(String atom) {
        double ffactor[][] = getFormFactor(atom);
        if (ffactor != null) {
            return (int) ffactor[0][0];
        }
        return -1;

    }

    /**
     * <p>
     * getFormFactorA</p>
     *
     * @param atom a {@link java.lang.String} object.
     * @return an array of double.
     */
    public static double[] getFormFactorA(String atom) {
        double ffactor[][] = getFormFactor(atom);
        if (ffactor != null) {
            return ffactor[1];
        }
        return null;
    }

    /**
     * <p>
     * getFormFactorB</p>
     *
     * @param atom a {@link java.lang.String} object.
     * @return an array of double.
     */
    public static double[] getFormFactorB(String atom) {
        double ffactor[][] = getFormFactor(atom);
        if (ffactor != null) {
            return ffactor[2];
        }
        return null;
    }

    /**
     * <p>
     * getFormFactor</p>
     *
     * @param atom a {@link java.lang.String} object.
     * @return an array of double.
     */
    public static double[][] getFormFactor(String atom) {
        double ffactor[][] = null;
        if (formfactors.containsKey(atom)) {
            ffactor = (double[][]) formfactors.get(atom);
        } else {
            String message = " Form factor for atom: " + atom
                    + " not found!\n";
            logger.severe(message);
        }
        return ffactor;
    }

    /**
     * <p>
     * f</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a double.
     */
    public double f(HKL hkl) {
        return fN(hkl, nGaussians);
    }

    /**
     * <p>
     * f_n</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @param nGaussians a int.
     * @return a double.
     */
    public double fN(HKL hkl, int nGaussians) {
        double sum = 0.0;

        for (int i = 0; i < nGaussians; i++) {
            sum += a[i] * exp(-twopi2 * Crystal.quad_form(hkl, u[i]));
        }
        return occupancy * sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double xyz[]) {
        return rhoN(f, lambda, xyz, nGaussians);
    }

    /**
     * <p>
     * rho_n</p>
     *
     * @param f a double.
     * @param lambda a double.
     * @param xyz an array of double.
     * @param nGaussians a int.
     * @return a double.
     */
    public double rhoN(double f, double lambda, double xyz[], int nGaussians) {
        assert (nGaussians > 0 && nGaussians <= nGaussians);
        diff(this.xyz, xyz, xyz);

        /**
         * Compare r^2 to form factor width^2 to avoid expensive sqrt.
         */
        if (rsq(xyz) > atom.getFormFactorWidth2()) {
            return f;
        }

        double sum = 0.0;
        for (int i = 0; i < nGaussians; i++) {
            sum += ainv[i] * exp(-0.5 * Crystal.quad_form(xyz, uinv[i]));
        }
        return f + (lambda * occupancy * twopi32 * sum);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void rhoGrad(double xyz[], double dfc, RefinementMode refinementMode) {
        rhoGradN(xyz, nGaussians, dfc, refinementMode);
    }

    /**
     * <p>
     * rho_grad_n</p>
     *
     * @param xyz an array of double.
     * @param nGaussians a int.
     * @param dfc a double.
     * @param refinementMode a
     * {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public void rhoGradN(double xyz[], int nGaussians, double dfc, RefinementMode refinementMode) {
        assert (nGaussians > 0 && nGaussians <= nGaussians);
        diff(this.xyz, xyz, dxyz);
        double r2 = rsq(dxyz);

        /**
         * Compare r^2 to form factor width^2 to avoid expensive sqrt.
         */
        if (r2 > atom.getFormFactorWidth2()) {
            return;
        }

        fill(gradp, 0.0);
        fill(gradu, 0.0);
        double aex;
        boolean refinexyz = false;
        boolean refineb = false;
        boolean refineanisou = false;
        boolean refineocc = false;
        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refinexyz = true;
        }
        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineb = true;
            if (hasAnisou) {
                refineanisou = true;
            }
        }
        if (refinementMode == RefinementMode.OCCUPANCIES
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineocc = true;
        }

        for (int i = 0; i < nGaussians; i++) {
            aex = ainv[i] * exp(-0.5 * Crystal.quad_form(dxyz, uinv[i]));

            if (refinexyz) {
                vec3Mat3(dxyz, uinv[i], resv);
                gradp[0] += aex * dot(resv, vx);
                gradp[1] += aex * dot(resv, vy);
                gradp[2] += aex * dot(resv, vz);
            }

            if (refineocc) {
                gradp[3] += aex;
            }

            if (refineb) {
                gradp[4] += aex * 0.5 * (r2 * binv[i] * binv[i] - 3.0 * binv[i]);

                if (refineanisou) {
                    scalarMat3Mat3(-1.0, uinv[i], u11, resm);
                    mat3Mat3(resm, uinv[i], jmat[0]);
                    scalarMat3Mat3(-1.0, uinv[i], u22, resm);
                    mat3Mat3(resm, uinv[i], jmat[1]);
                    scalarMat3Mat3(-1.0, uinv[i], u33, resm);
                    mat3Mat3(resm, uinv[i], jmat[2]);
                    scalarMat3Mat3(-1.0, uinv[i], u12, resm);
                    mat3Mat3(resm, uinv[i], jmat[3]);
                    scalarMat3Mat3(-1.0, uinv[i], u13, resm);
                    mat3Mat3(resm, uinv[i], jmat[4]);
                    scalarMat3Mat3(-1.0, uinv[i], u23, resm);
                    mat3Mat3(resm, uinv[i], jmat[5]);

                    gradu[0] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[0]) - uinv[i][0][0]);
                    gradu[1] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[1]) - uinv[i][1][1]);
                    gradu[2] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[2]) - uinv[i][2][2]);
                    gradu[3] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[3]) - uinv[i][0][1] * 2.0);
                    gradu[4] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[4]) - uinv[i][0][2] * 2.0);
                    gradu[5] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[5]) - uinv[i][1][2] * 2.0);
                }
            }
        }

        //double rho = occ * twopi32 * gradp[3];
        // x, y, z
        if (refinexyz) {
            atom.addToXYZGradient(
                    dfc * occupancy * -twopi32 * gradp[0],
                    dfc * occupancy * -twopi32 * gradp[1],
                    dfc * occupancy * -twopi32 * gradp[2]);
        }

        // occ
        if (refineocc) {
            atom.addToOccupancyGradient(dfc * twopi32 * gradp[3]);
        }

        // Biso
        if (refineb) {
            atom.addToTempFactorGradient(dfc * b2u(occupancy * twopi32 * gradp[4]));
            // Uaniso
            if (hasAnisou) {
                for (int i = 0; i < 6; i++) {
                    gradu[i] = dfc * occupancy * twopi32 * gradu[i];
                }
                atom.addToAnisouGradient(gradu);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double xyz[]) {
        update(xyz, u2b(uAdd));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double xyz[], double bAdd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
        bIso = atom.getTempFactor();
        uAdd = b2u(bAdd);
        occupancy = atom.getOccupancy();

        // check occ is valid
        if (occupancy < 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Negative occupancy for atom: " + atom.toString());
            sb.append(" \nResetting to 0.0\n");
            logger.warning(sb.toString());
            occupancy = 0.0;
            atom.setOccupancy(0.0);
        }

        // check if anisou changed
        if (atom.getAnisou() == null) {
            if (anisou == null) {
                anisou = new double[6];
            }
            hasAnisou = false;
        } else {
            hasAnisou = true;
        }

        if (hasAnisou) {
            // first check the ANISOU is valid
            anisou = atom.getAnisou();
            double det = determinant3(anisou);

            if (det <= 1e-14) {
                StringBuilder sb = new StringBuilder();
                sb.append(" Non-positive definite ANISOU for atom: " + atom.toString());
                sb.append("\n Resetting ANISOU based on isotropic B: (" + bIso + ")\n");
                logger.warning(sb.toString());

                anisou[0] = anisou[1] = anisou[2] = b2u(bIso);
                anisou[3] = anisou[4] = anisou[5] = 0.0;
            }
        } else {
            if (bIso < 0.0) {
                StringBuilder sb = new StringBuilder();
                sb.append(" Negative B factor for atom: " + atom.toString());
                sb.append("\n Resetting B to 0.01\n");
                logger.warning(sb.toString());
                bIso = 0.01;
                atom.setTempFactor(0.01);
            }
            anisou[0] = anisou[1] = anisou[2] = b2u(bIso);
            anisou[3] = anisou[4] = anisou[5] = 0.0;
        }

        for (int i = 0; i < nGaussians; i++) {
            u[i][0][0] = anisou[0] + b2u(b[i]) + uAdd;
            u[i][1][1] = anisou[1] + b2u(b[i]) + uAdd;
            u[i][2][2] = anisou[2] + b2u(b[i]) + uAdd;
            u[i][0][1] = u[i][1][0] = anisou[3];
            u[i][0][2] = u[i][2][0] = anisou[4];
            u[i][1][2] = u[i][2][1] = anisou[5];

            mat3Inverse(u[i], uinv[i]);

            double det = determinant3(u[i]);
            ainv[i] = a[i] / sqrt(det);
            // b[i] = pow(det, 0.33333333333);
            det = determinant3(uinv[i]);
            binv[i] = pow(det, oneThird);
        }
    }

    private static final double twopi2 = 2.0 * PI * PI;
    private static final double twopi32 = pow(2.0 * PI, -1.5);
    private final static double oneThird = 1.0 / 3.0;
    private static final double vx[] = {1.0, 0.0, 0.0};
    private static final double vy[] = {0.0, 1.0, 0.0};
    private static final double vz[] = {0.0, 0.0, 1.0};
    private static final double u11[][] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u22[][] = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u33[][] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    private static final double u12[][] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    private static final double u13[][] = {{0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    private static final double u23[][] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
    private static final HashMap formfactors = new HashMap();

    private static final String[] atoms = {"H", "He", "Li", "Be", "B", "C", "N", "O",
        "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
        "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
        "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Li+", "Be2+", "Cv", "O-", "F-", "Na+", "Mg2+", "Al3+", "Si0", "Si4+",
        "Cl-", "K+", "Ca2+", "Sc3+", "Ti2+", "Ti3+", "Ti4+", "V2+", "V3+",
        "V5+", "Cr2+", "Cr3+", "Mn2+", "Mn3+", "Mn4+", "Fe2+", "Fe3+", "Co2+",
        "Co3+", "Ni2+", "Ni3+", "Cu+", "Cu2+", "Zn2+", "Ga3+", "Ge4+", "Br-",
        "Rb+", "Sr2+", "Y3+", "Zr4+", "Nb3+", "Nb5+", "Mo3+", "Mo6+", "Ru3+",
        "Ru4+", "Rh3+", "Rh4+", "Pd2+", "Pd4+", "Ag+", "Ag2+", "Cd2+", "In3+",
        "Sn2+", "Sn4+", "Sb3+", "Sb5+", "I-",
        "Be2+_3g", "Mg2+_3g", "Al3+_3g", "Si4+_3g", "Ca2+_3g", "Sc3+_3g",
        "Ti2+_3g", "Ti3+_3g", "Ti4+_3g", "V2+_3g", "V3+_3g", "V5+_3g", "Cr2+_3g",
        "Cr3+_3g", "Mn2+_3g", "Mn3+_3g", "Mn4+_3g", "Fe2+_3g", "Fe3+_3g", "Co2+_3g",
        "Co3+_3g", "Ni2+_3g", "Ni3+_3g", "Cu2+_3g", "Zn2+_3g", "Ga3+_3g", "Ge4+_3g",
        "Sr2+_3g", "Y3+_3g", "Zr4+_3g", "Nb3+_3g", "Nb5+_3g", "Mo3+_3g", "Mo6+_3g",
        "Ru3+_3g", "Ru4+_3g", "Rh3+_3g", "Rh4+_3g", "Pd2+_3g", "Pd4+_3g", "Ag2+_3g",
        "Cd2+_3g", "In3+_3g", "Sn2+_3g", "Sn4+_3g", "Sb3+_3g", "Sb5+_3g", "Hg2+_3g",
        "H_3g", "He_3g", "Li_3g", "Be_3g", "B_3g", "C_3g", "N_3g", "O_3g",
        "F_3g", "Ne_3g", "Na_3g", "Mg_3g", "Al_3g", "Si_3g", "P_3g", "S_3g",
        "Cl_3g", "Ar_3g", "K_3g", "Ca_3g", "Sc_3g", "Ti_3g", "V_3g", "Cr_3g",
        "Mn_3g", "Fe_3g", "Co_3g", "Ni_3g", "Cu_3g", "Zn_3g", "Ga_3g", "Ge_3g",
        "As_3g", "Se_3g", "Br_3g", "Kr_3g", "Rb_3g", "Sr_3g", "Y_3g", "Zr_3g",
        "Nb_3g", "Mo_3g", "Tc_3g", "Ru_3g", "Rh_3g", "Pd_3g", "Ag_3g", "Cd_3g",
        "In_3g", "Sn_3g", "Sb_3g", "Te_3g", "I_3g", "Xe_3g", "Hg_3g"};
    private static final String[] atomsi = {"1", "2", "3", "4", "5", "6", "7", "8",
        "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
        "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44",
        "45", "46", "47", "48", "49", "50", "51", "52", "53", "54",
        "3_1", "4_2", "6_5", "8_-1", "9_-1", "11_1", "12_2", "13_3", "14_0", "14_4",
        "17_-1", "19_1", "20_2", "21_3", "22_2", "22_3", "22_4", "23_2", "23_3",
        "23_5", "24_2", "24_3", "25_2", "25_3", "25_4", "26_2", "26_3", "27_2",
        "27_3", "28_2", "28_3", "29_1", "29_2", "30_2", "31_3", "32_4", "35_-1",
        "37_1", "38_2", "39_3", "40_4", "41_3", "41_5", "42_3", "42_6", "44_3",
        "44_4", "45_3", "45_4", "46_2", "46_4", "47_1", "47_2", "48_2", "49_3",
        "50_2", "50_4", "51_3", "51_5", "53_-1",
        "4_2_3g", "12_2_3g", "13_3_3g", "14_4_3g", "20_2_3g", "21_3_3g",
        "22_2_3g", "22_3_3g", "22_4_3g", "23_2_3g", "23_3_3g", "23_5_3g", "24_2_3g",
        "24_3_3g", "25_2_3g", "25_3_3g", "25_4_3g", "26_2_3g", "26_3_3g", "27_2_3g",
        "27_3_3g", "28_2_3g", "28_3_3g", "29_2_3g", "30_2_3g", "31_3_3g", "32_4_3g",
        "38_2_3g", "39_3_3g", "40_4_3g", "41_3_3g", "41_5_3g", "42_3_3g", "42_6_3g",
        "44_3_3g", "44_4_3g", "45_3_3g", "45_4_3g", "46_2_3g", "46_4_3g", "47_2_3g",
        "48_2_3g", "49_3_3g", "50_2_3g", "50_4_3g", "51_3_3g", "51_5_3g", "80_2_3g",
        "1_3g", "2_3g", "3_3g", "4_3g", "5_3g", "6_3g", "7_3g", "8_3g",
        "9_3g", "10_3g", "11_3g", "12_3g", "13_3g", "14_3g", "15_3g", "16_3g", "17_3g", "18_3g", "19_3g", "20_3g",
        "21_3g", "22_3g", "23_3g", "24_3g", "25_3g", "26_3g", "27_3g", "28_3g", "29_3g", "30_3g", "31_3g", "32_3g",
        "33_3g", "34_3g", "35_3g", "36_3g", "37_3g", "38_3g", "39_3g", "40_3g", "41_3g", "42_3g", "43_3g", "44_3g",
        "45_3g", "46_3g", "47_3g", "48_3g", "49_3g", "50_3g", "51_3g", "52_3g", "53_3g", "54_3g", "80_3g"};

    private static final double[][][] ffactors = {
        /**
         * Su and Coppens data.
         */
        {{0},
        {0.43028, 0.28537, 0.17134, 0.09451, 0.01725, 0.00114},
        {23.02312, 10.20138, 51.25444, 4.13511, 1.35427, 0.24269}},
        {{1},
        {0.69475, 0.62068, 0.38661, 0.15223, 0.12661, 0.01907},
        {5.83366, 12.87682, 2.53296, 28.16171, 0.97507, 0.25308}},
        {{2},
        {0.84645, 0.81146, 0.81096, 0.26115, 0.26055, 0.00930},
        {4.63253, 1.71862, 97.87364, 0.50620, 200.00000, 0.00010}},
        {{3},
        {1.59261, 1.12768, 0.70296, 0.53815, 0.03863, 0.00010},
        {43.67397, 1.86275, 0.54243, 103.44910, 0.00010, 0.34975}},
        {{4},
        {2.07418, 1.20577, 1.07592, 0.52023, 0.12280, 0.00010},
        {23.39543, 1.07672, 60.93249, 0.27132, 0.27192, 0.11361}},
        {{5},
        {2.09921, 1.80832, 1.26159, 0.56775, 0.26303, 0.00010},
        {13.18997, 30.37956, 0.69255, 0.16381, 68.42774, 0.44083}},
        {{6},
        {2.45424, 2.15782, 1.05782, 0.57557, 0.44959, 0.30480},
        {18.66694, 8.31271, 0.46989, 42.44646, 0.08747, 0.47126}},
        {{7},
        {2.34752, 1.83006, 1.61538, 1.52402, 0.41423, 0.26867},
        {9.69710, 18.59876, 5.19879, 0.32408, 39.79099, 0.01150}},
        {{8},
        {2.96981, 2.04536, 1.78123, 1.52086, 0.42253, 0.26008},
        {7.52365, 15.41441, 3.79721, 0.25209, 33.76478, 0.00488}},
        {{9},
        {3.56413, 2.72559, 1.67359, 1.58884, 0.25468, 0.19320},
        {7.30559, 3.34491, 15.93226, 0.13859, 0.69111, 35.26368}},
        {{10},
        {4.16491, 2.38097, 1.70484, 1.59622, 0.66291, 0.48971},
        {4.23096, 9.48502, 0.12559, 1.98358, 172.13327, 82.23091}},
        {{11},
        {3.90882, 2.62159, 1.69157, 1.52610, 1.47907, 0.77262},
        {3.06041, 6.12146, 0.10357, 58.65022, 1.56940, 125.49980}},
        {{12},
        {4.25474, 3.58301, 2.37351, 1.72366, 0.99400, 0.07031},
        {3.76670, 1.69151, 45.27810, 0.09238, 113.96978, 17.47922}},
        {{13},
        {4.94976, 3.25403, 2.84957, 1.66053, 1.22949, 0.05611},
        {2.70254, 34.45314, 1.24059, 0.07201, 84.53648, 56.34208}},
        {{14},
        {6.48197, 4.31666, 1.73759, 1.35793, 1.10559, 0.00010},
        {1.89537, 27.61455, 0.50991, 66.28296, 0.00010, 12.05652}},
        {{15},
        {6.90565, 5.24410, 1.54516, 1.42922, 0.87564, 0.00010},
        {1.46764, 22.31576, 56.06328, 0.25588, 0.00010, 26.96892}},
        {{16},
        {7.13381, 6.26972, 1.82658, 1.62579, 0.14431, 0.00010},
        {1.17455, 18.57626, 0.07869, 48.08203, 0.07871, 23.23894}},
        {{17},
        {7.28551, 7.24549, 1.74775, 1.72174, 0.00010, 0.00010},
        {15.63295, 0.95562, 0.04456, 41.07550, 0.00617, 20.09628}},
        {{18},
        {8.13161, 7.43972, 1.42159, 1.12030, 0.88342, 0.00010},
        {12.73675, 0.77443, 0.00010, 200.00000, 36.18711, 82.98380}},
        {{19},
        {8.62965, 7.38765, 1.63044, 1.37681, 0.97538, 0.00010},
        {10.45238, 0.66036, 87.06258, 0.00010, 181.27760, 28.57890}},
        {{20},
        {9.18894, 7.36727, 1.60214, 1.33655, 0.78386, 0.72047},
        {9.02948, 0.57364, 137.40503, 0.00010, 51.53615, 53.74395}},
        {{21},
        {9.75861, 7.35354, 1.46842, 1.40591, 1.28669, 0.72609},
        {7.86172, 0.50107, 32.75146, 90.95131, 0.00010, 149.02872}},
        {{22},
        {10.25443, 7.34699, 1.84039, 1.72148, 1.22611, 0.61000},
        {6.86177, 0.43939, 23.70259, 79.72053, 0.00010, 149.36488}},
        {{23},
        {10.67225, 4.62093, 3.33159, 2.72784, 1.45281, 1.19090},
        {6.12143, 0.39293, 20.15470, 0.39293, 92.01317, 0.00010}},
        {{24},
        {10.98576, 7.35617, 2.92091, 1.65707, 1.08018, 0.99906},
        {5.27951, 0.34199, 14.55791, 54.87900, 0.00010, 118.26511}},
        {{25},
        {11.18858, 7.37206, 3.55141, 1.68125, 1.20893, 0.99652},
        {4.64599, 0.30327, 12.07655, 44.15316, 104.11866, 0.00010}},
        {{26},
        {11.41624, 7.38902, 4.21351, 1.80189, 1.26103, 0.91710},
        {4.12258, 0.27069, 10.36636, 38.32442, 97.14970, 0.00010}},
        {{27},
        {11.76300, 7.39888, 4.85491, 1.98079, 1.14857, 0.85325},
        {3.69729, 0.24374, 9.30593, 36.58880, 96.02875, 0.00010}},
        {{28},
        {11.87211, 7.37491, 6.08548, 1.94337, 0.86475, 0.85837},
        {3.34773, 0.22522, 8.46165, 27.95010, 98.02165, 0.00012}},
        {{29},
        {12.53020, 6.57092, 5.84880, 2.07610, 1.65893, 1.31388},
        {3.05828, 0.14326, 7.58930, 28.50706, 0.38369, 82.22092}},
        {{30},
        {10.69865, 7.89127, 4.74778, 3.83120, 2.59218, 1.23712},
        {3.44787, 0.15426, 2.07387, 8.38441, 34.93356, 99.34732}},
        {{31},
        {9.56335, 7.86994, 7.64215, 3.31296, 2.13351, 1.47704},
        {2.21494, 0.14284, 3.86490, 32.69417, 8.94286, 82.15827}},
        {{32},
        {10.86205, 7.83248, 5.48862, 4.21250, 2.56904, 2.03413},
        {2.10046, 0.13209, 3.33631, 26.38254, 5.81992, 63.74567}},
        {{33},
        {12.63843, 7.77316, 5.80645, 4.44296, 1.82898, 1.50938},
        {1.97006, 0.12167, 3.57609, 28.84348, 15.15766, 64.03025}},
        {{34},
        {12.56835, 7.70669, 5.76243, 4.78093, 2.48412, 1.69674},
        {1.79826, 0.11204, 2.98848, 25.62856, 14.95420, 55.44329}},
        {{35},
        {13.32373, 7.64645, 5.71351, 4.95009, 2.80427, 1.56038},
        {1.67399, 0.10346, 17.43646, 2.62566, 42.87908, 19.80281}},
        {{36},
        {17.73932, 7.70415, 5.33484, 4.92829, 1.28671, 0.00010},
        {1.68298, 0.09944, 12.80739, 23.59343, 200.00000, 77.16806}},
        {{37},
        {11.77920, 9.53489, 7.57120, 6.03047, 2.02653, 1.05652},
        {1.52266, 13.50271, 0.08995, 1.52251, 162.86971, 53.07068}},
        {{38},
        {17.89478, 9.91124, 7.40424, 2.14475, 1.64266, 0.00010},
        {1.37779, 12.18084, 0.08009, 137.73235, 49.81442, 0.42187}},
        {{39},
        {18.00877, 10.47108, 7.22234, 2.43263, 1.86405, 0.00010},
        {1.25042, 11.25972, 0.07050, 49.09408, 131.67513, 1.76480}},
        {{40},
        {18.18722, 11.07349, 7.02786, 3.35224, 1.35250, 0.00606},
        {1.13993, 10.82683, 0.06116, 38.71734, 115.18009, 1.19550}},
        {{41},
        {18.36000, 6.75320, 6.25470, 5.52972, 3.76774, 1.33338},
        {1.03291, 0.05000, 10.10135, 10.12179, 34.16693, 104.10497}},
        {{42},
        {18.53113, 12.72135, 6.39681, 2.88811, 1.72002, 0.74148},
        {0.93112, 9.26800, 0.03703, 31.91681, 110.11821, 44.07274}},
        {{43},
        {18.82022, 13.49636, 6.01136, 3.54102, 1.19962, 0.93207},
        {0.84363, 8.84277, 0.02355, 27.02179, 44.09284, 113.68484}},
        {{44},
        {19.15093, 14.43898, 4.66972, 4.66263, 1.22522, 0.85125},
        {0.75936, 8.27523, 26.67965, 0.00694, 97.04210, 0.00695}},
        {{45},
        {19.32300, 15.30162, 5.26970, 5.12338, 0.98021, 0.00010},
        {0.69750, 7.93132, 0.00010, 23.54133, 60.82499, 1.28291}},
        {{46},
        {19.28330, 16.71519, 5.18450, 4.77793, 1.03807, 0.00010},
        {0.64519, 7.48785, 0.00010, 24.79225, 100.31405, 2.33951}},
        {{47},
        {19.22320, 17.67107, 5.07851, 4.43017, 1.59588, 0.00010},
        {0.59542, 6.92490, 0.00010, 24.85505, 87.61222, 31.90172}},
        {{48},
        {19.16300, 18.59170, 4.95237, 4.27994, 2.00969, 0.00010},
        {0.54868, 6.39500, 0.00010, 26.18224, 93.70112, 8.23922}},
        {{49},
        {19.22704, 19.09869, 4.79841, 4.37320, 2.50037, 0.00010},
        {5.84698, 0.50421, 0.00010, 27.22571, 81.57248, 31.56814}},
        {{50},
        {19.04077, 13.05412, 6.63670, 4.95963, 4.60941, 2.69795},
        {0.46176, 5.31900, 5.31953, 28.54198, 0.00010, 72.65174}},
        {{51},
        {19.96327, 18.99686, 6.19148, 4.38583, 2.46194, 0.00010},
        {4.81879, 0.42169, 28.79858, 0.00010, 70.63864, 12.77096}},
        {{52},
        {18.97925, 15.69578, 7.06433, 4.42489, 4.10018, 2.73271},
        {0.38267, 4.34879, 26.93604, 4.35210, 0.00010, 61.59836}},
        {{53},
        {20.29787, 19.00556, 9.04165, 3.76022, 1.89561, 0.00010},
        {3.93838, 0.34588, 26.70066, 0.00010, 65.34476, 20.30305}},
        {{54},
        {0.79375, 0.54736, 0.46161, 0.13918, 0.05800, 0.00010},
        {2.88678, 1.16905, 6.18250, 0.31715, 12.60983, 28.15927}},
        {{55},
        {0.82577, 0.73691, 0.23557, 0.20135, 0.00034, 0.00010},
        {2.04212, 0.80252, 4.60157, 0.21162, 43.68258, 103.45510}},
        {{56},
        {2.03492, 1.64286, 0.68060, 0.67022, 0.51650, 0.45488},
        {25.99675, 11.77809, 0.51013, 0.97866, 0.16915, 57.91874}},
        {{57},
        {3.56378, 2.14950, 1.52760, 1.47980, 0.27065, 0.00010},
        {14.10561, 5.60491, 0.32801, 46.88862, 0.00980, 10.98084}},
        {{58},
        {3.22684, 2.47111, 1.59839, 1.28490, 1.11335, 0.30182},
        {4.95997, 14.45952, 0.17267, 11.39653, 43.30817, 0.96703}},
        {{59},
        {3.69529, 3.30459, 1.68333, 0.69149, 0.62431, 0.00088},
        {3.24183, 7.07179, 0.12279, 15.45334, 1.43664, 35.26383}},
        {{60},
        {4.30385, 2.58390, 1.71397, 1.39368, 0.00470, 0.00010},
        {4.02045, 1.85304, 0.10693, 8.78523, 58.58712, 125.50050}},
        {{61},
        {4.19367, 3.00032, 1.71590, 1.08840, 0.00167, 0.00010},
        {3.37134, 1.58637, 0.09158, 6.99679, 45.26456, 113.97270}},
        {{62},
        {5.49488, 3.33770, 2.38765, 1.59864, 1.17986, 0.00010},
        {2.60802, 37.46289, 1.09647, 0.06439, 80.52337, 56.27056}},
        {{63},
        {3.98392, 3.53675, 1.72808, 0.75103, 0.00013, 0.00010},
        {2.94648, 1.39488, 0.08069, 5.91604, 56.23176, 79.76580}},
        {{64},
        {7.13932, 6.34213, 2.29801, 1.97826, 0.22854, 0.00983},
        {1.18073, 19.52901, 61.04850, 0.08057, 23.18225, 0.09759}},
        {{65},
        {8.00372, 7.44077, 1.42217, 1.13491, 0.00010, 0.00010},
        {12.70476, 0.77473, 0.00010, 32.44270, 199.99900, 82.98298}},
        {{66},
        {8.66803, 7.39747, 1.38325, 0.55348, 0.00010, 0.00010},
        {10.62955, 0.66306, 0.00010, 30.98476, 199.99880, 82.97898}},
        {{67},
        {9.01395, 7.36477, 1.32160, 0.30179, 0.00010, 0.00010},
        {8.86658, 0.56771, 0.00010, 29.98133, 137.40030, 53.69811}},
        {{68},
        {9.67607, 7.35874, 1.66775, 1.29681, 0.00010, 0.00010},
        {7.92858, 0.50388, 23.88214, 0.00010, 92.10388, 145.58810}},
        {{69},
        {9.56376, 7.35320, 1.26997, 0.81496, 0.00010, 0.00010},
        {7.72729, 0.49604, 0.00010, 22.37931, 92.10560, 145.58920}},
        {{70},
        {9.22395, 7.35117, 1.23367, 0.19305, 0.00010, 0.00010},
        {7.44634, 0.48595, 0.00010, 28.20512, 92.10930, 145.59010}},
        {{71},
        {10.14209, 7.35015, 2.25361, 1.23887, 0.01533, 0.00010},
        {6.90615, 0.44224, 20.14575, 0.00010, 120.21700, 55.09812}},
        {{72},
        {10.05723, 7.34875, 1.38759, 1.20752, 0.00010, 0.00010},
        {6.75290, 0.43509, 18.25122, 0.00010, 120.22150, 55.11062}},
        {{73},
        {9.37695, 7.36389, 1.11621, 0.14450, 0.00010, 0.00010},
        {6.31625, 0.41568, 0.00010, 25.36044, 199.99870, 82.97847}},
        {{74},
        {10.54130, 4.41457, 2.93436, 2.87024, 1.17229, 0.06743},
        {6.04009, 0.38967, 0.38966, 16.94938, 0.00010, 59.98400}},
        {{75},
        {10.45597, 4.43683, 2.92505, 2.06149, 1.11981, 0.00120},
        {5.90641, 0.38863, 0.37041, 15.34221, 0.00010, 59.68271}},
        {{76},
        {10.86580, 7.35401, 3.49267, 1.09987, 0.18537, 0.00249},
        {5.30450, 0.34487, 14.15718, 0.00010, 38.60730, 100.13560}},
        {{77},
        {11.04414, 4.43611, 4.06737, 2.44502, 0.00559, 0.00189},
        {5.32462, 0.15971, 0.47488, 13.90108, 100.14020, 38.59723}},
        {{78},
        {10.80739, 7.37819, 1.80548, 1.00948, 0.00010, 0.00010},
        {5.12031, 0.33181, 12.46589, 0.00010, 100.14660, 38.60185}},
        {{79},
        {11.32394, 7.35828, 4.08542, 1.03934, 0.19438, 0.00010},
        {4.71611, 0.30793, 12.87900, 0.00024, 43.73118, 103.91920}},
        {{80},
        {11.27641, 7.37595, 3.32058, 0.98461, 0.04263, 0.00010},
        {4.63894, 0.30169, 11.63908, 0.00010, 44.10289, 103.92070}},
        {{81},
        {11.59539, 7.37601, 4.75131, 0.95818, 0.31843, 0.00010},
        {4.18474, 0.27510, 11.19206, 0.00010, 36.27509, 93.95933}},
        {{82},
        {11.58135, 7.38964, 4.01201, 0.91419, 0.10353, 0.00010},
        {4.13155, 0.27012, 10.32693, 0.00010, 35.20369, 93.95908}},
        {{83},
        {11.83838, 5.16446, 4.59215, 3.72826, 0.67719, 0.00010},
        {3.76040, 9.57707, 0.31557, 0.11646, 25.17286, 96.76703}},
        {{84},
        {12.08932, 7.37051, 4.53328, 0.89389, 0.11440, 0.00010},
        {3.73486, 0.24588, 9.52524, 0.00100, 36.54998, 96.77110}},
        {{85},
        {11.74994, 6.77249, 6.21229, 1.75552, 1.47560, 0.03461},
        {3.34714, 0.23831, 8.32820, 23.58346, 0.04331, 98.01738}},
        {{86},
        {11.83187, 5.78192, 5.77531, 2.46041, 1.14698, 0.00353},
        {3.33965, 0.25530, 8.03031, 0.08201, 19.99327, 98.02090}},
        {{87},
        {12.49609, 7.88148, 4.99190, 2.05602, 0.57505, 0.00010},
        {3.52509, 0.16619, 9.20541, 1.71372, 24.20427, 82.21923}},
        {{88},
        {10.80193, 7.89470, 5.30620, 3.91136, 0.08693, 0.00010},
        {3.67800, 0.15468, 2.08510, 9.11568, 34.76155, 99.34953}},
        {{89},
        {8.64238, 8.44015, 7.88210, 2.99985, 0.03590, 0.00010},
        {3.75852, 2.14595, 0.14366, 8.16207, 30.93576, 72.31449}},
        {{90},
        {14.72809, 7.73340, 4.08153, 3.89920, 2.84995, 2.70412},
        {1.87781, 0.11285, 23.45650, 3.65207, 21.50646, 68.50430}},
        {{91},
        {17.72736, 7.70846, 6.22707, 4.23320, 0.10456, 0.00010},
        {1.68258, 0.09962, 13.34713, 25.64859, 76.90928, 199.99860}},
        {{92},
        {13.56253, 9.15282, 7.57461, 4.23621, 1.47524, 0.00010},
        {1.52639, 13.37893, 0.09009, 1.50827, 28.97999, 162.86130}},
        {{93},
        {17.83594, 10.00061, 7.34299, 0.76995, 0.05161, 0.00010},
        {1.37290, 11.94201, 0.07979, 27.59179, 0.08311, 137.72530}},
        {{94},
        {17.88797, 10.57832, 7.18725, 0.34750, 0.00010, 0.00010},
        {1.24006, 10.60035, 0.06944, 29.00543, 131.45550, 1.67829}},
        {{95},
        {17.94269, 11.64938, 7.03542, 1.17571, 0.20353, 0.00010},
        {1.13911, 10.82291, 0.06147, 34.40293, 1.15832, 134.27490}},
        {{96},
        {17.35713, 10.99074, 7.04050, 0.57079, 0.04542, 0.00010},
        {1.13181, 9.52278, 0.06199, 1.11378, 134.27980, 38.40765}},
        {{97},
        {16.70847, 11.98967, 6.70451, 1.98553, 1.61267, 0.00010},
        {1.02628, 9.86398, 0.04848, 26.23584, 1.02613, 83.38388}},
        {{98},
        {16.84671, 11.18317, 6.67150, 1.21668, 0.08306, 0.00010},
        {1.01489, 8.31776, 0.04772, 1.01511, 36.37142, 83.39908}},
        {{99},
        {16.20121, 13.68489, 5.92693, 2.62037, 2.56751, 0.00010},
        {0.83651, 8.66621, 0.02083, 0.83653, 22.32915, 67.41669}},
        {{100},
        {15.97671, 13.58921, 5.91839, 2.79182, 1.72564, 0.00010},
        {0.83452, 8.38679, 0.02066, 0.83387, 21.20783, 67.42265}},
        {{101},
        {14.55243, 14.36520, 5.43109, 3.60085, 2.86567, 1.18601},
        {8.09600, 0.75250, 0.00422, 0.75381, 21.00325, 0.75895}},
        {{102},
        {14.57165, 14.10996, 5.40851, 3.65768, 1.90013, 1.35484},
        {7.90759, 0.75012, 0.00354, 0.75338, 19.97214, 0.75124}},
        {{103},
        {19.27390, 15.67787, 5.26036, 3.78685, 0.00010, 0.00010},
        {0.69511, 7.84482, 0.00010, 22.21775, 60.82368, 1.12994}},
        {{104},
        {19.16608, 15.58248, 5.24991, 1.97949, 0.02452, 0.00010},
        {0.69220, 7.50980, 0.00010, 19.35021, 0.69139, 60.83056}},
        {{105},
        {19.29333, 16.76786, 5.18419, 4.69146, 0.06334, 0.00010},
        {0.64534, 7.54710, 0.00010, 23.16034, 100.32570, 2.35114}},
        {{106},
        {19.26038, 16.76118, 5.17728, 3.80102, 0.00010, 0.00010},
        {0.64383, 7.44215, 0.00010, 21.24567, 100.31430, 2.43992}},
        {{107},
        {19.24328, 17.81622, 5.07556, 3.86538, 0.00010, 0.00010},
        {0.59548, 7.03822, 0.00010, 20.12238, 87.60555, 31.88584}},
        {{108},
        {19.15099, 19.02664, 5.11556, 1.72846, 1.00259, 0.00010},
        {0.55860, 6.79490, 0.00370, 25.60539, 8.23095, 93.69624}},
        {{109},
        {19.14517, 19.11002, 4.80720, 4.48861, 0.25075, 0.20103},
        {5.86776, 0.50516, 0.00010, 24.33452, 87.00222, 31.41846}},
        {{110},
        {19.71431, 19.14550, 4.79767, 2.34645, 0.00010, 0.00010},
        {6.04052, 0.50506, 0.00010, 16.17828, 87.05909, 31.49791}},
        {{111},
        {19.06093, 12.90928, 6.64901, 4.63278, 4.60732, 0.14140},
        {0.46390, 5.35884, 5.35853, 0.00010, 21.75129, 70.66362}},
        {{112},
        {19.55274, 19.11016, 4.62585, 1.75378, 0.96170, 0.00010},
        {5.57560, 0.46433, 0.00010, 15.08594, 5.57571, 70.66860}},
        {{113},
        {18.97534, 15.68841, 6.74714, 4.42194, 4.08431, 4.06854},
        {0.38165, 4.33217, 26.51128, 4.35007, 0.00013, 70.73529}},
        /**
         * cctbx - 3 Gaussians.
         */
        {{114},
        {1.05674999997, 0.644323530897, 0.298305313502},
        {1.10072695216, 3.30968118776, 0.263851354839}},
        {{115},
        {4.36297453919, 3.95402186199, 1.67261115359},
        {5.84416794424, 2.0522287083, 0.0990905852289}},
        {{116},
        {4.1982011524, 4.10980380127, 1.68447184661},
        {1.73308116771, 4.58815466909, 0.0858298380754}},
        {{117},
        {4.47168612436, 3.87243074662, 1.66071703302},
        {1.44515163245, 3.79224664386, 0.0714406145137}},
        {{118},
        {9.00614617514, 7.11474468263, 1.85112246997},
        {11.5043846011, 0.730949041574, 0.0473828343711}},
        {{119},
        {9.16075817591, 7.1050465418, 1.71459201986},
        {9.35476737709, 0.620771447674, 0.0342761622543}},
        {{120},
        {7.93431359573, 7.75702894924, 4.27186054645},
        {0.344341308343, 12.0850676027, 4.28886237826}},
        {{121},
        {10.0039669344, 6.48634784642, 2.47486846392},
        {8.68083236931, 0.63076461812, 0.080571202237}},
        {{122},
        {9.30135876449, 7.07501002111, 1.60674098372},
        {7.7294767073, 0.529160463879, 0.026439065502}},
        {{123},
        {7.97543708242, 7.63636302328, 5.34474255061},
        {0.309079820533, 11.6025702687, 4.25903863771}},
        {{124},
        {8.18999077376, 7.90783514463, 3.88300169133},
        {9.4272329602, 0.307235108471, 3.6928463073}},
        {{125},
        {9.35615562168, 7.03336806389, 1.60169919909},
        {6.58156213681, 0.47100169955, 0.0241326952344}},
        {{126},
        {7.99912950446, 7.38596545164, 6.58011100415},
        {0.283923739083, 4.4824019183, 12.0615843142}},
        {{127},
        {8.70462787589, 7.86457501857, 4.40400027591},
        {8.79040932546, 0.275593957027, 3.39537512496}},
        {{128},
        {9.56593817028, 8.07602885623, 5.33445562929},
        {4.56386883527, 0.259145646999, 12.8970300261}},
        {{129},
        {7.89836965773, 7.08191595383, 7.00167267967},
        {0.251091897431, 9.28582994856, 3.8036260062}},
        {{130},
        {8.78159491796, 7.82567400183, 4.37597058338},
        {7.25270024967, 0.247385674714, 3.06034043675}},
        {{131},
        {10.7556787804, 8.08468656259, 5.13923400598},
        {4.30170995299, 0.237147049617, 12.6864130192}},
        {{132},
        {11.1450231671, 7.57801275896, 4.22386820999},
        {7.10765696463, 0.214122838685, 2.25363704097}},
        {{133},
        {9.30087028578, 7.86055043642, 7.79103605212},
        {3.34124560949, 0.207379836561, 9.72882671072}},
        {{134},
        {9.35811523679, 7.88288171805, 6.7382513407},
        {3.45256987642, 0.209910830837, 8.73850168377}},
        {{135},
        {10.9574548068, 7.91364440258, 7.09037754465},
        {3.30397344741, 0.193332220057, 9.84576083535}},
        {{136},
        {9.03277595215, 8.16707734383, 7.76819167497},
        {2.96007681004, 7.67163956047, 0.188408691652}},
        {{137},
        {12.359628475, 7.87131386986, 6.73863127449},
        {3.18130979285, 0.178269967352, 9.74802018078}},
        {{138},
        {13.4378597184, 7.98154466967, 6.5512483394},
        {3.02261523908, 0.1687408754, 9.55493525266}},
        {{139},
        {11.4539, 8.90972752105, 7.60276831073},
        {2.37794666223, 6.59530146321, 0.142862536768}},
        {{140},
        {13.8480203642, 7.75041063864, 6.39192497212},
        {2.4567640923, 0.137907723422, 6.46766456568}},
        {{141},
        {18.0767914951, 10.2422601571, 7.64105919371},
        {1.56063518321, 15.475201363, 0.0908340434334}},
        {{142},
        {17.9433074992, 10.5236761903, 7.50305231095},
        {1.40475303467, 12.9752454537, 0.0824340965413}},
        {{143},
        {17.7734133861, 10.7877635507, 7.41457517106},
        {1.26744992396, 11.0724056318, 0.0764947481608}},
        {{144},
        {17.6183400651, 11.881198903, 8.43509868407},
        {1.31193510861, 12.7973587692, 0.105008616681}},
        {{145},
        {17.7398863368, 10.9188297535, 7.32567504128},
        {1.160593637, 9.69655399567, 0.0710750543894}},
        {{146},
        {18.3188870634, 13.0152411164, 7.58377586849},
        {1.14468107806, 12.0720126341, 0.072546006161}},
        {{147},
        {17.6781713727, 11.0390849019, 7.27524687075},
        {1.06872809116, 8.57978579159, 0.0667266185651}},
        {{148},
        {18.3425380073, 14.9935363225, 7.56968941267},
        {0.990303247192, 10.8757033155, 0.0658550377396}},
        {{149},
        {18.2136421714, 14.4015964784, 7.32335925891},
        {0.953792258618, 9.83764291281, 0.0610368319049}},
        {{150},
        {15.5169232348, 14.6084670706, 11.7836912415},
        {10.5233527721, 1.19805976888, 0.173293763471}},
        {{151},
        {18.2500753825, 15.4398316764, 7.24317937151},
        {0.883340061883, 9.34400156369, 0.0566949236825}},
        {{152},
        {21.4258951173, 12.2630134383, 10.2462505043},
        {0.382415912786, 13.6148868752, 3.88227600024}},
        {{153},
        {18.2501716117, 16.4875851039, 7.19092695047},
        {0.820905236411, 8.85744233571, 0.0533496503559}},
        {{154},
        {21.1212947884, 13.6396062605, 10.1654365258},
        {0.349668807221, 12.3702805418, 3.44568968188}},
        {{155},
        {20.8284525128, 15.0096107479, 10.0792911382},
        {0.320672873475, 11.3005350502, 3.06922248664}},
        {{156},
        {20.1774443662, 16.9816497138, 8.77051768304},
        {0.288373343583, 9.13033692607, 2.35002850502}},
        {{157},
        {22.1164438441, 18.4102698652, 7.39983890917},
        {0.312440051695, 4.73418318858, 19.5310272638}},
        {{158},
        {20.9547831947, 18.1555765089, 6.81840971331},
        {7.0152600456, 0.605887177524, 0.0401455609886}},
        {{159},
        {21.9870799745, 18.1894285843, 7.77246946156},
        {0.294069193351, 4.28026514278, 16.4477064302}},
        {{160},
        {21.1518685061, 17.9939494168, 6.78998333403},
        {6.32805369701, 0.566290989365, 0.0381542687902}},
        {{161},
        {30.112436978, 28.6235543927, 19.1340014396},
        {1.77440063269, 0.184015028333, 12.1530773066}},
        {{162},
        {0.502196691881, 0.373818972889, 0.123052629233},
        {13.478927669, 38.6665372038, 3.38304689207}},
        {{163},
        {1.08189530121, 0.467245206939, 0.4501601147},
        {6.70265036943, 19.957676141, 1.5821652553}},
        {{164},
        {1.07146599943, 1.05534384128, 0.863685495839},
        {4.30327380157, 114.129734127, 1.03599382371}},
        {{165},
        {2.06778108412, 1.03403401833, 0.878360479364},
        {54.3775616887, 2.21128244408, 0.548264395981}},
        {{166},
        {2.40390483689, 1.74082003938, 0.840968677045},
        {40.6492977422, 0.649340687236, 12.2463914516}},
        {{167},
        {2.51340127252, 1.74867019409, 1.72398202356},
        {31.8053433708, 0.445605499982, 10.5831679451}},
        {{168},
        {2.99954939487, 2.25583887941, 1.7278842283},
        {23.2726795155, 7.45433091596, 0.316224876669}},
        {{169},
        {3.21184129664, 3.04156392126, 1.73156010601},
        {18.8370006399, 5.90590162558, 0.241263012791}},
        {{170},
        {3.76051707547, 3.47766990973, 1.74594840518},
        {4.73185569767, 15.4384441173, 0.194238121265}},
        {{171},
        {4.38310831035, 3.83422263038, 1.76279016611},
        {3.82501909721, 12.6640899017, 0.161786329667}},
        {{172},
        {6.6351112713, 3.01293367286, 1.30238479723},
        {5.5442312674, 0.545797156151, 90.8590190323}},
        {{173},
        {6.67814146187, 2.96430131029, 2.32496434711},
        {4.14073262313, 0.428901927646, 71.0199150719}},
        {{174},
        {6.62238016224, 3.27781520856, 3.05972557508},
        {3.3038329346, 59.3145382378, 0.374151395498}},
        {{175},
        {6.08447684913, 4.2732659722, 3.59835544898},
        {2.83886951807, 46.6331390474, 0.418767784547}},
        {{176},
        {5.27941125504, 5.27354842119, 4.40634943423},
        {36.7775304861, 2.58877453733, 0.47041178021}},
        {{177},
        {6.83012748437, 6.13863224738, 2.99357763173},
        {0.664089368899, 30.1886951281, 3.52397406633}},
        {{178},
        {8.12176088935, 6.65345104038, 2.19043698007},
        {0.697928192858, 25.8943955401, 6.00738032174}},
        {{179},
        {8.1978704772, 7.25214271687, 2.51643840211},
        {0.603152037749, 22.1423332433, 6.28874083466}},
        {{180},
        {8.66788257982, 8.58340898529, 1.68051604224},
        {12.3159255661, 0.562749994595, 110.002273029}},
        {{181},
        {8.75937156177, 8.41257168569, 2.76798129934},
        {9.64475680185, 0.475138202259, 97.3905740758}},
        {{182},
        {9.51340433555, 8.43057223084, 3.00293226194},
        {8.72304034893, 0.421969259986, 85.1544003488}},
        {{183},
        {10.2799540619, 8.45913489426, 3.20116956042},
        {7.88841343088, 0.381270873952, 73.589906538}},
        {{184},
        {11.1458421487, 8.52770210726, 3.27333965837},
        {7.26900302269, 0.350572882934, 67.9662846278}},
        {{185},
        {12.1965729193, 8.6651380941, 3.06905493686},
        {6.89852281719, 0.333439737044, 51.4431406852}},
        {{186},
        {12.7428954458, 8.68059401604, 3.51227744547},
        {6.09717645893, 0.307362041711, 54.7923744615}},
        {{187},
        {13.5386808693, 8.76290322281, 3.6282889925},
        {5.5999743031, 0.290387429344, 49.6184457459}},
        {{188},
        {14.2511611792, 9.13683816117, 3.55259119666},
        {5.36759291075, 0.305437893913, 48.447700971}},
        {{189},
        {14.9767789806, 9.35030169796, 3.61038378024},
        {5.01191696697, 0.302607553312, 44.9232209459}},
        {{190},
        {15.9430456345, 9.6553328387, 3.33083077178},
        {4.82663324046, 0.306553792833, 35.4599869576}},
        {{191},
        {15.9676170894, 10.4463803921, 3.5271729716},
        {4.64159415074, 0.358107174158, 41.3067132366}},
        {{192},
        {16.4944618554, 10.3236971493, 4.10716029583},
        {4.12737191087, 0.318425619984, 44.5159627214}},
        {{193},
        {16.4458600761, 10.5599493419, 4.92533026545},
        {3.72882293769, 0.316174473194, 42.0934525241}},
        {{194},
        {17.1960878679, 9.69900138738, 6.04192031603},
        {3.12619083301, 0.234035807171, 36.5140401994}},
        {{195},
        {17.0504033455, 9.78301185102, 7.10176272738},
        {2.79743745985, 0.227209294309, 31.7723375223}},
        {{196},
        {17.5420861481, 9.08359865439, 8.30546540071},
        {2.3876252301, 0.174206423325, 27.4378361399}},
        {{197},
        {17.3989367325, 9.36450036612, 9.17346482754},
        {2.16761932439, 24.4439214341, 0.170251945105}},
        {{198},
        {23.910804374, 10.5360489582, 2.44931455737},
        {0.928223646688, 13.1425847975, 87.5230688031}},
        {{199},
        {23.3302207832, 11.1676164796, 3.40901694124},
        {0.805464539522, 10.7378777121, 93.295102231}},
        {{200},
        {23.5235691064, 11.4829330129, 3.916262423},
        {0.757601750963, 10.2442370891, 83.0146587908}},
        {{201},
        {23.3726972532, 12.2519628369, 4.31377561634},
        {0.684926546016, 9.56037750445, 74.3015056094}},
        {{202},
        {23.5031600172, 12.6386032107, 4.7681241134},
        {0.642934099047, 9.20719558971, 52.8701133916}},
        {{203},
        {23.1827814092, 13.2100575903, 5.50231022913},
        {0.577001089934, 8.23840117065, 43.9978888806}},
        {{204},
        {23.5346902744, 14.5111421298, 4.88999302843},
        {0.555290522238, 8.51654999704, 55.4141673419}},
        {{205},
        {23.122417304, 14.8125929435, 5.95928154858},
        {0.497739492816, 7.50901560896, 36.7823648354}},
        {{206},
        {23.3012159473, 15.9533741871, 5.6515846267},
        {0.473398855106, 7.49282025419, 36.2469928289}},
        {{207},
        {22.7629045939, 14.6226398022, 8.54242304795},
        {0.423182233992, 6.12063761773, 22.660232699}},
        {{208},
        {23.3122067553, 18.4080605903, 5.20642112685},
        {0.414662819069, 7.07241088531, 34.9460945594}},
        {{209},
        {22.9961521829, 19.1440914287, 5.75799203993},
        {0.378523645723, 6.39052944494, 35.7007532116}},
        {{210},
        {23.1150202427, 20.2733776757, 5.50800401336},
        {0.361318740795, 6.21302979882, 43.859239884}},
        {{211},
        {23.073115678, 20.721043971, 6.11619707401},
        {0.341229928758, 5.77195480985, 45.0386050033}},
        {{212},
        {22.7485769177, 21.0063112152, 7.16102625262},
        {0.312929289314, 5.18082235944, 41.8483442701}},
        {{213},
        {22.4342311303, 21.2105816232, 8.27410663564},
        {0.28840668083, 4.66744106421, 38.2957918128}},
        {{214},
        {22.1334760622, 21.3152722398, 9.47258237105},
        {0.265661657324, 4.19946212892, 34.5739885559}},
        {{215},
        {21.8411079371, 21.3868466973, 10.6938907155},
        {0.245775394234, 3.79473697198, 31.3171324821}},
        {{216},
        {42.250520114, 25.51363075, 12.0790174828},
        {0.382882567089, 4.03482946811, 23.3930785954}}
    };

    static {
        for (int i = 0; i < atoms.length; i++) {
            formfactors.put(atomsi[i], ffactors[i]);
        }
    }
}
