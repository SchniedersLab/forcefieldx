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
package ffx.xray;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static ffx.numerics.VectorMath.determinant3;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.mat3inverse;
import static ffx.numerics.VectorMath.mat3mat3;
import static ffx.numerics.VectorMath.scalarmat3mat3;
import static ffx.numerics.VectorMath.vec3mat3;
import static ffx.numerics.VectorMath.rsq;
import static ffx.numerics.VectorMath.b2u;
import static ffx.numerics.VectorMath.u2b;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.Arrays;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * This implementation uses the coefficients from International Tables, Vol. C,
 * chapter 4.4.4.
 *
 * @author Tim Fenn<br>
 *
 * @see <a href="http://dx.doi.org/10.1107/9780955360206000059410.1107/97809553602060000594" target="_blank">
 * V. F. Sears, Int. Tables Vol. C (2006). Table 4.4.4.1</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444909022707" target="_blank">
 * M. J. Schnieders, T. D. Fenn, V. S. Pande and A. T. Brunger,
 * Acta Cryst. (2009). D65 952-965.</a>
 */
public final class NeutronFormFactor implements FormFactor {

    private static final Logger logger = Logger.getLogger(ffx.xray.NeutronFormFactor.class.getName());
    private static final double twopi2 = 2.0 * PI * PI;
    private static final double twopi32 = pow(2.0 * PI, -1.5);
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
    private static final String[] atoms = {"H", "D", "He", "Li", "Be", "B", "C",
        "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
        "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Ru",
        "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
        "Pb", "Bi", "Th", "U"};
    private static final String[] atomsi = {"1", "1", "2", "3", "4", "5", "6",
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
        "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
        "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "44",
        "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56",
        "57", "58", "59", "60", "62", "63", "64", "65", "66", "67", "68", "69",
        "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81",
        "82", "83", "90", "92"};
    private static final double[][][] ffactors = {
        {{0}, {-3.7390, 0.0}},
        {{1}, {6.671, 0.0}},
        {{2}, {3.26, 0.0}},
        {{3}, {-1.90, 0.0}},
        {{4}, {7.79, 0.0}},
        {{5}, {5.30, 0.213}},
        {{6}, {6.6460, 0.0}},
        {{7}, {9.36, 0.0}},
        {{8}, {5.803, 0.0}},
        {{9}, {5.654, 0.0}},
        {{10}, {4.566, 0.0}},
        {{11}, {3.63, 0.0}},
        {{12}, {5.375, 0.0}},
        {{13}, {3.449, 0.0}},
        {{14}, {4.1491, 0.0}},
        {{15}, {5.13, 0.0}},
        {{16}, {2.847, 0.0}},
        {{17}, {9.5770, 0.0}},
        {{18}, {1.909, 0.0}},
        {{19}, {3.67, 0.0}},
        {{20}, {4.70, 0.0}},
        {{21}, {12.29, 0.0}},
        {{22}, {-3.370, 0.0}},
        {{23}, {-0.3824, 0.0}},
        {{24}, {3.635, 0.0}},
        {{25}, {-3.750, 0.0}},
        {{26}, {9.45, 0.0}},
        {{27}, {2.49, 0.0}},
        {{28}, {10.3, 0.0}},
        {{29}, {7.718, 0.0}},
        {{30}, {5.60, 0.0}},
        {{31}, {7.288, 0.0}},
        {{32}, {8.185, 0.0}},
        {{33}, {6.58, 0.0}},
        {{34}, {7.970, 0.0}},
        {{35}, {6.795, 0.0}},
        {{36}, {7.81, 0.0}},
        {{37}, {7.09, 0.0}},
        {{38}, {7.02, 0.0}},
        {{39}, {7.75, 0.0}},
        {{40}, {7.16, 0.0}},
        {{41}, {7.054, 0.0}},
        {{42}, {6.715, 0.0}},
        {{43}, {7.03, 0.0}},
        {{44}, {5.88, 0.0}},
        {{45}, {5.91, 0.0}},
        {{46}, {5.922, 0.0}},
        {{47}, {4.87, -0.70}},
        {{48}, {2.08, -0.0539}},
        {{49}, {6.225, 0.0}},
        {{50}, {5.57, 0.0}},
        {{51}, {5.80, 0.0}},
        {{52}, {5.28, 0.0}},
        {{53}, {4.92, 0.0}},
        {{54}, {5.42, 0.0}},
        {{55}, {5.07, 0.0}},
        {{56}, {8.24, 0.0}},
        {{57}, {4.84, 0.0}},
        {{58}, {4.58, 0.0}},
        {{59}, {7.69, 0.0}},
        {{60}, {0.80, -1.65}},
        {{61}, {7.22, -1.26}},
        {{62}, {6.5, -13.82}},
        {{63}, {7.38, 0.0}},
        {{64}, {16.9, -0.276}},
        {{65}, {8.01, 0.0}},
        {{66}, {7.79, 0.0}},
        {{67}, {7.07, 0.0}},
        {{68}, {12.43, 0.0}},
        {{69}, {7.21, 0.0}},
        {{70}, {7.77, 0.0}},
        {{71}, {6.91, 0.0}},
        {{72}, {4.86, 0.0}},
        {{73}, {9.2, 0.0}},
        {{74}, {10.7, 0.0}},
        {{75}, {10.6, 0.0}},
        {{76}, {9.60, 0.0}},
        {{77}, {7.63, 0.0}},
        {{78}, {12.692, 0.0}},
        {{79}, {8.776, 0.0}},
        {{80}, {9.405, 0.0}},
        {{81}, {8.532, 0.0}},
        {{82}, {10.31, 0.0}},
        {{83}, {8.417, 0.0}}
    };

    static {
        for (int i = 0; i < atoms.length; i++) {
            formfactors.put(atomsi[i], ffactors[i]);
        }
    }
    private final Atom atom;
    protected final int ffindex;
    private double xyz[] = new double[3];
    private double dxyz[] = new double[3];
    private double biso;
    private double uadd;
    private boolean hasanisou;
    private double uaniso[] = null;
    private double occ;
    private final double a[] = new double[2];
    private double ainv[] = new double[1];
    private double binv[] = new double[1];
    private double u[][][] = new double[1][3][3];
    private double uinv[][][] = new double[1][3][3];
    private double jmat[][][] = new double[6][3][3];
    private double gradp[] = new double[6];
    private double gradu[] = new double[6];
    private double resv[] = new double[3];
    private double resm[][] = new double[3][3];

    public NeutronFormFactor(Atom atom) {
        this(atom, 0.0, atom.getXYZ());
    }

    public NeutronFormFactor(Atom atom, double badd) {
        this(atom, badd, atom.getXYZ());
    }

    public NeutronFormFactor(Atom atom, double badd, double xyz[]) {
        this.atom = atom;
        this.uadd = b2u(badd);
        double ffactor[][];
        String key = "" + atom.getAtomicNumber();

        int atomindex = atom.getFormFactorIndex();
        if (atomindex < 0) {
            ffactor = getFormFactor(key);
            ffindex = (int) ffactor[0][0];
            atom.setFormFactorIndex(ffindex);
        } else {
            ffindex = atomindex;
            ffactor = ffactors[atomindex];
        }

        System.arraycopy(ffactor[1], 0, a, 0, ffactor[1].length);
        if (a[1] != 0.0) {
            logger.severe("complex neutron form factor method not supported");
        }
        occ = atom.getOccupancy();

        if (occ <= 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append("zero occ for atom: " + atom.toString() + "\n");
            sb.append("(atom will not contribute to electron density calculation)\n");
            logger.warning(sb.toString());
        }

        update(xyz, uadd);
    }

    public static int getFormFactorIndex(String atom) {
        double ffactor[][] = null;

        ffactor = getFormFactor(atom);
        if (ffactor != null) {
            return (int) ffactor[0][0];
        } else {
            return -1;
        }
    }

    public static double[] getFormFactorA(String atom) {
        double ffactor[][] = null;

        ffactor = getFormFactor(atom);
        if (ffactor != null) {
            return ffactor[1];
        } else {
            return null;
        }
    }

    public static double[][] getFormFactor(String atom) {
        double ffactor[][] = null;

        if (formfactors.containsKey(atom)) {
            ffactor = (double[][]) formfactors.get(atom);
        } else {
            String message = "Form factor for atom: " + atom
                    + " not found!";
            logger.severe(message);
        }

        return ffactor;
    }

    public double f(HKL hkl) {
        double sum = a[0] * exp(-twopi2 * Crystal.quad_form(hkl, u[0]));
        return occ * sum;
    }

    @Override
    public double rho(double f, double xyz[]) {
        diff(this.xyz, xyz, dxyz);
        double sum = ainv[0] * exp(-0.5 * Crystal.quad_form(dxyz, uinv[0]));
        return f + (occ * twopi32 * sum);
    }

    @Override
    public void rho_grad(double xyz[], double dfc, RefinementMode refinementmode) {
        diff(this.xyz, xyz, dxyz);
        double r2 = rsq(dxyz);
        double aex;
        double rho;
        Arrays.fill(gradp, 0.0);
        Arrays.fill(gradu, 0.0);

        boolean refinexyz = false;
        boolean refineb = false;
        boolean refineanisou = false;
        boolean refineocc = false;
        if (refinementmode == RefinementMode.COORDINATES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refinexyz = true;
        }
        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineb = true;
            if (hasanisou) {
                refineanisou = true;
            }
        }
        if (refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineocc = true;
        }

        aex = ainv[0] * exp(-0.5 * Crystal.quad_form(dxyz, uinv[0]));

        if (refinexyz) {
            vec3mat3(dxyz, uinv[0], resv);
            gradp[0] += aex * dot(resv, vx);
            gradp[1] += aex * dot(resv, vy);
            gradp[2] += aex * dot(resv, vz);
        }

        if (refineocc) {
            gradp[3] += aex;
        }

        if (refineb) {
            gradp[4] += aex * 0.5 * (r2 * binv[0] * binv[0] - 3.0 * binv[0]);

            if (refineanisou) {
                scalarmat3mat3(-1.0, uinv[0], u11, resm);
                mat3mat3(resm, uinv[0], jmat[0]);
                scalarmat3mat3(-1.0, uinv[0], u22, resm);
                mat3mat3(resm, uinv[0], jmat[1]);
                scalarmat3mat3(-1.0, uinv[0], u33, resm);
                mat3mat3(resm, uinv[0], jmat[2]);
                scalarmat3mat3(-1.0, uinv[0], u12, resm);
                mat3mat3(resm, uinv[0], jmat[3]);
                scalarmat3mat3(-1.0, uinv[0], u13, resm);
                mat3mat3(resm, uinv[0], jmat[4]);
                scalarmat3mat3(-1.0, uinv[0], u23, resm);
                mat3mat3(resm, uinv[0], jmat[5]);

                gradu[0] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[0]) - uinv[0][0][0]);
                gradu[1] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[1]) - uinv[0][1][1]);
                gradu[2] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[2]) - uinv[0][2][2]);
                gradu[3] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[3]) - uinv[0][0][1] * 2.0);
                gradu[4] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[4]) - uinv[0][0][2] * 2.0);
                gradu[5] += aex * 0.5 * (-Crystal.quad_form(dxyz, jmat[5]) - uinv[0][1][2] * 2.0);
            }
        }
        rho = occ * twopi32 * gradp[3];

        // x, y, z
        if (refinexyz) {
            atom.addToXYZGradient(
                    dfc * occ * -twopi32 * gradp[0],
                    dfc * occ * -twopi32 * gradp[1],
                    dfc * occ * -twopi32 * gradp[2]);
        }

        // occ
        if (refineocc) {
            atom.addToOccupancyGradient(dfc * twopi32 * gradp[3]);
        }

        // Biso
        if (refineb) {
            atom.addToTempFactorGradient(dfc * b2u(occ * twopi32 * gradp[4]));
            // Uaniso
            if (hasanisou) {
                for (int i = 0; i < 6; i++) {
                    gradu[i] = dfc * occ * twopi32 * gradu[i];
                }
                atom.addToAnisouGradient(gradu);
            }
        }
    }

    @Override
    public void update(double xyz[]) {
        update(xyz, u2b(uadd));
    }

    @Override
    public void update(double xyz[], double badd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
        biso = atom.getTempFactor();
        uadd = b2u(badd);
        occ = atom.getOccupancy();

        // check occ is valid
        if (occ < 0.0) {
            StringBuilder sb = new StringBuilder();
            sb.append("negative occupancy for atom: " + atom.toString() + "\n");
            sb.append("resetting to 0.0\n");
            logger.warning(sb.toString());
            occ = 0.0;
            atom.setOccupancy(0.0);
        }

        // check if anisou changed
        if (atom.getAnisou() == null) {
            if (uaniso == null) {
                uaniso = new double[6];
            }
            hasanisou = false;
        } else {
            hasanisou = true;
        }

        if (hasanisou) {
            // first check the ANISOU is valid
            uaniso = atom.getAnisou();
            double det = determinant3(uaniso);

            if (det <= 1e-10) {
                StringBuilder sb = new StringBuilder();
                sb.append("non-positive definite ANISOU for atom: " + atom.toString() + "\n");
                sb.append("resetting ANISOU based on isotropic B: (" + biso + ")\n");
                logger.warning(sb.toString());

                uaniso[0] = uaniso[1] = uaniso[2] = b2u(biso);
                uaniso[3] = uaniso[4] = uaniso[5] = 0.0;
            }
        } else {
            if (biso < 0.0) {
                StringBuilder sb = new StringBuilder();
                sb.append("negative B factor for atom: " + atom.toString() + "\n");
                sb.append("resetting B to 0.01\n");
                logger.warning(sb.toString());
                atom.setTempFactor(0.01);
                uaniso[0] = uaniso[1] = uaniso[2] = b2u(0.01);
            } else {
                uaniso[0] = uaniso[1] = uaniso[2] = b2u(biso);
            }
            uaniso[3] = uaniso[4] = uaniso[5] = 0.0;
        }

        u[0][0][0] = uaniso[0] + uadd;
        u[0][1][1] = uaniso[1] + uadd;
        u[0][2][2] = uaniso[2] + uadd;
        u[0][0][1] = u[0][1][0] = uaniso[3];
        u[0][0][2] = u[0][2][0] = uaniso[4];
        u[0][1][2] = u[0][2][1] = uaniso[5];

        mat3inverse(u[0], uinv[0]);

        double det = determinant3(u[0]);
        ainv[0] = a[0] / sqrt(det);
        // b[i] = pow(det, 0.33333333333);
        det = determinant3(uinv[0]);
        binv[0] = pow(det, 0.33333333333);
    }
}
