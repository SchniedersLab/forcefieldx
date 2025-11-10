// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import java.util.HashMap;
import java.util.logging.Logger;

import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.numerics.math.MatrixMath.determinant3;
import static ffx.numerics.math.MatrixMath.mat3Inverse;
import static ffx.numerics.math.MatrixMath.mat3Mat3;
import static ffx.numerics.math.MatrixMath.scalarMat3Mat3;
import static ffx.numerics.math.MatrixMath.vec3Mat3;
import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.numerics.math.ScalarMath.quadForm;
import static ffx.numerics.math.ScalarMath.u2b;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * This implementation uses the coefficients from International Tables, Vol. C, chapter 4.4.4.
 *
 * @author Timothy D. Fenn
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000594" target="_blank"> V. F. Sears, Int.
 *     Tables Vol. C (2006). Table 4.4.4.1</a>
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000600" target="_blank"> B. T. M. Willis,
 *     Int. Tables Vol. C (2006). Chapter 6.1.3</a>
 * @see <a href="http://dx.doi.org/10.1107/S0907444909022707" target="_blank"> M. J. Schnieders, T.
 *     D. Fenn, V. S. Pande and A. T. Brunger, Acta Cryst. (2009). D65 952-965.</a>
 * @since 1.0
 */
public final class NeutronFormFactor implements FormFactor {

  private static final Logger logger = Logger.getLogger(NeutronFormFactor.class.getName());

  private static final double twoPI2 = 2.0 * PI * PI;
  private static final double inverseTwoPI32 = pow(2.0 * PI, -1.5);
  private static final double oneThird = 1.0 / 3.0;
  private static final double[] vx = {1.0, 0.0, 0.0};
  private static final double[] vy = {0.0, 1.0, 0.0};
  private static final double[] vz = {0.0, 0.0, 1.0};
  private static final double[][] u11 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  private static final double[][] u22 = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
  private static final double[][] u33 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
  private static final double[][] u12 = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  private static final double[][] u13 = {{0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  private static final double[][] u23 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
  private static final HashMap<String, double[][]> formfactors = new HashMap<>();
  private static final String[] atomNames = {
      "H", "D", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
      "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
      "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
      "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
      "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
      "Bi", "Th", "U"
  };
  private static final String[] atomicNumbers = {
      "1_1", "1_2", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
      "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
      "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49",
      "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "62", "63", "64", "65", "66",
      "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82",
      "83", "90", "92"
  };
  private static final double[][][] neutronFormFactors = {
      {{0}, {-3.7390, 0.0}},  // H
      {{1}, {6.671, 0.0}},    // D
      {{2}, {3.26, 0.0}},     // He
      {{3}, {-1.90, 0.0}},    // Li
      {{4}, {7.79, 0.0}},     // Be
      {{5}, {5.30, 0.213}},   // B
      {{6}, {6.6460, 0.0}},   // C
      {{7}, {9.36, 0.0}},     // N
      {{8}, {5.803, 0.0}},    // O
      {{9}, {5.654, 0.0}},    // F
      {{10}, {4.566, 0.0}},   // Ne
      {{11}, {3.63, 0.0}},    // Na
      {{12}, {5.375, 0.0}},   // Mg
      {{13}, {3.449, 0.0}},   // Al
      {{14}, {4.1491, 0.0}},  // Si
      {{15}, {5.13, 0.0}},    // P
      {{16}, {2.847, 0.0}},   // S
      {{17}, {9.5770, 0.0}},  // Cl
      {{18}, {1.909, 0.0}},   // Ar
      {{19}, {3.67, 0.0}},    // K
      {{20}, {4.70, 0.0}},    // Ca
      {{21}, {12.29, 0.0}},   // Sc
      {{22}, {-3.370, 0.0}},  // Ti
      {{23}, {-0.3824, 0.0}}, // V
      {{24}, {3.635, 0.0}},   // Cr
      {{25}, {-3.750, 0.0}},  // Mn
      {{26}, {9.45, 0.0}},    // Fe
      {{27}, {2.49, 0.0}},    // Co
      {{28}, {10.3, 0.0}},    // Ni
      {{29}, {7.718, 0.0}},   // Cu
      {{30}, {5.60, 0.0}},    // Zn
      {{31}, {7.288, 0.0}},   // Ga
      {{32}, {8.185, 0.0}},   // Ge
      {{33}, {6.58, 0.0}},    // As
      {{34}, {7.970, 0.0}},   // Se
      {{35}, {6.795, 0.0}},   // Br
      {{36}, {7.81, 0.0}},    // Kr
      {{37}, {7.09, 0.0}},    // Rb
      {{38}, {7.02, 0.0}},    // Sr
      {{39}, {7.75, 0.0}},    // Y
      {{40}, {7.16, 0.0}},    // Zr
      {{41}, {7.054, 0.0}},   // Nb
      {{42}, {6.715, 0.0}},   // Mo
      {{43}, {6.80, 0.0}},    // Tc
      {{44}, {7.03, 0.0}},    // Ru
      {{45}, {5.88, 0.0}},    // Rh
      {{46}, {5.91, 0.0}},    // Pd
      {{47}, {5.922, 0.0}},   // Ag
      {{48}, {4.87, -0.70}},  // Cd
      {{49}, {2.08, -0.0539}},// In; nCNS 4.065
      {{50}, {6.225, 0.0}},   // Sn
      {{51}, {5.57, 0.0}},    // Sb
      {{52}, {5.80, 0.0}},    // Te; nCNS 5.68
      {{53}, {5.28, 0.0}},    // I
      {{54}, {4.92, 0.0}},    // Xe
      {{55}, {5.42, 0.0}},    // Cs
      {{56}, {5.07, 0.0}},    // Ba
      {{57}, {8.24, 0.0}},    // La
      {{58}, {4.84, 0.0}},    // Ce
      {{59}, {4.58, 0.0}},    // Pr
      {{60}, {7.69, 0.0}},    // Nd
                              // No PM
      {{61}, {0.80, -1.65}},  // Sm cCNS 0.0
      {{62}, {7.22, -1.26}},  // Eu nCNS 5.3
      {{63}, {6.5, -13.82}},  // Gd nCNS 9.5
      {{64}, {7.38, 0.0}},    // Tb
      {{65}, {16.9, -0.276}}, // Dy
      {{66}, {8.01, 0.0}},    // Ho
      {{67}, {7.79, 0.0}},    // Er
      {{68}, {7.07, 0.0}},    // Tm
      {{69}, {12.43, 0.0}},   // Yb
      {{70}, {7.21, 0.0}},    // Lu
      {{71}, {7.77, 0.0}},    // Hf
      {{72}, {6.91, 0.0}},    // Ta
      {{73}, {4.86, 0.0}},    // W
      {{74}, {9.2, 0.0}},     // Re
      {{75}, {10.7, 0.0}},    // Os
      {{76}, {10.6, 0.0}},    // Ir
      {{77}, {9.60, 0.0}},    // Pt
      {{78}, {7.63, 0.0}},    // Au
      {{79}, {12.692, 0.0}},  // Hg
      {{80}, {8.776, 0.0}},   // Tl
      {{81}, {9.405, 0.0}},   // Pb
      {{82}, {8.532, 0.0}},   // Bi
      {{83}, {10.31, 0.0}},   // Th
      {{84}, {8.417, 0.0}}    // U
  };

  static {
    for (int i = 0; i < atomNames.length; i++) {
      formfactors.put(atomicNumbers[i], neutronFormFactors[i]);
    }
  }

  private final Atom atom;
  private final double[] xyz = new double[3];
  private final double[] dxyz = new double[3];
  private final double[] a = new double[2];
  private final double[] ainv = new double[1];
  private final double[] binv = new double[1];
  private final double[][][] u = new double[1][3][3];
  private final double[][][] uinv = new double[1][3][3];
  private final double[][][] jmat = new double[6][3][3];
  private final double[] gradp = new double[6];
  private final double[] gradu = new double[6];
  private final double[] resv = new double[3];
  private final double[][] resm = new double[3][3];
  private double uadd;
  private double occ;
  private boolean hasAnisou;
  private double[] uaniso = null;

  /**
   * Constructor for NeutronFormFactor.
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   */
  public NeutronFormFactor(Atom atom) {
    this(atom, 0.0, atom.getXYZ(null));
  }

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
   * @param xyz an array of double.
   */
  public NeutronFormFactor(Atom atom, double badd, double[] xyz) {
    this.atom = atom;
    this.uadd = b2u(badd);
    double[][] ffactor;

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

    ffactor = getFormFactor(key);
    arraycopy(ffactor[1], 0, a, 0, ffactor[1].length);
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
   * getFormFactorA
   *
   * @param atom a {@link java.lang.String} object.
   * @return an array of double.
   */
  public static double[] getFormFactorA(String atom) {
    double[][] ffactor = getFormFactor(atom);
    if (ffactor != null) {
      return ffactor[1];
    } else {
      return null;
    }
  }

  /**
   * getFormFactorIndex
   *
   * @param atom a {@link java.lang.String} object.
   * @return a int.
   */
  public static int getFormFactorIndex(String atom) {
    double[][] ffactor = getFormFactor(atom);
    if (ffactor != null) {
      return (int) ffactor[0][0];
    } else {
      return -1;
    }
  }

  /**
   * getFormFactor
   *
   * @param atom a {@link java.lang.String} object.
   * @return an array of double.
   */
  private static double[][] getFormFactor(String atom) {
    double[][] ffactor = null;

    if (formfactors.containsKey(atom)) {
      ffactor = formfactors.get(atom);
    } else {
      String message = " Form factor for atom: " + atom + " not found!";
      logger.severe(message);
    }

    return ffactor;
  }

  /**
   * f
   *
   * @param hkl a {@link ffx.crystal.HKL} object.
   * @return a double.
   */
  public double f(HKL hkl) {
    double sum = a[0] * exp(-twoPI2 * hkl.quadForm(u[0]));
    return occ * sum;
  }

  /** {@inheritDoc} */
  @Override
  public double rho(double f, double lambda, double[] xyz) {
    sub(this.xyz, xyz, xyz);
    double r = length(xyz);
    if (r > atom.getFormFactorWidth()) {
      return f;
    }
    double sum = ainv[0] * exp(-0.5 * quadForm(xyz, uinv[0]));
    return f + (lambda * occ * inverseTwoPI32 * sum);
  }

  /** {@inheritDoc} */
  @Override
  public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementmode) {
    sub(this.xyz, xyz, dxyz);
    double r = length(dxyz);
    double r2 = r * r;
    fill(gradp, 0.0);
    fill(gradu, 0.0);

    if (r > atom.getFormFactorWidth()) {
      return;
    }

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
      if (hasAnisou) {
        refineanisou = true;
      }
    }
    if (refinementmode == RefinementMode.OCCUPANCIES
        || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
        || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
        || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
      refineocc = true;
    }

    double aex = ainv[0] * exp(-0.5 * quadForm(dxyz, uinv[0]));

    if (refinexyz) {
      vec3Mat3(dxyz, uinv[0], resv);
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
        scalarMat3Mat3(-1.0, uinv[0], u11, resm);
        mat3Mat3(resm, uinv[0], jmat[0]);
        scalarMat3Mat3(-1.0, uinv[0], u22, resm);
        mat3Mat3(resm, uinv[0], jmat[1]);
        scalarMat3Mat3(-1.0, uinv[0], u33, resm);
        mat3Mat3(resm, uinv[0], jmat[2]);
        scalarMat3Mat3(-1.0, uinv[0], u12, resm);
        mat3Mat3(resm, uinv[0], jmat[3]);
        scalarMat3Mat3(-1.0, uinv[0], u13, resm);
        mat3Mat3(resm, uinv[0], jmat[4]);
        scalarMat3Mat3(-1.0, uinv[0], u23, resm);
        mat3Mat3(resm, uinv[0], jmat[5]);

        gradu[0] += aex * 0.5 * (-quadForm(dxyz, jmat[0]) - uinv[0][0][0]);
        gradu[1] += aex * 0.5 * (-quadForm(dxyz, jmat[1]) - uinv[0][1][1]);
        gradu[2] += aex * 0.5 * (-quadForm(dxyz, jmat[2]) - uinv[0][2][2]);
        gradu[3] += aex * 0.5 * (-quadForm(dxyz, jmat[3]) - uinv[0][0][1] * 2.0);
        gradu[4] += aex * 0.5 * (-quadForm(dxyz, jmat[4]) - uinv[0][0][2] * 2.0);
        gradu[5] += aex * 0.5 * (-quadForm(dxyz, jmat[5]) - uinv[0][1][2] * 2.0);
      }
    }

    // double rho = occ * inverseTwoPI32 * gradp[3];
    // x, y, z
    if (refinexyz) {
      atom.addToXYZGradient(
          dfc * occ * -inverseTwoPI32 * gradp[0],
          dfc * occ * -inverseTwoPI32 * gradp[1],
          dfc * occ * -inverseTwoPI32 * gradp[2]);
    }

    // occ
    if (refineocc) {
      atom.addToOccupancyGradient(dfc * inverseTwoPI32 * gradp[3]);
    }

    // Biso
    if (refineb) {
      atom.addToTempFactorGradient(dfc * b2u(occ * inverseTwoPI32 * gradp[4]));
      // Uaniso
      if (hasAnisou) {
        for (int i = 0; i < 6; i++) {
          gradu[i] = dfc * occ * inverseTwoPI32 * gradu[i];
        }
        atom.addToAnisouGradient(gradu);
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public void update(double[] xyz) {
    update(xyz, u2b(uadd));
  }

  /** {@inheritDoc} */
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
      double det = determinant3(uaniso);
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
        sb.append("negative B factor for atom: ").append(atom).append("\n");
        sb.append("resetting B to 0.01\n");
        logger.warning(sb.toString());
        atom.setTempFactor(0.01);
        double val = b2u(0.01);
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

    u[0][0][0] = uaniso[0] + uadd;
    u[0][1][1] = uaniso[1] + uadd;
    u[0][2][2] = uaniso[2] + uadd;
    u[0][0][1] = uaniso[3];
    u[0][1][0] = uaniso[3];
    u[0][0][2] = uaniso[4];
    u[0][2][0] = uaniso[4];
    u[0][1][2] = uaniso[5];
    u[0][2][1] = uaniso[5];
    mat3Inverse(u[0], uinv[0]);

    double det = determinant3(u[0]);
    ainv[0] = a[0] / sqrt(det);
    det = determinant3(uinv[0]);
    binv[0] = pow(det,  oneThird);
  }
}
