// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.parameters;

import static ffx.numerics.math.DoubleMath.add;
import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.DoubleMath.normalize;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.potential.parameters.ForceField.ELEC_FORM.FIXED_CHARGE;
import static ffx.potential.parameters.ForceField.ForceFieldType.MULTIPOLE;
import static ffx.utilities.Constants.BOHR;
import static ffx.utilities.Constants.BOHR2;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField.ELEC_FORM;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The MultipoleType class defines a multipole in its local frame.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class MultipoleType extends BaseType implements Comparator<String> {

  /** Constant <code>zeroM</code> */
  public static final double[] zeroM =
      new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  /** Constant <code>zeroD</code> */
  public static final double[] zeroD = new double[] {0.0, 0.0, 0.0};
  /** Constant <code>chrg=t000</code> */
  public static final int t000 = 0;
  /** Constant <code>t100=1</code> */
  public static final int t100 = 1;
  /** Constant <code>t010=2</code> */
  public static final int t010 = 2;
  /** Constant <code>t001=3</code> */
  public static final int t001 = 3;
  /** Constant <code>t200=4</code> */
  public static final int t200 = 4;
  /** Constant <code>t020=5</code> */
  public static final int t020 = 5;
  /** Constant <code>t002=6</code> */
  public static final int t002 = 6;
  /** Constant <code>t110=7</code> */
  public static final int t110 = 7;
  /** Constant <code>t101=8</code> */
  public static final int t101 = 8;
  /** Constant <code>t011=9</code> */
  public static final int t011 = 9;
  /** Constant <code>t300=10</code> */
  public static final int t300 = 10;
  /** Constant <code>t030=11</code> */
  public static final int t030 = 11;
  /** Constant <code>t003=12</code> */
  public static final int t003 = 12;
  /** Constant <code>t210=13</code> */
  public static final int t210 = 13;
  /** Constant <code>t201=14</code> */
  public static final int t201 = 14;
  /** Constant <code>t120=15</code> */
  public static final int t120 = 15;
  /** Constant <code>t021=16</code> */
  public static final int t021 = 16;
  /** Constant <code>t102=17</code> */
  public static final int t102 = 17;
  /** Constant <code>t012=18</code> */
  public static final int t012 = 18;
  /** Constant <code>t111=19</code> */
  public static final int t111 = 19;

  private static final Logger logger = Logger.getLogger(MultipoleType.class.getName());
  /** Partial atomic charge (e). */
  public final double charge;
  /** Atomic dipole. 1 x 3 (e Angstroms). */
  public final double[] dipole;
  /** Atomic quadrupole. 3 x 3 (e Angstroms^2). */
  public final double[][] quadrupole;
  /** Local frame definition method. */
  public final MultipoleFrameDefinition frameDefinition;
  /** Atom types that define the local frame of this multipole. */
  public final int[] frameAtomTypes;
  /**
   * Charge, dipole, and quadrupole packed into tensor notation: c, dx, dy, dz, qxx, qyy, qzz, qxy,
   * qxz, qyz
   */
  private final double[] multipole;

  /**
   * Multipole Constructor. Conversion to electron Angstroms should be requested only when reading
   * multipole values from the force field file.
   *
   * @param multipole an array of {@link double} objects.
   * @param frameAtomTypes an array of {@link int} objects.
   * @param frameDefinition a {@link
   *     ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition} object.
   * @param convertFromBohr a boolean.
   */
  public MultipoleType(
      double[] multipole,
      int[] frameAtomTypes,
      MultipoleFrameDefinition frameDefinition,
      boolean convertFromBohr) {
    super(MULTIPOLE, frameAtomTypes);
    this.multipole = (convertFromBohr) ? bohrToElectronAngstroms(multipole) : multipole;
    this.frameAtomTypes = frameAtomTypes;
    this.frameDefinition = frameDefinition;
    charge = multipole[t000];
    dipole = unpackDipole(multipole);
    quadrupole = unpackQuad(multipole);
    checkMultipole();
  }

  /**
   * Constructor for MultipoleType.
   *
   * @param charge a double.
   * @param dipole an array of {@link double} objects.
   * @param quadrupole an array of {@link double} objects.
   * @param frameAtomTypes an array of {@link int} objects.
   * @param frameDefinition a {@link
   *     ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition} object.
   * @param convertFromBohr a boolean.
   */
  public MultipoleType(
      double charge,
      double[] dipole,
      double[][] quadrupole,
      int[] frameAtomTypes,
      MultipoleFrameDefinition frameDefinition,
      boolean convertFromBohr) {
    super(MULTIPOLE, frameAtomTypes);
    this.charge = charge;
    if (convertFromBohr) {
      this.multipole = bohrToElectronAngstroms(pack(charge, dipole, quadrupole));
      this.dipole = unpackDipole(multipole);
      this.quadrupole = unpackQuad(multipole);
    } else {
      this.multipole = pack(charge, dipole, quadrupole);
      this.dipole = dipole;
      this.quadrupole = quadrupole;
    }
    this.frameAtomTypes = frameAtomTypes;
    this.frameDefinition = frameDefinition;
    checkMultipole();
  }

  /**
   * assignMultipole.
   *
   * @param elecForm a {@link ffx.potential.parameters.ForceField.ELEC_FORM} object.
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param multipole an array of {@link double} objects.
   * @param i a int.
   * @param axisAtom an array of {@link int} objects.
   * @param frame an array of {@link
   *     ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition} objects.
   * @return a boolean.
   */
  public static boolean assignMultipole(
      ELEC_FORM elecForm,
      Atom atom,
      ForceField forceField,
      double[] multipole,
      int i,
      int[][] axisAtom,
      MultipoleFrameDefinition[] frame) {
    MultipoleType type = multipoleTypeFactory(elecForm, atom, forceField);
    if (type == null) {
      return false;
    }
    arraycopy(type.getMultipole(), 0, multipole, 0, 10);
    axisAtom[i] = atom.getAxisAtomIndices();
    frame[i] = atom.getMultipoleType().frameDefinition;
    return true;
  }

  /**
   * Average two MultipoleType instances. The atom types that define the frame of the new type must
   * be supplied.
   *
   * @param multipoleType1 a {@link ffx.potential.parameters.MultipoleType} object.
   * @param multipoleType2 a {@link ffx.potential.parameters.MultipoleType} object.
   * @param multipoleFrameTypes an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.MultipoleType} object.
   */
  public static MultipoleType averageTypes(
      MultipoleType multipoleType1, MultipoleType multipoleType2, int[] multipoleFrameTypes) {
    if (multipoleType1.frameDefinition != multipoleType2.frameDefinition) {
      return null;
    }
    MultipoleType[] types = {multipoleType1, multipoleType2};
    double[] weights = {0.5, 0.5};
    double[] averagedMultipole = weightMultipole(types, weights);
    if (averagedMultipole == null) {
      return null;
    }
    return new MultipoleType(
        averagedMultipole, multipoleFrameTypes, multipoleType1.frameDefinition, false);
  }

  /**
   * checkMultipoleChirality.
   *
   * @param multipole an array of {@link double} objects.
   * @param frame a {@link ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition} object.
   * @param frameCoords an array of {@link double} objects.
   * @param localOrigin an array of {@link double} objects.
   * @return Whether this multipole underwent chiral inversion.
   */
  public static boolean checkMultipoleChirality(
      double[] multipole,
      MultipoleFrameDefinition frame,
      double[] localOrigin,
      double[][] frameCoords) {
    if (frame != MultipoleFrameDefinition.ZTHENX) {
      return false;
    }
    double[] zAxis = new double[3];
    double[] xAxis = new double[3];
    double[] yAxis = new double[3];
    double[] yMinOrigin = new double[3];
    zAxis[0] = frameCoords[0][0];
    zAxis[1] = frameCoords[0][1];
    zAxis[2] = frameCoords[0][2];
    xAxis[0] = frameCoords[1][0];
    xAxis[1] = frameCoords[1][1];
    xAxis[2] = frameCoords[1][2];
    yAxis[0] = frameCoords[2][0];
    yAxis[1] = frameCoords[2][1];
    yAxis[2] = frameCoords[2][2];
    sub(localOrigin, yAxis, yMinOrigin);
    sub(zAxis, yAxis, zAxis);
    sub(xAxis, yAxis, xAxis);
    double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
    double c2 = xAxis[1] * yMinOrigin[2] - xAxis[2] * yMinOrigin[1];
    double c3 = yMinOrigin[1] * zAxis[2] - yMinOrigin[2] * zAxis[1];
    double vol = yMinOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
    return (vol < 0.0);
  }

  /**
   * Return the rotation matrix for the local to lab frame.
   *
   * @param frame the multipole frame definition
   * @param frameCoords the coordinates of the frame atoms
   * @param localOrigin the local origin of the frame
   * @return the rotation matrix
   */
  public static double[][] getRotationMatrix(
      MultipoleFrameDefinition frame, double[] localOrigin, double[][] frameCoords) {
    double[][] rotmat = new double[3][3];
    double[] zAxis = new double[3];
    double[] xAxis = new double[3];
    switch (frame) {
      case BISECTOR:
        zAxis[0] = frameCoords[0][0];
        zAxis[1] = frameCoords[0][1];
        zAxis[2] = frameCoords[0][2];
        xAxis[0] = frameCoords[1][0];
        xAxis[1] = frameCoords[1][1];
        xAxis[2] = frameCoords[1][2];
        sub(zAxis, localOrigin, zAxis);
        normalize(zAxis, zAxis);
        sub(xAxis, localOrigin, xAxis);
        normalize(xAxis, xAxis);
        add(xAxis, zAxis, zAxis);
        normalize(zAxis, zAxis);
        rotmat[0][2] = zAxis[0];
        rotmat[1][2] = zAxis[1];
        rotmat[2][2] = zAxis[2];
        double dot = dot(xAxis, zAxis);
        DoubleMath.scale(zAxis, dot, zAxis);
        sub(xAxis, zAxis, xAxis);
        normalize(xAxis, xAxis);
        break;
      case ZTHENBISECTOR:
        double[] yAxis = new double[3];
        zAxis[0] = frameCoords[0][0];
        zAxis[1] = frameCoords[0][1];
        zAxis[2] = frameCoords[0][2];
        xAxis[0] = frameCoords[1][0];
        xAxis[1] = frameCoords[1][1];
        xAxis[2] = frameCoords[1][2];
        yAxis[0] = frameCoords[2][0];
        yAxis[1] = frameCoords[2][1];
        yAxis[2] = frameCoords[2][2];
        sub(zAxis, localOrigin, zAxis);
        normalize(zAxis, zAxis);
        rotmat[0][2] = zAxis[0];
        rotmat[1][2] = zAxis[1];
        rotmat[2][2] = zAxis[2];
        sub(xAxis, localOrigin, xAxis);
        normalize(xAxis, xAxis);
        sub(yAxis, localOrigin, yAxis);
        normalize(yAxis, yAxis);
        add(xAxis, yAxis, xAxis);
        normalize(xAxis, xAxis);
        dot = dot(xAxis, zAxis);
        DoubleMath.scale(zAxis, dot, zAxis);
        sub(xAxis, zAxis, xAxis);
        normalize(xAxis, xAxis);
        break;
      case ZONLY:
        zAxis[0] = frameCoords[0][0];
        zAxis[1] = frameCoords[0][1];
        zAxis[2] = frameCoords[0][2];
        sub(zAxis, localOrigin, zAxis);
        normalize(zAxis, zAxis);
        rotmat[0][2] = zAxis[0];
        rotmat[1][2] = zAxis[1];
        rotmat[2][2] = zAxis[2];
        xAxis[0] = random();
        xAxis[1] = random();
        xAxis[2] = random();
        dot = dot(xAxis, zAxis);
        DoubleMath.scale(zAxis, dot, zAxis);
        sub(xAxis, zAxis, xAxis);
        normalize(xAxis, xAxis);
        break;
      case ZTHENX:
      default:
        zAxis[0] = frameCoords[0][0];
        zAxis[1] = frameCoords[0][1];
        zAxis[2] = frameCoords[0][2];
        xAxis[0] = frameCoords[1][0];
        xAxis[1] = frameCoords[1][1];
        xAxis[2] = frameCoords[1][2];
        sub(zAxis, localOrigin, zAxis);
        normalize(zAxis, zAxis);
        rotmat[0][2] = zAxis[0];
        rotmat[1][2] = zAxis[1];
        rotmat[2][2] = zAxis[2];
        sub(xAxis, localOrigin, xAxis);
        dot = dot(xAxis, zAxis);
        DoubleMath.scale(zAxis, dot, zAxis);
        sub(xAxis, zAxis, xAxis);
        normalize(xAxis, xAxis);
    }
    // Set the X elements.
    rotmat[0][0] = xAxis[0];
    rotmat[1][0] = xAxis[1];
    rotmat[2][0] = xAxis[2];
    // Set the Y elements.
    rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0] * rotmat[2][2];
    rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0] * rotmat[0][2];
    rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0] * rotmat[1][2];
    return rotmat;
  }

  /**
   * invertMultipoleChirality.
   *
   * @param mpole an array of {@link double} objects.
   */
  public static void invertMultipoleChirality(double[] mpole) {
    mpole[t010] = -mpole[t010];
    mpole[t110] = -mpole[t110];
    mpole[t011] = -mpole[t011];
  }

  /**
   * multipoleTypeFactory.
   *
   * @param elecForm   The electrostatics form being used.
   * @param atom       a {@link ffx.potential.bonded.Atom} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @return a {@link ffx.potential.parameters.MultipoleType} object.
   */
  public static MultipoleType multipoleTypeFactory(
      ELEC_FORM elecForm, Atom atom, ForceField forceField) {
    AtomType atomType = atom.getAtomType();
    if (atomType == null) {
      String message = " Multipoles can only be assigned to atoms that have been typed.";
      logger.severe(message);
      return null;
    }

    if (elecForm != FIXED_CHARGE) {
      PolarizeType polarizeType = forceField.getPolarizeType(atomType.getKey());
      if (polarizeType != null) {
        atom.setPolarizeType(polarizeType);
      } else {
        String message = " No polarization type was found for " + atom.toString();
        logger.info(message);
        double polarizability = 0.0;
        double thole = 0.0;
        int[] polarizationGroup = null;
        polarizeType = new PolarizeType(atomType.type, polarizability, thole, polarizationGroup);
        forceField.addForceFieldType(polarizeType);
        atom.setPolarizeType(polarizeType);
      }
    }

    // No reference atoms.
    String key1 = atomType.getKey();
    MultipoleType multipoleType = forceField.getMultipoleType(key1);
    if (multipoleType != null) {
      atom.setMultipoleType(multipoleType);
      atom.setAxisAtoms((Atom[]) null);
      return multipoleType;
    }

    // List sorted of 1-2 iteractions.
    List<Atom> n12 = atom.get12List();

    // No bonds.
    if (n12 == null || n12.size() < 1) {
      String message = "Multipoles can only be assigned after bonded relationships are defined.\n";
      logger.severe(message);
      return null;
    }

    // Only 1-2 connected atoms: 1 reference atom.
    for (Atom atom2 : n12) {
      String key = key1 + " " + atom2.getAtomType().getKey();
      multipoleType = forceField.getMultipoleType(key);
      if (multipoleType != null) {
        atom.setMultipoleType(multipoleType);
        atom.setAxisAtoms(atom2);
        return multipoleType;
      }
    }

    // Only 1-2 connected atoms: 2 reference atoms.
    for (Atom atom2 : n12) {
      String key2 = atom2.getAtomType().getKey();
      for (Atom atom3 : n12) {
        if (atom2 == atom3) {
          continue;
        }
        String key3 = atom3.getAtomType().getKey();
        String key = key1 + " " + key2 + " " + key3;
        multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
          atom.setMultipoleType(multipoleType);
          atom.setAxisAtoms(atom2, atom3);
          return multipoleType;
        }
      }
    }

    // Only 1-2 connected atoms: 3 reference atoms.
    for (Atom atom2 : n12) {
      String key2 = atom2.getAtomType().getKey();
      for (Atom atom3 : n12) {
        if (atom2 == atom3) {
          continue;
        }
        String key3 = atom3.getAtomType().getKey();
        for (Atom atom4 : n12) {
          if (atom2 == atom4 || atom3 == atom4) {
            continue;
          }
          String key4 = atom4.getAtomType().getKey();
          String key = key1 + " " + key2 + " " + key3 + " " + key4;
          multipoleType = forceField.getMultipoleType(key);
          if (multipoleType != null) {
            atom.setMultipoleType(multipoleType);
            atom.setAxisAtoms(atom2, atom3, atom4);
            return multipoleType;
          }
        }
      }
    }

    List<Atom> n13 = atom.get13List();

    // Revert to a reference atom definition that may include a 1-3 site. For example a hydrogen on
    // water.
    for (Atom atom2 : n12) {
      String key2 = atom2.getAtomType().getKey();
      for (Atom atom3 : n13) {
        String key3 = atom3.getAtomType().getKey();
        String key = key1 + " " + key2 + " " + key3;
        multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
          atom.setMultipoleType(multipoleType);
          atom.setAxisAtoms(atom2, atom3);
          return multipoleType;
        }
        for (Atom atom4 : n13) {
          if (atom4 != null && atom4 != atom3) {
            String key4 = atom4.getAtomType().getKey();
            key = key1 + " " + key2 + " " + key3 + " " + key4;
            multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
              atom.setMultipoleType(multipoleType);
              atom.setAxisAtoms(atom2, atom3, atom4);
              return multipoleType;
            }
          }
        }
      }
    }

    return null;
  }

  /**
   * Parse a single line multipole.
   *
   * @param input Input String.
   * @param tokens Input tokens.
   * @param br A BufferedReader instance.
   * @return a MultipoleType instance.
   * @since 1.0
   */
  public static MultipoleType parse(String input, String[] tokens, BufferedReader br) {
    if (tokens == null || tokens.length < 3 || tokens.length > 6) {
      logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
      return null;
    }

    try {
      int nTokens = tokens.length;
      ArrayList<Integer> frameAtoms = new ArrayList<>();
      // Loop over integer tokens (i.e. except for the initial 'multipole' keyword and final charge
      // value).
      for (int i = 1; i < nTokens - 1; i++) {
        int frameType = parseInt(tokens[i]);
        // Ignore atom types of '0'.
        if (frameType != 0) {
          frameAtoms.add(frameType);
        }
      }
      int nAtomTypes = frameAtoms.size();
      int[] atomTypes = new int[nAtomTypes];
      for (int i = 0; i < nAtomTypes; i++) {
        atomTypes[i] = frameAtoms.get(i);
      }
      // Last token is the monopole.
      double charge = parseDouble(tokens[nTokens - 1]);

      MultipoleFrameDefinition frameDefinition = null;
      if (nAtomTypes == 1) {
        frameDefinition = MultipoleFrameDefinition.NONE;
      } else if (nAtomTypes == 2) {
        frameDefinition = MultipoleFrameDefinition.ZONLY;
      } else if (nAtomTypes == 3) {
        // ZTHENX or BISECTOR
        if (atomTypes[1] < 0 || atomTypes[2] < 0) {
          frameDefinition = MultipoleFrameDefinition.BISECTOR;
        } else {
          frameDefinition = MultipoleFrameDefinition.ZTHENX;
        }
      } else if (nAtomTypes == 4) {
        // ZTHENBISECTOR or THREEFOLD
        if (atomTypes[2] < 0 && atomTypes[3] < 0) {
          frameDefinition = MultipoleFrameDefinition.ZTHENBISECTOR;
          if (atomTypes[1] < 0) {
            frameDefinition = MultipoleFrameDefinition.THREEFOLD;
          }
        }
      }

      // Warn if the frame definition is ambiguous.
      if (frameDefinition == null) {
        logger.log(Level.FINE, "Ambiguous MULTIPOLE type:\n{0}", input);
        frameDefinition = MultipoleFrameDefinition.ZTHENX;
      }

      for (int i = 0; i < nAtomTypes; i++) {
        atomTypes[i] = abs(atomTypes[i]);
      }

      input = br.readLine().split("#")[0];
      tokens = input.trim().split(" +");
      if (tokens.length != 3) {
        logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
        return null;
      }
      double[] dipole = new double[3];
      dipole[0] = parseDouble(tokens[0]);
      dipole[1] = parseDouble(tokens[1]);
      dipole[2] = parseDouble(tokens[2]);
      input = br.readLine().split("#")[0];
      tokens = input.trim().split(" +");
      if (tokens.length != 1) {
        logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
        return null;
      }
      double[][] quadrupole = new double[3][3];
      quadrupole[0][0] = parseDouble(tokens[0]);
      input = br.readLine().split("#")[0];
      tokens = input.trim().split(" +");
      if (tokens.length != 2) {
        logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
        return null;
      }
      quadrupole[1][0] = parseDouble(tokens[0]);
      quadrupole[1][1] = parseDouble(tokens[1]);
      input = br.readLine().split("#")[0];
      tokens = input.trim().split(" +");
      if (tokens.length != 3) {
        logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
        return null;
      }
      quadrupole[2][0] = parseDouble(tokens[0]);
      quadrupole[2][1] = parseDouble(tokens[1]);
      quadrupole[2][2] = parseDouble(tokens[2]);
      // Fill in symmetric components.
      quadrupole[0][1] = quadrupole[1][0];
      quadrupole[0][2] = quadrupole[2][0];
      quadrupole[1][2] = quadrupole[2][1];
      return new MultipoleType(charge, dipole, quadrupole, atomTypes, frameDefinition, true);
    } catch (Exception e) {
      String message = "Exception parsing MULTIPOLE type:\n" + input + "\n";
      logger.log(Level.SEVERE, message, e);
    }
    return null;
  }

  /**
   * Parse a single line multipole.
   *
   * @param input Input String.
   * @param tokens Input tokens.
   * @return a MultipoleType instance.
   * @since 1.0
   */
  public static MultipoleType parse(String input, String[] tokens) {
    if (tokens == null || tokens.length < 12 || tokens.length > 15) {
      logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
      return null;
    }

    try {
      int nTokens = tokens.length;
      ArrayList<Integer> frameAtoms = new ArrayList<>();
      // Loop over integer tokens (i.e. except for the initial 'multipole' keyword and final 10
      // multipole values).
      for (int i = 1; i < nTokens - 10; i++) {
        int frameType = parseInt(tokens[i]);
        // Ignore atom types of '0'.
        if (frameType != 0) {
          frameAtoms.add(frameType);
        }
      }
      int nAtomTypes = frameAtoms.size();
      int[] atomTypes = new int[nAtomTypes];
      for (int i = 0; i < nAtomTypes; i++) {
        atomTypes[i] = frameAtoms.get(i);
      }

      MultipoleFrameDefinition frameDefinition = null;
      if (nAtomTypes == 1) {
        frameDefinition = MultipoleFrameDefinition.NONE;
      } else if (nAtomTypes == 2) {
        frameDefinition = MultipoleFrameDefinition.ZONLY;
      } else if (nAtomTypes == 3) {
        // ZTHENX or BISECTOR
        if (atomTypes[1] < 0 || atomTypes[2] < 0) {
          frameDefinition = MultipoleFrameDefinition.BISECTOR;
        } else {
          frameDefinition = MultipoleFrameDefinition.ZTHENX;
        }
      } else if (nAtomTypes == 4) {
        // ZTHENBISECTOR or THREEFOLD
        if (atomTypes[2] < 0 && atomTypes[3] < 0) {
          frameDefinition = MultipoleFrameDefinition.ZTHENBISECTOR;
          if (atomTypes[1] < 0) {
            frameDefinition = MultipoleFrameDefinition.THREEFOLD;
          }
        }
      }

      // Warn if the frame definition is ambiguous.
      if (frameDefinition == null) {
        logger.log(Level.FINE, "Ambiguous MULTIPOLE type:\n{0}", input);
        frameDefinition = MultipoleFrameDefinition.ZTHENX;
      }

      for (int i = 0; i < nAtomTypes; i++) {
        atomTypes[i] = abs(atomTypes[i]);
      }

      double[] dipole = new double[3];
      double[][] quadrupole = new double[3][3];
      double charge = parseDouble(tokens[nTokens - 10]);
      dipole[0] = parseDouble(tokens[nTokens - 9]);
      dipole[1] = parseDouble(tokens[nTokens - 8]);
      dipole[2] = parseDouble(tokens[nTokens - 7]);
      quadrupole[0][0] = parseDouble(tokens[nTokens - 6]);
      quadrupole[1][0] = parseDouble(tokens[nTokens - 5]);
      quadrupole[1][1] = parseDouble(tokens[nTokens - 4]);
      quadrupole[2][0] = parseDouble(tokens[nTokens - 3]);
      quadrupole[2][1] = parseDouble(tokens[nTokens - 2]);
      quadrupole[2][2] = parseDouble(tokens[nTokens - 1]);
      // Fill in symmetric components.
      quadrupole[0][1] = quadrupole[1][0];
      quadrupole[0][2] = quadrupole[2][0];
      quadrupole[1][2] = quadrupole[2][1];
      return new MultipoleType(charge, dipole, quadrupole, atomTypes, frameDefinition, true);
    } catch (Exception e) {
      String message = "Exception parsing MULTIPOLE type:\n" + input + "\n";
      logger.log(Level.SEVERE, message, e);
    }
    return null;
  }

  /**
   * Map charge parameters to a Multipole instance.
   *
   * @param input Input string.
   * @param tokens Input string tokens.
   * @return a MultipoleType instance
   */
  public static MultipoleType parseChargeType(String input, String[] tokens) {
    if (tokens.length < 3) {
      logger.log(Level.WARNING, "Invalid CHARGE type:\n{0}", input);
    } else {
      try {
        int[] atomTypes = new int[] {parseInt(tokens[1]), 0, 0};
        double charge = parseDouble(tokens[2]);
        double[] dipole = new double[3];
        double[][] quadrupole = new double[3][3];
        MultipoleFrameDefinition frameDefinition = MultipoleFrameDefinition.NONE;
        return new MultipoleType(charge, dipole, quadrupole, atomTypes, frameDefinition, true);
      } catch (NumberFormatException e) {
        String message = "Exception parsing CHARGE type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * rotateDipole.
   *
   * @param rotmat an array of {@link double} objects.
   * @param dipole an array of {@link double} objects.
   * @param rotatedDipole an array of {@link double} objects.
   */
  public static void rotateDipole(double[][] rotmat, double[] dipole, double[] rotatedDipole) {
    for (int i = 0; i < 3; i++) {
      double[] rotmati = rotmat[i];
      for (int j = 0; j < 3; j++) {
        rotatedDipole[i] += rotmati[j] * dipole[j];
      }
    }
  }

  /**
   * rotateMultipole.
   *
   * @param rotmat an array of {@link double} objects.
   * @param dipole an array of {@link double} objects.
   * @param quadrupole an array of {@link double} objects.
   * @param rotatedDipole an array of {@link double} objects.
   * @param rotatedQuadrupole an array of {@link double} objects.
   */
  public static void rotateMultipole(
      double[][] rotmat,
      double[] dipole,
      double[][] quadrupole,
      double[] rotatedDipole,
      double[][] rotatedQuadrupole) {
    for (int i = 0; i < 3; i++) {
      double[] rotmati = rotmat[i];
      double[] quadrupolei = rotatedQuadrupole[i];
      for (int j = 0; j < 3; j++) {
        double[] rotmatj = rotmat[j];
        rotatedDipole[i] += rotmati[j] * dipole[j];
        if (j < i) {
          quadrupolei[j] = rotatedQuadrupole[j][i];
        } else {
          for (int k = 0; k < 3; k++) {
            double[] localQuadrupolek = quadrupole[k];
            quadrupolei[j] +=
                rotmati[k]
                    * (rotmatj[0] * localQuadrupolek[0]
                        + rotmatj[1] * localQuadrupolek[1]
                        + rotmatj[2] * localQuadrupolek[2]);
          }
        }
      }
    }
  }

  /**
   * scale.
   *
   * @param type a {@link ffx.potential.parameters.MultipoleType} object.
   * @param cdtScales an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  public static double[] scale(MultipoleType type, double[] cdtScales) {
    return scale(type.getMultipole(), cdtScales);
  }

  /**
   * scale.
   *
   * @param multipole an array of {@link double} objects.
   * @param cdtScales an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  public static double[] scale(double[] multipole, double[] cdtScales) {
    double chargeScale = cdtScales[0];
    double dipoleScale = cdtScales[1];
    double quadScale = cdtScales[2];
    return new double[] {
      multipole[t000] * chargeScale,
      multipole[t100] * dipoleScale,
      multipole[t010] * dipoleScale,
      multipole[t001] * dipoleScale,
      multipole[t200] * quadScale,
      multipole[t020] * quadScale,
      multipole[t002] * quadScale,
      multipole[t110] * quadScale,
      multipole[t101] * quadScale,
      multipole[t011] * quadScale
    };
  }

  /* Indices into a 1D tensor array based on compressed tensor notation. This makes multipole code easier to read. */

  /**
   * weightMultipoleTypes.
   *
   * @param types an array of {@link ffx.potential.parameters.MultipoleType} objects.
   * @param weights an array of {@link double} objects.
   * @param frameAtomTypes an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.MultipoleType} object.
   */
  public static MultipoleType weightMultipoleTypes(
      MultipoleType[] types, double[] weights, int[] frameAtomTypes) {
    double[] weightedMultipole = weightMultipole(types, weights);
    if (weightedMultipole == null) {
      return null;
    }
    return new MultipoleType(weightedMultipole, frameAtomTypes, types[0].frameDefinition, false);
  }

  private static double[] bohrToElectronAngstroms(double[] multipole) {
    return new double[] {
      multipole[t000],
      multipole[t100] *= BOHR,
      multipole[t010] *= BOHR,
      multipole[t001] *= BOHR,
      multipole[t200] *= BOHR2,
      multipole[t020] *= BOHR2,
      multipole[t002] *= BOHR2,
      multipole[t110] *= BOHR2,
      multipole[t101] *= BOHR2,
      multipole[t011] *= BOHR2
    };
  }

  /**
   * Pack charge, dipole, quad into 1d tensor array form.
   *
   * @param charge a double.
   * @param dipole an array of {@link double} objects.
   * @param quad an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  private static double[] pack(double charge, double[] dipole, double[][] quad) {
    return new double[] {
      charge,
      dipole[0],
      dipole[1],
      dipole[2],
      quad[0][0],
      quad[1][1],
      quad[2][2],
      quad[0][1],
      quad[0][2],
      quad[1][2]
    };
  }

  /** Unpack dipole from 1d tensor-form multipole. */
  private static double[] unpackDipole(double[] mpole) {
    return new double[] {mpole[t100], mpole[t010], mpole[t001]};
  }

  /** Unpack quadrupole from 1d tensor-form multipole. */
  private static double[][] unpackQuad(double[] mpole) {
    return new double[][] {
      {mpole[t200], mpole[t110], mpole[t101]},
      {mpole[t110], mpole[t020], mpole[t011]},
      {mpole[t101], mpole[t011], mpole[t002]}
    };
  }

  /**
   * Create a new multipole representing a weighted average.
   *
   * @param types an array of {@link ffx.potential.parameters.MultipoleType} objects.
   * @param weights an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  private static double[] weightMultipole(MultipoleType[] types, double[] weights) {
    if (types == null || weights == null || types.length != weights.length) {
      throw new IllegalArgumentException();
    }
    if (Arrays.asList(types).contains(null)) {
      // Multipoles have not yet been assigned.
      return null;
    }
    for (MultipoleType type : types) {
      if (type.frameDefinition != types[0].frameDefinition) {
        logger.warning(
            format(
                "Multipole frame definition mismatch during weighting:\n\t%s->%s,\n\t%s->%s",
                types[0].toString(),
                types[0].frameDefinition.toString(),
                type.toString(),
                type.frameDefinition.toString()));
        throw new IllegalArgumentException();
      }
    }
    double[] weightedMultipole = new double[10];
    fill(weightedMultipole, 0.0);
    for (int idx = 0; idx < types.length; idx++) {
      double[] multipole = types[idx].getMultipole();
      for (int comp = 0; comp < 10; comp++) {
        weightedMultipole[comp] += weights[idx] * multipole[comp];
      }
    }
    return weightedMultipole;
  }

  /** {@inheritDoc} */
  @Override
  public int compare(String s1, String s2) {
    String[] keys1 = s1.split(" ");
    String[] keys2 = s2.split(" ");

    int len = keys1.length;
    if (keys1.length > keys2.length) {
      len = keys2.length;
    }
    int[] c1 = new int[len];
    int[] c2 = new int[len];
    for (int i = 0; i < len; i++) {
      c1[i] = abs(parseInt(keys1[i]));
      c2[i] = abs(parseInt(keys2[i]));
      if (c1[i] < c2[i]) {
        return -1;
      } else if (c1[i] > c2[i]) {
        return 1;
      }
    }

    if (keys1.length < keys2.length) {
      return -1;
    } else if (keys1.length > keys2.length) {
      return 1;
    }

    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    MultipoleType multipoleType = (MultipoleType) o;
    return Arrays.equals(frameAtomTypes, multipoleType.frameAtomTypes);
  }

  /**
   * Getter for the field <code>charge</code>.
   *
   * @return An uneditable copy of this type's charge. To make changes, use getMultipoleReference().
   */
  public double getCharge() {
    return multipole[t000];
  }

  /**
   * Getter for the field <code>dipole</code>.
   *
   * @return An uneditable copy of this type's dipole. To make changes, use getMultipoleReference().
   */
  public double[] getDipole() {
    return new double[] {multipole[t100], multipole[t010], multipole[t001]};
  }

  /**
   * Getter for the field <code>multipole</code>.
   *
   * @return An uneditable copy of this type's multipole. To make changes, use
   *     getMultipoleReference().
   */
  public double[] getMultipole() {
    return new double[] {
      multipole[t000],
      multipole[t100],
      multipole[t010],
      multipole[t001],
      multipole[t200],
      multipole[t020],
      multipole[t002],
      multipole[t110],
      multipole[t101],
      multipole[t011]
    };
  }

  /**
   * Getter for the field <code>quadrupole</code>.
   *
   * @return An uneditable copy of this type's quadrupole. To make changes, use
   *     getMultipoleReference().
   */
  public double[][] getQuadrupole() {
    return new double[][] {
      {multipole[t200], multipole[t110], multipole[t101]},
      {multipole[t110], multipole[t020], multipole[t011]},
      {multipole[t101], multipole[t011], multipole[t002]}
    };
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Arrays.hashCode(frameAtomTypes);
  }

  /**
   * Nicely formatted multipole string. Dipole and qaudrupole are in electron-Bohrs and
   * electron-Bohrs^2, respectively.
   *
   * @return String
   */
  public String toCompactBohrString() {
    StringBuilder multipoleBuffer = new StringBuilder("mpol ");
    switch (frameDefinition) {
      case NONE:
        multipoleBuffer.append(format("(%3d):", frameAtomTypes[0]));
        break;
      case ZONLY:
        multipoleBuffer.append(format("(%3d,%3d):", frameAtomTypes[0], frameAtomTypes[1]));
        break;
      case ZTHENX:
        multipoleBuffer.append(
            format("(%3d,%3d,%3d):", frameAtomTypes[0], frameAtomTypes[1], frameAtomTypes[2]));
        break;
      case BISECTOR:
        multipoleBuffer.append(
            format("(%3d,%3d,%3d):", frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2]));
        break;
      case ZTHENBISECTOR:
        multipoleBuffer.append(
            format(
                "(%3d,%3d,%3d,%3d):",
                frameAtomTypes[0], frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
        break;
      case THREEFOLD:
        multipoleBuffer.append(
            format(
                "(%3d,%3d,%3d,%3d):",
                frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
    }
    multipoleBuffer.append(
        format(
            "[%6.3f / %6.3f %6.3f %6.3f / %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]",
            charge,
            dipole[0] / BOHR,
            dipole[1] / BOHR,
            dipole[2] / BOHR,
            quadrupole[0][0] / BOHR2,
            quadrupole[1][0] / BOHR2,
            quadrupole[1][1] / BOHR2,
            quadrupole[2][0] / BOHR2,
            quadrupole[2][1] / BOHR2,
            quadrupole[2][2] / BOHR2));
    return multipoleBuffer.toString();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return toBohrString();
  }

  /**
   * incrementType
   *
   * @param increment a int.
   */
  void incrementType(int increment) {
    for (int i = 0; i < frameAtomTypes.length; i++) {
      // Frame atom types of 0 are unchanged.
      if (frameAtomTypes[i] > 0) {
        frameAtomTypes[i] += increment;
      } else if (frameAtomTypes[i] < 0) {
        frameAtomTypes[i] -= increment;
      }
    }
    setKey(frameAtomTypes);
  }

  /**
   * Remap new atom types to known internal ones.
   *
   * @param typeMap a lookup between new atom types and known atom types.
   * @return a {@link ffx.potential.parameters.MultipoleType} object.
   */
  MultipoleType patchTypes(HashMap<AtomType, AtomType> typeMap) {
    int count = 0;
    int len = frameAtomTypes.length;

    // Look for a MultipoleType that contain a mapped atom class.
    for (AtomType newType : typeMap.keySet()) {
      for (int frameAtomType : frameAtomTypes) {
        if (frameAtomType == newType.type || frameAtomType == 0) {
          count++;
        }
      }
    }

    // If found, create a new MultipoleType that bridges to known classes.
    if (count > 0 && count < len) {
      int[] newFrame = Arrays.copyOf(frameAtomTypes, len);
      for (AtomType newType : typeMap.keySet()) {
        for (int i = 0; i < len; i++) {
          if (frameAtomTypes[i] == newType.type) {
            AtomType knownType = typeMap.get(newType);
            newFrame[i] = knownType.type;
          }
        }
      }
      return new MultipoleType(multipole, newFrame, frameDefinition, false);
    }
    return null;
  }

  private void checkMultipole() {
    double[][] quadrupole = unpackQuad(multipole);
    // Check symmetry.
    if (abs(quadrupole[0][1] - quadrupole[1][0]) > 1.0e-6) {
      logger.warning("Multipole component Qxy != Qyx");
      logger.info(this.toString());
    }
    if (abs(quadrupole[0][2] - quadrupole[2][0]) > 1.0e-6) {
      logger.warning("Multipole component Qxz != Qzx");
      logger.info(this.toString());
    }
    if (abs(quadrupole[1][2] - quadrupole[2][1]) > 1.0e-6) {
      logger.warning("Multipole component Qyz != Qzy");
      logger.info(this.toString());
    }
    // Warn if the multipole is not traceless.
    if (abs(quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2]) > 1.0e-5) {
      logger.log(
          Level.WARNING,
          format(
              "Multipole is not traceless: %12.8f",
              abs(quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2])));
      logger.info(this.toString());
    }
  }

  /**
   * Nicely formatted multipole string. Dipole and qaudrupole are in electron-Bohrs and
   * electron-Bohrs^2, respectively.
   *
   * @return String
   */
  private String toBohrString() {
    StringBuilder multipoleBuffer = new StringBuilder("multipole");
    switch (frameDefinition) {
      case NONE:
        multipoleBuffer.append(format("  %5d  %5s  %5s  %5s", frameAtomTypes[0], "", "", ""));
        break;
      case ZONLY:
        multipoleBuffer.append(
            format("  %5d  %5d  %5s  %5s", frameAtomTypes[0], frameAtomTypes[1], "", ""));
        break;
      case ZTHENX:
        multipoleBuffer.append(
            format(
                "  %5d  %5d  %5d  %5s",
                frameAtomTypes[0], frameAtomTypes[1], frameAtomTypes[2], ""));
        break;
      case BISECTOR:
        multipoleBuffer.append(
            format(
                "  %5d  %5d  %5d  %5s",
                frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2], ""));
        break;
      case ZTHENBISECTOR:
        multipoleBuffer.append(
            format(
                "  %5d  %5d  %5d  %5d",
                frameAtomTypes[0], frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
        break;
      case THREEFOLD:
        multipoleBuffer.append(
            format(
                "  %5d  %5d  %5d  %5d",
                frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
    }
    if (frameAtomTypes.length == 3) {
      multipoleBuffer.append("       ");
    }
    multipoleBuffer.append(
        format(
            "  % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f % 7.5f \\\n"
                + "%11$s % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f % 7.5f",
            multipole[t000],
            multipole[t100] / BOHR,
            multipole[t010] / BOHR,
            multipole[t001] / BOHR,
            multipole[t200] / BOHR2,
            multipole[t110] / BOHR2,
            multipole[t020] / BOHR2,
            multipole[t101] / BOHR2,
            multipole[t011] / BOHR2,
            multipole[t002] / BOHR2,
            "                                      "));
    return multipoleBuffer.toString();
  }

  /** The local multipole frame is defined by the Z-then-X or Bisector convention. */
  public enum MultipoleFrameDefinition {
    NONE,
    ZONLY,
    ZTHENX,
    BISECTOR,
    ZTHENBISECTOR,
    THREEFOLD,
  }
}
