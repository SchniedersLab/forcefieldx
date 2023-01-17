// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.potential.bonded;

import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.math.DoubleMath;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionType;
import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.logging.Logger;

/**
 * The Torsion class represents a torsional angle formed between four bonded atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Torsion extends BondedTerm implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(Torsion.class.getName());
  /** The force field Torsion type in use. */
  public TorsionType torsionType = null;
  /** Value of lambda. */
  private double lambda = 1.0;
  /** Value of dE/dL. */
  private double dEdL = 0.0;
  /** Flag to indicate lambda dependence. */
  private boolean lambdaTerm = false;
  /** Maps global lambda to either itself or 1 - global lambda. */
  private final DoubleUnaryOperator lambdaMapper;

  /**
   * Torsion constructor.
   *
   * @param an1 Angle that combines to form the Torsional Angle
   * @param an2 Angle that combines to form the Torsional Angle
   */
  public Torsion(Angle an1, Angle an2) {
    super();
    bonds = new Bond[3];
    bonds[1] = an1.getCommonBond(an2);
    bonds[0] = an1.getOtherBond(bonds[1]);
    bonds[2] = an2.getOtherBond(bonds[1]);
    lambdaMapper = (double d) -> d;
    initialize();
  }

  /**
   * Torsion constructor.
   *
   * @param a Angle that has one Atom in common with Bond b
   * @param b Bond that has one Atom in common with Angle A
   */
  public Torsion(Angle a, Bond b) {
    super();
    bonds = new Bond[3];
    bonds[0] = b;
    bonds[1] = a.getBond(0);
    bonds[2] = a.getBond(1);
    // See if bond 2 or bond 3 is the middle bond
    Atom atom = bonds[1].getCommonAtom(b);
    if (atom == null) {
      Bond temp = bonds[1];
      bonds[1] = bonds[2];
      bonds[2] = temp;
    }
    lambdaMapper = (double d) -> d;
    initialize();
  }

  /**
   * Create a Torsion from 3 connected bonds (no error checking)
   *
   * @param b1 Bond
   * @param b2 Bond
   * @param b3 Bond
   */
  public Torsion(Bond b1, Bond b2, Bond b3) {
    super();
    bonds = new Bond[3];
    bonds[0] = b1;
    bonds[1] = b2;
    bonds[2] = b3;
    lambdaMapper = (double d) -> d;
    initialize();
  }

  /**
   * Torsion Constructor.
   *
   * @param n Torsion id
   */
  public Torsion(String n) {
    super(n);
    lambdaMapper = (double d) -> d;
  }

  /**
   * Log that no TorsionType exists.
   *
   * @param a0 Atom 0.
   * @param a1 Atom 1.
   * @param a2 Atom 2.
   * @param a3 Atom 3.
   */
  public static void logNoTorsionType(Atom a0, Atom a1, Atom a2, Atom a3, ForceField forceField) {
    AtomType atomType0 = a0.getAtomType();
    AtomType atomType1 = a1.getAtomType();
    AtomType atomType2 = a2.getAtomType();
    AtomType atomType3 = a3.getAtomType();
    int[] c = {atomType0.atomClass, atomType1.atomClass, atomType2.atomClass, atomType3.atomClass};
    String key = TorsionType.sortKey(c);
    StringBuilder sb = new StringBuilder(
        format(
            "No TorsionType for key: %s\n %s -> %s\n %s -> %s\n %s -> %s\n %s -> %s",
            key, a0, atomType0, a1, atomType1, a2, atomType2, a3, atomType3));
    if (matchTorsions(a0, a1, a2, a3, forceField, sb, true) <= 0) {
      matchTorsions(a0, a1, a2, a3, forceField, sb, false);
    }
    logger.severe(sb.toString());
  }

  private static int matchTorsions(Atom a0, Atom a1, Atom a2, Atom a3,
      ForceField forceField, StringBuilder sb, boolean strict) {
    AtomType atomType0 = a0.getAtomType();
    AtomType atomType1 = a1.getAtomType();
    AtomType atomType2 = a2.getAtomType();
    AtomType atomType3 = a3.getAtomType();
    int c0 = atomType0.atomClass;
    int c1 = atomType1.atomClass;
    int c2 = atomType2.atomClass;
    int c3 = atomType3.atomClass;
    List<AtomType> types0 = forceField.getSimilarAtomTypes(atomType0);
    List<AtomType> types1 = forceField.getSimilarAtomTypes(atomType1);
    List<AtomType> types2 = forceField.getSimilarAtomTypes(atomType2);
    List<AtomType> types3 = forceField.getSimilarAtomTypes(atomType3);
    List<TorsionType> torsionTypes = new ArrayList<>();
    boolean match = false;
    for (AtomType type1 : types1) {
      for (AtomType type2 : types2) {
        // Match both inner atom classes.
        if ((type1.atomClass == c1 && type2.atomClass == c2) ||
            (type1.atomClass == c2 && type2.atomClass == c1)) {
          for (AtomType type0 : types0) {
            for (AtomType type3 : types3) {
              // Match one distal atom class.
              if (strict) {
                if ((type0.atomClass != c0) && (type0.atomClass != c3) &&
                    (type3.atomClass != c0) && (type3.atomClass != c3)) {
                  continue;
                }
              }
              TorsionType torsionType = forceField.getTorsionType(type0, type1, type2, type3);
              if (torsionType != null && !torsionTypes.contains(torsionType)) {
                if (!match) {
                  match = true;
                  sb.append("\n Similar Angle Types:");
                }
                torsionTypes.add(torsionType);
                sb.append(format("\n  %s", torsionType));
              }
            }
          }
        }
      }
    }
    return torsionTypes.size();
  }

  /**
   * Attempt to create a new Torsion based on the supplied bonds. There is no error checking to
   * enforce that the bonds make up a linear series of 4 bonded atoms.
   *
   * @param bond1 the first Bond.
   * @param middleBond the middle Bond.
   * @param bond3 the last Bond.
   * @param forceField the ForceField parameters to apply.
   * @return a new Torsion, or null.
   */
  static Torsion torsionFactory(Bond bond1, Bond middleBond, Bond bond3, ForceField forceField) {
    Atom a0 = bond1.getOtherAtom(middleBond);
    Atom a1 = middleBond.getAtom(0);
    Atom a2 = middleBond.getAtom(1);
    Atom a3 = bond3.getOtherAtom(middleBond);

    TorsionType torsionType = forceField.getTorsionType(a0.getAtomType(), a1.getAtomType(),
        a2.getAtomType(), a3.getAtomType());

    // No torsion type found.
    if (torsionType == null) {
      logNoTorsionType(a0, a1, a2, a3, forceField);
      return null;
    }

    Torsion torsion = new Torsion(bond1, middleBond, bond3);
    torsion.setTorsionType(torsionType);
    return torsion;
  }

  /**
   * Set the torsion type.
   *
   * @param torsionType The TorsionType.
   */
  public void setTorsionType(TorsionType torsionType) {
    this.torsionType = torsionType;
  }

  /**
   * compare
   *
   * @param a0 a {@link ffx.potential.bonded.Atom} object.
   * @param a1 a {@link ffx.potential.bonded.Atom} object.
   * @param a2 a {@link ffx.potential.bonded.Atom} object.
   * @param a3 a {@link ffx.potential.bonded.Atom} object.
   * @return a boolean.
   */
  public boolean compare(Atom a0, Atom a1, Atom a2, Atom a3) {
    if (a0 == atoms[0] && a1 == atoms[1] && a2 == atoms[2] && a3 == atoms[3]) {
      return true;
    }
    return (a0 == atoms[3] && a1 == atoms[2] && a2 == atoms[1] && a3 == atoms[0]);
  }

  /**
   * Compute the torsional angle in degrees.
   *
   * @return The torsion in degrees.
   */
  public double measure() {
    double angle = DoubleMath.dihedralAngle(
        atoms[0].getXYZ(null),
        atoms[1].getXYZ(null),
        atoms[2].getXYZ(null),
        atoms[3].getXYZ(null));
    return toDegrees(angle);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Evaluate the Torsional Angle energy.
   */
  @Override
  public double energy(
      boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    energy = 0.0;
    value = 0.0;
    dEdL = 0.0;
    var atomA = atoms[0];
    var atomB = atoms[1];
    var atomC = atoms[2];
    var atomD = atoms[3];
    var va = atomA.getXYZ();
    var vb = atomB.getXYZ();
    var vc = atomC.getXYZ();
    var vd = atomD.getXYZ();
    var vba = vb.sub(va);
    var vcb = vc.sub(vb);
    var vdc = vd.sub(vc);
    var vt = vba.X(vcb);
    var vu = vcb.X(vdc);
    var rt2 = vt.length2();
    var ru2 = vu.length2();
    var rtru2 = rt2 * ru2;
    if (rtru2 != 0.0) {
      var rr = sqrt(rtru2);
      var rcb = vcb.length();
      var cosine = vt.dot(vu) / rr;
      var sine = vcb.dot(vt.X(vu)) / (rcb * rr);
      value = toDegrees(acos(cosine));
      if (sine < 0.0) {
        value = -value;
      }
      var amp = torsionType.amplitude;
      var tsin = torsionType.sine;
      var tcos = torsionType.cosine;
      energy = amp[0] * (1.0 + cosine * tcos[0] + sine * tsin[0]);
      var dedphi = amp[0] * (cosine * tsin[0] - sine * tcos[0]);
      var cosprev = cosine;
      var sinprev = sine;
      var n = torsionType.terms;
      for (int i = 1; i < n; i++) {
        var cosn = cosine * cosprev - sine * sinprev;
        var sinn = sine * cosprev + cosine * sinprev;
        var phi = 1.0 + cosn * tcos[i] + sinn * tsin[i];
        var dphi = (1.0 + i) * (cosn * tsin[i] - sinn * tcos[i]);
        energy = energy + amp[i] * phi;
        dedphi = dedphi + amp[i] * dphi;
        cosprev = cosn;
        sinprev = sinn;
      }
      energy = torsionType.torsionUnit * energy * lambda;
      dEdL = torsionType.torsionUnit * energy;
      if (gradient || lambdaTerm) {
        dedphi = torsionType.torsionUnit * dedphi;
        var vca = vc.sub(va);
        var vdb = vd.sub(vb);
        var dedt = vt.X(vcb).scaleI(dedphi / (rt2 * rcb));
        var dedu = vu.X(vcb).scaleI(-dedphi / (ru2 * rcb));
        var ga = dedt.X(vcb);
        var gb = vca.X(dedt).addI(dedu.X(vdc));
        var gc = dedt.X(vba).addI(vdb.X(dedu));
        var gd = dedu.X(vcb);
        int ia = atomA.getIndex() - 1;
        int ib = atomB.getIndex() - 1;
        int ic = atomC.getIndex() - 1;
        int id = atomD.getIndex() - 1;
        if (lambdaTerm) {
          lambdaGrad.add(threadID, ia, ga);
          lambdaGrad.add(threadID, ib, gb);
          lambdaGrad.add(threadID, ic, gc);
          lambdaGrad.add(threadID, id, gd);
        }
        if (gradient) {
          grad.add(threadID, ia, ga.scaleI(lambda));
          grad.add(threadID, ib, gb.scaleI(lambda));
          grad.add(threadID, ic, gc.scaleI(lambda));
          grad.add(threadID, id, gd.scaleI(lambda));
        }
      }
    }

    return energy;
  }

  /**
   * If the specified atom is not a central atom of <b>this</b> torsion, the atom at the opposite end
   * is returned. These atoms are said to be 1-4 to each other.
   *
   * @param a Atom
   * @return Atom
   */
  public Atom get1_4(Atom a) {
    if (a == atoms[0]) {
      return atoms[3];
    }
    if (a == atoms[3]) {
      return atoms[0];
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public double getLambda() {
    return lambda;
  }

  /** {@inheritDoc} */
  @Override
  public void setLambda(double lambda) {
    if (applyAllLambda()) {
      lambdaTerm = true;
    }
    this.lambda = lambdaTerm ? lambdaMapper.applyAsDouble(lambda) : 1.0;
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    return 0.0;
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {
    if (lambdaTerm) {
      return dEdL;
    } else {
      return 0.0;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] gradient) {
    // The chain rule term is zero.
  }

  /** Log details for this Torsional Angle energy term. */
  public void log() {
    logger.info(
        format(
            " %-8s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f %10.4f",
            "Torsional-Angle",
            atoms[0].getIndex(),
            atoms[0].getAtomType().name,
            atoms[1].getIndex(),
            atoms[1].getAtomType().name,
            atoms[2].getIndex(),
            atoms[2].getAtomType().name,
            atoms[3].getIndex(),
            atoms[3].getAtomType().name,
            value,
            energy));
  }

  /**
   * {@inheritDoc}
   *
   * <p>Overidden toString Method returns the Term's id.
   */
  @Override
  public String toString() {
    return format("%s  (%7.1f,%7.2f)", id, value, energy);
  }

  /** Initialization */
  private void initialize() {
    atoms = new Atom[4];
    atoms[1] = bonds[0].getCommonAtom(bonds[1]);
    atoms[0] = bonds[0].get1_2(atoms[1]);
    atoms[2] = bonds[1].get1_2(atoms[1]);
    atoms[3] = bonds[2].get1_2(atoms[2]);
    atoms[0].setTorsion(this);
    atoms[1].setTorsion(this);
    atoms[2].setTorsion(this);
    atoms[3].setTorsion(this);
    setID_Key(false);
    value =
        dihedralAngle(
            atoms[0].getXYZ(null),
            atoms[1].getXYZ(null),
            atoms[2].getXYZ(null),
            atoms[3].getXYZ(null));
  }
}
