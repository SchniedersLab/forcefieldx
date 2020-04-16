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
package ffx.potential.bonded;

import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionType;
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
  /** Unit conversion. */
  public double units = 0.5;
  /** Value of lambda. */
  private double lambda = 1.0;
  /** Value of dE/dL. */
  private double dEdL = 0.0;
  /** Flag to indicate lambda dependence. */
  private boolean lambdaTerm = false;

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
    initialize();
  }

  /**
   * Torsion Constructor.
   *
   * @param n Torsion id
   */
  public Torsion(String n) {
    super(n);
  }

  /**
   * Log that no TorsionType exists.
   *
   * @param a0 Atom 0.
   * @param a1 Atom 1.
   * @param a2 Atom 2.
   * @param a3 Atom 3.
   * @param key The class key.
   */
  public static void logNoTorsionType(Atom a0, Atom a1, Atom a2, Atom a3, String key) {
    logger.severe(
        format(
            "No TorsionType for key: %s\n %s -> %s\n %s -> %s\n %s -> %s\n %s -> %s",
            key,
            a0.toString(),
            a0.getAtomType().toString(),
            a1.toString(),
            a1.getAtomType().toString(),
            a2.toString(),
            a2.getAtomType().toString(),
            a3.toString(),
            a3.getAtomType().toString()));
  }

  /**
   * Find a torsion based on the specified classes.
   *
   * @param c0 Atom class 0.
   * @param c1 Atom class 1.
   * @param c2 Atom class 2.
   * @param c3 Atom class 3.
   * @param forceField Force Field parameters to use.
   * @return A torsion type if it exists.
   */
  private static TorsionType getTorsionType(int c0, int c1, int c2, int c3, ForceField forceField) {
    int[] c = {c0, c1, c2, c3};
    String key = TorsionType.sortKey(c);
    return forceField.getTorsionType(key);
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

    int c0 = a0.getAtomType().atomClass;
    int c1 = a1.getAtomType().atomClass;
    int c2 = a2.getAtomType().atomClass;
    int c3 = a3.getAtomType().atomClass;

    TorsionType torsionType = getTorsionType(c0, c1, c2, c3, forceField);

    // Single wild card.
    if (torsionType == null) {
      if (c0 > c3) {
        torsionType = getTorsionType(c0, c1, c2, 0, forceField);
        if (torsionType == null) {
          torsionType = getTorsionType(0, c1, c2, c3, forceField);
        }
      } else {
        torsionType = getTorsionType(0, c1, c2, c3, forceField);
        if (torsionType == null) {
          torsionType = getTorsionType(c0, c1, c2, 0, forceField);
        }
      }
    }

    // Double wild card.
    if (torsionType == null) {
      torsionType = getTorsionType(0, c1, c2, 0, forceField);
    }

    // No torsion type found.
    if (torsionType == null) {
      int[] c = {c0, c1, c2, c3};
      String key = TorsionType.sortKey(c);
      logNoTorsionType(a0, a1, a2, a3, key);
      return null;
    }

    Torsion torsion = new Torsion(bond1, middleBond, bond3);
    torsion.torsionType = torsionType;
    torsion.units = forceField.getDouble("TORSIONUNIT", 1.0);

    return torsion;
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
      if (esvTerm) {
        esvDerivLocal = units * energy * dedesvChain * lambda;
      }
      energy = units * energy * esvLambda * lambda;
      dEdL = units * energy * esvLambda;
      if (gradient || lambdaTerm) {
        dedphi = units * dedphi * esvLambda;
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
   * If the specified atom is not a central atom of <b>this</b> torsion, the atom at the opposite
   * end is returned. These atoms are said to be 1-4 to each other.
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
      this.lambda = lambda;
      lambdaTerm = true;
    } else {
      this.lambda = 1.0;
    }
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
