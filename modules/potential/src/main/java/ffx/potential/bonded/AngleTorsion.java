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
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.AngleTorsionType;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionType;
import java.util.logging.Logger;

/**
 * The AngleTorsion class represents an angle torsion coupling between four bonded atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AngleTorsion extends BondedTerm implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(AngleTorsion.class.getName());
  /** Functional form for OpenMM. */
  private static final String mathForm;

  static {
    /*
     Defined constants:
     p1-p4 are particles 1-4.
     m is an angle number, from 1-2, representing angles p1-p2-p3, p2-p3-p4.
     n is a periodicity, from 1-3.

     k[m][n] is a set of 6 energy constants defined in the parameter file for this angle-torsion.

     aVal[m] is the current value for angle m.
     a[m] is the equilibrium value for angle m.

     tVal is the current value of the 1-2-3-4 dihedral angle.

     phi[m] is a phase offset constant; phi1 = phi3 = 0, phi2 = pi.
    */

    StringBuilder mathFormBuilder = new StringBuilder();
    for (int m = 1; m < 3; m++) {
      for (int n = 1; n < 4; n++) {
        // kmn * (am - am(equil)) * (1 + cos(n*tors + phi(n)))
        mathFormBuilder.append(
            format("k%d%d*(aVal%d-a%d)*(1+cos(%d*tVal+phi%d))+", m, n, m, m, n, n));
      }
    }
    int lenStr = mathFormBuilder.length();
    mathFormBuilder.replace(lenStr - 1, lenStr, ";");
    for (int m = 1; m < 3; m++) {
      mathFormBuilder.append(format("aVal%d=angle(p%d,p%d,p%d);", m, m, (m + 1), (m + 2)));
    }
    mathFormBuilder.append("tVal=dihedral(p1,p2,p3,p4)");
    mathForm = mathFormBuilder.toString();
  }

  /**
   * Angle Torsion force constants (may be reversed compared to storage in the AngleTorsionType
   * instance).
   */
  private final double[] constants = new double[6];
  /** The AngleTorsion may use more sine and cosine terms than are defined in the TorsionType. */
  private final double[] tsin = new double[] {0.0, 0.0, 0.0};

  private final double[] tcos = new double[] {1.0, 1.0, 1.0};
  /** First angle force field type. */
  public AngleType angleType1 = null;
  /** Second angle force field type. */
  public AngleType angleType2 = null;
  /** Angle Torsion force field type. */
  private AngleTorsionType angleTorsionType = null;
  /** Torsion force field type. */
  private TorsionType torsionType = null;
  /** Value of lambda. */
  private double lambda = 1.0;
  /** Value of dE/dL. */
  private double dEdL = 0.0;
  /** Flag to indicate lambda dependence. */
  private boolean lambdaTerm = false;

  /**
   * AngleTorsion constructor.
   *
   * @param an1 Angle that combines to form the Torsional Angle
   * @param an2 Angle that combines to form the Torsional Angle
   */
  public AngleTorsion(Angle an1, Angle an2) {
    super();
    bonds = new Bond[3];
    bonds[1] = an1.getCommonBond(an2);
    bonds[0] = an1.getOtherBond(bonds[1]);
    bonds[2] = an2.getOtherBond(bonds[1]);
    initialize();
  }

  /**
   * Create a AngleTorsion from 3 connected bonds (no error checking)
   *
   * @param b1 Bond
   * @param b2 Bond
   * @param b3 Bond
   */
  public AngleTorsion(Bond b1, Bond b2, Bond b3) {
    super();
    bonds = new Bond[3];
    bonds[0] = b1;
    bonds[1] = b2;
    bonds[2] = b3;
    initialize();
  }

  /**
   * AngleTorsion Constructor.
   *
   * @param n Torsion id
   */
  public AngleTorsion(String n) {
    super(n);
  }

  /**
   * Returns the mathematical form of an angle-torsion as an OpenMM-parsable String.
   *
   * @return Mathematical form of the angle-torsion coupling.
   */
  public static String angleTorsionForm() {
    return mathForm;
  }

  /**
   * Attempt to create a new AngleTorsion based on the supplied torsion.
   *
   * @param torsion the Torsion.
   * @param forceField the ForceField parameters to apply.
   * @return a new Torsion, or null.
   */
  static AngleTorsion angleTorsionFactory(Torsion torsion, ForceField forceField) {
    TorsionType torsionType = torsion.torsionType;
    String key = torsionType.getKey();
    AngleTorsionType angleTorsionType = forceField.getAngleTorsionType(key);

    if (angleTorsionType != null) {
      Bond bond1 = torsion.bonds[0];
      Bond middleBond = torsion.bonds[1];
      Bond bond3 = torsion.bonds[2];

      AngleTorsion angleTorsion = new AngleTorsion(bond1, middleBond, bond3);
      angleTorsion.angleTorsionType = angleTorsionType;
      angleTorsion.torsionType = torsion.torsionType;
      Atom atom1 = torsion.atoms[0];
      Atom atom2 = torsion.atoms[1];
      Atom atom3 = torsion.atoms[2];
      Atom atom4 = torsion.atoms[3];

      Angle angle1 = atom1.getAngle(atom2, atom3);
      Angle angle2 = atom2.getAngle(atom3, atom4);
      angleTorsion.angleType1 = angle1.angleType;
      angleTorsion.angleType2 = angle2.angleType;

      if (atom1.getAtomType().atomClass == angleTorsionType.atomClasses[0]
          && atom2.getAtomType().atomClass == angleTorsionType.atomClasses[1]
          && atom3.getAtomType().atomClass == angleTorsionType.atomClasses[2]
          && atom4.getAtomType().atomClass == angleTorsionType.atomClasses[3]) {
        angleTorsion.setFlipped(false);
      } else {
        angleTorsion.setFlipped(true);
      }
      return angleTorsion;
    }
    return null;
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
   * <p>Evaluate the Angle-Torsion energy.
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
    var rba2 = vba.length2();
    var rcb2 = vcb.length2();
    var rdc2 = vdc.length2();
    if (min(min(rba2, rcb2), rdc2) == 0.0) {
      return 0.0;
    }
    var rcb = sqrt(rcb2);
    var vt = vba.X(vcb);
    var vu = vcb.X(vdc);
    var rt2 = max(vt.length2(), 0.000001);
    var ru2 = max(vu.length2(), 0.000001);
    var rtru = sqrt(rt2 * ru2);
    var vca = vc.sub(va);
    var vdb = vd.sub(vb);
    var cosine = vt.dot(vu) / rtru;
    var sine = vcb.dot(vt.X(vu)) / (rcb * rtru);
    value = toDegrees(acos(cosine));
    if (sine < 0.0) {
      value = -value;
    }

    // Compute multiple angle trigonometry and phase terms
    var cosine2 = cosine * cosine - sine * sine;
    var sine2 = 2.0 * cosine * sine;
    var cosine3 = cosine * cosine2 - sine * sine2;
    var sine3 = cosine * sine2 + sine * cosine2;
    var phi1 = 1.0 + (cosine * tcos[0] + sine * tsin[0]);
    var dphi1 = (cosine * tsin[0] - sine * tcos[0]);
    var phi2 = 1.0 + (cosine2 * tcos[1] + sine2 * tsin[1]);
    var dphi2 = 2.0 * (cosine2 * tsin[1] - sine2 * tcos[1]);
    var phi3 = 1.0 + (cosine3 * tcos[2] + sine3 * tsin[2]);
    var dphi3 = 3.0 * (cosine3 * tsin[2] - sine3 * tcos[2]);

    // Set the angle-torsion parameters for the first angle
    var units = AngleTorsionType.units;
    var c1 = constants[0];
    var c2 = constants[1];
    var c3 = constants[2];
    var angle1 = toDegrees(acos(-vba.dot(vcb) / sqrt(rba2 * rcb2)));
    var dt1 = angle1 - angleType1.angle[0];
    var s1 = c1 * phi1 + c2 * phi2 + c3 * phi3;
    var e1 = units * dt1 * s1;

    // Set the angle-torsion values for the second angle
    var c4 = constants[3];
    var c5 = constants[4];
    var c6 = constants[5];
    var angle2 = toDegrees(acos(-vcb.dot(vdc) / sqrt(rcb2 * rdc2)));
    var dt2 = angle2 - angleType2.angle[0];
    var s2 = c4 * phi1 + c5 * phi2 + c6 * phi3;
    var e2 = units * dt2 * s2;
    energy = e1 + e2;
    if (esvTerm) {
      esvDerivLocal = energy * dedesvChain * lambda;
    }
    energy = energy * esvLambda * lambda;
    dEdL = energy * esvLambda;
    if (gradient || lambdaTerm) {
      // Compute derivative components for this interaction.
      var dedphi = units * dt1 * (c1 * dphi1 + c2 * dphi2 + c3 * dphi3);
      var ddt = units * toDegrees(s1);
      var vdt = vt.X(vcb).scaleI(dedphi / (rt2 * rcb));
      var vdu = vcb.X(vu).scaleI(dedphi / (ru2 * rcb));

      // Determine chain rule components for the first angle.
      var rt = sqrt(rt2);
      var sa = -ddt / (rba2 * rt);
      var sc = ddt / (rcb2 * rt);
      var ga = vt.X(vba).scaleI(sa).addI(vdt.X(vcb));
      var gb = vba.X(vt).scaleI(sa).addI(vt.X(vcb).scaleI(sc)).addI(vca.X(vdt)).addI(vdu.X(vdc));
      var gc = vcb.X(vt).scaleI(sc).addI(vdt.X(vba)).addI(vdb.X(vdu));
      var gd = vdu.X(vcb);

      // Compute derivative components for the 2nd angle.
      dedphi = units * dt2 * (c4 * dphi1 + c5 * dphi2 + c6 * dphi3);
      ddt = units * toDegrees(s2);
      vdt = vt.X(vcb).scaleI(dedphi / (rt2 * rcb));
      vdu = vcb.X(vu).scaleI(dedphi / (ru2 * rcb));

      // Increment chain rule components for the 2nd angle.
      var ur = sqrt(ru2);
      var sb = -ddt / (rcb2 * ur);
      var sd = ddt / (rdc2 * ur);
      ga.addI(vdt.X(vcb));
      gb.addI(vu.X(vcb).scaleI(sb)).addI(vca.X(vdt)).addI(vdu.X(vdc));
      gc.addI(vcb.X(vu).scaleI(sb)).addI(vu.X(vdc).scaleI(sd)).addI(vdt.X(vba)).addI(vdb.X(vdu));
      gd.addI(vdc.X(vu).scaleI(sd)).addI(vdu.X(vcb));
      ga.scaleI(esvLambda);
      gb.scaleI(esvLambda);
      gc.scaleI(esvLambda);
      gd.scaleI(esvLambda);
      var ia = atomA.getIndex() - 1;
      var ib = atomB.getIndex() - 1;
      var ic = atomC.getIndex() - 1;
      var id = atomD.getIndex() - 1;
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

    return energy;
  }

  /**
   * If the specified atom is not a central atom of <b>this</b> torsion, the atom at the opposite
   * end is returned. These atoms are said to be 1-4 to each other.
   *
   * @param a Atom
   * @return Atom
   */
  public Atom get1_4(ffx.potential.bonded.Atom a) {
    if (a == atoms[0]) {
      return atoms[3];
    }
    if (a == atoms[3]) {
      return atoms[0];
    }
    return null;
  }

  /**
   * Returns the array of stretch-torsion constants, in units of kcal/mol/degree.
   *
   * @return Stretch-torsion constants.
   */
  public double[] getConstants() {
    return copyOf(constants, constants.length);
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
    // The dEdXdL terms are zero.
  }

  /** Log details for this Angle-Torsion energy term. */
  public void log() {
    logger.info(
        format(
            " %-8s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f %10.4f",
            "Angle-Torsion",
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
    setID_Key(false);
    value =
        dihedralAngle(
            atoms[0].getXYZ(null),
            atoms[1].getXYZ(null),
            atoms[2].getXYZ(null),
            atoms[3].getXYZ(null));
  }

  /**
   * setFlipped.
   *
   * @param flipped a boolean.
   */
  private void setFlipped(boolean flipped) {
    if (flipped) {
      constants[0] = angleTorsionType.forceConstants[3];
      constants[1] = angleTorsionType.forceConstants[4];
      constants[2] = angleTorsionType.forceConstants[5];
      constants[3] = angleTorsionType.forceConstants[0];
      constants[4] = angleTorsionType.forceConstants[1];
      constants[5] = angleTorsionType.forceConstants[2];
    } else {
      arraycopy(angleTorsionType.forceConstants, 0, constants, 0, 6);
    }

    arraycopy(torsionType.sine, 0, tsin, 0, min(torsionType.sine.length, 3));
    arraycopy(torsionType.cosine, 0, tcos, 0, min(torsionType.cosine.length, 3));
  }
}
