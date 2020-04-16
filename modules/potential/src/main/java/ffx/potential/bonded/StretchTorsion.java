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
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.StretchTorsionType;
import ffx.potential.parameters.TorsionType;
import java.util.logging.Logger;

/**
 * The StretchTorsion class represents a coupling between a torsional angle and the three bonds
 * contained in the torsion, as defined in the 2017 AMOEBA nucleic acid force field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class StretchTorsion extends BondedTerm implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(StretchTorsion.class.getName());
  /** Functional form for OpenMM. */
  private static final String mathForm;

  static {
    /*
     Defined constants:
     p1-p4 are particles 1-4.
     m is a bond number, from 1-3, representing bonds p1-p2, p2-p3, p3-p4.
     n is a periodicity, from 1-3.

     k[m][n] is a set of 9 energy constants defined in the parameter file for this stretch-torsion.

     bVal[m] is the current bond distance for bond m.
     b[m] is the equilibrium distance for bond m.

     tVal is the current value of the 1-2-3-4 dihedral angle.

     phi[m] is a phase offset constant; phi1 = phi3 = 0, phi2 = pi.
    */

    StringBuilder mathFormBuilder = new StringBuilder();

    for (int m = 1; m < 4; m++) {
      for (int n = 1; n < 4; n++) {
        // kmn * (bm - bm(equil)) * (1 + cos(n*tors + phi(n)))
        mathFormBuilder.append(
            String.format("k%d%d*(bVal%d-b%d)*(1+cos(%d*tVal+phi%d))+", m, n, m, m, n, n));
      }
    }
    int lenStr = mathFormBuilder.length();
    mathFormBuilder.replace(lenStr - 1, lenStr, ";");

    for (int m = 1; m < 4; m++) {
      mathFormBuilder.append(String.format("bVal%d=distance(p%d,p%d);", m, m, (m + 1)));
    }

    mathFormBuilder.append("tVal=dihedral(p1,p2,p3,p4)");
    mathForm = mathFormBuilder.toString();
  }

  /**
   * Stretch Torsion force constants (may be reversed compared to storage in the StretchTorsionType
   * instance).
   */
  private final double[] constants = new double[9];
  /** First bond force field type. */
  public BondType bondType1 = null;
  /** Second bond force field type. */
  public BondType bondType2 = null;
  /** Third bond force field type. */
  public BondType bondType3 = null;
  /** Value of lambda. */
  private double lambda = 1.0;
  /** Value of dE/dL. */
  private double dEdL = 0.0;
  /** Flag to indicate lambda dependence. */
  private boolean lambdaTerm = false;
  /** Stretch Torsion force field type. */
  private StretchTorsionType stretchTorsionType = null;
  /** The StretchTorsion may use more sine and cosine terms than are defined in the TorsionType. */
  private double[] tsin = new double[] {0.0, 0.0, 0.0};

  private double[] tcos = new double[] {1.0, 1.0, 1.0};

  /** Torsion force field type. */
  private TorsionType torsionType = null;

  /**
   * Create a StretchTorsion from 3 connected bonds (no error checking)
   *
   * @param b1 Bond
   * @param b2 Bond
   * @param b3 Bond
   */
  private StretchTorsion(Bond b1, Bond b2, Bond b3) {
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
  private StretchTorsion(String n) {
    super(n);
  }

  /**
   * Attempt to create a new StretchTorsion based on the supplied torsion.
   *
   * @param torsion the Torsion.
   * @param forceField the ForceField parameters to apply.
   * @return a new StretchTorsion, or null.
   */
  public static StretchTorsion stretchTorsionFactory(Torsion torsion, ForceField forceField) {
    TorsionType torsionType = torsion.torsionType;
    String key = torsionType.getKey();
    StretchTorsionType stretchTorsionType = forceField.getStretchTorsionType(key);
    if (stretchTorsionType != null) {
      Bond bond1 = torsion.bonds[0];
      Bond middleBond = torsion.bonds[1];
      Bond bond3 = torsion.bonds[2];
      StretchTorsion stretchTorsion = new StretchTorsion(bond1, middleBond, bond3);
      stretchTorsion.stretchTorsionType = stretchTorsionType;
      stretchTorsion.torsionType = torsion.torsionType;
      stretchTorsion.bondType1 = bond1.bondType;
      stretchTorsion.bondType2 = middleBond.bondType;
      stretchTorsion.bondType3 = bond3.bondType;
      Atom atom1 = torsion.atoms[0];
      Atom atom2 = torsion.atoms[1];
      Atom atom3 = torsion.atoms[2];
      Atom atom4 = torsion.atoms[3];
      if (atom1.getAtomType().atomClass == stretchTorsionType.atomClasses[0]
          && atom2.getAtomType().atomClass == stretchTorsionType.atomClasses[1]
          && atom3.getAtomType().atomClass == stretchTorsionType.atomClasses[2]
          && atom4.getAtomType().atomClass == stretchTorsionType.atomClasses[3]) {
        stretchTorsion.setFlipped(false);
      } else {
        stretchTorsion.setFlipped(true);
      }
      return stretchTorsion;
    }
    return null;
  }

  /**
   * Returns the mathematical form of a stretch-torsion as an OpenMM-parsable String.
   *
   * @return Mathematical form of the stretch-torsion coupling.
   */
  public static String stretchTorsionForm() {
    return mathForm;
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
   * <p>Evaluate the Stretch-Torsion energy.
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
    var t = vba.X(vcb);
    var u = vcb.X(vdc);
    var rt2 = max(t.length2(), 0.000001);
    var ru2 = max(u.length2(), 0.000001);
    var rtru = sqrt(rt2 * ru2);
    var vca = vc.sub(va);
    var vdb = vd.sub(vb);
    var cosine = t.dot(u) / rtru;
    var sine = vcb.dot(t.X(u)) / (rcb * rtru);
    value = toDegrees(acos(cosine));
    if (sine < 0.0) {
      value = -value;
    }

    // Compute multiple angle trigonometry and phase terms
    var phi1 = 1.0 + (cosine * tcos[0] + sine * tsin[0]);
    var dphi1 = (cosine * tsin[0] - sine * tcos[0]);
    var cosine2 = cosine * cosine - sine * sine;
    var sine2 = 2.0 * cosine * sine;
    var phi2 = 1.0 + (cosine2 * tcos[1] + sine2 * tsin[1]);
    var dphi2 = 2.0 * (cosine2 * tsin[1] - sine2 * tcos[1]);
    var sine3 = cosine * sine2 + sine * cosine2;
    var cosine3 = cosine * cosine2 - sine * sine2;
    var phi3 = 1.0 + (cosine3 * tcos[2] + sine3 * tsin[2]);
    var dphi3 = 3.0 * (cosine3 * tsin[2] - sine3 * tcos[2]);

    // Get the stretch-torsion values for the first bond.
    var c1 = constants[0];
    var c2 = constants[1];
    var c3 = constants[2];
    var rba = sqrt(rba2);
    var dr1 = rba - bondType1.distance;
    var units = StretchTorsionType.units;
    var s1 = c1 * phi1 + c2 * phi2 + c3 * phi3;
    var e1 = units * dr1 * s1;

    // Get the stretch-torsion values for the second bond.
    var c4 = constants[3];
    var c5 = constants[4];
    var c6 = constants[5];
    var dr2 = rcb - bondType2.distance;
    var s2 = c4 * phi1 + c5 * phi2 + c6 * phi3;
    var e2 = units * dr2 * s2;

    // Get the stretch-torsion values for the third bond.
    var c7 = constants[6];
    var c8 = constants[7];
    var c9 = constants[8];
    var rdc = sqrt(rdc2);
    var dr3 = rdc - bondType3.distance;
    var s3 = c7 * phi1 + c8 * phi2 + c9 * phi3;
    var e3 = units * dr3 * s3;
    energy = e1 + e2 + e3;
    if (esvTerm) {
      esvDerivLocal = energy * dedesvChain * lambda;
    }
    energy = energy * esvLambda * lambda;
    dEdL = energy * esvLambda;
    if (gradient || lambdaTerm) {
      // Compute derivative components for the first bond.
      var dedphi = units * dr1 * (c1 * dphi1 + c2 * dphi2 + c3 * dphi3);
      var ddrd = vba.scale(units * s1 / rba);
      var dedt = t.X(vcb).scaleI(dedphi / (rt2 * rcb));
      var dedu = u.X(vcb).scaleI(-dedphi / (ru2 * rcb));
      // Determine chain rule components for the first bond.
      var ga = dedt.X(vcb).subI(ddrd);
      var gb = vca.X(dedt).addI(dedu.X(vdc)).addI(ddrd);
      var gc = dedt.X(vba).addI(vdb.X(dedu));
      var gd = dedu.X(vcb);

      // Compute derivative components for the 2nd bond.
      dedphi = units * dr2 * (c4 * dphi1 + c5 * dphi2 + c6 * dphi3);
      ddrd = vcb.scale(units * s2 / rcb);
      dedt = t.X(vcb).scaleI(dedphi / (rt2 * rcb));
      dedu = u.X(vcb).scaleI(-dedphi / (ru2 * rcb));
      // Accumulate derivative components.
      ga.addI(dedt.X(vcb));
      gb.addI(vca.X(dedt).addI(dedu.X(vdc)).subI(ddrd));
      gc.addI(dedt.X(vba).addI(vdb.X(dedu)).addI(ddrd));
      gd.addI(dedu.X(vcb));

      // Compute derivative components for the 3rd bond.
      dedphi = units * dr3 * (c7 * dphi1 + c8 * dphi2 + c9 * dphi3);
      ddrd = vdc.scale(units * s3 / rdc);
      dedt = t.X(vcb).scaleI(dedphi / (rt2 * rcb));
      dedu = u.X(vcb).scaleI(-dedphi / (ru2 * rcb));
      // Accumulate derivative components.
      ga.addI(dedt.X(vcb));
      gb.addI(vca.X(dedt).addI(dedu.X(vdc)));
      gc.addI(dedt.X(vba).addI(vdb.X(dedu)).subI(ddrd));
      gd.addI(dedu.X(vcb).addI(ddrd));

      // Apply ESV lambda
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
  public Atom get1_4(Atom a) {
    if (a == atoms[0]) {
      return atoms[3];
    }
    if (a == atoms[3]) {
      return atoms[0];
    }
    return null;
  }

  /**
   * Returns the array of stretch-torsion constants, in units of kcal/mol/A.
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
    // No contribution.
  }

  /** Log details for this Stretch-Torsional Angle energy term. */
  public void log() {
    logger.info(
        String.format(
            " %-8s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f %10.4f",
            "Stretch-Torsion",
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
    return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
  }

  /** Initialization */
  private void initialize() {
    atoms = new ffx.potential.bonded.Atom[4];
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
      constants[0] = stretchTorsionType.forceConstants[6];
      constants[1] = stretchTorsionType.forceConstants[7];
      constants[2] = stretchTorsionType.forceConstants[8];
      constants[3] = stretchTorsionType.forceConstants[3];
      constants[4] = stretchTorsionType.forceConstants[4];
      constants[5] = stretchTorsionType.forceConstants[5];
      constants[6] = stretchTorsionType.forceConstants[0];
      constants[7] = stretchTorsionType.forceConstants[1];
      constants[8] = stretchTorsionType.forceConstants[2];
    } else {
      arraycopy(stretchTorsionType.forceConstants, 0, constants, 0, 9);
    }

    tsin = new double[] {0.0, 0.0, 0.0};
    tcos = new double[] {1.0, 1.0, 1.0};
    arraycopy(torsionType.sine, 0, tsin, 0, min(torsionType.sine.length, 3));
    arraycopy(torsionType.cosine, 0, tcos, 0, min(torsionType.cosine.length, 3));
  }
}
