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

import static ffx.potential.parameters.StretchBendType.units;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.StretchBendType;
import java.util.logging.Logger;

/**
 * The StretchBend class represents a Stretch-Bend formed between three linearly bonded atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class StretchBend extends BondedTerm implements Comparable<BondedTerm> {

  private static final Logger logger = Logger.getLogger(StretchBend.class.getName());
  /** Equilibrium angle. */
  public final double angleEq;
  /** First equilibrium bond distance. */
  public final double bond0Eq;
  /** Second equilibrium bond distance. */
  public final double bond1Eq;
  /** Angle this Stretch-Bend is based on. */
  protected final Angle angle;
  /** Force constant. */
  public double force0, force1;
  /** Force field parameters to compute the Stretch-Bend energy. */
  private StretchBendType stretchBendType = null;
  /** Rigid scale factor to apply to the force constant. */
  private double rigidScale = 1.0;

  /**
   * Constructor for the Stretch-Bend class.
   *
   * @param a The Angle that this stretch-bend is based on.
   */
  public StretchBend(Angle a) {
    super();
    angle = a;
    atoms = a.atoms;
    bonds = a.bonds;
    angleEq = angle.angleType.angle[angle.nh];
    bond0Eq = bonds[0].bondType.distance;
    bond1Eq = bonds[1].bondType.distance;
    setID_Key(false);
  }

  /**
   * Attempt to create a new StretchBend if a StretchBendType exists for the specified Angle.
   *
   * @param angle the Angle to created the StrechBend around.
   * @param forceField the ForceField parameters to use.
   * @return a new StretchBend, or null.
   */
  static StretchBend stretchBendFactory(Angle angle, ForceField forceField) {
    StretchBendType stretchBendType = forceField.getStretchBendType(angle.getAngleType().getKey());
    if (stretchBendType == null) {
      return null;
    }
    StretchBend stretchBend = new StretchBend(angle);
    stretchBend.setStretchBendType(stretchBendType);
    return stretchBend;
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(BondedTerm sb) {
    if (!sb.getClass().isInstance(this)) {
      return super.compareTo(sb);
    }
    return angle.compareTo(((StretchBend) sb).angle);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Evaluate the Stretch-Bend energy.
   */
  @Override
  public double energy(
      boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    energy = 0.0;
    value = 0.0;
    var atomA = atoms[0];
    var atomB = atoms[1];
    var atomC = atoms[2];
    var ia = atomA.getIndex() - 1;
    var ib = atomB.getIndex() - 1;
    var ic = atomC.getIndex() - 1;
    var va = atomA.getXYZ();
    var vb = atomB.getXYZ();
    var vc = atomC.getXYZ();
    var vab = va.sub(vb);
    var vcb = vc.sub(vb);
    var rab2 = vab.length2();
    var rcb2 = vcb.length2();
    if (rab2 != 0.0 && rcb2 != 0.0) {
      var rab = sqrt(rab2);
      var rcb = sqrt(rcb2);
      var vp = vcb.X(vab);
      var rp = max(vp.length(), 0.000001);
      value = toDegrees(acos(min(1.0, max(-1.0, vab.dot(vcb) / (rab * rcb)))));
      var e0 = rab - bond0Eq;
      var e1 = rcb - bond1Eq;
      var dt = value - angleEq;
      var dr = force0 * e0 + force1 * e1;
      var prefactor = rigidScale * esvLambda;
      energy = prefactor * dr * dt;
      if (gradient) {
        var vdta = vab.X(vp).scaleI(-prefactor * dr * toDegrees(1.0 / (rab2 * rp)));
        var vdtc = vcb.X(vp).scaleI(prefactor * dr * toDegrees(1.0 / (rcb2 * rp)));
        var ga = vdta.addI(vab.scaleI(prefactor * force0 * dt / rab));
        var gc = vdtc.addI(vcb.scaleI(prefactor * force1 * dt / rcb));
        grad.add(threadID, ia, ga);
        grad.sub(threadID, ib, ga.add(gc));
        grad.add(threadID, ic, gc);
      }
    }
    if (esvTerm) {
      final var esvLambdaInv = (esvLambda != 0.0) ? 1 / esvLambda : 1.0;
      setEsvDeriv(energy * dedesvChain * esvLambdaInv);
    }
    return energy;
  }

  /** log */
  public void log() {
    logger.info(
        String.format(
            " %s %6d-%s %6d-%s %6d-%s" + "%7.4f %10.4f",
            "Stretch-Bend",
            atoms[0].getIndex(),
            atoms[0].getAtomType().name,
            atoms[1].getIndex(),
            atoms[1].getAtomType().name,
            atoms[2].getIndex(),
            atoms[2].getAtomType().name,
            value,
            energy));
  }

  /**
   * Setter for the field <code>rigidScale</code>.
   *
   * @param rigidScale a double.
   */
  public void setRigidScale(double rigidScale) {
    this.rigidScale = rigidScale;
  }

  /**
   * Setter for the field <code>stretchBendType</code>.
   *
   * @param stretchBendType a {@link ffx.potential.parameters.StretchBendType} object.
   */
  public void setStretchBendType(StretchBendType stretchBendType) {
    this.stretchBendType = stretchBendType;
    // Match the atom class of the angle to the atom class of the stretch-bend type.
    if (atoms[0].getAtomType().atomClass == stretchBendType.atomClasses[0]) {
      force0 = units * stretchBendType.forceConstants[0];
      force1 = units * stretchBendType.forceConstants[1];
    } else {
      force0 = units * stretchBendType.forceConstants[1];
      force1 = units * stretchBendType.forceConstants[0];
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Overidden toString Method returns the Term's id.
   */
  @Override
  public String toString() {
    return String.format(
        "%s  (%7.2f,%7.2f,%7.1f,%7.2f)", id, bonds[0].value, bonds[1].value, angle.value, energy);
  }
}
