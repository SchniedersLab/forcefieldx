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

import static ffx.potential.parameters.UreyBradleyType.cubic;
import static ffx.potential.parameters.UreyBradleyType.quartic;
import static ffx.potential.parameters.UreyBradleyType.units;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.UreyBradleyType;
import java.util.logging.Logger;

/**
 * The UreyBradley class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class UreyBradley extends BondedTerm {

  private static final Logger logger = Logger.getLogger(UreyBradley.class.getName());

  /** Force field parameters to compute the Stretch-Bend energy. */
  public final UreyBradleyType ureyBradleyType;
  /** The Angle this UreyBradley term is based on. */
  protected final Angle angle;
  /** Scale factor to apply to Urey-Bradley term. */
  private double rigidScale = 1.0;

  /**
   * Constructor for the UreyBradley class.
   *
   * @param a a {@link ffx.potential.bonded.Angle} object.
   * @param ureyBradleyType The Urey-Bradley type to apply.
   */
  public UreyBradley(Angle a, UreyBradleyType ureyBradleyType) {
    super();
    angle = a;
    bonds = a.bonds;
    atoms = a.atoms;
    this.ureyBradleyType = ureyBradleyType;
    setID_Key(false);
  }

  /**
   * Attempt to create a new UreyBradley for the specified Angle.
   *
   * @param angle the Angle to create the UreyBradley from.
   * @param forceField the ForceField parameters to apply.
   * @return a new UreyBradley, or null.
   */
  public static UreyBradley ureyBradlyFactory(Angle angle, ForceField forceField) {
    if (angle == null) {
      return null;
    }
    UreyBradleyType ureyBradleyType = forceField.getUreyBradleyType(angle.angleType.getKey());
    if (ureyBradleyType == null) {
      return null;
    }
    return new UreyBradley(angle, ureyBradleyType);
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(BondedTerm ub) {
    if (!ub.getClass().isInstance(this)) {
      return super.compareTo(ub);
    }
    return angle.compareTo(((UreyBradley) ub).angle);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Evaluate the Urey-Bradley energy.
   */
  @Override
  public double energy(
      boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    var atomA = atoms[0];
    var atomC = atoms[2];
    var va = atomA.getXYZ();
    var vc = atomC.getXYZ();
    var vac = va.sub(vc);
    value = vac.length();
    var dv = value - ureyBradleyType.distance;
    var dv2 = dv * dv;
    energy =
        units
            * rigidScale
            * ureyBradleyType.forceConstant
            * dv2
            * (1.0 + cubic * dv + quartic * dv2)
            * esvLambda;
    if (gradient) {
      var deddt =
          2.0
              * units
              * rigidScale
              * ureyBradleyType.forceConstant
              * dv
              * (1.0 + 1.5 * cubic * dv + 2.0 * quartic * dv2)
              * esvLambda;
      var de = 0.0;
      if (value > 0.0) {
        de = deddt / value;
      }
      var ia = atomA.getIndex() - 1;
      var ic = atomC.getIndex() - 1;
      grad.add(threadID, ia, vac.scaleI(de));
      grad.sub(threadID, ic, vac);
    }
    if (esvTerm) {
      var esvLambdaInv = (esvLambda != 0.0) ? 1 / esvLambda : 1.0;
      setEsvDeriv(energy * dedesvChain * esvLambdaInv);
    }
    return energy;
  }

  /** log */
  public void log() {
    logger.info(
        String.format(
            " %s %6d-%s %6d-%s %6.4f  %6.4f  %10.4f",
            "Urey-Bradley",
            atoms[0].getIndex(),
            atoms[0].getAtomType().name,
            atoms[2].getIndex(),
            atoms[2].getAtomType().name,
            ureyBradleyType.distance,
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
}
