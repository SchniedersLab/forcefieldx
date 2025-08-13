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
package ffx.potential.terms;

import ffx.potential.bonded.Angle;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.parameters.AngleType;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;

/**
 * Angle potential energy term using {@link ffx.potential.bonded.Angle} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AnglePotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(AnglePotentialEnergy.class.getName());

  /**
   * Internal list of Angle instances.
   */
  private final List<Angle> angles = new ArrayList<>();

  /**
   * Create an AnglePotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public AnglePotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create an AnglePotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public AnglePotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create an AnglePotentialEnergy initialized with a list of angles and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Force group identifier.
   * @param angles     List of Angle instances to add (null-safe).
   */
  public AnglePotentialEnergy(String name, int forceGroup, List<Angle> angles) {
    super(name, forceGroup);
    if (angles != null) {
      Collections.sort(angles);
      this.angles.addAll(angles);
      logger.info(format("  Angles:                            %10d", getNumberOfAngles()));
    }
  }

  /**
   * Get the number of terms in this potential energy term.
   *
   * @return The number of terms, which is the same as the number of angles.
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfAngles();
  }

  /**
   * Get an array of BondedTerms in this term.
   *
   * @return Array of BondedTerms, which are actually Angles in this case.
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getAngleArray();
  }

  /**
   * Create an AnglePotentialEnergy initialized with a collection of angles.
   *
   * @param name   Name for this term (may be null).
   * @param angles Collection of Angle instances to add (null-safe).
   */
  public AnglePotentialEnergy(String name, Collection<Angle> angles) {
    super(name);
    if (angles != null) {
      this.angles.addAll(angles);
    }
  }

  /**
   * Add an Angle to this term.
   *
   * @param angle Angle to add (ignored if null).
   * @return true if the angle was added.
   */
  public boolean addAngle(Angle angle) {
    if (angle == null) {
      return false;
    }
    return angles.add(angle);
  }

  /**
   * Add an array of Angles to this term.
   *
   * @param angles Array of Angle instances to add.
   * @return true if the angles were added.
   */
  public boolean addAngles(Angle[] angles) {
    if (angles == null) {
      return false;
    }
    Collections.addAll(this.angles, angles);
    return true;
  }

  /**
   * Add a list of Angles to this term.
   *
   * @param angles List of Angle instances to add.
   * @return true if the angles were added.
   */
  public boolean addAngles(List<Angle> angles) {
    if (angles == null) {
      return false;
    }
    this.angles.addAll(angles);
    return true;
  }

  /**
   * Remove an Angle from this term.
   *
   * @param angle Angle to remove (ignored if null).
   * @return true if the angle was present and removed.
   */
  public boolean removeAngle(Angle angle) {
    if (angle == null) {
      return false;
    }
    return angles.remove(angle);
  }

  /**
   * Get the Angle at a given index.
   *
   * @param index Index in the internal list.
   * @return Angle at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public Angle getAngle(int index) {
    return angles.get(index);
  }

  /**
   * Get an unmodifiable view of the Angles in this term.
   *
   * @return Unmodifiable List of Angles.
   */
  public List<Angle> getAngles() {
    return Collections.unmodifiableList(angles);
  }

  /**
   * Get an array of Angles in this term.
   *
   * @return Array of Angles.
   */
  public Angle[] getAngleArray() {
    return angles.toArray(new Angle[0]);
  }

  /**
   * Get the number of Angles in this term.
   *
   * @return The number of Angles.
   */
  public int getNumberOfAngles() {
    return angles.size();
  }

  /**
   * Get the String used for OpenMM angle energy expressions.
   * @return String representing the angle energy expression.
   */
  public String getAngleEnergyString() {
    AngleType angleType = angles.getFirst().angleType;
    String energy;
    if (angleType.angleFunction == AngleType.AngleFunction.SEXTIC) {
      energy = format("""
              k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6);
              d=%.15g*theta-theta0;
              """,
          angleType.cubic, angleType.quartic, angleType.pentic, angleType.sextic, 180.0 / PI);
    } else {
      energy = format("""
              k*(d^2);
              d=%.15g*theta-theta0;
              """,
          180.0 / PI);
    }
    return energy;
  }

  public String getInPlaneAngleEnergyString() {
    AngleType angleType = angles.getFirst().angleType;
    String energy = format("""
            k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6);
            d=theta-theta0;
            theta = %.15g*pointangle(x1, y1, z1, projx, projy, projz, x3, y3, z3);
            projx = x2-nx*dot;
            projy = y2-ny*dot;
            projz = z2-nz*dot;
            dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
            nx = px/norm;
            ny = py/norm;
            nz = pz/norm;
            norm = sqrt(px*px + py*py + pz*pz);
            px = (d1y*d2z-d1z*d2y);
            py = (d1z*d2x-d1x*d2z);
            pz = (d1x*d2y-d1y*d2x);
            d1x = x1-x4;
            d1y = y1-y4;
            d1z = z1-z4;
            d2x = x3-x4;
            d2y = y3-y4;
            d2z = z3-z4;
            """,
        angleType.cubic, angleType.quartic, angleType.pentic, angleType.sextic, 180.0 / PI);
    return energy;
  }

  /**
   * Log the details of Angle interactions.
   */
  @Override
  public void log() {
    if (getNumberOfAngles() <= 0) {
      return;
    }
    logger.info("\n Angle Bending Interactions:");
    for (Angle angle : getAngles()) {
      logger.info(" Angle \t" + angle.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfAngles() <= 0) {
      return "";
    }
    StringBuilder sb = new StringBuilder();
    sb.append(format("REMARK   3   %s %g (%d)\n", "ANGLE BENDING              : ", getEnergy(), getNumberOfAngles()));
    sb.append(format("REMARK   3   %s %g\n", "ANGLE RMSD                 : ", getRMSD()));
    return sb.toString();
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f (%8.5f)\n", "Angle Bending     ",
        getEnergy(), getNumberOfAngles(), getTime(), getRMSD());
  }

}