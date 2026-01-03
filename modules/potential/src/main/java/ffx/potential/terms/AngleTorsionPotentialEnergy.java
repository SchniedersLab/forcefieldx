// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.BondedTerm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Angle-Torsion potential energy term using {@link ffx.potential.bonded.AngleTorsion} instances.
 */
public class AngleTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(AngleTorsionPotentialEnergy.class.getName());


  /**
   * Internal list of AngleTorsion instances.
   */
  private final List<AngleTorsion> angleTorsions = new ArrayList<>();

  /**
   * Create an AngleTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public AngleTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create an AngleTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public AngleTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create an AngleTorsionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name          Name for this term.
   * @param forceGroup    Force group identifier.
   * @param angleTorsions List of AngleTorsion instances to add (null-safe).
   */
  public AngleTorsionPotentialEnergy(String name, int forceGroup, List<AngleTorsion> angleTorsions) {
    super(name, forceGroup);
    if (angleTorsions != null) {
      Collections.sort(angleTorsions);
      this.angleTorsions.addAll(angleTorsions);
      logger.info(format("  Angle-Torsions:                    %10d", getNumberOfAngleTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfAngleTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getAngleTorsionArray();
  }

  /**
   * Create an AngleTorsionPotentialEnergy initialized with a collection of terms.
   *
   * @param name          Name for this term (may be null).
   * @param angleTorsions Collection of AngleTorsion instances to add (null-safe).
   */
  public AngleTorsionPotentialEnergy(String name, Collection<AngleTorsion> angleTorsions) {
    super(name);
    if (angleTorsions != null) {
      this.angleTorsions.addAll(angleTorsions);
    }
  }

  /**
   * Add an AngleTorsion to this term.
   *
   * @param angleTorsion AngleTorsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addAngleTorsion(AngleTorsion angleTorsion) {
    if (angleTorsion == null) {
      return false;
    }
    return angleTorsions.add(angleTorsion);
  }

  /**
   * Add an array of AngleTorsions to this term.
   *
   * @param angleTorsions Array of AngleTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addAngleTorsions(AngleTorsion[] angleTorsions) {
    if (angleTorsions == null) {
      return false;
    }
    Collections.addAll(this.angleTorsions, angleTorsions);
    return true;
  }

  /**
   * Add a list of AngleTorsions to this term.
   *
   * @param angleTorsions List of AngleTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addAngleTorsions(List<AngleTorsion> angleTorsions) {
    if (angleTorsions == null) {
      return false;
    }
    this.angleTorsions.addAll(angleTorsions);
    return true;
  }

  /**
   * Remove an AngleTorsion from this term.
   *
   * @param angleTorsion AngleTorsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeAngleTorsion(AngleTorsion angleTorsion) {
    if (angleTorsion == null) {
      return false;
    }
    return angleTorsions.remove(angleTorsion);
  }

  /**
   * Get the AngleTorsion at a given index.
   *
   * @param index Index in the internal list.
   * @return AngleTorsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public AngleTorsion getAngleTorsion(int index) {
    return angleTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of the AngleTorsions in this term.
   *
   * @return Unmodifiable List of AngleTorsions.
   */
  public List<AngleTorsion> getAngleTorsions() {
    return Collections.unmodifiableList(angleTorsions);
  }

  /**
   * Get an array of AngleTorsions in this term.
   *
   * @return Array of AngleTorsions.
   */
  public AngleTorsion[] getAngleTorsionArray() {
    return angleTorsions.toArray(new AngleTorsion[0]);
  }

  /**
   * Get the number of AngleTorsions in this term.
   *
   * @return The number of AngleTorsions.
   */
  public int getNumberOfAngleTorsions() {
    return angleTorsions.size();
  }

  /**
   * Log the details of Angle-Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfAngleTorsions() <= 0) {
      return;
    }
    logger.info("\n Angle-Torsion Interactions:");
    for (AngleTorsion angleTorsion : getAngleTorsions()) {
      logger.info(" Angle-Torsion \t" + angleTorsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfAngleTorsions() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "ANGLE-TORSION              : ", getEnergy(), getNumberOfAngleTorsions());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Angle-Torsion     ",
        getEnergy(), getNumberOfAngleTorsions(), getTime());
  }
}
