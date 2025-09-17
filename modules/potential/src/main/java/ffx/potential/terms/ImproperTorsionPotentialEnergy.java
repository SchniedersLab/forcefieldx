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

import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.ImproperTorsion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 * Improper Torsion potential energy term using {@link ffx.potential.bonded.ImproperTorsion} instances.
 */
public class ImproperTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(ImproperTorsionPotentialEnergy.class.getName());


  /**
   * Internal list of ImproperTorsion instances.
   */
  private final List<ImproperTorsion> improperTorsions = new ArrayList<>();

  /**
   * Create an ImproperTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public ImproperTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create an ImproperTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public ImproperTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create an ImproperTorsionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name             Name for this term.
   * @param forceGroup       Force group identifier.
   * @param improperTorsions List of ImproperTorsion instances to add (null-safe).
   */
  public ImproperTorsionPotentialEnergy(String name, int forceGroup, List<ImproperTorsion> improperTorsions) {
    super(name, forceGroup);
    if (improperTorsions != null) {
      Collections.sort(improperTorsions);
      this.improperTorsions.addAll(improperTorsions);
      logger.info(format("  Improper Torsions:                 %10d", getNumberOfImproperTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfImproperTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getImproperTorsionArray();
  }

  /**
   * Create an ImproperTorsionPotentialEnergy initialized with a collection of terms.
   *
   * @param name             Name for this term (may be null).
   * @param improperTorsions Collection of ImproperTorsion instances to add (null-safe).
   */
  public ImproperTorsionPotentialEnergy(String name, Collection<ImproperTorsion> improperTorsions) {
    super(name);
    if (improperTorsions != null) {
      this.improperTorsions.addAll(improperTorsions);
    }
  }

  /**
   * Add an ImproperTorsion to this term.
   *
   * @param improperTorsion ImproperTorsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addImproperTorsion(ImproperTorsion improperTorsion) {
    if (improperTorsion == null) {
      return false;
    }
    return improperTorsions.add(improperTorsion);
  }

  /**
   * Add an array of ImproperTorsions to this term.
   *
   * @param improperTorsions Array of ImproperTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addImproperTorsions(ImproperTorsion[] improperTorsions) {
    if (improperTorsions == null) {
      return false;
    }
    Collections.addAll(this.improperTorsions, improperTorsions);
    return true;
  }

  /**
   * Add a list of ImproperTorsions to this term.
   *
   * @param improperTorsions List of ImproperTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addImproperTorsions(List<ImproperTorsion> improperTorsions) {
    if (improperTorsions == null) {
      return false;
    }
    this.improperTorsions.addAll(improperTorsions);
    return true;
  }

  /**
   * Remove an ImproperTorsion from this term.
   *
   * @param improperTorsion ImproperTorsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeImproperTorsion(ImproperTorsion improperTorsion) {
    if (improperTorsion == null) {
      return false;
    }
    return improperTorsions.remove(improperTorsion);
  }

  /**
   * Get the ImproperTorsion at a given index.
   *
   * @param index Index in the internal list.
   * @return ImproperTorsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public ImproperTorsion getImproperTorsion(int index) {
    return improperTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of the ImproperTorsions in this term.
   *
   * @return Unmodifiable List of ImproperTorsions.
   */
  public List<ImproperTorsion> getImproperTorsions() {
    return Collections.unmodifiableList(improperTorsions);
  }

  /**
   * Get an array of ImproperTorsions in this term.
   *
   * @return Array of ImproperTorsions.
   */
  public ImproperTorsion[] getImproperTorsionArray() {
    return improperTorsions.toArray(new ImproperTorsion[0]);
  }

  /**
   * Get the number of ImproperTorsions in this term.
   *
   * @return The number of ImproperTorsions.
   */
  public int getNumberOfImproperTorsions() {
    return improperTorsions.size();
  }

  /**
   * Log the details of Improper interactions.
   */
  @Override
  public void log() {
    if (getNumberOfImproperTorsions() <= 0) {
      return;
    }
    logger.info("\n Improper Interactions:");
    for (ImproperTorsion improperTorsion : getImproperTorsions()) {
      logger.info(" Improper \t" + improperTorsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfImproperTorsions() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "IMPROPER TORSION           : ", getEnergy(), getNumberOfImproperTorsions());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Improper Torsion  ",
        getEnergy(), getNumberOfImproperTorsions(), getTime());
  }
}
