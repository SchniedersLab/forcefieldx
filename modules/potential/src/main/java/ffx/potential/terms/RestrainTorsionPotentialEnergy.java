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
import ffx.potential.bonded.Torsion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * Restrain-Torsion potential energy term using {@link ffx.potential.bonded.Torsion} instances.
 * Method names are specific to restrain torsions for clarity.
 */
public class RestrainTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(RestrainTorsionPotentialEnergy.class.getName());

  /**
   * Internal list of Torsion instances for restrain torsions.
   */
  private final List<Torsion> restrainTorsions = new ArrayList<>();

  /**
   * Create a RestrainTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public RestrainTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a RestrainTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public RestrainTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a RestrainTorsionPotentialEnergy initialized with a list of torsions and force group.
   *
   * @param name             Name for this term.
   * @param forceGroup       Force group identifier.
   * @param restrainTorsions List of Torsion instances to add (null-safe).
   */
  public RestrainTorsionPotentialEnergy(String name, int forceGroup, List<Torsion> restrainTorsions) {
    super(name, forceGroup);
    if (restrainTorsions != null) {
      Collections.sort(restrainTorsions);
      this.restrainTorsions.addAll(restrainTorsions);
      logger.info(String.format("  Restrain Torsions:                   %10d", getNumberOfRestrainTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfRestrainTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getRestrainTorsionArray();
  }

  /**
   * Create a RestrainTorsionPotentialEnergy initialized with a collection of torsions.
   *
   * @param name             Name for this term (may be null).
   * @param restrainTorsions Collection of Torsion instances to add (null-safe).
   */
  public RestrainTorsionPotentialEnergy(String name, Collection<Torsion> restrainTorsions) {
    super(name);
    if (restrainTorsions != null) {
      this.restrainTorsions.addAll(restrainTorsions);
    }
  }

  /**
   * Add a single restrain torsion.
   *
   * @param torsion Torsion to add (ignored if null).
   * @return true if added.
   */
  public boolean addRestrainTorsion(Torsion torsion) {
    if (torsion == null) {
      return false;
    }
    return restrainTorsions.add(torsion);

  }

  /**
   * Add an array of restrain torsions.
   *
   * @param torsions Array of Torsion instances to add.
   * @return true if added.
   */
  public boolean addRestrainTorsions(Torsion[] torsions) {
    if (torsions == null) {
      return false;
    }
    Collections.addAll(this.restrainTorsions, torsions);
    return true;
  }

  /**
   * Add a list of restrain torsions.
   *
   * @param torsions List of Torsion instances to add.
   * @return true if added.
   */
  public boolean addRestrainTorsions(List<Torsion> torsions) {
    if (torsions == null) {
      return false;
    }
    this.restrainTorsions.addAll(torsions);
    return true;
  }

  /**
   * Remove a restrain torsion.
   *
   * @param torsion Torsion to remove (ignored if null).
   * @return true if removed.
   */
  public boolean removeRestrainTorsion(Torsion torsion) {
    if (torsion == null) {
      return false;
    }
    return restrainTorsions.remove(torsion);
  }

  /**
   * Get a restrain torsion at an index.
   *
   * @param index Index.
   * @return Torsion.
   * @throws IndexOutOfBoundsException if invalid index.
   */
  public Torsion getRestrainTorsion(int index) {
    return restrainTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of restrain torsions.
   *
   * @return List of Torsion.
   */
  public List<Torsion> getRestrainTorsions() {
    return Collections.unmodifiableList(restrainTorsions);
  }

  /**
   * Get an array of restrain torsions.
   *
   * @return Array of Torsion.
   */
  public Torsion[] getRestrainTorsionArray() {
    return restrainTorsions.toArray(new Torsion[0]);
  }

  /**
   * Get the number of restrain torsions.
   *
   * @return Number of restrain torsions.
   */
  public int getNumberOfRestrainTorsions() {
    return restrainTorsions.size();
  }

  @Override
  public String toPDBString() {
    if (getNumberOfRestrainTorsions() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "RESTRAIN TORSION           : ", getEnergy(), getNumberOfRestrainTorsions());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Restrain Torsion  ",
        getEnergy(), getNumberOfRestrainTorsions(), getTime());
  }

  /**
   * Log the details of Restrain Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfRestrainTorsions() <= 0) {
      return;
    }
    logger.info("\n Restrain Torsion Interactions:");
    for (Torsion torsion : getRestrainTorsions()) {
      logger.info(" Restrain Torsion \t" + torsion.toString());
    }
  }
}
