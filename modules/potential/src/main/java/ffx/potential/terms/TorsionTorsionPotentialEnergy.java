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
import ffx.potential.bonded.TorsionTorsion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * Torsion-Torsion potential energy term using {@link ffx.potential.bonded.TorsionTorsion} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TorsionTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(TorsionTorsionPotentialEnergy.class.getName());


  /**
   * Internal list of TorsionTorsion instances.
   */
  private final List<TorsionTorsion> torsionTorsions = new ArrayList<>();

  /**
   * Create a TorsionTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public TorsionTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a TorsionTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public TorsionTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a TorsionTorsionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name            Name for this term.
   * @param forceGroup      Force group identifier.
   * @param torsionTorsions List of TorsionTorsion instances to add (null-safe).
   */
  public TorsionTorsionPotentialEnergy(String name, int forceGroup, List<TorsionTorsion> torsionTorsions) {
    super(name, forceGroup);
    if (torsionTorsions != null) {
      Collections.sort(torsionTorsions);
      this.torsionTorsions.addAll(torsionTorsions);
      logger.info(String.format("  Torsion-Torsions:                  %10d", getNumberOfTorsionTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfTorsionTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getTorsionTorsionArray();
  }

  /**
   * Create a TorsionTorsionPotentialEnergy initialized with a collection of terms.
   *
   * @param name            Name for this term (may be null).
   * @param torsionTorsions Collection of TorsionTorsion instances to add (null-safe).
   */
  public TorsionTorsionPotentialEnergy(String name, Collection<TorsionTorsion> torsionTorsions) {
    super(name);
    if (torsionTorsions != null) {
      this.torsionTorsions.addAll(torsionTorsions);
    }
  }

  /**
   * Add a TorsionTorsion to this term.
   *
   * @param torsionTorsion TorsionTorsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addTorsionTorsion(TorsionTorsion torsionTorsion) {
    if (torsionTorsion == null) {
      return false;
    }
    return torsionTorsions.add(torsionTorsion);
  }

  /**
   * Add an array of TorsionTorsions to this term.
   *
   * @param torsionTorsions Array of TorsionTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addTorsionTorsions(TorsionTorsion[] torsionTorsions) {
    if (torsionTorsions == null) {
      return false;
    }
    Collections.addAll(this.torsionTorsions, torsionTorsions);
    return true;
  }

  /**
   * Add a list of TorsionTorsions to this term.
   *
   * @param torsionTorsions List of TorsionTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addTorsionTorsions(List<TorsionTorsion> torsionTorsions) {
    if (torsionTorsions == null) {
      return false;
    }
    this.torsionTorsions.addAll(torsionTorsions);
    return true;
  }

  /**
   * Remove a TorsionTorsion from this term.
   *
   * @param torsionTorsion TorsionTorsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeTorsionTorsion(TorsionTorsion torsionTorsion) {
    if (torsionTorsion == null) {
      return false;
    }
    return torsionTorsions.remove(torsionTorsion);
  }

  /**
   * Get the TorsionTorsion at a given index.
   *
   * @param index Index in the internal list.
   * @return TorsionTorsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public TorsionTorsion getTorsionTorsion(int index) {
    return torsionTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of the TorsionTorsions in this term.
   *
   * @return Unmodifiable List of TorsionTorsions.
   */
  public List<TorsionTorsion> getTorsionTorsions() {
    return Collections.unmodifiableList(torsionTorsions);
  }

  /**
   * Get an array of TorsionTorsions in this term.
   *
   * @return Array of TorsionTorsions.
   */
  public TorsionTorsion[] getTorsionTorsionArray() {
    return torsionTorsions.toArray(new TorsionTorsion[0]);
  }

  /**
   * Get the number of TorsionTorsions in this term.
   *
   * @return The number of TorsionTorsions.
   */
  public int getNumberOfTorsionTorsions() {
    return torsionTorsions.size();
  }

  /**
   * Log the details of Torsion-Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfTorsionTorsions() <= 0) {
      return;
    }
    logger.info("\n Torsion-Torsion Interactions:");
    for (TorsionTorsion torsionTorsion : getTorsionTorsions()) {
      logger.info(" Torsion-Torsion \t" + torsionTorsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfTorsionTorsions() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "TORSION-TORSION            : ", getEnergy(), getNumberOfTorsionTorsions());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Torsion-Torsion   ",
        getEnergy(), getNumberOfTorsionTorsions(), getTime());
  }
}
