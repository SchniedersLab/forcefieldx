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
import ffx.potential.bonded.StretchTorsion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 * Stretch-Torsion potential energy term using {@link ffx.potential.bonded.StretchTorsion} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class StretchTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(StretchTorsionPotentialEnergy.class.getName());


  /**
   * Internal list of StretchTorsion instances.
   */
  private final List<StretchTorsion> stretchTorsions = new ArrayList<>();

  /**
   * Create a StretchTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public StretchTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a StretchTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public StretchTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a StretchTorsionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name            Name for this term.
   * @param forceGroup      Force group identifier.
   * @param stretchTorsions List of StretchTorsion instances to add (null-safe).
   */
  public StretchTorsionPotentialEnergy(String name, int forceGroup, List<StretchTorsion> stretchTorsions) {
    super(name, forceGroup);
    if (stretchTorsions != null) {
      Collections.sort(stretchTorsions);
      this.stretchTorsions.addAll(stretchTorsions);
      logger.info(format("  Stretch-Torsions:                  %10d", getNumberOfStretchTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfStretchTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getStretchTorsionArray();
  }

  /**
   * Create a StretchTorsionPotentialEnergy initialized with a collection of terms.
   *
   * @param name            Name for this term (may be null).
   * @param stretchTorsions Collection of StretchTorsion instances to add (null-safe).
   */
  public StretchTorsionPotentialEnergy(String name, Collection<StretchTorsion> stretchTorsions) {
    super(name);
    if (stretchTorsions != null) {
      this.stretchTorsions.addAll(stretchTorsions);
    }
  }

  /**
   * Add a StretchTorsion to this term.
   *
   * @param stretchTorsion StretchTorsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addStretchTorsion(StretchTorsion stretchTorsion) {
    if (stretchTorsion == null) {
      return false;
    }
    return stretchTorsions.add(stretchTorsion);
  }

  /**
   * Add an array of StretchTorsions to this term.
   *
   * @param stretchTorsions Array of StretchTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addStretchTorsions(StretchTorsion[] stretchTorsions) {
    if (stretchTorsions == null) {
      return false;
    }
    Collections.addAll(this.stretchTorsions, stretchTorsions);
    return true;
  }

  /**
   * Add a list of StretchTorsions to this term.
   *
   * @param stretchTorsions List of StretchTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addStretchTorsions(List<StretchTorsion> stretchTorsions) {
    if (stretchTorsions == null) {
      return false;
    }
    this.stretchTorsions.addAll(stretchTorsions);
    return true;
  }

  /**
   * Remove a StretchTorsion from this term.
   *
   * @param stretchTorsion StretchTorsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeStretchTorsion(StretchTorsion stretchTorsion) {
    if (stretchTorsion == null) {
      return false;
    }
    return stretchTorsions.remove(stretchTorsion);
  }

  /**
   * Get the StretchTorsion at a given index.
   *
   * @param index Index in the internal list.
   * @return StretchTorsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public StretchTorsion getStretchTorsion(int index) {
    return stretchTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of the StretchTorsions in this term.
   *
   * @return Unmodifiable List of StretchTorsions.
   */
  public List<StretchTorsion> getStretchTorsions() {
    return Collections.unmodifiableList(stretchTorsions);
  }

  /**
   * Get an array of StretchTorsions in this term.
   *
   * @return Array of StretchTorsions.
   */
  public StretchTorsion[] getStretchTorsionArray() {
    return stretchTorsions.toArray(new StretchTorsion[0]);
  }

  /**
   * Get the number of StretchTorsions in this term.
   *
   * @return The number of StretchTorsions.
   */
  public int getNumberOfStretchTorsions() {
    return stretchTorsions.size();
  }

  /**
   * Log the details of Stretch-Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfStretchTorsions() <= 0) {
      return;
    }
    logger.info("\n Stretch-Torsion Interactions:");
    for (StretchTorsion stretchTorsion : getStretchTorsions()) {
      logger.info(" Stretch-Torsion \t" + stretchTorsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfStretchTorsions() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "STRETCH-TORSION            : ", getEnergy(), getNumberOfStretchTorsions());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Stretch-Torsion   ",
        getEnergy(), getNumberOfStretchTorsions(), getTime());
  }
}
