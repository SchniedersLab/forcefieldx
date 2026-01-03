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

import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.RestrainPosition;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * Restrain-Position potential energy term using {@link ffx.potential.bonded.RestrainPosition} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainPositionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(RestrainPositionPotentialEnergy.class.getName());


  /**
   * Internal list of RestrainPosition instances.
   */
  private final List<RestrainPosition> restrainPositions = new ArrayList<>();

  /**
   * Create a RestrainPositionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public RestrainPositionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a RestrainPositionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public RestrainPositionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a RestrainPositionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name              Name for this term.
   * @param forceGroup        Force group identifier.
   * @param restrainPositions List of RestrainPosition instances (null-safe).
   */
  public RestrainPositionPotentialEnergy(String name, int forceGroup, List<RestrainPosition> restrainPositions) {
    super(name, forceGroup);
    if (restrainPositions != null) {
      Collections.sort(restrainPositions);
      this.restrainPositions.addAll(restrainPositions);
      logger.info(String.format("  Restrain Positions:                 %10d", getNumberOfRestrainPositions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfRestrainPositions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getRestrainPositionArray();
  }

  /**
   * Create a RestrainPositionPotentialEnergy initialized with a collection of terms.
   *
   * @param name              Name for this term (may be null).
   * @param restrainPositions Collection of RestrainPosition instances to add (null-safe).
   */
  public RestrainPositionPotentialEnergy(String name, Collection<RestrainPosition> restrainPositions) {
    super(name);
    if (restrainPositions != null) {
      this.restrainPositions.addAll(restrainPositions);
    }
  }

  /**
   * Add a RestrainPosition to this term.
   *
   * @param restrainPosition RestrainPosition to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addRestrainPosition(RestrainPosition restrainPosition) {
    if (restrainPosition == null) {
      return false;
    }
    return restrainPositions.add(restrainPosition);
  }

  /**
   * Add an array of RestrainPositions to this term.
   *
   * @param restrainPositions Array of RestrainPosition instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainPositions(RestrainPosition[] restrainPositions) {
    if (restrainPositions == null) {
      return false;
    }
    Collections.addAll(this.restrainPositions, restrainPositions);
    return true;
  }

  /**
   * Add a list of RestrainPositions to this term.
   *
   * @param restrainPositions List of RestrainPosition instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainPositions(List<RestrainPosition> restrainPositions) {
    if (restrainPositions == null) {
      return false;
    }
    this.restrainPositions.addAll(restrainPositions);
    return true;
  }

  /**
   * Remove a RestrainPosition from this term.
   *
   * @param restrainPosition RestrainPosition to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeRestrainPosition(RestrainPosition restrainPosition) {
    if (restrainPosition == null) {
      return false;
    }
    return restrainPositions.remove(restrainPosition);
  }

  /**
   * Get the RestrainPosition at a given index.
   *
   * @param index Index in the internal list.
   * @return RestrainPosition at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public RestrainPosition getRestrainPosition(int index) {
    return restrainPositions.get(index);
  }

  /**
   * Get an unmodifiable view of the RestrainPositions in this term.
   *
   * @return Unmodifiable List of RestrainPositions.
   */
  public List<RestrainPosition> getRestrainPositions() {
    return Collections.unmodifiableList(restrainPositions);
  }

  /**
   * Get an array of RestrainPositions in this term.
   *
   * @return Array of RestrainPositions.
   */
  public RestrainPosition[] getRestrainPositionArray() {
    return restrainPositions.toArray(new RestrainPosition[0]);
  }

  /**
   * Get the number of RestrainPositions in this term.
   *
   * @return The number of RestrainPositions.
   */
  public int getNumberOfRestrainPositions() {
    return restrainPositions.size();
  }

  /**
   * Get the mathematical form of the Restrain Position interaction.
   * @return The mathematical form of the Restrain Position interaction.
   */
  public static String getRestrainPositionEnergyString() {
    return "k0*periodicdistance(x,y,z,x0,y0,z0)^2";
  }

  /**
   * Log the details of Restrain Position interactions.
   */
  @Override
  public void log() {
    if (getNumberOfRestrainPositions() <= 0) {
      return;
    }
    logger.info("\n Restrain Position Interactions:");
    for (RestrainPosition restrainPosition : getRestrainPositions()) {
      logger.info(" Restrain Position \t" + restrainPosition.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfRestrainPositions() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "COORDINATE RESTRAINTS      : ", getEnergy(), getNumberOfRestrainPositions());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Restrain Position ",
        getEnergy(), getNumberOfRestrainPositions(), getTime());
  }
}
