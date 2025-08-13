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
import ffx.potential.bonded.RestrainDistance;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * Restrain-Distance potential energy term using {@link ffx.potential.bonded.RestrainDistance} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainDistancePotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(RestrainDistancePotentialEnergy.class.getName());

  /**
   * Log the details of Restrain Distance interactions.
   */
  @Override
  public void log() {
    if (getNumberOfRestrainDistances() <= 0) {
      return;
    }
    logger.info("\n Restrain Distance Interactions:");
    for (RestrainDistance restrainDistance : getRestrainDistances()) {
      logger.info(" Restrain Distance \t" + restrainDistance.toString());
    }
  }

  /**
   * Internal list of RestrainDistance instances.
   */
  private final List<RestrainDistance> restrainDistances = new ArrayList<>();

  /**
   * Create a RestrainDistancePotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public RestrainDistancePotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a RestrainDistancePotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public RestrainDistancePotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a RestrainDistancePotentialEnergy initialized with a list of terms and force group.
   *
   * @param name              Name for this term.
   * @param forceGroup        Force group identifier.
   * @param restrainDistances List of RestrainDistance instances (null-safe).
   */
  public RestrainDistancePotentialEnergy(String name, int forceGroup, List<RestrainDistance> restrainDistances) {
    super(name, forceGroup);
    if (restrainDistances != null) {
      Collections.sort(restrainDistances);
      this.restrainDistances.addAll(restrainDistances);
      logger.info(String.format("  Restrain Distances:                 %10d", getNumberOfRestrainDistances()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfRestrainDistances();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getRestrainDistanceArray();
  }

  /**
   * Create a RestrainDistancePotentialEnergy initialized with a collection of terms.
   *
   * @param name              Name for this term (may be null).
   * @param restrainDistances Collection of RestrainDistance instances to add (null-safe).
   */
  public RestrainDistancePotentialEnergy(String name, Collection<RestrainDistance> restrainDistances) {
    super(name);
    if (restrainDistances != null) {
      this.restrainDistances.addAll(restrainDistances);
    }
  }

  /**
   * Add a RestrainDistance to this term.
   *
   * @param restrainDistance RestrainDistance to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addRestrainDistance(RestrainDistance restrainDistance) {
    if (restrainDistance == null) {
      return false;
    }
    return restrainDistances.add(restrainDistance);
  }

  /**
   * Add an array of RestrainDistances to this term.
   *
   * @param restrainDistances Array of RestrainDistance instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainDistances(RestrainDistance[] restrainDistances) {
    if (restrainDistances == null) {
      return false;
    }
    Collections.addAll(this.restrainDistances, restrainDistances);
    return true;
  }

  /**
   * Add a list of RestrainDistances to this term.
   *
   * @param restrainDistances List of RestrainDistance instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainDistances(List<RestrainDistance> restrainDistances) {
    if (restrainDistances == null) {
      return false;
    }
    this.restrainDistances.addAll(restrainDistances);
    return true;
  }

  /**
   * Remove a RestrainDistance from this term.
   *
   * @param restrainDistance RestrainDistance to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeRestrainDistance(RestrainDistance restrainDistance) {
    if (restrainDistance == null) {
      return false;
    }
    return restrainDistances.remove(restrainDistance);
  }

  /**
   * Get the RestrainDistance at a given index.
   *
   * @param index Index in the internal list.
   * @return RestrainDistance at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public RestrainDistance getRestrainDistance(int index) {
    return restrainDistances.get(index);
  }

  /**
   * Get an unmodifiable view of the RestrainDistances in this term.
   *
   * @return Unmodifiable List of RestrainDistances.
   */
  public List<RestrainDistance> getRestrainDistances() {
    return Collections.unmodifiableList(restrainDistances);
  }

  /**
   * Get an array of RestrainDistances in this term.
   *
   * @return Array of RestrainDistances.
   */
  public RestrainDistance[] getRestrainDistanceArray() {
    return restrainDistances.toArray(new RestrainDistance[0]);
  }

  /**
   * Get the number of RestrainDistances in this term.
   *
   * @return The number of RestrainDistances.
   */
  public int getNumberOfRestrainDistances() {
    return restrainDistances.size();
  }

  @Override
  public String toPDBString() {
    if (getNumberOfRestrainDistances() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "RESTRAIN DISTANCE          : ", getEnergy(), getNumberOfRestrainDistances());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Restrain Distance ",
        getEnergy(), getNumberOfRestrainDistances(), getTime());
  }
}
