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
import ffx.potential.bonded.UreyBradley;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 * Urey-Bradley potential energy term using {@link ffx.potential.bonded.UreyBradley} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class UreyBradleyPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(UreyBradleyPotentialEnergy.class.getName());


  /**
   * Internal list of UreyBradley instances.
   */
  private final List<UreyBradley> ureyBradleys = new ArrayList<>();

  /**
   * Create a UreyBradleyPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public UreyBradleyPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a UreyBradleyPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public UreyBradleyPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a UreyBradleyPotentialEnergy initialized with a list of UreyBradleys and force group.
   *
   * @param name         Name for this term.
   * @param forceGroup   Force group identifier.
   * @param ureyBradleys List of UreyBradley instances to add (null-safe).
   */
  public UreyBradleyPotentialEnergy(String name, int forceGroup, List<UreyBradley> ureyBradleys) {
    super(name, forceGroup);
    if (ureyBradleys != null) {
      Collections.sort(ureyBradleys);
      this.ureyBradleys.addAll(ureyBradleys);
      logger.info(format("  Urey-Bradleys:                     %10d", getNumberOfUreyBradleys()));
    }
  }

  /**
   * Get the number of terms in this potential energy term.
   *
   * @return The number of terms, which is the same as the number of Urey-Bradley interactions.
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfUreyBradleys();
  }

  /**
   * Get an array of BondedTerms in this term.
   *
   * @return Array of BondedTerms, which are actually UreyBradleys in this case.
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getUreyBradleyArray();
  }

  /**
   * Create a UreyBradleyPotentialEnergy initialized with a collection of UreyBradleys.
   *
   * @param name         Name for this term (may be null).
   * @param ureyBradleys Collection of UreyBradley instances to add (null-safe).
   */
  public UreyBradleyPotentialEnergy(String name, Collection<UreyBradley> ureyBradleys) {
    super(name);
    if (ureyBradleys != null) {
      this.ureyBradleys.addAll(ureyBradleys);
    }
  }

  /**
   * Add a UreyBradley to this term.
   *
   * @param ureyBradley UreyBradley to add (ignored if null).
   * @return true if the UreyBradley was added.
   */
  public boolean addUreyBradley(UreyBradley ureyBradley) {
    if (ureyBradley == null) {
      return false;
    }
    return ureyBradleys.add(ureyBradley);
  }

  /**
   * Add an array of UreyBradleys to this term.
   *
   * @param ureyBradleys Array of UreyBradley instances to add.
   * @return true if the UreyBradleys were added.
   */
  public boolean addUreyBradleys(UreyBradley[] ureyBradleys) {
    if (ureyBradleys == null) {
      return false;
    }
    Collections.addAll(this.ureyBradleys, ureyBradleys);
    return true;
  }

  /**
   * Add a list of UreyBradleys to this term.
   *
   * @param ureyBradleys List of UreyBradley instances to add.
   * @return true if the UreyBradleys were added.
   */
  public boolean addUreyBradleys(List<UreyBradley> ureyBradleys) {
    if (ureyBradleys == null) {
      return false;
    }
    this.ureyBradleys.addAll(ureyBradleys);
    return true;
  }

  /**
   * Remove a UreyBradley from this term.
   *
   * @param ureyBradley UreyBradley to remove (ignored if null).
   * @return true if the UreyBradley was present and removed.
   */
  public boolean removeUreyBradley(UreyBradley ureyBradley) {
    if (ureyBradley == null) {
      return false;
    }
    return ureyBradleys.remove(ureyBradley);
  }

  /**
   * Get the UreyBradley at a given index.
   *
   * @param index Index in the internal list.
   * @return UreyBradley at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public UreyBradley getUreyBradley(int index) {
    return ureyBradleys.get(index);
  }

  /**
   * Get an unmodifiable view of the UreyBradleys in this term.
   *
   * @return Unmodifiable List of UreyBradleys.
   */
  public List<UreyBradley> getUreyBradleys() {
    return Collections.unmodifiableList(ureyBradleys);
  }

  /**
   * Get an array of UreyBradleys in this term.
   *
   * @return Array of UreyBradleys.
   */
  public UreyBradley[] getUreyBradleyArray() {
    return ureyBradleys.toArray(new UreyBradley[0]);
  }

  /**
   * Get the number of UreyBradleys in this term.
   *
   * @return The number of UreyBradleys.
   */
  public int getNumberOfUreyBradleys() {
    return ureyBradleys.size();
  }

  /**
   * Log the details of Urey-Bradley interactions.
   */
  @Override
  public void log() {
    if (getNumberOfUreyBradleys() <= 0) {
      return;
    }
    logger.info("\n Urey-Bradley Interactions:");
    for (UreyBradley ureyBradley : getUreyBradleys()) {
      logger.info("Urey-Bradley \t" + ureyBradley.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfUreyBradleys() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "UREY-BRADLEY               : ", getEnergy(), getNumberOfUreyBradleys());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Urey-Bradley      ",
        getEnergy(), getNumberOfUreyBradleys(), getTime());
  }

}
