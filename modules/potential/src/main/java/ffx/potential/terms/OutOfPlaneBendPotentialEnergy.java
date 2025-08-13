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
import ffx.potential.bonded.OutOfPlaneBend;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 * Out-of-Plane Bend potential energy term using {@link ffx.potential.bonded.OutOfPlaneBend} instances.
 */
public class OutOfPlaneBendPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(OutOfPlaneBendPotentialEnergy.class.getName());


  /**
   * Internal list of OutOfPlaneBend instances.
   */
  private final List<OutOfPlaneBend> outOfPlaneBends = new ArrayList<>();

  /**
   * Create an OutOfPlaneBendPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public OutOfPlaneBendPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create an OutOfPlaneBendPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public OutOfPlaneBendPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create an OutOfPlaneBendPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name            Name for this term.
   * @param forceGroup      Force group identifier.
   * @param outOfPlaneBends List of OutOfPlaneBend instances to add (null-safe).
   */
  public OutOfPlaneBendPotentialEnergy(String name, int forceGroup, List<OutOfPlaneBend> outOfPlaneBends) {
    super(name, forceGroup);
    if (outOfPlaneBends != null) {
      Collections.sort(outOfPlaneBends);
      this.outOfPlaneBends.addAll(outOfPlaneBends);
      logger.info(format("  Out-of-Plane Bends:                %10d", getNumberOfOutOfPlaneBends()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfOutOfPlaneBends();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getOutOfPlaneBendArray();
  }

  /**
   * Create an OutOfPlaneBendPotentialEnergy initialized with a collection of terms.
   *
   * @param name            Name for this term (may be null).
   * @param outOfPlaneBends Collection of OutOfPlaneBend instances to add (null-safe).
   */
  public OutOfPlaneBendPotentialEnergy(String name, Collection<OutOfPlaneBend> outOfPlaneBends) {
    super(name);
    if (outOfPlaneBends != null) {
      this.outOfPlaneBends.addAll(outOfPlaneBends);
    }
  }

  /**
   * Add an OutOfPlaneBend to this term.
   *
   * @param outOfPlaneBend OutOfPlaneBend to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addOutOfPlaneBend(OutOfPlaneBend outOfPlaneBend) {
    if (outOfPlaneBend == null) {
      return false;
    }
    return outOfPlaneBends.add(outOfPlaneBend);
  }

  /**
   * Add an array of OutOfPlaneBends to this term.
   *
   * @param outOfPlaneBends Array of OutOfPlaneBend instances to add.
   * @return true if they were added.
   */
  public boolean addOutOfPlaneBends(OutOfPlaneBend[] outOfPlaneBends) {
    if (outOfPlaneBends == null) {
      return false;
    }
    Collections.addAll(this.outOfPlaneBends, outOfPlaneBends);
    return true;
  }

  /**
   * Add a list of OutOfPlaneBends to this term.
   *
   * @param outOfPlaneBends List of OutOfPlaneBend instances to add.
   * @return true if they were added.
   */
  public boolean addOutOfPlaneBends(List<OutOfPlaneBend> outOfPlaneBends) {
    if (outOfPlaneBends == null) {
      return false;
    }
    this.outOfPlaneBends.addAll(outOfPlaneBends);
    return true;
  }

  /**
   * Remove an OutOfPlaneBend from this term.
   *
   * @param outOfPlaneBend OutOfPlaneBend to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeOutOfPlaneBend(OutOfPlaneBend outOfPlaneBend) {
    if (outOfPlaneBend == null) {
      return false;
    }
    return outOfPlaneBends.remove(outOfPlaneBend);
  }

  /**
   * Get the OutOfPlaneBend at a given index.
   *
   * @param index Index in the internal list.
   * @return OutOfPlaneBend at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public OutOfPlaneBend getOutOfPlaneBend(int index) {
    return outOfPlaneBends.get(index);
  }

  /**
   * Get an unmodifiable view of the OutOfPlaneBends in this term.
   *
   * @return Unmodifiable List of OutOfPlaneBends.
   */
  public List<OutOfPlaneBend> getOutOfPlaneBends() {
    return Collections.unmodifiableList(outOfPlaneBends);
  }

  /**
   * Get an array of OutOfPlaneBends in this term.
   *
   * @return Array of OutOfPlaneBends.
   */
  public OutOfPlaneBend[] getOutOfPlaneBendArray() {
    return outOfPlaneBends.toArray(new OutOfPlaneBend[0]);
  }

  /**
   * Get the number of OutOfPlaneBends in this term.
   *
   * @return The number of OutOfPlaneBends.
   */
  public int getNumberOfOutOfPlaneBends() {
    return outOfPlaneBends.size();
  }

  /**
   * Log the details of Out-of-Plane Bend interactions.
   */
  @Override
  public void log() {
    if (getNumberOfOutOfPlaneBends() <= 0) {
      return;
    }
    logger.info("\n Out-of-Plane Bend Interactions:");
    for (OutOfPlaneBend outOfPlaneBend : getOutOfPlaneBends()) {
      logger.info(" Out-of-Plane Bend \t" + outOfPlaneBend.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfOutOfPlaneBends() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "OUT-OF-PLANE BEND          : ", getEnergy(), getNumberOfOutOfPlaneBends());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Out-of-Plane Bend ",
        getEnergy(), getNumberOfOutOfPlaneBends(), getTime());
  }
}
