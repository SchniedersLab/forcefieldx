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
import ffx.potential.bonded.PiOrbitalTorsion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 * Pi-Orbital Torsion potential energy term using {@link ffx.potential.bonded.PiOrbitalTorsion} instances.
 */
public class PiOrbitalTorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(PiOrbitalTorsionPotentialEnergy.class.getName());


  /**
   * Internal list of PiOrbitalTorsion instances.
   */
  private final List<PiOrbitalTorsion> piOrbitalTorsions = new ArrayList<>();

  /**
   * Create a PiOrbitalTorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public PiOrbitalTorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a PiOrbitalTorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public PiOrbitalTorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a PiOrbitalTorsionPotentialEnergy initialized with a list of terms and force group.
   *
   * @param name              Name for this term.
   * @param forceGroup        Force group identifier.
   * @param piOrbitalTorsions List of PiOrbitalTorsion instances to add (null-safe).
   */
  public PiOrbitalTorsionPotentialEnergy(String name, int forceGroup, List<PiOrbitalTorsion> piOrbitalTorsions) {
    super(name, forceGroup);
    if (piOrbitalTorsions != null) {
      Collections.sort(piOrbitalTorsions);
      this.piOrbitalTorsions.addAll(piOrbitalTorsions);
      logger.info(String.format("  Pi-Orbital Torsions:               %10d", getNumberOfPiOrbitalTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfPiOrbitalTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getPiOrbitalTorsionArray();
  }

  /**
   * Create a PiOrbitalTorsionPotentialEnergy initialized with a collection of terms.
   *
   * @param name              Name for this term (may be null).
   * @param piOrbitalTorsions Collection of PiOrbitalTorsion instances to add (null-safe).
   */
  public PiOrbitalTorsionPotentialEnergy(String name, Collection<PiOrbitalTorsion> piOrbitalTorsions) {
    super(name);
    if (piOrbitalTorsions != null) {
      this.piOrbitalTorsions.addAll(piOrbitalTorsions);
    }
  }

  /**
   * Add a PiOrbitalTorsion to this term.
   *
   * @param piOrbitalTorsion PiOrbitalTorsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addPiOrbitalTorsion(PiOrbitalTorsion piOrbitalTorsion) {
    if (piOrbitalTorsion == null) {
      return false;
    }
    return piOrbitalTorsions.add(piOrbitalTorsion);
  }

  /**
   * Add an array of PiOrbitalTorsions to this term.
   *
   * @param piOrbitalTorsions Array of PiOrbitalTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addPiOrbitalTorsions(PiOrbitalTorsion[] piOrbitalTorsions) {
    if (piOrbitalTorsions == null) {
      return false;
    }
    Collections.addAll(this.piOrbitalTorsions, piOrbitalTorsions);
    return true;
  }

  /**
   * Add a list of PiOrbitalTorsions to this term.
   *
   * @param piOrbitalTorsions List of PiOrbitalTorsion instances to add.
   * @return true if they were added.
   */
  public boolean addPiOrbitalTorsions(List<PiOrbitalTorsion> piOrbitalTorsions) {
    if (piOrbitalTorsions == null) {
      return false;
    }
    this.piOrbitalTorsions.addAll(piOrbitalTorsions);
    return true;
  }

  /**
   * Remove a PiOrbitalTorsion from this term.
   *
   * @param piOrbitalTorsion PiOrbitalTorsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removePiOrbitalTorsion(PiOrbitalTorsion piOrbitalTorsion) {
    if (piOrbitalTorsion == null) {
      return false;
    }
    return piOrbitalTorsions.remove(piOrbitalTorsion);
  }

  /**
   * Get the PiOrbitalTorsion at a given index.
   *
   * @param index Index in the internal list.
   * @return PiOrbitalTorsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public PiOrbitalTorsion getPiOrbitalTorsion(int index) {
    return piOrbitalTorsions.get(index);
  }

  /**
   * Get an unmodifiable view of the PiOrbitalTorsions in this term.
   *
   * @return Unmodifiable List of PiOrbitalTorsions.
   */
  public List<PiOrbitalTorsion> getPiOrbitalTorsions() {
    return Collections.unmodifiableList(piOrbitalTorsions);
  }

  /**
   * Get an array of PiOrbitalTorsions in this term.
   *
   * @return Array of PiOrbitalTorsions.
   */
  public PiOrbitalTorsion[] getPiOrbitalTorsionArray() {
    return piOrbitalTorsions.toArray(new PiOrbitalTorsion[0]);
  }

  /**
   * Get the number of PiOrbitalTorsions in this term.
   *
   * @return The number of PiOrbitalTorsions.
   */
  public int getNumberOfPiOrbitalTorsions() {
    return piOrbitalTorsions.size();
  }

  /**
   * Log the details of Pi-Orbital Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfPiOrbitalTorsions() <= 0) {
      return;
    }
    logger.info("\n Pi-Orbital Torsion Interactions:");
    for (PiOrbitalTorsion piOrbitalTorsion : getPiOrbitalTorsions()) {
      logger.info(" Pi-Torsion \t" + piOrbitalTorsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfPiOrbitalTorsions() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "PI-ORBITAL TORSION         : ", getEnergy(), getNumberOfPiOrbitalTorsions());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Pi-Orbital Torsion",
        getEnergy(), getNumberOfPiOrbitalTorsions(), getTime());
  }
}
