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
import static java.lang.String.format;

/**
 * Torsion potential energy term using {@link ffx.potential.bonded.Torsion} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TorsionPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(TorsionPotentialEnergy.class.getName());


  /**
   * Internal list of Torsion instances.
   */
  private final List<Torsion> torsions = new ArrayList<>();

  /**
   * Create a TorsionPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public TorsionPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a TorsionPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public TorsionPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a TorsionPotentialEnergy initialized with a list of torsions and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Force group identifier.
   * @param torsions   List of Torsion instances to add (null-safe).
   */
  public TorsionPotentialEnergy(String name, int forceGroup, List<Torsion> torsions) {
    super(name, forceGroup);
    if (torsions != null) {
      Collections.sort(torsions);
      this.torsions.addAll(torsions);
      logger.info(format("  Torsions:                          %10d", getNumberOfTorsions()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfTorsions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getTorsionArray();
  }

  /**
   * Create a TorsionPotentialEnergy initialized with a collection of torsions.
   *
   * @param name     Name for this term (may be null).
   * @param torsions Collection of Torsion instances to add (null-safe).
   */
  public TorsionPotentialEnergy(String name, Collection<Torsion> torsions) {
    super(name);
    if (torsions != null) {
      this.torsions.addAll(torsions);
    }
  }

  /**
   * Add a Torsion to this term.
   *
   * @param torsion Torsion to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addTorsion(Torsion torsion) {
    if (torsion == null) {
      return false;
    }
    return torsions.add(torsion);
  }

  /**
   * Add an array of Torsions to this term.
   *
   * @param torsions Array of Torsion instances to add.
   * @return true if they were added.
   */
  public boolean addTorsions(Torsion[] torsions) {
    if (torsions == null) {
      return false;
    }
    Collections.addAll(this.torsions, torsions);
    return true;
  }

  /**
   * Add a list of Torsions to this term.
   *
   * @param torsions List of Torsion instances to add.
   * @return true if they were added.
   */
  public boolean addTorsions(List<Torsion> torsions) {
    if (torsions == null) {
      return false;
    }
    this.torsions.addAll(torsions);
    return true;
  }

  /**
   * Remove a Torsion from this term.
   *
   * @param torsion Torsion to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeTorsion(Torsion torsion) {
    if (torsion == null) {
      return false;
    }
    return torsions.remove(torsion);
  }

  /**
   * Get the Torsion at a given index.
   *
   * @param index Index in the internal list.
   * @return Torsion at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public Torsion getTorsion(int index) {
    return torsions.get(index);
  }

  /**
   * Get an unmodifiable view of the Torsions in this term.
   *
   * @return Unmodifiable List of Torsions.
   */
  public List<Torsion> getTorsions() {
    return Collections.unmodifiableList(torsions);
  }

  /**
   * Get an array of Torsions in this term.
   *
   * @return Array of Torsions.
   */
  public Torsion[] getTorsionArray() {
    return torsions.toArray(new Torsion[0]);
  }

  /**
   * Get the number of Torsions in this term.
   *
   * @return The number of Torsions.
   */
  public int getNumberOfTorsions() {
    return torsions.size();
  }

  /**
   * Set the lambda value for all Torsions in this term.
   * @param lambda Lambda value to set for all Torsions.
   */
  public void setLambda(double lambda) {
    for (Torsion torsion : torsions) {
      torsion.setLambda(lambda);
    }
  }

  /**
   * Get the energy contribution from all Torsions in this term.
   * @return Total energy from all Torsions.
   */
  public double getdEdL() {
    double dEdL = 0.0;
    for (Torsion torsion : getTorsions()) {
      dEdL += torsion.getdEdL();
    }
    return dEdL;
  }

  /**
   * Get the energy contribution from all Torsions in this term.
   * @return Total energy from all Torsions.
   */
  public double getd2EdL2() {
    double d2EdLambda2 = 0.0;
    for (Torsion torsion : getTorsions()) {
      d2EdLambda2 += torsion.getd2EdL2();
    }
    return d2EdLambda2;
  }

  /**
   * Log the details of Torsion interactions.
   */
  @Override
  public void log() {
    if (getNumberOfTorsions() <= 0) {
      return;
    }
    logger.info("\n Torsion Angle Interactions:");
    for (Torsion torsion : getTorsions()) {
      logger.info(" Torsion \t" + torsion.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfTorsions() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "TORSIONAL ANGLE            : ", getEnergy(), getNumberOfTorsions());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Torsional Angle   ",
        getEnergy(), getNumberOfTorsions(), getTime());
  }
}
