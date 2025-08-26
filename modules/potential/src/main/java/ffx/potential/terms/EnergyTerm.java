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

import edu.rit.pj.reduction.SharedDouble;
import ffx.potential.bonded.BondedTerm;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Base class for potential energy terms used by ForceFieldEnergy.
 * <p>
 * This abstract class holds minimal, common metadata shared by energy terms:
 * a human-readable name and an integer force-group identifier. Specific energy
 * terms (bonded, nonbonded, restraints, etc.) should extend this class.
 */
public abstract class EnergyTerm {

  private static final Logger logger = Logger.getLogger(EnergyTerm.class.getName());

  /**
   * Optional name for this energy term (for reporting/logging).
   */
  protected String name;

  /**
   * Integer force group this term belongs to. Defaults to 0.
   */
  protected int forceGroup = 0;

  /**
   * Total potential energy for this term.
   */
  protected SharedDouble energy = new SharedDouble();

  /**
   * Shared double for RMSD calculations.
   */
  private final SharedDouble sharedRMSD = new SharedDouble();

  /**
   * The time taken to compute the energy for this term in nanoseconds.
   */
  private long computeTime = 0L;

  /**
   * Create an unnamed EnergyTerm with default force group 0.
   */
  protected EnergyTerm() {
  }

  /**
   * Create a named EnergyTerm with default force group 0.
   *
   * @param name Name for this energy term.
   */
  protected EnergyTerm(String name) {
    this.name = name;
  }

  /**
   * Create a named EnergyTerm with a specified force group.
   *
   * @param name       Name for this energy term.
   * @param forceGroup Integer force group identifier.
   */
  protected EnergyTerm(String name, int forceGroup) {
    this.name = name;
    this.forceGroup = forceGroup;
  }

  /**
   * Get the number of BondedTerms in this term.
   *
   * @return The number of BondedTerms.
   */
  abstract public int getNumberOfTerms();

  /**
   * Get an array of BondedTerms in this term.
   *
   * @return Array of BondedTerms.
   */
  abstract public BondedTerm[] getBondedTermsArray();

  /**
   * Log the details of this energy term.
   */
  abstract public void log();

  /**
   * Get a string representation of this energy term.
   *
   * @return A string representation of the energy term.
   */
  abstract public String toString();

  /**
   * Get a PDB-style REMARK representation of this energy term.
   * @return A PDB REMARK string for this energy term.
   */
  abstract public String toPDBString();

  /**
   * Get the current force group identifier.
   *
   * @return force group as an int
   */
  public int getForceGroup() {
    return forceGroup;
  }

  /**
   * Set the force group identifier.
   *
   * @param forceGroup integer force group
   */
  public void setForceGroup(int forceGroup) {
    this.forceGroup = forceGroup;
  }

  /**
   * Get the name of this energy term.
   *
   * @return name (may be null)
   */
  public String getName() {
    return name;
  }

  /**
   * Set the name of this energy term.
   *
   * @param name A human-readable name
   */
  public void setName(String name) {
    this.name = name;
  }

  /**
   * Get the total potential energy of this term.
   *
   * @return Energy value as a double.
   */
  public double getEnergy() {
    return energy.get();
  }

  /**
   * Set the total potential energy of this term.
   *
   * @param e Energy value to set.
   */
  protected void setEnergy(double e) {
    energy.set(e);
  }

  /**
   * Add energy to the total energy for this term and return the new total.
   *
   * @param e Energy value to add to the total.
   * @return The new total energy after adding the provided energy.
   */
  protected double addAndGetEnergy(double e) {
    return energy.addAndGet(e);
  }

  /**
   * Set the RMSD accumulator value for this term.
   *
   * @param rmsd RMSD value to set for this term.
   */
  protected void setRMSD(double rmsd) {
    sharedRMSD.set(rmsd);
  }

  /**
   * Get the RMSD for this term.
   */
  public double getRMSD() {
    return sqrt(sharedRMSD.get() / getNumberOfTerms());
  }

  /**
   * Add RMSD to the shared RMSD and return the new total.
   *
   * @param rmsd RMSD value to add to the shared RMSD.
   * @return The new total RMSD after adding the provided value.
   */
  protected double addAndGetRMSD(double rmsd) {
    return sharedRMSD.addAndGet(rmsd);
  }

  /**
   * Set the starting time for the energy computation.
   */
  protected void startTime() {
    computeTime = -System.nanoTime();
  }

  /**
   * Stop the timer for the energy computation.
   */
  public void stopTime() {
    computeTime += System.nanoTime();
  }

  /**
   * Get the time taken to compute the energy for this term in seconds.
   *
   * @return The time in seconds.
   */
  public double getTime() {
    return computeTime * 1.0e-9;
  }


}
