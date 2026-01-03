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
package ffx.openmm.drude;

import ffx.openmm.NoseHooverIntegrator;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance;

/**
 * This Integrator simulates systems that include Drude particles.  It applies two different Nose-Hoover
 * chain thermostats to the different parts of the system.  The first is applied to ordinary particles (ones
 * that are not part of a Drude particle pair), as well as to the center of mass of each Drude particle pair.
 * A second thermostat, typically with a much lower temperature, is applied to the relative internal
 * displacement of each pair.
 * <p>
 * This integrator can optionally set an upper limit on how far any Drude particle is ever allowed to
 * get from its parent particle.  This can sometimes help to improve stability.  The limit is enforced
 * with a hard wall constraint.  By default the limit is set to 0.02 nm.
 * <p>
 * This Integrator requires the System to include a DrudeForce, which it uses to identify the Drude
 * particles.
 */
public class DrudeNoseHooverIntegrator extends NoseHooverIntegrator {

  /**
   * Create a DrudeNoseHooverIntegrator.
   *
   * @param stepSize         the step size with which to integrator the system (in picoseconds)
   * @param temperature      the target temperature for the system (in Kelvin).
   * @param drudeTemperature the target temperature for the Drude particles, relative to their parent atom (in Kelvin).
   * @param frequency        the frequency of the system's interaction with the heat bath (in inverse picoseconds).
   * @param drudeFrequency   the frequency of the drude particles' interaction with the heat bath (in inverse picoseconds).
   * @param chainLength      the number of beads in the Nose-Hoover chain.
   * @param drudeChainLength the number of beads in the Nose-Hoover chain for Drude particles.
   * @param numMTS           the number of step in the  multiple time step chain propagation algorithm.
   */
  public DrudeNoseHooverIntegrator(double stepSize, double temperature, double drudeTemperature,
                                   double frequency, double drudeFrequency, int chainLength,
                                   int drudeChainLength, int numMTS) {
    super(OpenMM_DrudeNoseHooverIntegrator_create(stepSize, temperature, drudeTemperature,
        frequency, drudeFrequency, chainLength, drudeChainLength, numMTS));
  }

  /**
   * Compute the kinetic energy of the Drude particles.
   * <p>
   * This method calculates the total kinetic energy associated with the Drude particles
   * in the system. This can be used to monitor the energy distribution between nuclear
   * and electronic degrees of freedom during the simulation.
   *
   * @return The Drude kinetic energy in kJ/mol.
   */
  public double computeDrudeKineticEnergy() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(pointer);
  }

  /**
   * Compute the instantaneous temperature of the Drude particles.
   * <p>
   * This method calculates the current kinetic temperature of the Drude particles
   * based on their velocities. This provides a measure of the thermal state
   * of the electronic degrees of freedom in the system.
   *
   * @return The instantaneous Drude temperature in Kelvin.
   */
  public double computeDrudeTemperature() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(pointer);
  }

  /**
   * Compute the instantaneous temperature of the real atoms.
   * <p>
   * This method calculates the current kinetic temperature of the real atoms
   * based on their velocities. This provides a measure of the thermal state
   * of the nuclear degrees of freedom in the system.
   *
   * @return The instantaneous system temperature in Kelvin.
   */
  public double computeSystemTemperature() {
    return OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(pointer);
  }

  /**
   * Compute the total kinetic energy of the system.
   * <p>
   * This method calculates the total kinetic energy of all particles in the system,
   * including both real atoms and Drude particles. This is useful for monitoring
   * energy conservation and the overall thermal state of the system.
   *
   * @return The total kinetic energy in kJ/mol.
   */
  public double computeTotalKineticEnergy() {
    return OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(pointer);
  }

  /**
   * Destroy the integrator.
   * <p>
   * This method releases the memory associated with the DrudeNoseHooverIntegrator object.
   * After calling this method, the integrator should not be used.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_DrudeNoseHooverIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the maximum allowed distance between Drude particles and their parent atoms.
   * <p>
   * This method returns the constraint distance that limits how far Drude particles
   * can move from their parent atoms. This constraint prevents the polarization
   * from becoming unphysically large and maintains numerical stability.
   *
   * @return The maximum Drude distance in nanometers.
   */
  public double getMaxDrudeDistance() {
    return OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(pointer);
  }

  /**
   * Set the maximum allowed distance between Drude particles and their parent atoms.
   * <p>
   * This method sets the constraint distance that limits how far Drude particles
   * can move from their parent atoms. This constraint is essential for maintaining
   * numerical stability and preventing unphysical polarization. A typical value
   * is around 0.02 nm.
   *
   * @param distance The maximum Drude distance in nanometers.
   */
  public void setMaxDrudeDistance(double distance) {
    OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(pointer, distance);
  }
}