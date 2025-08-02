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
package ffx.openmm.drude;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance;

/**
 * DrudeNoseHooverIntegrator implements a Nosé-Hoover integrator for systems with Drude oscillators.
 * <p>
 * This integrator extends the basic DrudeIntegrator by implementing Nosé-Hoover dynamics, which
 * provides deterministic temperature control through extended system variables. The Nosé-Hoover
 * approach maintains canonical ensemble sampling while preserving the time-reversible nature
 * of the dynamics, making it ideal for equilibrium simulations of Drude polarizable systems.
 * <p>
 * The integrator uses separate thermostat chains for real atoms and Drude particles:
 * <ul>
 * <li>Real atoms are coupled to a Nosé-Hoover chain at the system temperature</li>
 * <li>Drude particles are coupled to a separate chain at the Drude temperature</li>
 * <li>Multiple chain lengths can be specified for improved temperature control</li>
 * </ul>
 * <p>
 * The Nosé-Hoover method provides excellent temperature control without the stochastic
 * noise inherent in Langevin dynamics, making it particularly suitable for studies
 * requiring precise temperature control and long-time dynamical properties.
 */
public class DrudeNoseHooverIntegrator extends DrudeIntegrator {

  /**
   * Create a new DrudeNoseHooverIntegrator.
   * <p>
   * This constructor initializes a Drude Nosé-Hoover integrator with the specified parameters.
   * The integrator will use Nosé-Hoover dynamics to maintain temperature control for both
   * real atoms and Drude particles using separate thermostat chains.
   *
   * @param stepSize         The integration step size in picoseconds.
   * @param temperature      The target temperature for real atoms in Kelvin.
   * @param drudeTemperature The target temperature for Drude particles in Kelvin.
   * @param frequency        The thermostat frequency for real atoms in 1/picoseconds.
   * @param drudeFrequency   The thermostat frequency for Drude particles in 1/picoseconds.
   * @param chainLength      The length of the Nosé-Hoover chain for real atoms.
   * @param drudeChainLength The length of the Nosé-Hoover chain for Drude particles.
   * @param numMTS           The number of multiple time step levels.
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
  @Override
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
  @Override
  public void setMaxDrudeDistance(double distance) {
    OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(pointer, distance);
  }
}