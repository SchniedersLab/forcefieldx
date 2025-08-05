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

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_computeSystemTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_getDrudeFriction;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_getFriction;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_getTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_setDrudeFriction;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_setFriction;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_setTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeLangevinIntegrator_step;

/**
 * DrudeLangevinIntegrator implements a Langevin integrator for systems with Drude oscillators.
 * <p>
 * This integrator extends the basic DrudeIntegrator by adding Langevin dynamics, which provides
 * stochastic temperature control through friction and random forces. The Langevin approach is
 * particularly effective for maintaining proper thermal equilibrium in Drude polarizable systems.
 * <p>
 * The integrator uses separate friction coefficients for real atoms and Drude particles:
 * <ul>
 * <li>Real atoms experience friction at the system temperature</li>
 * <li>Drude particles experience separate friction at the Drude temperature</li>
 * <li>Both friction coefficients can be independently controlled</li>
 * </ul>
 * <p>
 * The Langevin dynamics help maintain proper temperature distribution while allowing
 * the Drude oscillators to respond appropriately to local electric fields. This approach
 * provides better temperature control compared to simple velocity rescaling methods.
 */
public class DrudeLangevinIntegrator extends DrudeIntegrator {

  /**
   * Create a new DrudeLangevinIntegrator.
   * <p>
   * This constructor initializes a Drude Langevin integrator with the specified parameters.
   * The integrator will use Langevin dynamics to maintain temperature control for both
   * real atoms and Drude particles with independent friction coefficients.
   *
   * @param stepSize         The integration step size in picoseconds.
   * @param temperature      The target temperature for real atoms in Kelvin.
   * @param friction         The friction coefficient for real atoms in 1/picoseconds.
   * @param drudeTemperature The target temperature for Drude particles in Kelvin.
   * @param drudeFriction    The friction coefficient for Drude particles in 1/picoseconds.
   */
  public DrudeLangevinIntegrator(double stepSize, double temperature, double friction,
                                 double drudeTemperature, double drudeFriction) {
    super(OpenMM_DrudeLangevinIntegrator_create(stepSize, temperature, friction,
        drudeTemperature, drudeFriction));
  }

  /**
   * Compute the instantaneous temperature of the Drude particles.
   * <p>
   * This method calculates the current kinetic temperature of the Drude particles
   * based on their velocities. This can be used to monitor the thermal equilibration
   * of the electronic degrees of freedom during the simulation.
   *
   * @return The instantaneous Drude temperature in Kelvin.
   */
  public double computeDrudeTemperature() {
    return OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(pointer);
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
    return OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(pointer);
  }

  /**
   * Destroy the integrator.
   * <p>
   * This method releases the memory associated with the DrudeLangevinIntegrator object.
   * After calling this method, the integrator should not be used.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_DrudeLangevinIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the friction coefficient for Drude particles.
   * <p>
   * This method returns the friction coefficient applied to Drude particles
   * in the Langevin dynamics. Higher friction values lead to stronger coupling
   * to the heat bath and faster temperature equilibration.
   *
   * @return The Drude friction coefficient in 1/picoseconds.
   */
  public double getDrudeFriction() {
    return OpenMM_DrudeLangevinIntegrator_getDrudeFriction(pointer);
  }

  /**
   * Get the friction coefficient for real atoms.
   * <p>
   * This method returns the friction coefficient applied to real atoms
   * in the Langevin dynamics. This parameter controls the strength of
   * coupling between the system and the thermal reservoir.
   *
   * @return The friction coefficient in 1/picoseconds.
   */
  public double getFriction() {
    return OpenMM_DrudeLangevinIntegrator_getFriction(pointer);
  }

  /**
   * Get the target temperature for real atoms.
   * <p>
   * This method returns the temperature at which the real atoms are
   * thermostated using Langevin dynamics. This is typically the desired
   * simulation temperature for the nuclear degrees of freedom.
   *
   * @return The target temperature in Kelvin.
   */
  public double getTemperature() {
    return OpenMM_DrudeLangevinIntegrator_getTemperature(pointer);
  }

  /**
   * Set the friction coefficient for Drude particles.
   * <p>
   * This method sets the friction coefficient applied to Drude particles
   * in the Langevin dynamics. Typical values range from 10-100 1/ps,
   * with higher values providing stronger temperature control.
   *
   * @param friction The Drude friction coefficient in 1/picoseconds.
   */
  public void setDrudeFriction(double friction) {
    OpenMM_DrudeLangevinIntegrator_setDrudeFriction(pointer, friction);
  }

  /**
   * Set the friction coefficient for real atoms.
   * <p>
   * This method sets the friction coefficient applied to real atoms
   * in the Langevin dynamics. Typical values range from 1-10 1/ps,
   * balancing temperature control with dynamic behavior preservation.
   *
   * @param friction The friction coefficient in 1/picoseconds.
   */
  public void setFriction(double friction) {
    OpenMM_DrudeLangevinIntegrator_setFriction(pointer, friction);
  }

  /**
   * Set the target temperature for real atoms.
   * <p>
   * This method sets the temperature at which the real atoms are
   * thermostated using Langevin dynamics. This should be set to the
   * desired simulation temperature for the system.
   *
   * @param temperature The target temperature in Kelvin.
   */
  public void setTemperature(double temperature) {
    OpenMM_DrudeLangevinIntegrator_setTemperature(pointer, temperature);
  }

  /**
   * Integrate the system forward in time by the specified number of time steps.
   * <p>
   * This method advances the simulation using Langevin dynamics for both
   * real atoms and Drude particles. Each step applies the appropriate
   * friction and random forces to maintain the target temperatures.
   *
   * @param steps The number of steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeLangevinIntegrator_step(pointer, steps);
  }
}