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

import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.Integrator;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_getDrudeTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_getMaxDrudeDistance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_setDrudeTemperature;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_setMaxDrudeDistance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeIntegrator_step;

/**
 * DrudeIntegrator implements a specialized integrator for systems with Drude oscillators.
 * <p>
 * This integrator is designed to handle the dual-temperature dynamics required for
 * Drude polarizable force fields. It maintains separate temperatures for the real
 * atoms and the Drude particles, allowing for proper equilibration of both the
 * nuclear and electronic degrees of freedom.
 * <p>
 * The Drude integrator uses a dual-thermostat approach where:
 * <ul>
 * <li>Real atoms are coupled to a thermostat at the system temperature</li>
 * <li>Drude particles are coupled to a separate thermostat at a lower temperature</li>
 * <li>The relative motion between Drude particles and their parent atoms is constrained</li>
 * </ul>
 * <p>
 * This approach ensures that the Drude oscillators remain close to their equilibrium
 * positions while allowing for proper thermal fluctuations of the induced dipoles.
 */
public class DrudeIntegrator extends Integrator {

  /**
   * Create a new DrudeIntegrator.
   * <p>
   * This constructor initializes a Drude integrator with the specified step size.
   * The integrator will use default values for the Drude temperature and maximum
   * Drude distance, which can be adjusted using the appropriate setter methods.
   *
   * @param pointer A pointer to the native OpenMM Drude integrator object.
   */
  public DrudeIntegrator(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Create a new DrudeIntegrator.
   * <p>
   * This constructor initializes a Drude integrator with the specified step size.
   * The integrator will use default values for the Drude temperature and maximum
   * Drude distance, which can be adjusted using the appropriate setter methods.
   *
   * @param stepSize The integration step size in picoseconds.
   */
  public DrudeIntegrator(double stepSize) {
    pointer = OpenMM_DrudeIntegrator_create(stepSize);
  }

  /**
   * Destroy the integrator.
   * <p>
   * This method releases the memory associated with the DrudeIntegrator object.
   * After calling this method, the integrator should not be used.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_DrudeIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the temperature of the Drude particles.
   * <p>
   * This method returns the temperature at which the Drude particles are
   * thermostated. This temperature is typically much lower than the system
   * temperature to keep the Drude oscillators close to their equilibrium positions.
   *
   * @return The Drude temperature in Kelvin.
   */
  public double getDrudeTemperature() {
    return OpenMM_DrudeIntegrator_getDrudeTemperature(pointer);
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
    return OpenMM_DrudeIntegrator_getMaxDrudeDistance(pointer);
  }

  /**
   * Get the random number seed used by this integrator.
   * <p>
   * This method returns the seed value used to initialize the random number
   * generator for the stochastic components of the integration algorithm.
   * The same seed will produce reproducible trajectories.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_DrudeIntegrator_getRandomNumberSeed(pointer);
  }

  /**
   * Set the temperature of the Drude particles.
   * <p>
   * This method sets the temperature at which the Drude particles are
   * thermostated. The Drude temperature should typically be much lower
   * than the system temperature (e.g., 1 K) to maintain the Drude oscillators
   * close to their equilibrium positions while allowing thermal fluctuations.
   *
   * @param temperature The Drude temperature in Kelvin.
   */
  public void setDrudeTemperature(double temperature) {
    OpenMM_DrudeIntegrator_setDrudeTemperature(pointer, temperature);
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
    OpenMM_DrudeIntegrator_setMaxDrudeDistance(pointer, distance);
  }

  /**
   * Set the random number seed used by this integrator.
   * <p>
   * This method sets the seed value for the random number generator used
   * in the stochastic components of the integration algorithm. Setting the
   * same seed will produce reproducible trajectories, which is useful for
   * debugging and validation purposes.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_DrudeIntegrator_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Integrate the system forward in time by the specified number of time steps.
   * <p>
   * This method advances the simulation by the given number of integration steps.
   * Each step advances the system by the step size specified in the constructor.
   * The method handles both the real atoms and Drude particles according to
   * their respective thermostats and constraints.
   *
   * @param steps The number of steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeIntegrator_step(pointer, steps);
  }
}