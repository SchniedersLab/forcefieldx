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
package ffx.openmm;

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_addSubsystemThermostat;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_addThermostat;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_computeHeatBathEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getMaximumPairDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getNumThermostats;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getRelativeTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_getThermostat;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_hasSubsystemThermostats;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_setCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_setMaximumPairDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_setRelativeTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_setTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NoseHooverIntegrator_step;

/**
 * This class implements the Nosé-Hoover integrator for constant temperature molecular dynamics.
 * The Nosé-Hoover method uses extended system dynamics to maintain constant temperature by
 * coupling the system to a heat bath through additional degrees of freedom (thermostats).
 * <p>
 * The integrator can handle multiple thermostats, each controlling different parts of the system
 * or different degrees of freedom. This allows for more sophisticated temperature control
 * compared to simple velocity rescaling methods.
 */
public class NoseHooverIntegrator extends Integrator {

  /**
   * Create a NoseHooverIntegrator with a single thermostat.
   *
   * @param stepSize The step size with which to integrate the system (in ps).
   */
  public NoseHooverIntegrator(double stepSize) {
    pointer = OpenMM_NoseHooverIntegrator_create(stepSize);
  }

  /**
   * Create a NoseHooverIntegrator with detailed thermostat parameters.
   *
   * @param stepSize           The step size with which to integrate the system (in ps).
   * @param temperature        The temperature of the heat bath (in Kelvin).
   * @param collisionFrequency The collision frequency (in 1/ps).
   * @param numMTS             The number of multiple time step levels.
   * @param numYoshidaSuzuki   The number of Yoshida-Suzuki steps.
   * @param numNoseHoover      The number of Nosé-Hoover chain thermostats.
   */
  public NoseHooverIntegrator(double stepSize, double temperature, double collisionFrequency,
                              int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    pointer = OpenMM_NoseHooverIntegrator_create_2(stepSize, temperature, collisionFrequency,
        numMTS, numYoshidaSuzuki, numNoseHoover);
  }

  /**
   * Add a subsystem thermostat to the integrator.
   *
   * @param particles                  The particles to be controlled by this thermostat.
   * @param chainWeights               The weights for the thermostat chain.
   * @param temperature                The temperature for this thermostat (in Kelvin).
   * @param collisionFrequency         The collision frequency for this thermostat (in 1/ps).
   * @param relativeTemperature        The relative temperature scaling factor.
   * @param relativeCollisionFrequency The relative collision frequency scaling factor.
   * @param numMTS                     The number of multiple time step levels.
   * @param numYoshidaSuzuki           The number of Yoshida-Suzuki steps.
   * @param numNoseHoover              The number of Nosé-Hoover chain thermostats.
   * @return The index of the thermostat that was added.
   */
  public int addSubsystemThermostat(PointerByReference particles, PointerByReference chainWeights,
                                    double temperature, double collisionFrequency,
                                    double relativeTemperature, double relativeCollisionFrequency,
                                    int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    return OpenMM_NoseHooverIntegrator_addSubsystemThermostat(pointer, particles, chainWeights,
        temperature, collisionFrequency,
        relativeTemperature, relativeCollisionFrequency,
        numMTS, numYoshidaSuzuki, numNoseHoover);
  }

  /**
   * Add a thermostat to the integrator.
   *
   * @param temperature        The temperature for this thermostat (in Kelvin).
   * @param collisionFrequency The collision frequency for this thermostat (in 1/ps).
   * @param numMTS             The number of multiple time step levels.
   * @param numYoshidaSuzuki   The number of Yoshida-Suzuki steps.
   * @param numNoseHoover      The number of Nosé-Hoover chain thermostats.
   * @return The index of the thermostat that was added.
   */
  public int addThermostat(double temperature, double collisionFrequency,
                           int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    return OpenMM_NoseHooverIntegrator_addThermostat(pointer, temperature, collisionFrequency,
        numMTS, numYoshidaSuzuki, numNoseHoover);
  }

  /**
   * Compute the total energy of all heat baths.
   *
   * @return The total heat bath energy.
   */
  public double computeHeatBathEnergy() {
    return OpenMM_NoseHooverIntegrator_computeHeatBathEnergy(pointer);
  }

  /**
   * Destroy the integrator.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_NoseHooverIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the collision frequency for a thermostat.
   *
   * @param thermostat The index of the thermostat.
   * @return The collision frequency (in 1/ps).
   */
  public double getCollisionFrequency(int thermostat) {
    return OpenMM_NoseHooverIntegrator_getCollisionFrequency(pointer, thermostat);
  }

  /**
   * Get the maximum pair distance for neighbor list updates.
   *
   * @return The maximum pair distance (in nm).
   */
  public double getMaximumPairDistance() {
    return OpenMM_NoseHooverIntegrator_getMaximumPairDistance(pointer);
  }

  /**
   * Get the number of thermostats.
   *
   * @return The number of thermostats.
   */
  public int getNumThermostats() {
    return OpenMM_NoseHooverIntegrator_getNumThermostats(pointer);
  }

  /**
   * Get the relative collision frequency for a thermostat.
   *
   * @param thermostat The index of the thermostat.
   * @return The relative collision frequency scaling factor.
   */
  public double getRelativeCollisionFrequency(int thermostat) {
    return OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency(pointer, thermostat);
  }

  /**
   * Get the relative temperature for a thermostat.
   *
   * @param thermostat The index of the thermostat.
   * @return The relative temperature scaling factor.
   */
  public double getRelativeTemperature(int thermostat) {
    return OpenMM_NoseHooverIntegrator_getRelativeTemperature(pointer, thermostat);
  }

  /**
   * Get the temperature for a thermostat.
   *
   * @param thermostat The index of the thermostat.
   * @return The temperature (in Kelvin).
   */
  public double getTemperature(int thermostat) {
    return OpenMM_NoseHooverIntegrator_getTemperature(pointer, thermostat);
  }

  /**
   * Get a reference to a thermostat.
   *
   * @param thermostat The index of the thermostat.
   * @return A reference to the thermostat object.
   */
  public PointerByReference getThermostat(int thermostat) {
    return OpenMM_NoseHooverIntegrator_getThermostat(pointer, thermostat);
  }

  /**
   * Check if the integrator has subsystem thermostats.
   *
   * @return 1 if subsystem thermostats are present, 0 otherwise.
   */
  public int hasSubsystemThermostats() {
    return OpenMM_NoseHooverIntegrator_hasSubsystemThermostats(pointer);
  }

  /**
   * Set the collision frequency for a thermostat.
   *
   * @param collisionFrequency The collision frequency (in 1/ps).
   * @param thermostat         The index of the thermostat.
   */
  public void setCollisionFrequency(double collisionFrequency, int thermostat) {
    OpenMM_NoseHooverIntegrator_setCollisionFrequency(pointer, collisionFrequency, thermostat);
  }

  /**
   * Set the maximum pair distance for neighbor list updates.
   *
   * @param distance The maximum pair distance (in nm).
   */
  public void setMaximumPairDistance(double distance) {
    OpenMM_NoseHooverIntegrator_setMaximumPairDistance(pointer, distance);
  }

  /**
   * Set the relative collision frequency for a thermostat.
   *
   * @param relativeCollisionFrequency The relative collision frequency scaling factor.
   * @param thermostat                 The index of the thermostat.
   */
  public void setRelativeCollisionFrequency(double relativeCollisionFrequency, int thermostat) {
    OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency(pointer, relativeCollisionFrequency, thermostat);
  }

  /**
   * Set the relative temperature for a thermostat.
   *
   * @param relativeTemperature The relative temperature scaling factor.
   * @param thermostat          The index of the thermostat.
   */
  public void setRelativeTemperature(double relativeTemperature, int thermostat) {
    OpenMM_NoseHooverIntegrator_setRelativeTemperature(pointer, relativeTemperature, thermostat);
  }

  /**
   * Set the temperature for a thermostat.
   *
   * @param temperature The temperature (in Kelvin).
   * @param thermostat  The index of the thermostat.
   */
  public void setTemperature(double temperature, int thermostat) {
    OpenMM_NoseHooverIntegrator_setTemperature(pointer, temperature, thermostat);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps The number of time steps to take.
   */
  public void step(int steps) {
    OpenMM_NoseHooverIntegrator_step(pointer, steps);
  }
}