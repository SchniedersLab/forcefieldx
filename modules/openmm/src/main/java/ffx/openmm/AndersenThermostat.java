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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_getDefaultCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_getDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;

/**
 * This class uses the Andersen method to maintain constant temperature.
 */
public class AndersenThermostat extends Force {

  /**
   * OpenMM AndersenThermostat constructor.
   *
   * @param temperature the default temperature of the heat bath (in Kelvin).
   * @param frequency   the default collision frequency (in 1/ps)
   */
  public AndersenThermostat(double temperature, double frequency) {
    pointer = OpenMM_AndersenThermostat_create(temperature, frequency);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AndersenThermostat_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the default collision frequency (in 1/ps).
   *
   * @return the default collision frequency (in 1/ps).
   */
  public double getDefaultCollisionFrequency() {
    return OpenMM_AndersenThermostat_getDefaultCollisionFrequency(pointer);
  }

  /**
   * Get the default temperature of the heat bath (in Kelvin).
   *
   * @return the default temperature of the heat bath, measured in Kelvin.
   */
  public double getDefaultTemperature() {
    return OpenMM_AndersenThermostat_getDefaultTemperature(pointer);
  }

  /**
   * Get the random number seed. See setRandomNumberSeed() for details.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_AndersenThermostat_getRandomNumberSeed(pointer);
  }

  /**
   * Set the default collision frequency. This will affect any new Contexts you create,
   * but not ones that already exist.
   *
   * @param frequency the default collision frequency (in 1/ps).
   */
  public void setDefaultCollisionFrequency(double frequency) {
    OpenMM_AndersenThermostat_setDefaultCollisionFrequency(pointer, frequency);
  }

  /**
   * Set the default temperature of the heat bath. This will affect any new Contexts
   * you create, but not ones that already exist.
   *
   * @param temperature the default temperature of the heat bath (in Kelvin).
   */
  public void setDefaultTemperature(double temperature) {
    OpenMM_AndersenThermostat_setDefaultTemperature(pointer, temperature);
  }

  /**
   * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of collisions will be different.  On
   * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
   * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
   * results on successive runs, even if those runs were initialized identically.
   * <p>
   * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
   * is created from this Force. This is done to ensure that each Context receives unique random seeds
   * without you needing to set them explicitly.
   *
   * @param seed the random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_AndersenThermostat_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Returns whether this force makes use of periodic boundary conditions.
   *
   * @return the Andersen Thermostat always returns false.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions(pointer);
    return pbc != OpenMM_False;
  }
}