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

import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getScaleX;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getScaleY;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_getScaleZ;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_setFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions;

/**
 * This class uses a Monte Carlo algorithm to adjust the size of the periodic box,
 * simulating the effect of constant pressure. It assumes the simulation is running
 * at constant temperature, and the box size is adjusted to maintain constant pressure.
 * Unlike MonteCarloBarostat, this version allows independent scaling of the three
 * box dimensions, making it suitable for anisotropic systems.
 * <p>
 * This class is most useful for simulating a system at constant pressure when
 * anisotropic scaling is desired, such as for layered materials or systems with
 * preferred orientations.
 */
public class MonteCarloAnisotropicBarostat extends Force {

  /**
   * Create a MonteCarloAnisotropicBarostat.
   *
   * @param defaultPressure    The default pressure acting on each axis (in bar).
   * @param defaultTemperature The default temperature at which the system is being maintained (in Kelvin).
   * @param scaleX             Whether to scale the X dimension of the periodic box.
   * @param scaleY             Whether to scale the Y dimension of the periodic box.
   * @param scaleZ             Whether to scale the Z dimension of the periodic box.
   * @param frequency          The frequency at which Monte Carlo pressure changes should be attempted (in time steps).
   */
  public MonteCarloAnisotropicBarostat(OpenMM_Vec3 defaultPressure, double defaultTemperature,
                                       int scaleX, int scaleY, int scaleZ, int frequency) {
    pointer = OpenMM_MonteCarloAnisotropicBarostat_create(defaultPressure, defaultTemperature,
        scaleX, scaleY, scaleZ, frequency);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_MonteCarloAnisotropicBarostat_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the default pressure (in bar).
   *
   * @return The default pressure acting on each axis.
   */
  public OpenMM_Vec3 getDefaultPressure() {
    return OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure(pointer);
  }

  /**
   * Get the default temperature at which the system is being maintained (in Kelvin).
   *
   * @return The default temperature.
   */
  public double getDefaultTemperature() {
    return OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature(pointer);
  }

  /**
   * Get the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.
   *
   * @return The frequency of pressure change attempts.
   */
  public int getFrequency() {
    return OpenMM_MonteCarloAnisotropicBarostat_getFrequency(pointer);
  }

  /**
   * Get the random number seed. See setRandomNumberSeed() for details.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed(pointer);
  }

  /**
   * Get whether to scale the X dimension of the periodic box.
   *
   * @return 1 if the X dimension should be scaled, 0 otherwise.
   */
  public int getScaleX() {
    return OpenMM_MonteCarloAnisotropicBarostat_getScaleX(pointer);
  }

  /**
   * Get whether to scale the Y dimension of the periodic box.
   *
   * @return 1 if the Y dimension should be scaled, 0 otherwise.
   */
  public int getScaleY() {
    return OpenMM_MonteCarloAnisotropicBarostat_getScaleY(pointer);
  }

  /**
   * Get whether to scale the Z dimension of the periodic box.
   *
   * @return 1 if the Z dimension should be scaled, 0 otherwise.
   */
  public int getScaleZ() {
    return OpenMM_MonteCarloAnisotropicBarostat_getScaleZ(pointer);
  }

  /**
   * Set the default pressure acting on each axis (in bar).
   *
   * @param pressure The default pressure acting on each axis.
   */
  public void setDefaultPressure(OpenMM_Vec3 pressure) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure(pointer, pressure);
  }

  /**
   * Set the default temperature at which the system is being maintained (in Kelvin).
   *
   * @param temperature The default temperature.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature(pointer, temperature);
  }

  /**
   * Set the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.
   *
   * @param frequency The frequency of pressure change attempts.
   */
  public void setFrequency(int frequency) {
    OpenMM_MonteCarloAnisotropicBarostat_setFrequency(pointer, frequency);
  }

  /**
   * Set the random number seed. The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way. It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of random numbers will be different.
   * On the other hand, no guarantees are made about the behavior of simulations that use the same seed.
   * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
   * results on successive runs, even if those runs were initialized identically.
   * <p>
   * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
   * is created from this Force. This is done to ensure that each Context receives unique random seeds
   * without you needing to set them explicitly.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Returns whether this force makes use of periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}