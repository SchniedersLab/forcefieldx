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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getXYMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_getZMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setXYMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_setZMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions;

/**
 * This class uses a Monte Carlo algorithm to adjust the size of the periodic box,
 * simulating the effect of constant pressure and surface tension on a membrane system.
 * It assumes the simulation is running at constant temperature, and the box size is
 * adjusted to maintain constant pressure in the XY plane and constant surface tension
 * in the Z direction.
 * <p>
 * This class is most useful for simulating membrane systems at constant pressure and
 * surface tension, where the membrane is oriented in the XY plane and the normal
 * direction is along the Z axis.
 */
public class MonteCarloMembraneBarostat extends Force {

  /**
   * Create a MonteCarloMembraneBarostat.
   *
   * @param defaultPressure       The default pressure acting on the system (in bar).
   * @param defaultSurfaceTension The default surface tension acting on the membrane (in bar*nm).
   * @param defaultTemperature    The default temperature at which the system is being maintained (in Kelvin).
   * @param xymode                The mode for scaling the XY dimensions (0 = isotropic, 1 = anisotropic).
   * @param zmode                 The mode for scaling the Z dimension (0 = constant volume, 1 = constant pressure).
   * @param frequency             The frequency at which Monte Carlo pressure changes should be attempted (in time steps).
   */
  public MonteCarloMembraneBarostat(double defaultPressure, double defaultSurfaceTension,
                                    double defaultTemperature, int xymode, int zmode, int frequency) {
    pointer = OpenMM_MonteCarloMembraneBarostat_create(defaultPressure, defaultSurfaceTension,
        defaultTemperature, xymode, zmode, frequency);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_MonteCarloMembraneBarostat_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the default pressure (in bar).
   *
   * @return The default pressure acting on the system.
   */
  public double getDefaultPressure() {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultPressure(pointer);
  }

  /**
   * Get the default surface tension (in bar*nm).
   *
   * @return The default surface tension acting on the membrane.
   */
  public double getDefaultSurfaceTension() {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension(pointer);
  }

  /**
   * Get the default temperature at which the system is being maintained (in Kelvin).
   *
   * @return The default temperature.
   */
  public double getDefaultTemperature() {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature(pointer);
  }

  /**
   * Get the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.
   *
   * @return The frequency of pressure change attempts.
   */
  public int getFrequency() {
    return OpenMM_MonteCarloMembraneBarostat_getFrequency(pointer);
  }

  /**
   * Get the random number seed. See setRandomNumberSeed() for details.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed(pointer);
  }

  /**
   * Get the mode for scaling the XY dimensions.
   *
   * @return The XY scaling mode (0 = isotropic, 1 = anisotropic).
   */
  public int getXYMode() {
    return OpenMM_MonteCarloMembraneBarostat_getXYMode(pointer);
  }

  /**
   * Get the mode for scaling the Z dimension.
   *
   * @return The Z scaling mode (0 = constant volume, 1 = constant pressure).
   */
  public int getZMode() {
    return OpenMM_MonteCarloMembraneBarostat_getZMode(pointer);
  }

  /**
   * Set the default pressure acting on the system (in bar).
   *
   * @param pressure The default pressure acting on the system.
   */
  public void setDefaultPressure(double pressure) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultPressure(pointer, pressure);
  }

  /**
   * Set the default surface tension acting on the membrane (in bar*nm).
   *
   * @param surfaceTension The default surface tension acting on the membrane.
   */
  public void setDefaultSurfaceTension(double surfaceTension) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension(pointer, surfaceTension);
  }

  /**
   * Set the default temperature at which the system is being maintained (in Kelvin).
   *
   * @param temperature The default temperature.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature(pointer, temperature);
  }

  /**
   * Set the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.
   *
   * @param frequency The frequency of pressure change attempts.
   */
  public void setFrequency(int frequency) {
    OpenMM_MonteCarloMembraneBarostat_setFrequency(pointer, frequency);
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
    OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Set the mode for scaling the XY dimensions.
   *
   * @param xymode The XY scaling mode (0 = isotropic, 1 = anisotropic).
   */
  public void setXYMode(int xymode) {
    OpenMM_MonteCarloMembraneBarostat_setXYMode(pointer, xymode);
  }

  /**
   * Set the mode for scaling the Z dimension.
   *
   * @param zmode The Z scaling mode (0 = constant volume, 1 = constant pressure).
   */
  public void setZMode(int zmode) {
    OpenMM_MonteCarloMembraneBarostat_setZMode(pointer, zmode);
  }

  /**
   * Returns whether this force makes use of periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}