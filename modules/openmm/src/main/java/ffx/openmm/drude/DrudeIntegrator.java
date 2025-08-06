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
 * A base class to encapsulate features common to Drude integrators.
 */
public class DrudeIntegrator extends Integrator {

  /**
   * Create a DrudeIntegrator.
   *
   * @param pointer A pointer to the native OpenMM Drude integrator object.
   */
  public DrudeIntegrator(PointerByReference pointer) {
    super(pointer);
  }

  /**
   * Create a DrudeIntegrator.
   *
   * @param stepSize the step size with which to integrator the system (in picoseconds)
   */
  public DrudeIntegrator(double stepSize) {
    super(OpenMM_DrudeIntegrator_create(stepSize));
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
   * Get the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
   *
   * @return the temperature of the heat bath, measured in Kelvin
   */
  public double getDrudeTemperature() {
    return OpenMM_DrudeIntegrator_getDrudeTemperature(pointer);
  }

  /**
   * Get the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
   * with a hard wall constraint.  The default value is 0.02.  If this distance is set to 0, the hard wall constraint is omitted.
   */
  public double getMaxDrudeDistance() {
    return OpenMM_DrudeIntegrator_getMaxDrudeDistance(pointer);
  }

  /**
   * Get the random number seed.  See setRandomNumberSeed() for details.
   */
  public int getRandomNumberSeed() {
    return OpenMM_DrudeIntegrator_getRandomNumberSeed(pointer);
  }

  /**
   * Set the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
   *
   * @param temperature the temperature of the heat bath, measured in Kelvin
   */
  public void setDrudeTemperature(double temperature) {
    OpenMM_DrudeIntegrator_setDrudeTemperature(pointer, temperature);
  }

  /**
   * Set the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
   * with a hard wall constraint.  The default value is 0.02.  If this distance is set to 0, the hard wall constraint is omitted.
   */
  public void setMaxDrudeDistance(double distance) {
    OpenMM_DrudeIntegrator_setMaxDrudeDistance(pointer, distance);
  }

  /**
   * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of random forces will be different.  On
   * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
   * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
   * results on successive runs, even if those runs were initialized identically.
   * <p>
   * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
   * is created from this Force. This is done to ensure that each Context receives unique random seeds
   * without you needing to set them explicitly.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_DrudeIntegrator_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps the number of time steps to take
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeIntegrator_step(pointer, steps);
  }
}