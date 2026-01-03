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
package ffx.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_computeCurrentPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_getDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_getDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_getFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions;

/**
 * This class uses a Monte Carlo algorithm to adjust the size of the periodic box, simulating the
 * effect of constant pressure.
 * <p>
 * This class assumes the simulation is also being run at constant temperature, and requires you
 * to specify the system temperature (since it affects the acceptance probability for Monte Carlo
 * moves).  It does not actually perform temperature regulation, however.  You must use another
 * mechanism along with it to maintain the temperature, such as LangevinIntegrator or AndersenThermostat.
 */
public class MonteCarloBarostat extends Force {

  /**
   * Create a MonteCarloBarostat.
   *
   * @param pressure    the default pressure acting on the system (in bar)
   * @param temperature the default temperature at which the system is being maintained (in Kelvin)
   * @param frequency   the frequency at which Monte Carlo pressure changes should be attempted (in time steps)
   */
  public MonteCarloBarostat(double pressure, double temperature, int frequency) {
    super(OpenMM_MonteCarloBarostat_create(pressure, temperature, frequency));
  }

  /**
   * Compute the instantaneous pressure of a system to which this barostat is applied.
   * <p>
   * The pressure is computed from the molecular virial, using a finite difference to
   * calculate the derivative of potential energy with respect to volume.  For most systems
   * in equilibrium, the time average of the instantaneous pressure should equal the
   * pressure applied by the barostat.  Fluctuations around the average value can be
   * extremely large, however, and it may take a very long simulation to accurately
   * compute the average.
   *
   * @param context the Context for which to compute the current pressure
   * @return the instantaneous pressure
   */
  public double computeCurrentPressure(Context context) {
    return OpenMM_MonteCarloBarostat_computeCurrentPressure(pointer, context.getPointer());
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_MonteCarloBarostat_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the default pressure acting on the system (in bar).
   *
   * @return the default pressure acting on the system, measured in bar.
   */
  public double getDefaultPressure() {
    return OpenMM_MonteCarloBarostat_getDefaultPressure(pointer);
  }

  /**
   * Get the default temperature at which the system is being maintained, measured in Kelvin.
   *
   * @return the default temperature at which the system is being maintained, measured in Kelvin.
   */
  public double getDefaultTemperature() {
    return OpenMM_MonteCarloBarostat_getDefaultTemperature(pointer);
  }

  /**
   * Get the frequency.
   *
   * @return The frequency.
   */
  public int getFrequency() {
    return OpenMM_MonteCarloBarostat_getFrequency(pointer);
  }

  /**
   * Get the random number seed.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_MonteCarloBarostat_getRandomNumberSeed(pointer);
  }

  /**
   * Set the default pressure acting on the system.  This will affect any new Contexts you create,
   * but not ones that already exist.
   *
   * @param pressure the default pressure acting on the system, measured in bar.
   */
  public void setDefaultPressure(double pressure) {
    OpenMM_MonteCarloBarostat_setDefaultPressure(pointer, pressure);
  }

  /**
   * Set the default temperature at which the system is being maintained.  This will affect any new Contexts you create,
   * but not ones that already exist.
   *
   * @param temperature the system temperature, measured in Kelvin.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMM_MonteCarloBarostat_setDefaultTemperature(pointer, temperature);
  }

  /**
   * Set the frequency.
   *
   * @param frequency The frequency.
   */
  public void setFrequency(int frequency) {
    OpenMM_MonteCarloBarostat_setFrequency(pointer, frequency);
  }

  /**
   * Set the random number seed.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_MonteCarloBarostat_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Does the force use periodic boundary conditions?
   *
   * @return True if the force uses periodic boundary conditions.
   */
  public boolean usesPeriodicBoundaryConditions() {
    int periodic = OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions(pointer);
    return periodic == 1;
  }

}