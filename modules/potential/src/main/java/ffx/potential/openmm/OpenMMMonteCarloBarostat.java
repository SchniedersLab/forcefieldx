// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setDefaultPressure;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setRandomNumberSeed;

/**
 * OpenMM MonteCarloBarostat.
 */
public class OpenMMMonteCarloBarostat extends OpenMMForce {

  /**
   * OpenMM MonteCarloBarostat constructor.
   *
   * @param pressure    The pressure.
   * @param temperature The temperature.
   * @param frequency   The frequency to apply the barostat.
   */
  public OpenMMMonteCarloBarostat(double pressure, double temperature, int frequency) {
    forcePointer = OpenMM_MonteCarloBarostat_create(pressure, temperature, frequency);
  }

  /**
   * Set the random number seed.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_MonteCarloBarostat_setRandomNumberSeed(forcePointer, seed);
  }

  /**
   * Set the frequency.
   *
   * @param frequency The frequency.
   */
  public void setFrequency(int frequency) {
    OpenMM_MonteCarloBarostat_setFrequency(forcePointer, frequency);
  }

  /**
   * Set the default temperature.
   *
   * @param temperature The temperature.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMM_MonteCarloBarostat_setDefaultTemperature(forcePointer, temperature);
  }

  /**
   * Set the default pressure.
   *
   * @param pressure The pressure.
   */
  public void setDefaultPressure(double pressure) {
    OpenMM_MonteCarloBarostat_setDefaultPressure(forcePointer, pressure);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (forcePointer != null) {
      OpenMM_MonteCarloBarostat_destroy(forcePointer);
      forcePointer = null;
    }
  }

}
