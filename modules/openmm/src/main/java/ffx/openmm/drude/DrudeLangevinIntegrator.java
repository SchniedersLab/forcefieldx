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
 * This Integrator simulates systems that include Drude particles.  It applies two different Langevin
 * thermostats to different parts of the system.  The first is applied to ordinary particles (ones that
 * are not part of a Drude particle pair), as well as to the center of mass of each Drude particle pair.
 * A second thermostat, typically with a much lower temperature, is applied to the relative internal
 * displacement of each pair.
 * <p>
 * This integrator can optionally set an upper limit on how far any Drude particle is ever allowed to
 * get from its parent particle.  This can sometimes help to improve stability.  The limit is enforced
 * with a hard wall constraint.  By default the limit is set to 0.02 nm.
 * <p>
 * This Integrator requires the System to include a DrudeForce, which it uses to identify the Drude
 * particles.
 */
public class DrudeLangevinIntegrator extends DrudeIntegrator {

  /**
   * Create a DrudeLangevinIntegrator.
   *
   * @param stepSize         the step size with which to integrator the system (in picoseconds)
   * @param temperature      the temperature of the main heat bath (in Kelvin)
   * @param friction         the friction coefficient which couples the system to the main heat bath (in inverse picoseconds)
   * @param drudeTemperature the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin)
   * @param drudeFriction    the friction coefficient which couples the system to the heat bath applied to internal coordinates of Drude particles (in inverse picoseconds)
   */
  public DrudeLangevinIntegrator(double stepSize, double temperature, double friction,
                                 double drudeTemperature, double drudeFriction) {
    super(OpenMM_DrudeLangevinIntegrator_create(stepSize, temperature, friction,
        drudeTemperature, drudeFriction));
  }

  /**
   * Compute the instantaneous temperature of the Drude system, measured in Kelvin.
   * This is calculated based on the kinetic energy of the internal motion of Drude pairs
   * and should remain close to the prescribed Drude temperature.
   */
  public double computeDrudeTemperature() {
    return OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(pointer);
  }

  /**
   * Compute the instantaneous temperature of the System, measured in Kelvin.
   * This is calculated based on the kinetic energy of the ordinary particles (ones
   * not attached to a Drude particle), as well as the center of mass motion of the
   * Drude particle pairs.  It does not include the internal motion of the pairs.
   * On average, this should be approximately equal to the value returned by
   * getTemperature().
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
   * Get the friction coefficient which determines how strongly the internal coordinates of Drude particles
   * are coupled to the heat bath (in inverse ps).
   *
   * @return the friction coefficient, measured in 1/ps
   */
  public double getDrudeFriction() {
    return OpenMM_DrudeLangevinIntegrator_getDrudeFriction(pointer);
  }

  /**
   * Get the friction coefficient which determines how strongly the system is coupled to
   * the main heat bath (in inverse ps).
   *
   * @return the friction coefficient, measured in 1/ps
   */
  public double getFriction() {
    return OpenMM_DrudeLangevinIntegrator_getFriction(pointer);
  }

  /**
   * Get the temperature of the main heat bath (in Kelvin).
   *
   * @return the temperature of the heat bath, measured in Kelvin
   */
  public double getTemperature() {
    return OpenMM_DrudeLangevinIntegrator_getTemperature(pointer);
  }

  /**
   * Set the friction coefficient which determines how strongly the internal coordinates of Drude particles
   * are coupled to the heat bath (in inverse ps).
   *
   * @param friction the friction coefficient, measured in 1/ps
   */
  public void setDrudeFriction(double friction) {
    OpenMM_DrudeLangevinIntegrator_setDrudeFriction(pointer, friction);
  }

  /**
   * Set the friction coefficient which determines how strongly the system is coupled to
   * the main heat bath (in inverse ps).
   *
   * @param friction the friction coefficient, measured in 1/ps
   */
  public void setFriction(double friction) {
    OpenMM_DrudeLangevinIntegrator_setFriction(pointer, friction);
  }

  /**
   * Set the temperature of the main heat bath (in Kelvin).
   *
   * @param temperature the temperature of the heat bath, measured in Kelvin
   */
  public void setTemperature(double temperature) {
    OpenMM_DrudeLangevinIntegrator_setTemperature(pointer, temperature);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps the number of time steps to take
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeLangevinIntegrator_step(pointer, steps);
  }
}