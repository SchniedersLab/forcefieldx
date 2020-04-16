// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.dynamics.thermostats;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Constraint;
import ffx.numerics.Potential.VARIABLE_TYPE;
import java.util.Collections;
import java.util.List;

/**
 * Thermostat a molecular dynamics trajectory to an external bath using the Berendsen weak-coupling
 * thermostat.
 *
 * @author Michael J. Schnieders derived from TINKER temperature control by Alan Grossfield and Jay
 *     Ponder
 * @see <a href="http://link.aip.org/link/?JCP/81/3684">H. J. C. Berendsen, J. P. M. Postma, W. F.
 *     van Gunsteren, A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling to an External
 *     Bath", Journal of Chemical Physics, 81, 3684-3690 (1984)</a>
 */
public class Berendsen extends Thermostat {

  /** Berendsen time constant (psec). */
  private double tau;

  /**
   * Constructor for Berendsen.
   *
   * @param n Number of degrees of freedom.
   * @param x Atomic coordinates.
   * @param v Velocities.
   * @param mass Mass of each degrees of freedom.
   * @param type The VARIABLE_TYPE of each variable.
   * @param targetTemperature The target temperatures.
   * @param tau Berendsen thermostat time constant (psec).
   */
  public Berendsen(
      int n,
      double[] x,
      double[] v,
      double[] mass,
      VARIABLE_TYPE[] type,
      double targetTemperature,
      double tau) {
    this(n, x, v, mass, type, targetTemperature, tau, Collections.emptyList());
  }

  public Berendsen(
      int n,
      double[] x,
      double[] v,
      double[] mass,
      VARIABLE_TYPE[] type,
      double targetTemperature,
      double tau,
      List<Constraint> constraints) {
    super(n, x, v, mass, type, targetTemperature, constraints);
    this.name = ThermostatEnum.BERENDSEN;
    this.tau = tau;
  }

  /**
   * Constructor for Berendsen.
   *
   * @param n Number of degrees of freedom.
   * @param x Atomic coordinates.
   * @param v Velocities.
   * @param mass Mass of each degrees of freedom.
   * @param type The VARIABLE_TYPE of each variable.
   * @param targetTemperature The target temperatures.
   */
  public Berendsen(
      int n,
      double[] x,
      double[] v,
      double[] mass,
      VARIABLE_TYPE[] type,
      double targetTemperature) {
    this(n, x, v, mass, type, targetTemperature, 0.2e0);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Full step velocity modification.
   */
  @Override
  public void fullStep(double dt) {
    double ratio = targetTemperature / currentTemperature;
    double scale = sqrt(1.0 + (dt / tau) * (ratio - 1.0));
    for (int i = 0; i < nVariables; i++) {
      v[i] *= scale;
    }
  }

  /**
   * Getter for the field <code>tau</code>.
   *
   * @return a double.
   */
  public double getTau() {
    return tau;
  }

  /**
   * Setter for the field <code>tau</code>.
   *
   * @param tau a double.
   */
  public void setTau(double tau) {
    this.tau = tau;
  }

  /**
   * {@inheritDoc}
   *
   * <p>No velocity modifications are made by the Berendsen method at the half-step.
   */
  @Override
  public void halfStep(double dt) {}

  /**
   * Add Thermostat details to the kinetic energy and temperature details.
   *
   * @return Description of the thermostat, kinetic energy and temperature.
   */
  public String toThermostatString() {
    return format("\n Berendsen Thermostat (tau = %8.3f psec)\n  %s", tau, super.toString());
  }
}
