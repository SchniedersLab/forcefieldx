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
package ffx.algorithms.dynamics.thermostats;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.SystemState;
import ffx.numerics.Constraint;
import ffx.numerics.Potential.VARIABLE_TYPE;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Thermostat a molecular dynamics trajectory to an external bath using the Bussi, Donadio, and
 * Parrinello method. This method is similar to Berendsen thermostat, but generates a canonical
 * distribution.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://dx.doi.org/10.1016/j.cpc.2008.01.006">G. Bussi and M. Parrinello,
 *     "Stochastic Thermostats: Comparison of Local and Global Schemes", Computer Physics
 *     Communications, 179, 26-29 (2008)</a>
 * @since 1.0
 */
public class Bussi extends Thermostat {

  /** The random number generator used to perturb velocities. */
  private final Random bussiRandom;
  /** Bussi thermostat time constant (psec). */
  private double tau;

  /**
   * Constructor for Bussi.
   *
   * @param state The MDState to operate on.
   * @param type the VARIABLE_TYPE of each variable.
   * @param targetTemperature The target temperature.
   * @param tau Bussi thermostat time constant (psec).
   */
  public Bussi(SystemState state, VARIABLE_TYPE[] type, double targetTemperature, double tau) {
    this(state, type, targetTemperature, tau, Collections.emptyList());
  }

  public Bussi(SystemState state, VARIABLE_TYPE[] type, double targetTemperature, double tau,
      List<Constraint> constraints) {
    super(state, type, targetTemperature, constraints);
    this.name = ThermostatEnum.BUSSI;
    this.tau = tau;
    this.bussiRandom = new Random();
  }

  /**
   * Constructor for Bussi.
   *
   * @param state The MDState to operate on.
   * @param type the VARIABLE_TYPE of each variable.
   * @param targetTemperature a double.
   */
  public Bussi(SystemState state, VARIABLE_TYPE[] type, double targetTemperature) {
    this(state, type, targetTemperature, 0.2e0);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Full step velocity modification.
   */
  @Override
  public void fullStep(double dt) {
    double expTau = exp(-dt / tau);
    double tempRatio = targetTemperature / state.getTemperature();
    double rate = (1.0 - expTau) * tempRatio / degreesOfFreedom;
    double r = bussiRandom.nextGaussian();
    double s = 0.0;
    for (int i = 0; i < degreesOfFreedom - 1; i++) {
      double si = bussiRandom.nextGaussian();
      s += si * si;
    }
    double scale = expTau + (s + r * r) * rate + 2.0 * r * sqrt(expTau * rate);
    scale = sqrt(scale);
    if (r + sqrt(expTau / rate) < 0.0) {
      scale = -scale;
    }
    double[] v = state.v();
    double[] mass = state.getMass();
    for (int i = 0; i < state.getNumberOfVariables(); i++) {
      if (mass[i] > 0.0) {
        v[i] *= scale;
      }
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
   * <p>No velocity modifications are made by the Bussi method at the half-step.
   */
  @Override
  public void halfStep(double dt) {
  }

  /**
   * {@inheritDoc}
   *
   * <p>Initialize the Random number generator used to apply random forces to the particles.
   */
  public void setRandomSeed(long seed) {
    bussiRandom.setSeed(seed);
  }

  /**
   * Add Thermostat details to the kinetic energy and temperature details.
   *
   * @return Description of the thermostat, kinetic energy and temperature.
   */
  public String toThermostatString() {
    return format("\n Bussi Thermostat (tau = %8.3f psec)\n%s", tau, super.toString());
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return "Bussi";
  }
}
