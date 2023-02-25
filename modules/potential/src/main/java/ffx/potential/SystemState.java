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
package ffx.potential;

import java.util.Arrays;

/**
 * The current state of the molecular dynamics simulation.
 */
public class SystemState {

  /**
   * Number of dynamics variables. The length of x, v, a, aPrevious, gradient, and mass.
   */
  protected final int numberOfVariables;
  /** Coordinates. */
  protected final double[] x;
  /** Velocities. */
  protected final double[] v;
  /** Accelerations. */
  protected final double[] a;
  /** Previous accelerations. */
  protected final double[] aPrevious;
  /** The gradient. */
  protected final double[] gradient;
  /** Mass for each degree of freedom. */
  protected final double[] mass;
  /** Current temperature. */
  double temperature;
  /** Current kinetic energy. */
  double kineticEnergy;
  /** Current potential energy. */
  double potentialEnergy;

  /**
   * Constructor for MDState.
   *
   * @param numberOfVariables The number of variables.
   */
  public SystemState(int numberOfVariables) {
    this.numberOfVariables = numberOfVariables;
    x = new double[numberOfVariables];
    v = new double[numberOfVariables];
    a = new double[numberOfVariables];
    aPrevious = new double[numberOfVariables];
    gradient = new double[numberOfVariables];
    mass = new double[numberOfVariables];
  }

  /**
   * Get an unmodifiable view of the current state.
   */
  public UnmodifiableState getUnmodifiableState() {
    return new UnmodifiableState(x, v, a, aPrevious, mass, gradient, kineticEnergy,
        potentialEnergy, temperature);
  }

  /**
   * Get the number of variables.
   *
   * @return The number of variables.
   */
  public int getNumberOfVariables() {
    return numberOfVariables;
  }

  /**
   * Revert the current state to the passed UnmodifiableMDState.
   *
   * @param state The state to revert to.
   */
  public void revertState(UnmodifiableState state) {
    assert (state.x().length == numberOfVariables);
    System.arraycopy(state.x(), 0, x, 0, numberOfVariables);
    System.arraycopy(state.v(), 0, v, 0, numberOfVariables);
    System.arraycopy(state.a(), 0, a, 0, numberOfVariables);
    System.arraycopy(state.aPrevious(), 0, aPrevious, 0, numberOfVariables);
    System.arraycopy(state.mass(), 0, mass, 0, numberOfVariables);
    System.arraycopy(state.gradient(), 0, gradient, 0, numberOfVariables);
    kineticEnergy = state.kineticEnergy();
    potentialEnergy = state.potentialEnergy();
    temperature = state.temperature();
  }

  /**
   * Set the mass of each degree of freedom.
   *
   * @param mass The mass of each degree of freedom.
   */
  public void setMass(double[] mass) {
    assert (mass.length == numberOfVariables);
    System.arraycopy(mass, 0, this.mass, 0, numberOfVariables);
  }

  /**
   * Set the coordinates via a copy of the passed array into the internal array.
   *
   * @param x The coordinates.
   */
  public void setCoordinates(double[] x) {
    assert (x.length == numberOfVariables);
    System.arraycopy(x, 0, this.x, 0, numberOfVariables);
  }

  /**
   * Set the velocities via a copy of the passed array into the internal array.
   *
   * @param v The velocities.
   */
  public void setVelocities(double[] v) {
    assert (v.length == numberOfVariables);
    System.arraycopy(v, 0, this.v, 0, numberOfVariables);
  }

  /**
   * Set the accelerations via a copy of the passed array into the internal array.
   *
   * @param a The accelerations.
   */
  public void setAccelerations(double[] a) {
    assert (a.length == numberOfVariables);
    System.arraycopy(a, 0, this.a, 0, numberOfVariables);
  }

  /**
   * Set the previous accelerations via a copy of the passed array into the internal array.
   *
   * @param aPrevious The previous accelerations.
   */
  public void setPreviousAccelerations(double[] aPrevious) {
    assert (aPrevious.length == numberOfVariables);
    System.arraycopy(aPrevious, 0, this.aPrevious, 0, numberOfVariables);
  }

  /**
   * Get a reference to the internal coordinates array.
   *
   * @return The coordinates.
   */
  public double[] x() {
    return x;
  }

  /**
   * Get a reference to the internal velocities array.
   *
   * @return The velocities.
   */
  public double[] v() {
    return v;
  }

  /**
   * Get a reference to the internal accelerations array.
   *
   * @return The accelerations.
   */
  public double[] a() {
    return a;
  }

  /**
   * Get a reference to the internal previous accelerations array.
   *
   * @return The previous accelerations.
   */
  public double[] aPrevious() {
    return aPrevious;
  }

  /**
   * Get a reference to the internal mass array.
   *
   * @return The mass.
   */
  public double[] mass() {
    return mass;
  }

  /**
   * Get a reference to the internal gradient array.
   *
   * @return The gradient.
   */
  public double[] gradient() {
    return gradient;
  }

  /**
   * Get a copy of the internal coordinate array.
   *
   * @return The coordinates.
   */
  public double[] getCoordinatesCopy() {
    return Arrays.copyOf(x, numberOfVariables);
  }

  /**
   * Copy the current accelerations to the previous accelerations.
   */
  public void copyAccelerationsToPrevious() {
    System.arraycopy(a, 0, aPrevious, 0, numberOfVariables);
  }

  /**
   * Set the temperature.
   *
   * @param temperature The temperature.
   */
  public void setTemperature(double temperature) {
    this.temperature = temperature;
  }

  /**
   * Set the kinetic energy.
   *
   * @param kineticEnergy The kinetic energy.
   */
  public void setKineticEnergy(double kineticEnergy) {
    this.kineticEnergy = kineticEnergy;
  }

  /**
   * Set the potential energy.
   *
   * @param potentialEnergy The potential energy.
   */
  public void setPotentialEnergy(double potentialEnergy) {
    this.potentialEnergy = potentialEnergy;
  }

  /**
   * Get the temperature.
   */
  public double getTemperature() {
    return temperature;
  }

  /**
   * Get the kinetic energy.
   */
  public double getKineticEnergy() {
    return kineticEnergy;
  }

  /**
   * Get the potential energy.
   */
  public double getPotentialEnergy() {
    return potentialEnergy;
  }

  /**
   * Get the total energy as the sum of the kinetic and potential energies.
   */
  public double getTotalEnergy() {
    return kineticEnergy + potentialEnergy;
  }

}
