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
package ffx.openmm.ffm;

import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Simulates a {@link System} using one or more Nose-Hoover chain thermostats with the "middle" leapfrog
 * propagation algorithm (J. Phys. Chem. A 2019, 123, 6056-6079).
 *
 * <p>Thermostats are added with {@link #addThermostat(double, double, int, int, int)} or {@link
 * #addSubsystemThermostat(IntArray, BondArray, double, double, double, double, int, int, int)}; temperatures, collision
 * frequencies and the maximum pair distance can then be adjusted. The OpenMM header does not state whether such
 * changes reach an existing {@link Context}.</p>
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>This class is not final so that {@code ffx.openmm.ffm.drude.DrudeNoseHooverIntegrator} can extend it, as the JNA
 * Drude class extends the JNA {@link ffx.openmm.NoseHooverIntegrator}.</p>
 *
 * <p>FFM/JNA differences: the JNA class uses {@code PointerByReference} for native handles and containers, while this
 * class uses {@link MemorySegment}, {@link IntArray} and {@link BondArray}. The JNA six-argument constructor passes its
 * arguments positionally to {@code OpenMM_NoseHooverIntegrator_create_2}, whose native order is (temperature,
 * collisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki); this class reorders them. The JNA single-argument
 * documentation describes a single thermostat, but the OpenMM header states that the step-size constructor creates a bare
 * integrator with no thermostats.</p>
 */
public class NoseHooverIntegrator extends Integrator {

  /**
   * Create a bare leapfrog integrator with no thermostats; add them with {@link #addThermostat(double, double,
   * int, int, int)}. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The OpenMM header documents this constructor as creating no thermostats. The JNA and earlier FFM text that
   * described "a single thermostat" disagrees with the header.</p>
   *
   * @param stepSize step size with which to integrate the system, in ps.
   */
  public NoseHooverIntegrator(double stepSize) {
    super(create(stepSize));
  }

  /**
   * Wrap an existing native Nose-Hoover integrator. No native object is created or runtime initialized here.
   *
   * <p>Ownership of the handle is transferred to this wrapper, which destroys it in {@link #destroy()}; this constructor
   * exists mainly for subclasses (such as the Drude Nose-Hoover wrapper) that create the native object themselves. The JNA
   * counterpart accepts a {@code PointerByReference}.</p>
   *
   * @param pointer non-null native integrator handle whose ownership is transferred to this wrapper.
   * @throws IllegalArgumentException if the handle address is zero.
   */
  public NoseHooverIntegrator(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Create an integrator with explicit thermostat parameters. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * <p>The Java parameter order (stepSize first) is kept from the JNA constructor but values are passed to {@code
   * OpenMM_NoseHooverIntegrator_create_2} in its native order (temperature, collisionFrequency, stepSize, chain length,
   * numMTS, numYoshidaSuzuki). The JNA constructor forwards its arguments positionally without reordering, so it does not
   * match the native order.</p>
   *
   * <p>Parameter-name mapping: the Java parameter {@code numNoseHoover} is OpenMM's {@code chainLength} (the number of
   * beads in each Nose-Hoover chain), not a number of thermostats. {@code numMTS} is the number of steps in the
   * multiple-time-step chain propagation, and {@code numYoshidaSuzuki} is the number of terms in the Yoshida-Suzuki
   * decomposition (the header requires 1, 3, 5 or 7). The values are passed to the native function in the order
   * chain length, numMTS, numYoshidaSuzuki.</p>
   *
   * @param stepSize           step size with which to integrate the system, in ps.
   * @param temperature        target temperature of the system, in K.
   * @param collisionFrequency frequency of interaction with the heat bath, in 1/ps.
   * @param numMTS             number of steps in the multiple-time-step chain propagation.
   * @param numYoshidaSuzuki   number of terms in the Yoshida-Suzuki decomposition (1, 3, 5 or 7 per the header).
   * @param numNoseHoover      number of beads in the Nose-Hoover chain (OpenMM's chain length).
   */
  public NoseHooverIntegrator(double stepSize, double temperature, double collisionFrequency,
                              int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    super(create(stepSize, temperature, collisionFrequency, numMTS, numYoshidaSuzuki,
        numNoseHoover));
  }

  /**
   * Add a Nose-Hoover chain thermostat controlling a subset of particles and/or connected particle pairs.
   *
   * <p>Only the listed particles are thermostated. For listed pairs both the absolute center-of-mass motion and the
   * relative motion are thermostated independently. If both containers are empty, all particles are thermostated. The
   * containers are only read and remain owned by the caller; the JNA counterpart passes {@code PointerByReference}
   * arguments (named {@code particles} and {@code chainWeights}) instead.</p>
   *
   * <p>Parameter-name mapping: the Java parameter {@code numNoseHoover} is OpenMM's {@code chainLength} (the number of
   * beads in each Nose-Hoover chain), not a number of thermostats. {@code numMTS} is the number of steps in the
   * multiple-time-step chain propagation, and {@code numYoshidaSuzuki} is the number of terms in the Yoshida-Suzuki
   * decomposition (the header requires 1, 3, 5 or 7). The values are passed to the native function in the order
   * chain length, numMTS, numYoshidaSuzuki.</p>
   *
   * @param particles                  {@link IntArray} of particle indices to thermostat; must not be null.
   * @param particlePairs              {@link BondArray} of connected particle pairs whose center-of-mass and relative
   *     motion are thermostated; must not be null.
   * @param temperature                target temperature for each pair's absolute center-of-mass motion, in K.
   * @param collisionFrequency         frequency of interaction with the heat bath for the pairs' center-of-mass motion,
   *     in 1/ps.
   * @param relativeTemperature        target temperature for each pair's relative motion, in K.
   * @param relativeCollisionFrequency frequency of interaction with the heat bath for the pairs' relative motion, in 1/ps.
   * @param numMTS                     number of steps in the multiple-time-step chain propagation.
   * @param numYoshidaSuzuki           number of terms in the Yoshida-Suzuki decomposition (1, 3, 5 or 7 per the header).
   * @param numNoseHoover              number of beads in each Nose-Hoover chain (OpenMM's chain length).
   * @return index of the thermostat that was added.
   * @throws NullPointerException if a container is null.
   */
  public int addSubsystemThermostat(IntArray particles, BondArray particlePairs,
                                    double temperature, double collisionFrequency, double relativeTemperature,
                                    double relativeCollisionFrequency, int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_addSubsystemThermostat(
        getPointer(), particles.getPointer(), particlePairs.getPointer(), temperature,
        collisionFrequency, relativeTemperature, relativeCollisionFrequency, numNoseHoover,
        numMTS, numYoshidaSuzuki);
  }

  /**
   * Add a simple Nose-Hoover chain thermostat controlling the temperature of the full system.
   *
   * <p>Parameter-name mapping: the Java parameter {@code numNoseHoover} is OpenMM's {@code chainLength} (the number of
   * beads in each Nose-Hoover chain), not a number of thermostats. {@code numMTS} is the number of steps in the
   * multiple-time-step chain propagation, and {@code numYoshidaSuzuki} is the number of terms in the Yoshida-Suzuki
   * decomposition (the header requires 1, 3, 5 or 7). The values are passed to the native function in the order
   * chain length, numMTS, numYoshidaSuzuki.</p>
   *
   * @param temperature        target temperature of the system, in K.
   * @param collisionFrequency frequency of interaction with the heat bath, in 1/ps.
   * @param numMTS             number of steps in the multiple-time-step chain propagation.
   * @param numYoshidaSuzuki   number of terms in the Yoshida-Suzuki decomposition (1, 3, 5 or 7 per the header).
   * @param numNoseHoover      number of beads in the Nose-Hoover chain (OpenMM's chain length).
   * @return index of the thermostat that was added.
   */
  public int addThermostat(double temperature, double collisionFrequency, int numMTS,
                           int numYoshidaSuzuki, int numNoseHoover) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_addThermostat(
        getPointer(), temperature, collisionFrequency, numNoseHoover, numMTS, numYoshidaSuzuki);
  }

  /**
   * Compute the total (potential plus kinetic) heat-bath energy of all heat baths of this integrator at the
   * current time.
   *
   * <p>Call only after the integrator is bound to a live {@link Context}; the header does not say what happens
   * otherwise, and a related Drude integrator query made without a context terminated the JVM during testing of
   * this port. The header states no unit; OpenMM energies are in kJ/mol.</p>
   *
   * @return total heat-bath energy.
   */
  public double computeHeatBathEnergy() {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_computeHeatBathEnergy(getPointer());
  }

  /**
   * Destroy the native integrator.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. If the
   * integrator was passed to a {@link Context}, that context's destruction also destroys the integrator, so do not call
   * this for an integrator that a live context owns.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_NoseHooverIntegrator_destroy);
  }

  /**
   * Get the collision frequency for the absolute motion of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param thermostat index of the thermostat.
   * @return collision frequency, in 1/ps.
   */
  public double getCollisionFrequency(int thermostat) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getCollisionFrequency(
        getPointer(), thermostat);
  }

  /**
   * Get the maximum distance a connected pair may stray from each other. Zero means there is no constraint on the
   * intra-pair separation.
   *
   * @return maximum pair distance, in nm.
   */
  public double getMaximumPairDistance() {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getMaximumPairDistance(getPointer());
  }

  /**
   * Get the number of Nose-Hoover chains registered with this integrator.
   *
   * @return number of thermostats.
   */
  public int getNumThermostats() {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getNumThermostats(getPointer());
  }

  /**
   * Get the collision frequency for the relative motion of the pairs of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param thermostat index of the thermostat.
   * @return relative collision frequency, in 1/ps.
   */
  public double getRelativeCollisionFrequency(int thermostat) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency(
        getPointer(), thermostat);
  }

  /**
   * Get the temperature for the relative motion of the pairs of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param thermostat index of the thermostat.
   * @return relative temperature, in K.
   */
  public double getRelativeTemperature(int thermostat) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getRelativeTemperature(
        getPointer(), thermostat);
  }

  /**
   * Get the temperature controlling the absolute particle motion of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param thermostat index of the thermostat.
   * @return temperature, in K.
   */
  public double getTemperature(int thermostat) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getTemperature(getPointer(), thermostat);
  }

  /**
   * Get the native Nose-Hoover chain object of a thermostat.
   *
   * <p>The returned segment is a borrowed, non-owning handle to a native {@code OpenMM_NoseHooverChain} that is valid
   * only while this integrator is alive; it must not be destroyed by the caller. It is not wrapped in a Java class. The
   * JNA counterpart returns a {@code PointerByReference}.</p>
   *
   * @param thermostat index of the thermostat.
   * @return borrowed native chain handle.
   */
  public MemorySegment getThermostat(int thermostat) {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_getThermostat(getPointer(), thermostat);
  }

  /**
   * Report whether subsystem thermostats are present.
   *
   * <p>The OpenMM header documents the native result as false if the integrator was set up with the constructor that
   * thermostats the whole system, and true otherwise. This method returns the native integer unchanged (nonzero for
   * true); it is not converted to a boolean, as in the JNA class.</p>
   *
   * @return nonzero if subsystem thermostats are present, zero otherwise.
   */
  public int hasSubsystemThermostats() {
    return OpenMMNative.OpenMM_NoseHooverIntegrator_hasSubsystemThermostats(getPointer());
  }

  /**
   * Set the collision frequency for the absolute motion of a thermostat. The header describes the value as "in
   * picosecond", which is a header typo for 1/ps; this class documents 1/ps, consistent with the getter.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param frequency  collision frequency, in 1/ps.
   * @param thermostat index of the thermostat.
   */
  public void setCollisionFrequency(double frequency, int thermostat) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_setCollisionFrequency(
        getPointer(), frequency, thermostat);
  }

  /**
   * Set the maximum distance a connected pair may stray from each other, implemented as a hard wall. Zero omits
   * the constraint so pairs may separate by any distance.
   *
   * @param distance maximum pair distance, in nm.
   */
  public void setMaximumPairDistance(double distance) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_setMaximumPairDistance(getPointer(), distance);
  }

  /**
   * Set the collision frequency for the relative motion of the pairs of a thermostat. The header describes the
   * value as "in picosecond", a typo for 1/ps.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param frequency  relative collision frequency, in 1/ps.
   * @param thermostat index of the thermostat.
   */
  public void setRelativeCollisionFrequency(double frequency, int thermostat) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency(
        getPointer(), frequency, thermostat);
  }

  /**
   * Set the temperature for the relative motion of the pairs of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param temperature relative temperature, in K.
   * @param thermostat  index of the thermostat.
   */
  public void setRelativeTemperature(double temperature, int thermostat) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_setRelativeTemperature(
        getPointer(), temperature, thermostat);
  }

  /**
   * Set the temperature controlling the absolute particle motion of a thermostat.
   *
   * <p>The index of a thermostat is the value OpenMM names {@code chainID}; thermostats are numbered from 0 in the
   * order they were added.</p>
   *
   * @param temperature temperature, in K.
   * @param thermostat  index of the thermostat.
   */
  public void setTemperature(double temperature, int thermostat) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_setTemperature(
        getPointer(), temperature, thermostat);
  }

  /**
   * Advance the simulation by a series of fixed time steps. This overrides {@link Integrator#step(int)} with the
   * Nose-Hoover native step call.
   *
   * @param steps number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_NoseHooverIntegrator_step(getPointer(), steps);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param stepSize step size, in ps.
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create(double stepSize) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_NoseHooverIntegrator_create(stepSize);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param stepSize           step size, in ps.
   * @param temperature        target temperature, in K.
   * @param collisionFrequency collision frequency, in 1/ps.
   * @param numMTS             multiple-time-step count.
   * @param numYoshidaSuzuki   Yoshida-Suzuki term count.
   * @param numNoseHoover      chain length.
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create(double stepSize, double temperature,
                                      double collisionFrequency, int numMTS, int numYoshidaSuzuki, int numNoseHoover) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_NoseHooverIntegrator_create_2(temperature, collisionFrequency,
        stepSize, numNoseHoover, numMTS, numYoshidaSuzuki);
  }
}
