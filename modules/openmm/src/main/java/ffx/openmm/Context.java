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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_applyConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_applyVelocityConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_computeVirtualSites;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getStepCount;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getSystem;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getTime;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_reinitialize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setStepCount;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setTime;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocitiesToTemperature;
import static ffx.openmm.Vec3Array.toVec3Array;

/**
 * A Context stores the complete state of a simulation.  More specifically, it includes:
 *
 * <ul>
 * <li>The current time</li>
 * <li>The position of each particle</li>
 * <li>The velocity of each particle</li>
 * <li>The values of configurable parameters defined by Force objects in the System</li>
 * </ul>
 * <p>
 * You can retrieve a snapshot of the current state at any time by calling getState().  This
 * allows you to record the state of the simulation at various points, either for analysis
 * or for checkpointing.  getState() can also be used to retrieve the current forces on each
 * particle and the current energy of the System.
 */
public class Context {

  /**
   * Context pointer.
   */
  private PointerByReference pointer;

  /**
   * The integrator used for this context.
   */
  protected Integrator integrator;

  /**
   * The platform used for this context.
   */
  protected Platform platform;

  /**
   * Constructor.
   */
  public Context() {
    pointer = null;
    integrator = null;
    platform = null;
  }

  /**
   * Construct a new Context in which to run a simulation, explicitly specifying what Platform should be used
   * to perform calculations.
   *
   * @param system     the System which will be simulated
   * @param integrator the Integrator which will be used to simulate the System
   * @param platform   the Platform to use for calculations
   */
  public Context(System system, Integrator integrator, Platform platform) {
    pointer = OpenMM_Context_create_2(system.getPointer(), integrator.getPointer(), platform.getPointer());
    this.integrator = integrator;
    this.platform = platform;
  }

  /**
   * Update the positions of particles so that all distance constraints are satisfied.  This also recomputes
   * the locations of all virtual sites.
   *
   * @param tol the distance tolerance within which constraints must be satisfied.
   */
  public void applyConstraints(double tol) {
    OpenMM_Context_applyConstraints(pointer, tol);
  }

  /**
   * Update the velocities of particles so the net velocity of each constrained distance is zero.
   *
   * @param tol the velocity tolerance within which constraints must be satisfied.
   */
  public void applyVelocityConstraints(double tol) {
    OpenMM_Context_applyVelocityConstraints(pointer, tol);
  }

  /**
   * Recompute the locations of all virtual sites.  There is rarely a reason to call
   * this, since virtual sites are also updated by applyConstraints().  This is only
   * for the rare situations when you want to enforce virtual sites but <i>not</i>
   * constraints.
   */
  public void computeVirtualSites() {
    OpenMM_Context_computeVirtualSites(pointer);
  }

  /**
   * Destroy the context.
   */
  public void destroy() {
    if (integrator != null) {
      integrator.destroy();
      integrator = null;
    }
    if (pointer != null) {
      OpenMM_Context_destroy(pointer);
      pointer = null;
      // The platform is handled by the Context destroy method.
      platform = null;
    }
  }

  /**
   * Get the value of an adjustable parameter defined by a Force object in the System.
   *
   * @param name the name of the parameter to get
   * @return the value of the parameter
   */
  public double getParameter(String name) {
    return OpenMM_Context_getParameter(pointer, name);
  }

  /**
   * Get the value of an adjustable parameter defined by a Force object in the System.
   *
   * @param name the name of the parameter to get
   * @return the value of the parameter
   */
  public double getParameter(Pointer name) {
    return OpenMM_Context_getParameter(pointer, name);
  }

  /**
   * Get all adjustable parameters that have been defined by Force objects in the System, along
   * with their current values.
   *
   * @return a PointerByReference to a map containing the values of all parameters
   */
  public PointerByReference getParameters() {
    return OpenMM_Context_getParameters(pointer);
  }

  /**
   * Get the Platform being used for calculations.
   *
   * @return the Platform being used for calculations
   */
  public Platform getPlatform() {
    return platform;
  }

  /**
   * Get Integrator being used by this context.
   *
   * @return the Integrator being used by this context
   */
  public Integrator getIntegrator() {
    return integrator;
  }

  /**
   * Get the pointer to the context.
   *
   * @return The pointer to the context.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get a State object recording the current state information stored in this context.
   *
   * @param types              the set of data types which should be stored in the State object.  This
   *                           should be a union of DataType values, e.g. (State::Positions | State::Velocities).
   * @param enforcePeriodicBox if false, the position of each particle will be whatever position
   *                           is stored in the Context, regardless of periodic boundary conditions.  If true, particle
   *                           positions will be translated so the center of every molecule lies in the same periodic box.
   * @return a State object recording the current state information
   */
  public State getState(int types, int enforcePeriodicBox) {
    return new State(OpenMM_Context_getState(pointer, types, enforcePeriodicBox));
  }

  /**
   * Get a State object recording the current state information stored in this context.
   *
   * @param types              the set of data types which should be stored in the State object.  This
   *                           should be a union of DataType values, e.g. (State::Positions | State::Velocities).
   * @param enforcePeriodicBox if false, the position of each particle will be whatever position
   *                           is stored in the Context, regardless of periodic boundary conditions.  If true, particle
   *                           positions will be translated so the center of every molecule lies in the same periodic box.
   * @param groups             a set of bit flags for which force groups to include when computing forces
   *                           and energies.  Group i will be included if (groups&amp;(1&lt;&lt;i)) != 0.  The default value includes all groups.
   * @return a State object recording the current state information
   */
  public State getState(int types, int enforcePeriodicBox, int groups) {
    return new State(OpenMM_Context_getState_2(pointer, types, enforcePeriodicBox, groups));
  }

  /**
   * Get the current step count.
   *
   * @return the current step count
   */
  public long getStepCount() {
    return OpenMM_Context_getStepCount(pointer);
  }

  /**
   * Get System being simulated in this context.
   *
   * @return the System being simulated in this context
   */
  public System getSystem() {
    PointerByReference systemPointer = OpenMM_Context_getSystem(pointer);
    return new System(systemPointer);
  }

  /**
   * Get the current time of the simulation (in picoseconds).
   *
   * @return the current time of the simulation (in picoseconds)
   */
  public double getTime() {
    return OpenMM_Context_getTime(pointer);
  }

  /**
   * Does the context have a pointer?
   *
   * @return True if the context pointer is not null.
   */
  public boolean hasContextPointer() {
    return pointer != null;
  }

  /**
   * When a Context is created, it caches information about the System being simulated
   * and the Force objects contained in it.  This means that, if the System or Forces are then
   * modified, the Context does not see the changes.  Call reinitialize() to force
   * the Context to rebuild its internal representation of the System and pick up any changes
   * that have been made.
   * <p>
   * This is an expensive operation, so you should try to avoid calling it too frequently.
   * Most Force classes have an updateParametersInContext() method that provides a less expensive
   * way of updating certain types of information.  However, this method is the only way to
   * make some types of changes, so it is sometimes necessary to call it.
   * <p>
   * By default, reinitializing a Context causes all state information (positions, velocities,
   * etc.) to be discarded.  You can optionally tell it to try to preserve state information.
   * It does this by internally creating a checkpoint, then reinitializing the Context, then
   * loading the checkpoint.  Be aware that if the System has changed in a way that prevents
   * the checkpoint from being loaded (such as changing the number of particles), this will
   * throw an exception and the state information will be lost.
   *
   * @param preserveState if true, try to preserve state information; if false, discard all state information
   */
  public void reinitialize(int preserveState) {
    if (pointer != null) {
      OpenMM_Context_reinitialize(pointer, preserveState);
    }
  }

  /**
   * Set the value of an adjustable parameter defined by a Force object in the System.
   *
   * @param name  the name of the parameter to set
   * @param value the value of the parameter
   */
  public void setParameter(String name, double value) {
    OpenMM_Context_setParameter(pointer, name, value);
  }

  /**
   * Set the value of an adjustable parameter defined by a Force object in the System.
   *
   * @param name  the name of the parameter to set
   * @param value the value of the parameter
   */
  public void setParameter(Pointer name, double value) {
    OpenMM_Context_setParameter(pointer, name, value);
  }

  /**
   * Set the vectors defining the axes of the periodic box (measured in nm).  They will affect
   * any Force that uses periodic boundary conditions.
   * <p>
   * Triclinic boxes are supported, but the vectors must satisfy certain requirements.  In particular,
   * a must point in the x direction, b must point "mostly" in the y direction, and c must point "mostly"
   * in the z direction.  See the documentation for details.
   *
   * @param a the vector defining the first edge of the periodic box
   * @param b the vector defining the second edge of the periodic box
   * @param c the vector defining the third edge of the periodic box
   */
  public void setPeriodicBoxVectors(OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c) {
    OpenMM_Context_setPeriodicBoxVectors(pointer, a, b, c);
  }

  /**
   * Set the positions of all particles in the System (measured in nm).  This method simply sets the positions
   * without checking to see whether they satisfy distance constraints.  If you want constraints to be
   * enforced, call applyConstraints() after setting the positions.
   *
   * @param positions a vector whose length equals the number of particles in the System.  The i'th element
   *                  contains the position of the i'th particle.
   */
  public void setPositions(double[] positions) {
    Vec3Array vec3Array = toVec3Array(positions);
    OpenMM_Context_setPositions(pointer, vec3Array.getPointer());
    vec3Array.destroy();
  }

  /**
   * Copy information from a State object into this Context.  This restores the Context to
   * approximately the same state it was in when the State was created.  If the State does not include
   * a piece of information (e.g. positions or velocities), that aspect of the Context is
   * left unchanged.
   * <p>
   * Even when all possible information is included in the State, the effect of calling this method
   * is still less complete than loadCheckpoint().  For example, it does not restore the internal
   * states of random number generators.  On the other hand, it has the advantage of not being hardware
   * specific.
   *
   * @param state the State object to copy information from
   */
  public void setState(State state) {
    OpenMM_Context_setState(pointer, state.getPointer());
  }

  /**
   * Set the current step count.
   *
   * @param steps the current step count
   */
  public void setStepCount(long steps) {
    OpenMM_Context_setStepCount(pointer, steps);
  }

  /**
   * Set the current time of the simulation (in picoseconds).
   *
   * @param time the current time of the simulation (in picoseconds)
   */
  public void setTime(double time) {
    OpenMM_Context_setTime(pointer, time);
  }

  /**
   * Set the velocities of all particles in the System (measured in nm/picosecond).
   *
   * @param velocities a vector whose length equals the number of particles in the System.  The i'th element
   *                   contains the velocity of the i'th particle.
   */
  public void setVelocities(double[] velocities) {
    Vec3Array velArray = toVec3Array(velocities);
    OpenMM_Context_setVelocities(pointer, velArray.getPointer());
    velArray.destroy();
  }

  /**
   * Set the velocities of all particles in the System to random values chosen from a Boltzmann
   * distribution at a given temperature.
   *
   * @param temperature the temperature for which to select the velocities (measured in Kelvin)
   * @param randomSeed  the random number seed to use when selecting velocities
   */
  public void setVelocitiesToTemperature(double temperature, int randomSeed) {
    OpenMM_Context_setVelocitiesToTemperature(pointer, temperature, randomSeed);
  }

  /**
   * Update the context.
   *
   * @param system     The system to simulate.
   * @param integrator The integrator to use for simulating the system.
   * @param platform   The platform to use for performing computations.
   */
  public void updateContext(System system, Integrator integrator, Platform platform) {
    // Destroy the old context.
    destroy();
    pointer = OpenMM_Context_create_2(system.getPointer(), integrator.getPointer(), platform.getPointer());
    this.integrator = integrator;
    this.platform = platform;
  }
}