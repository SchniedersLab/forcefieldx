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

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_addIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getCurrentIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getNumIntegrators;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setCurrentIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_step;

/**
 * This class allows you to use multiple integration algorithms within a single simulation,
 * switching back and forth between them.  To use it, create whatever other Integrators
 * you need, then add all of them to a CompoundIntegrator:
 * <pre>
 *    CompoundIntegrator compoundIntegrator = new CompoundIntegrator();
 *    compoundIntegrator.addIntegrator(new VerletIntegrator(0.001));
 *    compoundIntegrator.addIntegrator(new LangevinIntegrator(300.0, 1.0, 0.001));
 * </pre>
 * <p>
 * Next create a Context, specifying the CompoundIntegrator as the Integrator to use for
 * the Context:
 * <pre>
 *     Context context = new Context(system, compoundIntegrator);
 * </pre>
 * <p>
 * Finally, call setCurrentIntegrator() to set which Integrator is active.  That one will
 * be used for all calls to step() until the next time you change it.
 * <pre>
 *     compoundIntegrator.setCurrentIntegrator(0);
 *     compoundIntegrator.step(1000); // Take 1000 steps of Verlet dynamics
 *     compoundIntegrator.setCurrentIntegrator(1);
 *     compoundIntegrator.step(1000); // Take 1000 steps of Langevin dynamics
 * </pre>
 * <p>
 * When switching between integrators, it is important to make sure they are compatible with
 * each other, and that they will interpret the positions and velocities in the same way.
 * Remember that leapfrog style integrators assume the positions and velocities are offset
 * from each other by half a time step.  When switching between a leapfrog and non-leapfrog
 * integrator, you must first adjust the velocities to avoid introducing error.  This is also
 * true when switching between two leapfrog integrators that use different step sizes,
 * since they will interpret the velocities as corresponding to different times.
 */
public class CompoundIntegrator extends Integrator {

  /**
   * Create a CompoundIntegrator.
   */
  public CompoundIntegrator() {
    super(OpenMM_CompoundIntegrator_create());
  }

  /**
   * Add an Integrator to this CompoundIntegrator.  The Integrator object should have
   * been created on the heap with the "new" operator.  The CompoundIntegrator takes over
   * ownership of it, and deletes it when the CompoundIntegrator itself is deleted.
   * All Integrators must be added before the Context is created.
   *
   * @param integrator the Integrator to add
   * @return the index of the Integrator that was added
   */
  public int addIntegrator(Integrator integrator) {
    return OpenMM_CompoundIntegrator_addIntegrator(pointer, integrator.getPointer());
  }

  /**
   * Destroy the integrator.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CompoundIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
   * This method calls getConstraintTolerance() on whichever Integrator has been set as current.
   */
  @Override
  public double getConstraintTolerance() {
    return OpenMM_CompoundIntegrator_getConstraintTolerance(pointer);
  }

  /**
   * Get the index of the current Integrator.
   */
  public int getCurrentIntegrator() {
    return OpenMM_CompoundIntegrator_getCurrentIntegrator(pointer);
  }

  /**
   * Get which force groups to use for integration.  By default, all force groups
   * are included.  This is interpreted as a set of bit flags: the forces from group i
   * will be included if (groups&amp;(1&lt;&lt;i)) != 0.
   * <p>
   * This method returns the integration force groups for the current Integrator.
   */
  @Override
  public int getIntegrationForceGroups() {
    return OpenMM_CompoundIntegrator_getIntegrationForceGroups(pointer);
  }

  /**
   * Get a reference to one of the Integrators that have been added to this CompoundIntegrator.
   *
   * @param index the index of the Integrator to get
   */
  public PointerByReference getIntegrator(int index) {
    return OpenMM_CompoundIntegrator_getIntegrator(pointer, index);
  }

  /**
   * Get the number of Integrators that have been added to this CompoundIntegrator.
   */
  public int getNumIntegrators() {
    return OpenMM_CompoundIntegrator_getNumIntegrators(pointer);
  }

  /**
   * Get the size of each time step, in picoseconds.  This method calls getStepSize() on
   * whichever Integrator has been set as current.
   *
   * @return the step size, measured in ps
   */
  @Override
  public double getStepSize() {
    return OpenMM_CompoundIntegrator_getStepSize(pointer);
  }

  /**
   * Set the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
   * This method calls setConstraintTolerance() on whichever Integrator has been set as current.
   */
  @Override
  public void setConstraintTolerance(double tol) {
    OpenMM_CompoundIntegrator_setConstraintTolerance(pointer, tol);
  }

  /**
   * Set the current Integrator.
   *
   * @param index the index of the Integrator to use
   */
  public void setCurrentIntegrator(int index) {
    OpenMM_CompoundIntegrator_setCurrentIntegrator(pointer, index);
  }

  /**
   * Set which force groups to use for integration.  By default, all force groups
   * are included.  This is interpreted as a set of bit flags: the forces from group i
   * will be included if (groups&amp;(1&lt;&lt;i)) != 0.
   * <p>
   * Calling this method sets the integration force groups for all Integrators
   * contained in this CompoundIntegrator.
   */
  @Override
  public void setIntegrationForceGroups(int groups) {
    OpenMM_CompoundIntegrator_setIntegrationForceGroups(pointer, groups);
  }

  /**
   * Set the size of each time step, in picoseconds.  This method calls setStepSize() on
   * whichever Integrator has been set as current.
   *
   * @param size the step size, measured in ps
   */
  @Override
  public void setStepSize(double size) {
    OpenMM_CompoundIntegrator_setStepSize(pointer, size);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.  This method
   * calls step() on whichever Integrator has been set as current.
   *
   * @param steps the number of time steps to take
   */
  @Override
  public void step(int steps) {
    OpenMM_CompoundIntegrator_step(pointer, steps);
  }
}