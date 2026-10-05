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
 * Base class for FFM-backed OpenMM integrators.
 *
 * <p>An integrator is associated with at most one active OpenMM context. In this FFM façade,
 * {@link Context} destroys its associated integrator when that context is destroyed.</p>
 */
public abstract class Integrator extends OpenMMHandle {

  /**
   * Wrap a native integrator handle.
   *
   * @param pointer native OpenMM integrator handle.
   */
  public Integrator(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * @return constraint tolerance as a fraction of each constrained distance.
   */
  public double getConstraintTolerance() {
    return OpenMMNative.OpenMM_Integrator_getConstraintTolerance(getPointer());
  }

  /**
   * @return bit mask of force groups evaluated by this integrator; group {@code i} is included
   *         when bit {@code i} is set.
   */
  public int getIntegrationForceGroups() {
    return OpenMMNative.OpenMM_Integrator_getIntegrationForceGroups(getPointer());
  }

  /**
   * @return integration time step in picoseconds.
   */
  public double getStepSize() {
    return OpenMMNative.OpenMM_Integrator_getStepSize(getPointer());
  }

  /**
   * Set the distance tolerance within which constraints must be maintained, as a fraction of the
   * constrained distance.
   *
   * @param tolerance constraint tolerance as a fraction of the constrained distance.
   */
  public void setConstraintTolerance(double tolerance) {
    OpenMMNative.OpenMM_Integrator_setConstraintTolerance(getPointer(), tolerance);
  }

  /**
   * Set the force groups evaluated by this integrator. This is a bit mask; group {@code i} is
   * included when bit {@code i} is set.
   *
   * @param groups force-group bit mask.
   */
  public void setIntegrationForceGroups(int groups) {
    OpenMMNative.OpenMM_Integrator_setIntegrationForceGroups(getPointer(), groups);
  }

  /**
   * Rebind this wrapper to a native integrator handle without releasing the previous handle.
   *
   * <p>The previously referenced handle is not destroyed. The caller remains responsible for it,
   * matching the pointer-replacement behavior of the JNA façade.</p>
   *
   * @param pointer non-null replacement native integrator handle with a nonzero address.
   */
  public void setPointer(MemorySegment pointer) {
    rebindPointer(pointer);
  }

  /**
   * Set the integration time step. The native API leaves the effect undefined for integrators
   * that use variable time steps.
   *
   * @param stepSize time step in picoseconds.
   */
  public void setStepSize(double stepSize) {
    OpenMMNative.OpenMM_Integrator_setStepSize(getPointer(), stepSize);
  }

  /**
   * Advance the simulation in the context associated with this integrator.
   *
   * @param steps number of integration time steps to take.
   */
  public void step(int steps) {
    OpenMMNative.OpenMM_Integrator_step(getPointer(), steps);
  }
}
