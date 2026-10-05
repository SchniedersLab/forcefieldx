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

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Objects;

/**
 * An interaction between groups of four particles that varies periodically with the torsion angle between
 * them.
 *
 * <p>Call {@link #addTorsion(int, int, int, int, int, double, double)} once for each torsion, then add the force to a
 * {@link System}. After a torsion has been added, its parameters can be modified with {@link
 * #setTorsionParameters(int, int, int, int, int, int, double, double)}; this has no effect on contexts that already
 * exist unless {@link #updateParametersInContext(Context)} is called. The OpenMM header does not state the
 * energy expression.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.PeriodicTorsionForce} is not;
 * {@link #getTorsionParameters(int)} returns a {@link TorsionParameters} record, whereas the JNA methods return void
 * and fill {@code IntByReference}/{@code IntBuffer} (and double) output arguments.</p>
 */
public class PeriodicTorsionForce extends Force {
  /**
   * Create a native force. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The JNA counterpart calls the native create function directly without this step.</p>
   */
  public PeriodicTorsionForce() {
    super(create());
  }

  /**
   * Add a periodic torsion term to the force field.
   *
   * <p>The OpenMM header documents the periodic torsion force constant without a unit; kJ/mol is the unit of OpenMM
   * energy and is the unit documented by the previous FFM text, but it is not stated by the header or the JNA wrapper.</p>
   *
   * @param particle1   index of the first particle forming the torsion.
   * @param particle2   index of the second particle forming the torsion.
   * @param particle3   index of the third particle forming the torsion.
   * @param particle4   index of the fourth particle forming the torsion.
   * @param periodicity periodicity of the torsion (an integer).
   * @param phase       phase offset of the torsion, in radians.
   * @param k           force constant of the torsion (see note above on units).
   * @return index of the torsion that was added.
   */
  public int addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    return OpenMMNative.OpenMM_PeriodicTorsionForce_addTorsion(getPointer(), particle1, particle2, particle3, particle4, periodicity, phase, k);
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_PeriodicTorsionForce_destroy);
  }

  /**
   * Get the force field parameters for a periodic torsion term.
   *
   * <p>Native out-parameters are copied into the returned record; no native memory is retained. The JNA
   * counterpart returns void and fills caller-supplied reference/buffer arguments.</p>
   *
   * @param index index of the torsion for which to get parameters.
   * @return a {@link TorsionParameters} record holding the four particle indices, periodicity, phase in radians and
   *     force constant.
   */
  public TorsionParameters getTorsionParameters(int index) {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment p1 = a.allocate(ValueLayout.JAVA_INT), p2 = a.allocate(ValueLayout.JAVA_INT), p3 = a.allocate(ValueLayout.JAVA_INT), p4 = a.allocate(ValueLayout.JAVA_INT), n = a.allocate(ValueLayout.JAVA_INT), phase = a.allocate(ValueLayout.JAVA_DOUBLE), k = a.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_PeriodicTorsionForce_getTorsionParameters(getPointer(), index, p1, p2, p3, p4, n, phase, k);
      return new TorsionParameters(p1.get(ValueLayout.JAVA_INT, 0), p2.get(ValueLayout.JAVA_INT, 0), p3.get(ValueLayout.JAVA_INT, 0), p4.get(ValueLayout.JAVA_INT, 0), n.get(ValueLayout.JAVA_INT, 0), phase.get(ValueLayout.JAVA_DOUBLE, 0), k.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the number of periodic torsion terms in the potential function.
   *
   * @return number of torsion terms.
   */
  public int getNumTorsions() {
    return OpenMMNative.OpenMM_PeriodicTorsionForce_getNumTorsions(getPointer());
  }

  /**
   * Set the force field parameters for a periodic torsion term.
   *
   * <p>Existing contexts are not affected until {@link #updateParametersInContext(Context)} is called, and that
   * call cannot apply a change of particles.</p>
   *
   * @param index       index of the torsion for which to set parameters.
   * @param particle1   index of the first particle forming the torsion.
   * @param particle2   index of the second particle forming the torsion.
   * @param particle3   index of the third particle forming the torsion.
   * @param particle4   index of the fourth particle forming the torsion.
   * @param periodicity periodicity of the torsion (an integer).
   * @param phase       phase offset of the torsion, in radians.
   * @param k           force constant of the torsion (units as for {@link #addTorsion(int, int, int, int, int, double, double)}).
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    OpenMMNative.OpenMM_PeriodicTorsionForce_setTorsionParameters(getPointer(), index, particle1, particle2, particle3, particle4, periodicity, phase, k);
  }

  /**
   * Set whether this force applies periodic boundary conditions when calculating displacements.
   *
   * <p>The OpenMM header notes this is usually not appropriate for bonded forces, but can be useful in some
   * situations. The Java boolean is converted with {@link OpenMMBooleans#toNative(boolean)}. This setting is
   * not covered by {@link #updateParametersInContext(Context)}; the OpenMM header documents that method as
   * updating only per-torsion parameter values, so apply it before creating a {@link Context}.</p>
   *
   * @param periodic true to apply periodic boundary conditions to displacements, false otherwise.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions(getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Copy the per-torsion parameters stored in this force into an existing {@link Context}.
   *
   * <p>Call the parameter setters first, then this method; changes made to this force have no effect on contexts
   * that already exist until it is called. Only the values of per-torsion parameters (periodicity, phase and force constant) are updated. The
   * set of particles involved in a torsion cannot be changed and new torsions cannot be added.</p>
   *
   * @param context live context created from a system containing this force; must not be null.
   * @throws NullPointerException if {@code context} is null.
   * @throws IllegalStateException if this force or the context has been destroyed.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_PeriodicTorsionForce_updateParametersInContext(getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * <p>Overrides {@link Force#usesPeriodicBoundaryConditions()} with the PeriodicTorsionForce-specific native call. The native
   * flag is the value set by {@link #setUsesPeriodicBoundaryConditions(boolean)}.</p>
   *
   * @return true when periodic boundary conditions are applied to displacements, false otherwise.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Immutable copy of the parameters of one periodic torsion term, as returned by {@link #getTorsionParameters(int)}.
   *
   * @param particle1   index of the first particle forming the torsion.
   * @param particle2   index of the second particle forming the torsion.
   * @param particle3   index of the third particle forming the torsion.
   * @param particle4   index of the fourth particle forming the torsion.
   * @param periodicity periodicity of the torsion.
   * @param phase       phase offset of the torsion, in radians.
   * @param k           force constant of the torsion; the OpenMM header states no unit (see {@link PeriodicTorsionForce}).
   */
  public record TorsionParameters(int particle1, int particle2, int particle3, int particle4, int periodicity,
                                  double phase, double k) {
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_PeriodicTorsionForce_create.makeInvoker().apply();
  }
}
