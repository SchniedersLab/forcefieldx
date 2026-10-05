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
 * An interaction between particle pairs that varies harmonically with their distance, using the potential
 * {@code E = 1/2 k (r - length)^2}.
 *
 * <p>Create the force, call {@link #addBond(int, int, double, double)} once for each bond, then add the force to a
 * {@link System}. After a bond has been added, its parameters can be modified with {@link
 * #setBondParameters(int, int, int, double, double)}; this has no effect on contexts that already exist unless
 * {@link #updateParametersInContext(Context)} is called.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.HarmonicBondForce} is not;
 * {@link #getBondParameters(int)} returns a {@link BondParameters} record, whereas the JNA method returns void and
 * fills {@code IntByReference}/{@code IntBuffer}/{@code DoubleByReference} output arguments. The OpenMM header
 * does not document the energy expression; it is stated here from the standard harmonic form.</p>
 */
public class HarmonicBondForce extends Force {

  /**
   * Create a native force. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The JNA counterpart calls the native create function directly without this step.</p>
   */
  public HarmonicBondForce() {
    super(create());
  }

  /**
   * Add a bond term to the force field.
   *
   * @param particle1 index of the first particle connected by the bond.
   * @param particle2 index of the second particle connected by the bond.
   * @param length    equilibrium length of the bond, in nm.
   * @param k         harmonic force constant of the bond, in kJ/(mol nm^2).
   * @return index of the bond that was added.
   */
  public int addBond(int particle1, int particle2, double length, double k) {
    return OpenMMNative.OpenMM_HarmonicBondForce_addBond(
        getPointer(), particle1, particle2, length, k);
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_HarmonicBondForce_destroy);
  }

  /**
   * Get the force field parameters for a bond term.
   *
   * <p>Native out-parameters are read into confined temporary memory and copied into the returned record; no
   * native memory is retained. The JNA counterpart instead returns void and fills caller-supplied reference/buffer
   * arguments.</p>
   *
   * @param index index of the bond for which to get parameters.
   * @return a {@link BondParameters} record holding the particle indices, equilibrium length in nm and force
   *     constant in kJ/(mol nm^2).
   */
  public BondParameters getBondParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment particle1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment particle2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment length = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment k = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_HarmonicBondForce_getBondParameters(
          getPointer(), index, particle1, particle2, length, k);
      return new BondParameters(
          particle1.get(ValueLayout.JAVA_INT, 0),
          particle2.get(ValueLayout.JAVA_INT, 0),
          length.get(ValueLayout.JAVA_DOUBLE, 0),
          k.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the number of harmonic bond stretch terms in the potential function.
   *
   * @return number of bond terms.
   */
  public int getNumBonds() {
    return OpenMMNative.OpenMM_HarmonicBondForce_getNumBonds(getPointer());
  }

  /**
   * Set the force field parameters for a bond term.
   *
   * <p>Existing contexts are not affected until {@link #updateParametersInContext(Context)} is called, and that
   * call cannot apply a change of particles.</p>
   *
   * @param index     index of the bond for which to set parameters.
   * @param particle1 index of the first particle connected by the bond.
   * @param particle2 index of the second particle connected by the bond.
   * @param length    equilibrium length of the bond, in nm.
   * @param k         harmonic force constant of the bond, in kJ/(mol nm^2).
   */
  public void setBondParameters(int index, int particle1, int particle2, double length, double k) {
    OpenMMNative.OpenMM_HarmonicBondForce_setBondParameters(
        getPointer(), index, particle1, particle2, length, k);
  }

  /**
   * Set whether this force applies periodic boundary conditions when calculating displacements.
   *
   * <p>The OpenMM header notes this is usually not appropriate for bonded forces, but can be useful in some
   * situations. The Java boolean is converted with {@link OpenMMBooleans#toNative(boolean)}. This setting is
   * not covered by {@link #updateParametersInContext(Context)}; the OpenMM header documents that method as
   * updating only per-bond parameter values, so apply it before creating a {@link Context}.</p>
   *
   * @param periodic true to apply periodic boundary conditions to displacements, false otherwise.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Copy the per-bond parameters stored in this force into an existing {@link Context}.
   *
   * <p>Call the parameter setters first, then this method; changes made to this force have no effect on contexts
   * that already exist until it is called. Only the values of per-bond parameters (equilibrium length and force constant) are updated. The
   * set of particles involved in a bond cannot be changed and new bonds cannot be added.</p>
   *
   * @param context live context created from a system containing this force; must not be null.
   * @throws NullPointerException if {@code context} is null.
   * @throws IllegalStateException if this force or the context has been destroyed.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_HarmonicBondForce_updateParametersInContext(
        getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * <p>Overrides {@link Force#usesPeriodicBoundaryConditions()} with the HarmonicBondForce-specific native call. The native
   * flag is the value set by {@link #setUsesPeriodicBoundaryConditions(boolean)}.</p>
   *
   * @return true when periodic boundary conditions are applied to displacements, false otherwise.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Immutable copy of the parameters of one harmonic bond term, as returned by {@link #getBondParameters(int)}.
   *
   * @param particle1 index of the first particle connected by the bond.
   * @param particle2 index of the second particle connected by the bond.
   * @param length    equilibrium length of the bond, in nm.
   * @param k         harmonic force constant of the bond, in kJ/(mol nm^2).
   */
  public record BondParameters(int particle1, int particle2, double length, double k) {
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_HarmonicBondForce_create.makeInvoker().apply();
  }
}
