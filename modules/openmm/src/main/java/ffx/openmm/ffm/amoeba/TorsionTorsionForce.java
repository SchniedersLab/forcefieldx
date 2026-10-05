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
package ffx.openmm.ffm.amoeba;

import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * AMOEBA torsion-torsion interaction force. Each term couples two adjacent torsions through a
 * two-dimensional energy grid. Define grids before terms that reference them; grid values contain
 * torsion angles, energies, and optionally derivatives. This is normally a bonded force without
 * periodic displacements unless explicitly enabled. This façade exposes no
 * {@code updateParametersInContext} operation; changes to terms, grids, or periodic-boundary
 * settings after Context creation require Context reinitialization.
 */
public class TorsionTorsionForce extends Force {

  /**
   * Particle indices, chirality-check atom, and energy-grid index for one term.
   *
   * @param particle1 first particle in the torsion-torsion chain.
   * @param particle2 second particle in the chain.
   * @param particle3 central particle in the chain.
   * @param particle4 fourth particle in the chain.
   * @param particle5 fifth particle in the chain.
   * @param chiralCheckAtomIndex particle connected to particle3 but not particle2 or particle4,
   *     used for the chirality check.
   * @param gridIndex index of the associated energy grid.
   */
  public record TorsionParameters(
      int particle1, int particle2, int particle3, int particle4, int particle5,
      int chiralCheckAtomIndex, int gridIndex) {}

  /**
   * Create an empty AMOEBA torsion-torsion force.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public TorsionTorsionForce() {
    super(create());
  }

  /**
   * Add a torsion-torsion term to the force field.
   *
   * @param particle1 index of the first particle connected by the torsion-torsion.
   * @param particle2 index of the second particle connected by the torsion-torsion.
   * @param particle3 index of the third particle connected by the torsion-torsion.
   * @param particle4 index of the fourth particle connected by the torsion-torsion.
   * @param particle5 index of the fifth particle connected by the torsion-torsion.
   * @param chiralCheckAtomIndex index of the particle connected to particle3, but not particle2 or
   *     particle4, used in the chirality check.
   * @param gridIndex index of the grid to use.
   * @return index of the torsion-torsion term added.
   */
  public int addTorsionTorsion(
      int particle1, int particle2, int particle3, int particle4, int particle5,
      int chiralCheckAtomIndex, int gridIndex) {
    return OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(
        getPointer(), particle1, particle2, particle3, particle4, particle5,
        chiralCheckAtomIndex, gridIndex);
  }

  /** Release the native force owned by this façade. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AmoebaTorsionTorsionForce_destroy);
  }

  /**
   * Get the number of defined torsion-torsion grids.
   *
   * @return grid count.
   */
  public int getNumTorsionTorsionGrids() {
    return OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(getPointer());
  }

  /**
   * Get the number of defined torsion-torsion terms.
   *
   * @return term count.
   */
  public int getNumTorsionTorsions() {
    return OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(getPointer());
  }

  /**
   * Get a term's particle indices and grid index as copied Java values.
   *
   * @param index torsion-torsion term index.
   * @return copied particle, chirality-check, and grid indices.
   */
  public TorsionParameters getTorsionTorsionParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment[] values = new MemorySegment[7];
      for (int i = 0; i < values.length; i++) {
        values[i] = arena.allocate(ValueLayout.JAVA_INT);
      }
      OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(
          getPointer(), index, values[0], values[1], values[2], values[3],
          values[4], values[5], values[6]);
      return new TorsionParameters(
          values[0].get(ValueLayout.JAVA_INT, 0),
          values[1].get(ValueLayout.JAVA_INT, 0),
          values[2].get(ValueLayout.JAVA_INT, 0),
          values[3].get(ValueLayout.JAVA_INT, 0),
          values[4].get(ValueLayout.JAVA_INT, 0),
          values[5].get(ValueLayout.JAVA_INT, 0),
          values[6].get(ValueLayout.JAVA_INT, 0));
    }
  }

  /**
   * Get the native grid handle at the specified index.
   *
   * <p>The returned grid is borrowed from this force. Do not destroy it; it is valid only while
   * this force remains alive and until that grid is replaced.</p>
   *
   * @param index grid index.
   * @return borrowed native grid handle.
   */
  public MemorySegment getTorsionTorsionGrid(int index) {
    return OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(
        getPointer(), index);
  }

  /**
   * Replace a torsion-torsion term's particle indices and grid index.
   *
   * @param index torsion-torsion term index.
   * @param particle1 first particle in the torsion-torsion chain.
   * @param particle2 second particle in the chain.
   * @param particle3 central particle in the chain.
   * @param particle4 fourth particle in the chain.
   * @param particle5 fifth particle in the chain.
   * @param chiralCheckAtomIndex particle connected to particle3 but not particle2 or particle4,
   *     used for the chirality check.
   * @param gridIndex index of the energy grid used by this term.
   */
  public void setTorsionTorsionParameters(
      int index, int particle1, int particle2, int particle3, int particle4, int particle5,
      int chiralCheckAtomIndex, int gridIndex) {
    OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(
        getPointer(), index, particle1, particle2, particle3, particle4, particle5,
        chiralCheckAtomIndex, gridIndex);
  }

  /**
   * Set a torsion-torsion grid. Each point supplies either three values (the two torsion
   * coordinates and energy) or six values (coordinates, energy, two first derivatives, and the
   * mixed derivative). If derivatives are omitted, OpenMM fits a two-dimensional spline to
   * calculate them. The header describes these as torsion-angle coordinates, energy, and energy
   * derivatives but does not specify the angle unit; use a consistent unit for grid coordinates
   * and derivatives. Energies are in OpenMM energy units.
   *
   * @param gridIndex index of the grid to replace.
   * @param grid caller-owned native grid handle; its contents are copied into the force.
   * @throws NullPointerException if {@code grid} is null.
   */
  public void setTorsionTorsionGrid(int gridIndex, DoubleArray3D grid) {
    OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(
        getPointer(), gridIndex, grid.getPointer());
  }

  /**
   * Set a torsion-torsion grid using its native handle. The native force copies the grid; the
   * caller retains ownership of the supplied handle.
   *
   * @param gridIndex index of the grid to replace.
   * @param grid valid caller-owned native {@code OpenMM_3D_DoubleArray} handle.
   */
  public void setTorsionTorsionGrid(int gridIndex, MemorySegment grid) {
    OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(
        getPointer(), gridIndex, grid);
  }

  /**
   * Set whether this bonded force applies periodic boundary conditions when calculating
   * displacements. This is usually unnecessary for bonded terms but can be useful in some systems.
   *
   * @param periodic true to use periodic boundary conditions.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Report whether this force uses periodic boundary conditions for displacements.
   *
   * @return true if periodic displacements are enabled.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(
            getPointer()));
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_AmoebaTorsionTorsionForce_create.makeInvoker().apply();
  }
}
