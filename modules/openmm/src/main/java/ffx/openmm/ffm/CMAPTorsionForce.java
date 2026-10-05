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

/**
 * Adds an energy correction map (CMAP) interaction between pairs of dihedral angles.
 *
 * <p>Each map is a {@code size} by {@code size} grid of energies, interpolated with a natural
 * cubic spline at arbitrary angle pairs. Map values use OpenMM energy units (kJ/mol), and angles
 * are in radians. Create maps before torsion terms that refer to their indices. Each torsion
 * specifies four particle indices for the first dihedral followed by four for the second.
 */
public class CMAPTorsionForce extends Force {

  /**
   * A CMAP map's size and copied energy values.
   *
   * @param size number of grid points along each angle dimension
   * @param energy Java-owned copy of {@code size*size} energies in kJ/mol, laid out as
   *     {@code energy[i + size*j]}; {@code i} indexes the first angle and {@code j} the second
   *     angle, at grid angles {@code i*2*pi/size} and {@code j*2*pi/size}
   */
  public record MapParameters(int size, double[] energy) {}
  /**
   * A torsion term's map and particle indices.
   *
   * @param map index of the energy map used by this torsion term
   * @param a1 first particle of the first dihedral
   * @param a2 second particle of the first dihedral
   * @param a3 third particle of the first dihedral
   * @param a4 fourth particle of the first dihedral
   * @param b1 first particle of the second dihedral
   * @param b2 second particle of the second dihedral
   * @param b3 third particle of the second dihedral
   * @param b4 fourth particle of the second dihedral
   */
  public record TorsionParameters(
      int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {}

  /** Create an empty CMAP torsion force with no maps or torsion terms. */
  public CMAPTorsionForce() {
    super(create());
  }

  /**
   * Add a square energy map.
   *
   * @param size number of grid points along each angle dimension
   * @param energy native array containing {@code size*size} energies in kJ/mol, laid out as
   *     {@code energy[i + size*j]}; {@code i} indexes the first angle and {@code j} the second,
   *     sampled at {@code i*2*pi/size} and {@code j*2*pi/size}
   * @return index of the newly added map
   */
  public int addMap(int size, DoubleArray energy) {
    return OpenMMNative.OpenMM_CMAPTorsionForce_addMap(
        getPointer(), size, energy.getPointer());
  }

  /**
   * Add an energy map from a caller-provided native double-array handle.
   *
   * @param size number of grid points along each angle dimension
   * @param energy native handle to {@code size*size} kJ/mol values in
   *     {@code i + size*j} order; the handle remains caller-owned
   * @return index of the newly added map
   */
  public int addMap(int size, MemorySegment energy) {
    return OpenMMNative.OpenMM_CMAPTorsionForce_addMap(getPointer(), size, energy);
  }

  /**
   * Add an energy map from Java values; the values are copied into temporary native storage and
   * OpenMM retains its own map data.
   *
   * @param size number of grid points along each angle dimension
   * @param energy {@code size*size} energies in kJ/mol, in {@code i + size*j} order
   * @return index of the newly added map
   */
  public int addMap(int size, double[] energy) {
    try (DoubleArray values = toNative(energy)) {
      return addMap(size, values);
    }
  }

  /**
   * Add a torsion term using the specified map and ordered particle quartets.
   *
   * @param map index of a previously added energy map
   * @param a1 first particle of the first dihedral
   * @param a2 second particle of the first dihedral
   * @param a3 third particle of the first dihedral
   * @param a4 fourth particle of the first dihedral
   * @param b1 first particle of the second dihedral
   * @param b2 second particle of the second dihedral
   * @param b3 third particle of the second dihedral
   * @param b4 fourth particle of the second dihedral
   * @return index of the newly added torsion term
   */
  public int addTorsion(
      int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    return OpenMMNative.OpenMM_CMAPTorsionForce_addTorsion(
        getPointer(), map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /** Destroy the native CMAP force. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CMAPTorsionForce_destroy);
  }

  /**
   * Get a map's size and energy values as a Java-owned copy.
   *
   * @param index index of the map
   * @return map size and copied {@code size*size} energies in kJ/mol and {@code i + size*j}
   *     order
   */
  public MapParameters getMapParameters(int index) {
    try (Arena arena = Arena.ofConfined(); DoubleArray energy = new DoubleArray(0)) {
      MemorySegment size = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CMAPTorsionForce_getMapParameters(
          getPointer(), index, size, energy.getPointer());
      return new MapParameters(size.get(ValueLayout.JAVA_INT, 0), copy(energy));
    }
  }

  /**
   * Copy map values into a caller-owned native array and return the map size.
   *
   * @param index index of the map
   * @param energy destination native array, which must have capacity for the map's
   *     {@code size*size} energies; the caller retains ownership
   * @return number of grid points along each map dimension
   */
  public int getMapParameters(int index, DoubleArray energy) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment size = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CMAPTorsionForce_getMapParameters(
          getPointer(), index, size, energy.getPointer());
      return size.get(ValueLayout.JAVA_INT, 0);
    }
  }

  /**
   * Copy map values into a caller-owned native array handle and return the map size.
   *
   * @param index index of the map
   * @param energy destination native array handle, which must have capacity for the map's
   *     {@code size*size} energies; the caller retains ownership
   * @return number of grid points along each map dimension
   */
  public int getMapParameters(int index, MemorySegment energy) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment size = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CMAPTorsionForce_getMapParameters(
          getPointer(), index, size, energy);
      return size.get(ValueLayout.JAVA_INT, 0);
    }
  }

  /** @return number of energy maps currently defined */
  public int getNumMaps() {
    return OpenMMNative.OpenMM_CMAPTorsionForce_getNumMaps(getPointer());
  }

  /** @return number of CMAP torsion terms currently defined */
  public int getNumTorsions() {
    return OpenMMNative.OpenMM_CMAPTorsionForce_getNumTorsions(getPointer());
  }

  /**
   * Get a torsion's map and particle indices as copied scalar values.
   *
   * @param index index of the torsion term
   * @return its map index and first-then-second dihedral particle indices
   */
  public TorsionParameters getTorsionParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment[] values = new MemorySegment[9];
      for (int i = 0; i < values.length; i++) {
        values[i] = arena.allocate(ValueLayout.JAVA_INT);
      }
      OpenMMNative.OpenMM_CMAPTorsionForce_getTorsionParameters(
          getPointer(), index, values[0], values[1], values[2], values[3], values[4],
          values[5], values[6], values[7], values[8]);
      int[] ints = new int[values.length];
      for (int i = 0; i < ints.length; i++) {
        ints[i] = values[i].get(ValueLayout.JAVA_INT, 0);
      }
      return new TorsionParameters(
          ints[0], ints[1], ints[2], ints[3], ints[4],
          ints[5], ints[6], ints[7], ints[8]);
    }
  }

  /**
   * Replace the values and size of an existing map.
   *
   * @param index index of the map
   * @param size number of grid points along each angle dimension
   * @param energy native array containing {@code size*size} kJ/mol energies in
   *     {@code i + size*j} order
   */
  public void setMapParameters(int index, int size, DoubleArray energy) {
    OpenMMNative.OpenMM_CMAPTorsionForce_setMapParameters(
        getPointer(), index, size, energy.getPointer());
  }

  /**
   * Replace an existing map using a caller-provided native array handle.
   *
   * @param index index of the map
   * @param size number of grid points along each angle dimension
   * @param energy native handle to {@code size*size} kJ/mol energies in {@code i + size*j} order;
   *     the handle remains caller-owned
   */
  public void setMapParameters(int index, int size, MemorySegment energy) {
    OpenMMNative.OpenMM_CMAPTorsionForce_setMapParameters(
        getPointer(), index, size, energy);
  }

  /**
   * Replace an existing map from Java values, copied through temporary native storage.
   *
   * @param index index of the map
   * @param size number of grid points along each angle dimension
   * @param energy {@code size*size} kJ/mol energies in {@code i + size*j} order
   */
  public void setMapParameters(int index, int size, double[] energy) {
    try (DoubleArray values = toNative(energy)) {
      setMapParameters(index, size, values);
    }
  }

  /**
   * Replace an existing torsion's map and ordered particle indices.
   *
   * @param index index of the torsion term to replace
   * @param map index of the energy map to use
   * @param a1 first particle of the first dihedral
   * @param a2 second particle of the first dihedral
   * @param a3 third particle of the first dihedral
   * @param a4 fourth particle of the first dihedral
   * @param b1 first particle of the second dihedral
   * @param b2 second particle of the second dihedral
   * @param b3 third particle of the second dihedral
   * @param b4 fourth particle of the second dihedral
   */
  public void setTorsionParameters(
      int index, int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    OpenMMNative.OpenMM_CMAPTorsionForce_setTorsionParameters(
        getPointer(), index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Set whether periodic boundary conditions are used when calculating particle displacements.
   * This is usually inappropriate for a bonded force, but can be useful in some systems.
   *
   * @param periodic {@code true} to use periodic boundary conditions
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Copy supported parameter changes into an existing context without reinitializing it.
   * Update the force's map values and/or torsion map indices before calling this method. OpenMM
   * does not allow this operation to change a map's size, the particles in a torsion, or the
   * number of maps or torsions.
   *
   * @param context context to update; this FFM wrapper performs no update if it has no native
   *     context handle
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_CMAPTorsionForce_updateParametersInContext(
                getPointer(), pointer));
  }

  /**
   * Check whether periodic boundary conditions are used for displacement calculations.
   *
   * @return {@code true} if the force uses periodic boundary conditions
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_CMAPTorsionForce_create.makeInvoker().apply();
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray result = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) {
      result.set(i, values[i]);
    }
    return result;
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int i = 0; i < result.length; i++) {
      result[i] = values.get(i);
    }
    return result;
  }
}
