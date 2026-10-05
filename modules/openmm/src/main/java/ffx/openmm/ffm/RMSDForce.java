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
 * Computes the aligned root-mean-square deviation between current coordinates and a reference
 * structure. The selected particles are aligned to the reference before RMSD is evaluated; an
 * empty particle selection means all system particles. The reference-position array must contain
 * one position for every system particle, including unselected particles. Positions and the
 * resulting RMSD are in nm. This force is intended particularly as a collective variable for use
 * with {@code CustomCVForce}.
 */
public class RMSDForce extends Force {

  /**
   * Create an RMSD force from FFM array wrappers.
   *
   * @param referencePositions reference position for every system particle, in system particle
   *     order and nm
   * @param particles indices of particles included in the alignment and RMSD; an empty array
   *     selects all particles
   */
  public RMSDForce(Vec3Array referencePositions, IntArray particles) {
    super(create(referencePositions, particles));
  }

  /**
   * Create from FFM array wrappers in the legacy JNA façade's particles-first argument order.
   *
   * @param particles indices of particles included in the alignment and RMSD; an empty array
   *     selects all particles
   * @param referencePositions reference position for every system particle, in system particle
   *     order and nm
   */
  public RMSDForce(IntArray particles, Vec3Array referencePositions) {
    this(referencePositions, particles);
  }

  /**
   * Create from raw native handles using the JNA façade's particles-first argument order.
   *
   * @param particles native integer-array handle containing selected particle indices; an empty
   *     array selects all particles
   * @param referencePositions native vector-array handle with one nm reference position per
   *     system particle
   */
  public RMSDForce(MemorySegment particles, MemorySegment referencePositions) {
    this(create(particles, referencePositions));
  }

  /**
   * Create an RMSD force from Java vectors.
   *
   * @param referencePositions reference position for every system particle, in system particle
   *     order and nm; the values are copied to temporary native storage
   * @param particles indices of particles included in the alignment and RMSD; an empty array
   *     selects all particles
   */
  public RMSDForce(Vec3[] referencePositions, int[] particles) {
    this(create(referencePositions, particles));
  }

  /**
   * Create an RMSD force from packed Java coordinates and particle indices.
   *
   * @param referencePositions reference coordinates packed as {@code x,y,z} triples, one nm
   *     position per system particle in system order
   * @param particles indices of particles included in the alignment and RMSD; an empty array
   *     selects all particles
   */
  public RMSDForce(double[] referencePositions, int[] particles) {
    this(create(referencePositions, particles));
  }

  private RMSDForce(MemorySegment pointer) {
    super(pointer);
  }

  /** Destroy the native RMSD force. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_RMSDForce_destroy);
  }

  /**
   * Get selected particle indices as a Java-owned copy. An empty array means all particles.
   *
   * @return copied selected particle indices
   */
  public int[] getParticles() {
    try (IntArray particles = new IntArray(0)) {
      copy(OpenMMNative.OpenMM_RMSDForce_getParticles(getPointer()), particles);
      return copy(particles);
    }
  }

  /**
   * Get a borrowed native handle to the particle-index array.
   *
   * <p>The array is owned by this force. Do not free or destroy the returned handle; it becomes
   * invalid when the force is destroyed.
   *
   * @return borrowed native array handle
   */
  public MemorySegment getParticlesHandle() {
    return OpenMMNative.OpenMM_RMSDForce_getParticles(getPointer());
  }

  /**
   * Get reference positions as Java-owned {@link Vec3} values, copied from the native array.
   *
   * @return copied reference positions in system particle order and nm
   */
  public Vec3[] getReferencePositions() {
    MemorySegment nativeValues =
        OpenMMNative.OpenMM_RMSDForce_getReferencePositions(getPointer());
    Vec3Array borrowed = new Vec3Array(nativeValues);
    try {
      Vec3[] values = new Vec3[borrowed.getSize()];
      for (int i = 0; i < values.length; i++) {
        values[i] = borrowed.get(i);
      }
      return values;
    } finally {
      borrowed.invalidate();
    }
  }

  /**
   * Get a borrowed native handle to the reference-position array.
   *
   * <p>The array is owned by this force. Do not free or destroy the returned handle; it becomes
   * invalid when the force is destroyed.
   *
   * @return borrowed native array handle
   */
  public MemorySegment getReferencePositionsHandle() {
    return OpenMMNative.OpenMM_RMSDForce_getReferencePositions(getPointer());
  }

  /**
   * Set particle indices from an FFM integer array; OpenMM copies the values.
   *
   * @param particles selected particle indices; empty selects all system particles
   */
  public void setParticles(IntArray particles) {
    OpenMMNative.OpenMM_RMSDForce_setParticles(getPointer(), particles.getPointer());
  }

  /**
   * Set particle indices from a native integer-array handle; OpenMM copies the values during the
   * call and the caller retains ownership of the handle.
   *
   * @param particles native array of selected particle indices; empty selects all system
   *     particles
   */
  public void setParticles(MemorySegment particles) {
    OpenMMNative.OpenMM_RMSDForce_setParticles(getPointer(), particles);
  }

  /**
   * Set particle indices from Java values, copied through temporary native storage.
   *
   * @param particles selected particle indices; empty selects all system particles
   */
  public void setParticles(int[] particles) {
    try (IntArray values = toIntegers(particles)) {
      setParticles(values);
    }
  }

  /**
   * Set reference positions from an FFM vector array; OpenMM copies the values.
   *
   * @param positions one reference position in nm for every system particle, in system order
   */
  public void setReferencePositions(Vec3Array positions) {
    OpenMMNative.OpenMM_RMSDForce_setReferencePositions(getPointer(), positions.getPointer());
  }

  /**
   * Set reference positions from a native vector-array handle; OpenMM copies the values during
   * the call and the caller retains ownership of the handle.
   *
   * @param positions native array containing one nm reference position per system particle
   */
  public void setReferencePositions(MemorySegment positions) {
    OpenMMNative.OpenMM_RMSDForce_setReferencePositions(getPointer(), positions);
  }

  /**
   * Set reference positions from Java vectors, copied through temporary native storage.
   *
   * @param positions one reference position in nm for every system particle, in system order
   */
  public void setReferencePositions(Vec3[] positions) {
    try (Vec3Array values = toVectors(positions)) {
      setReferencePositions(values);
    }
  }

  /**
   * Set reference positions from packed Java coordinates, copied through temporary native
   * storage.
   *
   * @param positions coordinates packed as {@code x,y,z} triples, one nm position per system
   *     particle in system order
   */
  public void setReferencePositions(double[] positions) {
    try (Vec3Array values = Vec3Array.toVec3Array(positions)) {
      setReferencePositions(values);
    }
  }

  /**
   * Copy the current reference positions and particle selection to an existing context without
   * reinitializing it.
   *
   * @param context context to update; this FFM wrapper performs no update if it has no native
   *     context handle
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_RMSDForce_updateParametersInContext(getPointer(), pointer));
  }

  /** RMSDForce does not use periodic boundary conditions. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_RMSDForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(Vec3Array positions, IntArray particles) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_RMSDForce_create(
        positions.getPointer(), particles.getPointer());
  }

  private static MemorySegment create(
      MemorySegment particles, MemorySegment referencePositions) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_RMSDForce_create(referencePositions, particles);
  }

  private static MemorySegment create(Vec3[] positions, int[] particles) {
    try (Vec3Array nativePositions = toVectors(positions);
         IntArray nativeParticles = toIntegers(particles)) {
      return create(nativePositions, nativeParticles);
    }
  }

  private static MemorySegment create(double[] positions, int[] particles) {
    try (Vec3Array nativePositions = Vec3Array.toVec3Array(positions);
         IntArray nativeParticles = toIntegers(particles)) {
      return create(nativePositions, nativeParticles);
    }
  }

  private static Vec3Array toVectors(Vec3[] positions) {
    Vec3Array result = new Vec3Array(0);
    for (Vec3 position : positions) {
      result.append(position);
    }
    return result;
  }

  private static IntArray toIntegers(int[] values) {
    IntArray result = new IntArray(values.length);
    for (int i = 0; i < values.length; i++) {
      result.set(i, values[i]);
    }
    return result;
  }

  private static void copy(MemorySegment nativeArray, IntArray destination) {
    IntArray source = new IntArray(nativeArray);
    try {
      destination.resize(source.getSize());
      for (int i = 0; i < source.getSize(); i++) {
        destination.set(i, source.get(i));
      }
    } finally {
      source.invalidate();
    }
  }

  private static int[] copy(IntArray values) {
    int[] result = new int[values.getSize()];
    for (int i = 0; i < result.length; i++) {
      result[i] = values.get(i);
    }
    return result;
  }
}
