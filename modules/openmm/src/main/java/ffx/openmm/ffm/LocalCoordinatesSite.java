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
import java.util.Objects;

/**
 * A virtual site placed at a fixed location in a local coordinate system computed from other particles. The origin
 * and the directions of the x and y axes are each a weighted sum of the particle positions:
 * {@code origin = sum(wo_i r_i)}, {@code xdir = sum(wx_i r_i)}, {@code ydir = sum(wy_i r_i)}. The origin weights must add
 * to one (weights of {@code (1, 0, 0)} put the origin at particle 1) and the x and y weights must each add to zero (for
 * example {@code (-1, 0.5, 0.5)} points the x axis from particle 1 to the midpoint of particles 2 and 3).
 *
 * <p>OpenMM computes {@code zdir = xdir x ydir}, recomputes {@code ydir = zdir x xdir} to make the axes orthogonal,
 * normalizes all three axes, and places the site at {@code origin + x*xdir + y*ydir + z*zdir} using the local position.</p>
 *
 * <p>Wrapper/header notes: inputs are copied into temporary native memory or read from the supplied arrays, which the caller
 * keeps owning. The header's {@code Vec3}-returning {@code getOriginWeights()}, {@code getXWeights()} and {@code
 * getYWeights()} overloads (which throw unless exactly three particles define the site) are not mapped; this class returns
 * {@code double[]} of one value per defining particle, which works for any particle count. The header gives no unit for the
 * local position; as with OpenMM positions it is a distance (nm).</p>
 *
 * <p>Ownership: this wrapper owns the native virtual site until it is passed to {@link System#setVirtualSite(int,
 * VirtualSite)}. OpenMM then assumes native ownership, and the system invalidates this wrapper when it is destroyed, so a
 * site given to a system should not also be destroyed directly. A virtual site has no {@code updateParametersInContext}
 * operation and is not a {@link Force}.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.LocalCoordinatesSite} is not. The JNA
 * first constructor takes {@code (originWeights, xWeights, yWeights, zWeights, localPosition)} as {@code PointerByReference}
 * values, with an extra {@code zWeights} and no particle list, which does not match the header's {@code (particles,
 * originWeights, xWeights, yWeights, localPosition)}; this class follows the header. The JNA three-particle constructor
 * names its direction arguments {@code xdir} and {@code ydir}, where the header calls them {@code xWeights} and {@code
 * yWeights}. The JNA getters fill a {@code PointerByReference}.</p>
 */
public class LocalCoordinatesSite extends VirtualSite {

  /**
   * Create a local-coordinate site from variable-length particle and weight arrays. The FFM runtime is initialized
   * first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The origin weights must sum to one; the x and y direction weights must each sum to zero. The arrays are read during
   * construction and remain owned by the caller. The header does not state that the arrays must have equal lengths.</p>
   *
   * @param particles     indices of the particles defining the coordinate system; must not be null.
   * @param originWeights weight factors used to compute the origin; must not be null.
   * @param xWeights      weight factors used to compute the x direction; must not be null.
   * @param yWeights      weight factors used to compute the y direction; must not be null.
   * @param localPosition position of the site in the local coordinate system; must not be null.
   * @throws NullPointerException if an argument is null.
   */
  public LocalCoordinatesSite(IntArray particles, DoubleArray originWeights, DoubleArray xWeights,
                              DoubleArray yWeights, Vec3 localPosition) {
    super(create(particles, originWeights, xWeights, yWeights, localPosition));
  }

  /**
   * Create a local-coordinate site that depends on exactly three particles. The FFM runtime is initialized first
   * through {@link OpenMMRuntime#initialize()}.
   *
   * @param particle1     index of the first particle.
   * @param particle2     index of the second particle.
   * @param particle3     index of the third particle.
   * @param originWeights weights for the three particles when computing the origin; their components sum to one.
   * @param xWeights      weights for the three particles when computing the x direction; their components sum to zero.
   * @param yWeights      weights for the three particles when computing the y direction; their components sum to zero.
   * @param localPosition position of the site in the local coordinate system.
   * @throws NullPointerException if a {@link Vec3} argument is null.
   */
  public LocalCoordinatesSite(int particle1, int particle2, int particle3, Vec3 originWeights,
                              Vec3 xWeights, Vec3 yWeights, Vec3 localPosition) {
    super(create(particle1, particle2, particle3, originWeights, xWeights, yWeights, localPosition));
  }

  /**
   * Destroy the native virtual site.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. Do not call this
   * for a site whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_LocalCoordinatesSite_destroy);
  }

  /**
   * Get the position of the site in the local coordinate system.
   *
   * <p>The native vector is copied into the returned record.</p>
   *
   * @return copied local position.
   */
  public Vec3 getLocalPosition() {
    return Vec3.fromNative(OpenMMNative.OpenMM_LocalCoordinatesSite_getLocalPosition(getPointer()));
  }

  /**
   * Get the weight factors used to compute the origin.
   *
   * @return new array with one value per defining particle.
   */
  public double[] getOriginWeights() {
    return getWeights(OpenMMNative::OpenMM_LocalCoordinatesSite_getOriginWeights);
  }

  /**
   * Fill an existing OpenMM double array with the origin weights.
   *
   * @param weights destination array, resized by OpenMM; must not be null and stays owned by the caller.
   * @throws NullPointerException if {@code weights} is null.
   */
  public void getOriginWeights(DoubleArray weights) {
    getWeights(weights, OpenMMNative::OpenMM_LocalCoordinatesSite_getOriginWeights);
  }

  /**
   * Get the weight factors used to compute the x direction.
   *
   * @return new array with one value per defining particle.
   */
  public double[] getXWeights() {
    return getWeights(OpenMMNative::OpenMM_LocalCoordinatesSite_getXWeights);
  }

  /**
   * Fill an existing OpenMM double array with the x-direction weights.
   *
   * @param weights destination array, resized by OpenMM; must not be null and stays owned by the caller.
   * @throws NullPointerException if {@code weights} is null.
   */
  public void getXWeights(DoubleArray weights) {
    getWeights(weights, OpenMMNative::OpenMM_LocalCoordinatesSite_getXWeights);
  }

  /**
   * Get the weight factors used to compute the y direction.
   *
   * @return new array with one value per defining particle.
   */
  public double[] getYWeights() {
    return getWeights(OpenMMNative::OpenMM_LocalCoordinatesSite_getYWeights);
  }

  /**
   * Fill an existing OpenMM double array with the y-direction weights.
   *
   * @param weights destination array, resized by OpenMM; must not be null and stays owned by the caller.
   * @throws NullPointerException if {@code weights} is null.
   */
  public void getYWeights(DoubleArray weights) {
    getWeights(weights, OpenMMNative::OpenMM_LocalCoordinatesSite_getYWeights);
  }

  /**
   * Create the native virtual site after loading the FFM runtime.
   *
   * @param particles     defining particle indices.
   * @param originWeights origin weights.
   * @param xWeights      x-direction weights.
   * @param yWeights      y-direction weights.
   * @param localPosition local position.
   * @return native virtual-site handle owned by the new wrapper.
   */
  private static MemorySegment create(IntArray particles, DoubleArray originWeights,
                                      DoubleArray xWeights, DoubleArray yWeights, Vec3 localPosition) {
    Objects.requireNonNull(particles, "Particles cannot be null.");
    Objects.requireNonNull(originWeights, "Origin weights cannot be null.");
    Objects.requireNonNull(xWeights, "X weights cannot be null.");
    Objects.requireNonNull(yWeights, "Y weights cannot be null.");
    Objects.requireNonNull(localPosition, "Local position cannot be null.");
    OpenMMRuntime.initialize();
    try (var arena = java.lang.foreign.Arena.ofConfined()) {
      return OpenMMNative.OpenMM_LocalCoordinatesSite_create(
          particles.getPointer(), originWeights.getPointer(), xWeights.getPointer(),
          yWeights.getPointer(), localPosition.toNative(arena));
    }
  }

  /**
   * Create the native virtual site after loading the FFM runtime.
   *
   * @param particle1     first particle index.
   * @param particle2     second particle index.
   * @param particle3     third particle index.
   * @param originWeights origin weights.
   * @param xWeights      x-direction weights.
   * @param yWeights      y-direction weights.
   * @param localPosition local position.
   * @return native virtual-site handle owned by the new wrapper.
   */
  private static MemorySegment create(int particle1, int particle2, int particle3,
                                      Vec3 originWeights, Vec3 xWeights, Vec3 yWeights, Vec3 localPosition) {
    Objects.requireNonNull(originWeights, "Origin weights cannot be null.");
    Objects.requireNonNull(xWeights, "X weights cannot be null.");
    Objects.requireNonNull(yWeights, "Y weights cannot be null.");
    Objects.requireNonNull(localPosition, "Local position cannot be null.");
    OpenMMRuntime.initialize();
    try (var arena = java.lang.foreign.Arena.ofConfined()) {
      return OpenMMNative.OpenMM_LocalCoordinatesSite_create_2(
          particle1, particle2, particle3, originWeights.toNative(arena), xWeights.toNative(arena),
          yWeights.toNative(arena), localPosition.toNative(arena));
    }
  }

  private double[] getWeights(WeightGetter getter) {
    try (DoubleArray weights = new DoubleArray(getNumParticles())) {
      getter.get(getPointer(), weights.getPointer());
      double[] copied = new double[weights.getSize()];
      for (int i = 0; i < copied.length; i++) {
        copied[i] = weights.get(i);
      }
      return copied;
    }
  }

  private void getWeights(DoubleArray weights, WeightGetter getter) {
    Objects.requireNonNull(weights, "Weights destination cannot be null.");
    getter.get(getPointer(), weights.getPointer());
  }

  @FunctionalInterface
  private interface WeightGetter {
    void get(MemorySegment site, MemorySegment weights);
  }
}
