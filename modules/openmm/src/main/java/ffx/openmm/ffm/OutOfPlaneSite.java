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
 * A virtual site computed from three other particles. With {@code r1} the position of particle 1, {@code r12} the
 * vector from particle 1 to particle 2 and {@code r13} the vector from particle 1 to particle 3, the site is at
 * {@code r1 + w12*r12 + w13*r13 + wcross*(r12 x r13)}. The three weights are user-specified, which allows the site to lie
 * out of the plane of the particles. The weights {@code w12} and {@code w13} are unitless, while {@code wcross} has
 * units of inverse distance.
 *
 * <p>Ownership: this wrapper owns the native virtual site until it is passed to {@link System#setVirtualSite(int,
 * VirtualSite)}. OpenMM then assumes native ownership, and the system invalidates this wrapper when it is destroyed, so a
 * site given to a system should not also be destroyed directly. A virtual site has no {@code updateParametersInContext}
 * operation and is not a {@link Force}.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.OutOfPlaneSite} is not.</p>
 */
public class OutOfPlaneSite extends VirtualSite {

  /**
   * Create an out-of-plane virtual site. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * @param particle1   index of the first particle.
   * @param particle2   index of the second particle.
   * @param particle3   index of the third particle.
   * @param weight12    weight factor for the vector from particle 1 to particle 2; dimensionless.
   * @param weight13    weight factor for the vector from particle 1 to particle 3; dimensionless.
   * @param weightCross weight factor for the cross product of those vectors, in inverse distance units.
   */
  public OutOfPlaneSite(int particle1, int particle2, int particle3, double weight12,
                        double weight13, double weightCross) {
    super(create(particle1, particle2, particle3, weight12, weight13, weightCross));
  }

  /**
   * Destroy the native virtual site.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. Do not call this
   * for a site whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_OutOfPlaneSite_destroy);
  }

  /**
   * Get the weight factor for the vector from particle 1 to particle 2.
   *
   * @return dimensionless weight factor.
   */
  public double getWeight12() {
    return OpenMMNative.OpenMM_OutOfPlaneSite_getWeight12(getPointer());
  }

  /**
   * Get the weight factor for the vector from particle 1 to particle 3.
   *
   * @return dimensionless weight factor.
   */
  public double getWeight13() {
    return OpenMMNative.OpenMM_OutOfPlaneSite_getWeight13(getPointer());
  }

  /**
   * Get the weight factor for the cross product of the two vectors.
   *
   * @return weight factor, in inverse distance units.
   */
  public double getWeightCross() {
    return OpenMMNative.OpenMM_OutOfPlaneSite_getWeightCross(getPointer());
  }

  /**
   * Create the native virtual site after loading the FFM runtime.
   *
   * @param particle1   first particle index.
   * @param particle2   second particle index.
   * @param particle3   third particle index.
   * @param weight12    weight for the vector 1 to 2.
   * @param weight13    weight for the vector 1 to 3.
   * @param weightCross weight for the cross product.
   * @return native virtual-site handle owned by the new wrapper.
   */
  private static MemorySegment create(int particle1, int particle2, int particle3, double weight12,
                                      double weight13, double weightCross) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_OutOfPlaneSite_create(
        particle1, particle2, particle3, weight12, weight13, weightCross);
  }
}
