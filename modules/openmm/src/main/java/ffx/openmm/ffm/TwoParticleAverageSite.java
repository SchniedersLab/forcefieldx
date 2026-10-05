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
 * A virtual site whose position is a weighted average of the positions of two other particles, so it lies on the
 * line through them. The weights normally sum to 1, although OpenMM does not strictly require it.
 *
 * <p>Ownership: this wrapper owns the native virtual site until it is passed to {@link System#setVirtualSite(int,
 * VirtualSite)}. OpenMM then assumes native ownership, and the system invalidates this wrapper when it is destroyed, so a
 * site given to a system should not also be destroyed directly. A virtual site has no {@code updateParametersInContext}
 * operation and is not a {@link Force}.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.TwoParticleAverageSite} is not.</p>
 */
public class TwoParticleAverageSite extends VirtualSite {

  /**
   * Create a two-particle average virtual site. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * @param particle1 index of the first particle.
   * @param particle2 index of the second particle.
   * @param weight1   weight factor (typically between 0 and 1) for the first particle; dimensionless.
   * @param weight2   weight factor (typically between 0 and 1) for the second particle; dimensionless.
   */
  public TwoParticleAverageSite(int particle1, int particle2, double weight1, double weight2) {
    super(create(particle1, particle2, weight1, weight2));
  }

  /**
   * Destroy the native virtual site.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. Do not call this
   * for a site whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_TwoParticleAverageSite_destroy);
  }

  /**
   * Get the weight factor used for a particle this virtual site depends on.
   *
   * @param particle dependency index, from 0 up to but not including {@link #getNumParticles()} (the header says
   *     "between 0 and getNumParticles()").
   * @return weight factor used for that particle.
   */
  public double getWeight(int particle) {
    return OpenMMNative.OpenMM_TwoParticleAverageSite_getWeight(getPointer(), particle);
  }

  /**
   * Create the native virtual site after loading the FFM runtime.
   *
   * @param particle1 first particle index.
   * @param particle2 second particle index.
   * @param weight1   weight for the first particle.
   * @param weight2   weight for the second particle.
   * @return native virtual-site handle owned by the new wrapper.
   */
  private static MemorySegment create(int particle1, int particle2, double weight1, double weight2) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_TwoParticleAverageSite_create(particle1, particle2, weight1, weight2);
  }
}
