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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ThreeParticleAverageSite_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ThreeParticleAverageSite_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ThreeParticleAverageSite_getWeight;

/**
 * This is a VirtualSite that computes the particle location as a weighted average of three
 * other particles. The virtual site is positioned at:
 * <p>
 * r = weight1*r1 + weight2*r2 + weight3*r3
 * <p>
 * where r1, r2, and r3 are the positions of the three particles, and weight1, weight2, and
 * weight3 are the corresponding weights. This is useful for representing sites like the
 * center of mass of a group of atoms or other geometrically defined positions.
 */
public class ThreeParticleAverageSite extends VirtualSite {

  /**
   * Create a ThreeParticleAverageSite.
   *
   * @param particle1 The index of the first particle.
   * @param particle2 The index of the second particle.
   * @param particle3 The index of the third particle.
   * @param weight1   The weight for the first particle.
   * @param weight2   The weight for the second particle.
   * @param weight3   The weight for the third particle.
   */
  public ThreeParticleAverageSite(int particle1, int particle2, int particle3,
                                  double weight1, double weight2, double weight3) {
    super(OpenMM_ThreeParticleAverageSite_create(particle1, particle2, particle3, weight1, weight2, weight3));
  }

  /**
   * Destroy the virtual site.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_ThreeParticleAverageSite_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the weight for one of the particles.
   *
   * @param particle The index of the particle (0, 1, or 2) for which to get the weight.
   * @return The weight for the specified particle.
   */
  public double getWeight(int particle) {
    return OpenMM_ThreeParticleAverageSite_getWeight(pointer, particle);
  }
}