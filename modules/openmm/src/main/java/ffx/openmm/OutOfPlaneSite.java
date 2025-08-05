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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_OutOfPlaneSite_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_OutOfPlaneSite_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_OutOfPlaneSite_getWeight12;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_OutOfPlaneSite_getWeight13;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_OutOfPlaneSite_getWeightCross;

/**
 * This is a VirtualSite that computes the particle location based on three other particles.
 * The virtual site is placed out of the plane defined by the three particles, at a position
 * determined by the cross product of two vectors formed by the particles.
 * <p>
 * If r1, r2, and r3 are the positions of the three particles, then the virtual site location is:
 * <p>
 * r = r1 + weight12*(r2-r1) + weight13*(r3-r1) + weightCross*((r2-r1) x (r3-r1))
 * <p>
 * This is useful for representing sites like lone pairs on atoms in trigonal planar geometries,
 * where the virtual site is positioned above or below the plane of the three defining particles.
 */
public class OutOfPlaneSite extends VirtualSite {

  /**
   * Create an OutOfPlaneSite.
   *
   * @param particle1   The index of the first particle.
   * @param particle2   The index of the second particle.
   * @param particle3   The index of the third particle.
   * @param weight12    The weight for the vector from particle1 to particle2.
   * @param weight13    The weight for the vector from particle1 to particle3.
   * @param weightCross The weight for the cross product of the two vectors.
   */
  public OutOfPlaneSite(int particle1, int particle2, int particle3,
                        double weight12, double weight13, double weightCross) {
    super(OpenMM_OutOfPlaneSite_create(particle1, particle2, particle3, weight12, weight13, weightCross));
  }

  /**
   * Destroy the virtual site.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_OutOfPlaneSite_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the weight for the vector from particle1 to particle2.
   *
   * @return The weight12 parameter.
   */
  public double getWeight12() {
    return OpenMM_OutOfPlaneSite_getWeight12(pointer);
  }

  /**
   * Get the weight for the vector from particle1 to particle3.
   *
   * @return The weight13 parameter.
   */
  public double getWeight13() {
    return OpenMM_OutOfPlaneSite_getWeight13(pointer);
  }

  /**
   * Get the weight for the cross product of the two vectors.
   *
   * @return The weightCross parameter.
   */
  public double getWeightCross() {
    return OpenMM_OutOfPlaneSite_getWeightCross(pointer);
  }
}