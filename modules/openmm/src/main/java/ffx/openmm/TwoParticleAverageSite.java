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
package ffx.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_TwoParticleAverageSite_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_TwoParticleAverageSite_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_TwoParticleAverageSite_getWeight;

/**
 * This is a VirtualSite that computes the particle location as a weighted average
 * of two other particle's locations. This means the virtual site is on the
 * line passing through the two particles.
 */
public class TwoParticleAverageSite extends VirtualSite {

  /**
   * Create a new TwoParticleAverageSite virtual site. Normally weight1 and weight2
   * should add up to 1, although this is not strictly required.
   *
   * @param particle1 the index of the first particle
   * @param particle2 the index of the second particle
   * @param weight1   the weight factor (typically between 0 and 1) for the first particle
   * @param weight2   the weight factor (typically between 0 and 1) for the second particle
   */
  public TwoParticleAverageSite(int particle1, int particle2, double weight1, double weight2) {
    super(OpenMM_TwoParticleAverageSite_create(particle1, particle2, weight1, weight2));
  }

  /**
   * Destroy the virtual site.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_TwoParticleAverageSite_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the weight factor used for a particle this virtual site depends on.
   *
   * @param particle the particle to get (between 0 and getNumParticles())
   * @return the weight factor used for that particle
   */
  public double getWeight(int particle) {
    return OpenMM_TwoParticleAverageSite_getWeight(pointer, particle);
  }
}