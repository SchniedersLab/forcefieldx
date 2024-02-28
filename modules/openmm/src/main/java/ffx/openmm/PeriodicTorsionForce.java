// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_updateParametersInContext;

/**
 * Periodic Torsion Force.
 */
public class PeriodicTorsionForce extends Force {

  public PeriodicTorsionForce() {
    pointer = OpenMM_PeriodicTorsionForce_create();
  }

  /**
   * Add a torsion to the PeriodicTorsionForce.
   *
   * @param particle1   Index of the first atom.
   * @param particle2   Index of the second atom.
   * @param particle3   Index of the third atom.
   * @param particle4   Index of the fourth atom.
   * @param periodicity The periodicity of the torsion.
   * @param phase       The phase of the torsion.
   * @param k           The force constant for the torsion.
   */
  public void addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    OpenMM_PeriodicTorsionForce_addTorsion(pointer, particle1, particle2, particle3, particle4, periodicity, phase, k);
  }

  /**
   * Set the parameters for a torsion.
   *
   * @param index       The index of the torsion for which to set parameters.
   * @param particle1   The index of the first atom in the torsion.
   * @param particle2   The index of the second atom in the torsion.
   * @param particle3   The index of the third atom in the torsion.
   * @param particle4   The index of the fourth atom in the torsion.
   * @param periodicity The periodicity of the torsion.
   * @param phase       The phase of the torsion.
   * @param k           The force constant for the torsion.
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(pointer, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
  }

  /**
   * Update the parameters for a torsion in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_PeriodicTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  public void destroy() {
    if (pointer != null) {
      OpenMM_PeriodicTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

}
