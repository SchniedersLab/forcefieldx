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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_updateParametersInContext;

/**
 * Harmonic Bond Force.
 */
public class HarmonicBondForce extends Force {

  public HarmonicBondForce() {
    pointer = OpenMM_HarmonicBondForce_create();
  }

  /**
   * Add a Harmonic Bond.
   *
   * @param i1     Index of the first atom.
   * @param i2     Index of the second atom.
   * @param length The equilibrium bond length.
   * @param k      The force constant.
   */
  public void addBond(int i1, int i2, double length, double k) {
    OpenMM_HarmonicBondForce_addBond(pointer, i1, i2, length, k);
  }

  /**
   * Set the bond parameters.
   *
   * @param i      The bond index.
   * @param i1     Index of the first atom.
   * @param i2     Index of the second atom.
   * @param length The equilibrium bond length.
   * @param k      The force constant.
   */
  public void setBondParameters(int i, int i1, int i2, double length, double k) {
    OpenMM_HarmonicBondForce_addBond(pointer, i1, i2, length, k);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_HarmonicBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

}
