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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addAngle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addPerAngleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_updateParametersInContext;

/**
 * Custom Angle Force.
 */
public class CustomAngleForce extends Force {

  public CustomAngleForce(String energy) {
    pointer = OpenMM_CustomAngleForce_create(energy);
  }

  /**
   * Add a per-angle parameter to the OpenMM System.
   *
   * @param name The name of the parameter.
   */
  public void addPerAngleParameter(String name) {
    OpenMM_CustomAngleForce_addPerAngleParameter(pointer, name);
  }

  /**
   * Add an angle force to the OpenMM System.
   *
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param i3         The index of the third atom.
   * @param parameters The parameters for the angle.
   */
  public void addAngle(int i1, int i2, int i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_addAngle(pointer, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Set the parameters for one angle in the OpenMM System.
   *
   * @param index      The index of the angle.
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param i3         The index of the third atom.
   * @param parameters The angle parameters.
   */
  public void setAngleParameters(int index, int i1, int i2, int i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_setAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomAngleForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

}
