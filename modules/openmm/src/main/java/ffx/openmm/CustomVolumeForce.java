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

import com.sun.jna.Pointer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_usesPeriodicBoundaryConditions;

/**
 * This class computes an energy that depends only on the volume of the periodic box, or more generally
 * on the box shape as specified by the elements of the box vectors.  Because the energy does not
 * depend on particle positions, it does not apply any forces to particles.  It is primarily useful
 * for constant pressure simulations, where the volume-dependent energy can influence the behavior
 * of the barostat.  Energy terms of this sort are often used for pressure matching in coarse grained
 * force fields.
 *
 * <p>To use this class, create a CustomVolumeForce object, passing an algebraic expression to the constructor
 * that defines the energy.  The expression may depend on the following variables.
 *
 * <ul>
 * <li>v: The volume of the periodic box in nm^3.</li>
 * <li>ax: The x component of the first box vector in nm.  (The y and z components are always zero.)</li>
 * <li>bx, by: The x and y components of the second box vector in nm.  (The z component is always zero.)</li>
 * <li>cx, cy, cz: The x, y and z components of the third box vector in nm.</li>
 * <li>Global parameters that you define by calling addGlobalParameter().</li>
 * </ul>
 *
 * <p>The initial value of a global parameter is specified in the call to addGlobalParameter().  Their values
 * can be modified during a simulation by calling Context::setParameter().
 *
 * <p>Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x &lt; 0, 1 otherwise.  delta(x) = 1 if x = 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */
public class CustomVolumeForce extends Force {

  /**
   * Create a CustomVolumeForce.
   *
   * @param energy an algebraic expression giving the energy as a function of the box shape
   */
  public CustomVolumeForce(String energy) {
    super(OpenMM_CustomVolumeForce_create(energy));
  }

  /**
   * Add a new global parameter that the interaction may depend on.  The default value provided to
   * this method is the initial value of the parameter in newly created Contexts.  You can change
   * the value at any time by calling setParameter() on the Context.
   *
   * @param name         the name of the parameter
   * @param defaultValue the default value of the parameter
   * @return the index of the parameter that was added
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomVolumeForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a new global parameter that the energy may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomVolumeForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomVolumeForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the algebraic expression that defines the energy.
   *
   * @return the energy expression
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomVolumeForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index the index of the parameter for which to get the default value
   * @return the parameter default value
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomVolumeForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index the index of the parameter for which to get the name
   * @return the parameter name
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomVolumeForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of global parameters that the energy depends on.
   *
   * @return The number of parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomVolumeForce_getNumGlobalParameters(pointer);
  }

  /**
   * Set the algebraic expression that defines the energy.
   *
   * @param energy the energy expression
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomVolumeForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the algebraic expression that defines the energy.
   *
   * @param energy the energy expression
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomVolumeForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        the index of the parameter for which to set the default value
   * @param defaultValue the default value of the parameter
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomVolumeForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index the index of the parameter for which to set the name
   * @param name  the name of the parameter
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomVolumeForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index the index of the parameter for which to set the name
   * @param name  the name of the parameter
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomVolumeForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Returns whether or not this force makes use of periodic boundary conditions.  Because this
   * class is only applicable to periodic systems, this always returns true.
   *
   * @return true if the force uses periodic boundary conditions
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomVolumeForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}