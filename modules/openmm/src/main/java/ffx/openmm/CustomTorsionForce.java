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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_addPerTorsionParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getNumPerTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getNumTorsions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getPerTorsionParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_getTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setPerTorsionParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions;

/**
 * This class implements torsions whose interaction energy is defined by a user-supplied
 * algebraic expression. The expression may depend on the torsion angle, as well as on
 * global and per-torsion parameters.
 * <p>
 * To use this class, create a CustomTorsionForce object, passing an algebraic expression to the
 * constructor that defines the interaction energy of each torsion. Then call addPerTorsionParameter()
 * to define per-torsion parameters, addGlobalParameter() to define global parameters, and addTorsion()
 * to define torsions and specify their parameter values.
 */
public class CustomTorsionForce extends Force {

  /**
   * Create a CustomTorsionForce.
   *
   * @param energy The algebraic expression that gives the interaction energy of each torsion as a function of theta, the torsion angle.
   */
  public CustomTorsionForce(String energy) {
    pointer = OpenMM_CustomTorsionForce_create(energy);
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   *
   * @param name The name of the parameter.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomTorsionForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   *
   * @param name The name of the parameter.
   */
  public void addEnergyParameterDerivative(Pointer name) {
    OpenMM_CustomTorsionForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a new global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomTorsionForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a new global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomTorsionForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a new per-torsion parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerTorsionParameter(String name) {
    return OpenMM_CustomTorsionForce_addPerTorsionParameter(pointer, name);
  }

  /**
   * Add a new per-torsion parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerTorsionParameter(Pointer name) {
    return OpenMM_CustomTorsionForce_addPerTorsionParameter(pointer, name);
  }

  /**
   * Add a torsion to the Force.
   *
   * @param particle1  The index of the first particle forming the torsion.
   * @param particle2  The index of the second particle forming the torsion.
   * @param particle3  The index of the third particle forming the torsion.
   * @param particle4  The index of the fourth particle forming the torsion.
   * @param parameters The list of parameters for the new torsion.
   * @return The index of the torsion that was added.
   */
  public int addTorsion(int particle1, int particle2, int particle3, int particle4, PointerByReference parameters) {
    return OpenMM_CustomTorsionForce_addTorsion(pointer, particle1, particle2, particle3, particle4, parameters);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the algebraic expression that gives the interaction energy of each torsion.
   *
   * @return The energy expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomTorsionForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of a parameter with respect to which the derivative of the energy should be computed.
   *
   * @param index The index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives().
   * @return The parameter name.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter for which to get the default value.
   * @return The parameter default value.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter for which to get the name.
   * @return The parameter name.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomTorsionForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters that the interaction depends on.
   *
   * @return The number of parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomTorsionForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-torsion parameters that the interaction depends on.
   *
   * @return The number of parameters.
   */
  public int getNumPerTorsionParameters() {
    return OpenMM_CustomTorsionForce_getNumPerTorsionParameters(pointer);
  }

  /**
   * Get the number of torsions for which force field parameters have been defined.
   *
   * @return The number of torsions.
   */
  public int getNumTorsions() {
    return OpenMM_CustomTorsionForce_getNumTorsions(pointer);
  }

  /**
   * Get the name of a per-torsion parameter.
   *
   * @param index The index of the parameter for which to get the name.
   * @return The parameter name.
   */
  public String getPerTorsionParameterName(int index) {
    Pointer p = OpenMM_CustomTorsionForce_getPerTorsionParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the force field parameters for a torsion.
   *
   * @param index      The index of the torsion for which to get parameters.
   * @param particle1  The index of the first particle forming the torsion (output).
   * @param particle2  The index of the second particle forming the torsion (output).
   * @param particle3  The index of the third particle forming the torsion (output).
   * @param particle4  The index of the fourth particle forming the torsion (output).
   * @param parameters The list of parameters (output).
   */
  public void getTorsionParameters(int index, IntByReference particle1, IntByReference particle2,
                                   IntByReference particle3, IntByReference particle4, PointerByReference parameters) {
    OpenMM_CustomTorsionForce_getTorsionParameters(pointer, index, particle1, particle2, particle3, particle4, parameters);
  }

  /**
   * Get the force field parameters for a torsion.
   *
   * @param index      The index of the torsion for which to get parameters.
   * @param particle1  The index of the first particle forming the torsion (output).
   * @param particle2  The index of the second particle forming the torsion (output).
   * @param particle3  The index of the third particle forming the torsion (output).
   * @param particle4  The index of the fourth particle forming the torsion (output).
   * @param parameters The list of parameters (output).
   */
  public void getTorsionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                   IntBuffer particle3, IntBuffer particle4, PointerByReference parameters) {
    OpenMM_CustomTorsionForce_getTorsionParameters(pointer, index, particle1, particle2, particle3, particle4, parameters);
  }

  /**
   * Set the algebraic expression that gives the interaction energy of each torsion.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomTorsionForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the algebraic expression that gives the interaction energy of each torsion.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomTorsionForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter for which to set the default value.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomTorsionForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomTorsionForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-torsion parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setPerTorsionParameterName(int index, String name) {
    OpenMM_CustomTorsionForce_setPerTorsionParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-torsion parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setPerTorsionParameterName(int index, Pointer name) {
    OpenMM_CustomTorsionForce_setPerTorsionParameterName(pointer, index, name);
  }

  /**
   * Set the force field parameters for a torsion.
   *
   * @param index      The index of the torsion for which to set parameters.
   * @param particle1  The index of the first particle forming the torsion.
   * @param particle2  The index of the second particle forming the torsion.
   * @param particle3  The index of the third particle forming the torsion.
   * @param particle4  The index of the fourth particle forming the torsion.
   * @param parameters The list of parameters for the torsion.
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, PointerByReference parameters) {
    OpenMM_CustomTorsionForce_setTorsionParameters(pointer, index, particle1, particle2, particle3, particle4, parameters);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the per-torsion parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}