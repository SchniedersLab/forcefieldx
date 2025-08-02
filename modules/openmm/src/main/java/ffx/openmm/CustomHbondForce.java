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
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addAcceptor;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addDonor;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addPerAcceptorParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addPerDonorParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getAcceptorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getDonorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumAcceptors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumDonors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumPerAcceptorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumPerDonorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getPerAcceptorParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getPerDonorParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setAcceptorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setDonorParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setPerAcceptorParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_setPerDonorParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions;

/**
 * This class implements hydrogen bond interactions between donors and acceptors.
 * It allows for custom energy expressions to define the hydrogen bond potential.
 * The force is computed between hydrogen bond donors (typically a hydrogen atom
 * bonded to an electronegative atom) and acceptors (typically lone pairs on
 * electronegative atoms).
 * <p>
 * To use this class, create a CustomHbondForce object, then call addDonor() and
 * addAcceptor() to define the donors and acceptors. The energy expression can
 * include built-in variables and custom parameters.
 */
public class CustomHbondForce extends Force {

  /**
   * Create a new CustomHbondForce.
   *
   * @param energy The energy expression for the hydrogen bond interaction.
   */
  public CustomHbondForce(String energy) {
    pointer = OpenMM_CustomHbondForce_create(energy);
  }

  /**
   * Add an acceptor to the force.
   *
   * @param a1         The index of the first atom that defines the acceptor.
   * @param a2         The index of the second atom that defines the acceptor.
   * @param a3         The index of the third atom that defines the acceptor.
   * @param parameters The parameters for the acceptor.
   * @return The index of the acceptor that was added.
   */
  public int addAcceptor(int a1, int a2, int a3, PointerByReference parameters) {
    return OpenMM_CustomHbondForce_addAcceptor(pointer, a1, a2, a3, parameters);
  }

  /**
   * Add a donor to the force.
   *
   * @param d1         The index of the first atom that defines the donor.
   * @param d2         The index of the second atom that defines the donor.
   * @param d3         The index of the third atom that defines the donor.
   * @param parameters The parameters for the donor.
   * @return The index of the donor that was added.
   */
  public int addDonor(int d1, int d2, int d3, PointerByReference parameters) {
    return OpenMM_CustomHbondForce_addDonor(pointer, d1, d2, d3, parameters);
  }

  /**
   * Add an exclusion to the force.
   *
   * @param donor    The index of the donor.
   * @param acceptor The index of the acceptor.
   * @return The index of the exclusion that was added.
   */
  public int addExclusion(int donor, int acceptor) {
    return OpenMM_CustomHbondForce_addExclusion(pointer, donor, acceptor);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @param min      The minimum value of the independent variable for which the function is defined.
   * @param max      The maximum value of the independent variable for which the function is defined.
   * @return The index of the function that was added.
   */
  public int addFunction(String name, PointerByReference function, double min, double max) {
    return OpenMM_CustomHbondForce_addFunction(pointer, name, function, min, max);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @param min      The minimum value of the independent variable for which the function is defined.
   * @param max      The maximum value of the independent variable for which the function is defined.
   * @return The index of the function that was added.
   */
  public int addFunction(Pointer name, PointerByReference function, double min, double max) {
    return OpenMM_CustomHbondForce_addFunction(pointer, name, function, min, max);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomHbondForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomHbondForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a per-acceptor parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerAcceptorParameter(String name) {
    return OpenMM_CustomHbondForce_addPerAcceptorParameter(pointer, name);
  }

  /**
   * Add a per-acceptor parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerAcceptorParameter(Pointer name) {
    return OpenMM_CustomHbondForce_addPerAcceptorParameter(pointer, name);
  }

  /**
   * Add a per-donor parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerDonorParameter(String name) {
    return OpenMM_CustomHbondForce_addPerDonorParameter(pointer, name);
  }

  /**
   * Add a per-donor parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerDonorParameter(Pointer name) {
    return OpenMM_CustomHbondForce_addPerDonorParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomHbondForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(Pointer name, PointerByReference function) {
    return OpenMM_CustomHbondForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomHbondForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for an acceptor.
   *
   * @param index      The index of the acceptor.
   * @param a1         The index of the first atom that defines the acceptor (output).
   * @param a2         The index of the second atom that defines the acceptor (output).
   * @param a3         The index of the third atom that defines the acceptor (output).
   * @param parameters The parameters for the acceptor (output).
   */
  public void getAcceptorParameters(int index, IntByReference a1, IntByReference a2, IntByReference a3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_getAcceptorParameters(pointer, index, a1, a2, a3, parameters);
  }

  /**
   * Get the parameters for an acceptor.
   *
   * @param index      The index of the acceptor.
   * @param a1         The index of the first atom that defines the acceptor (output).
   * @param a2         The index of the second atom that defines the acceptor (output).
   * @param a3         The index of the third atom that defines the acceptor (output).
   * @param parameters The parameters for the acceptor (output).
   */
  public void getAcceptorParameters(int index, IntBuffer a1, IntBuffer a2, IntBuffer a3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_getAcceptorParameters(pointer, index, a1, a2, a3, parameters);
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_CustomHbondForce_getCutoffDistance(pointer);
  }

  /**
   * Get the parameters for a donor.
   *
   * @param index      The index of the donor.
   * @param d1         The index of the first atom that defines the donor (output).
   * @param d2         The index of the second atom that defines the donor (output).
   * @param d3         The index of the third atom that defines the donor (output).
   * @param parameters The parameters for the donor (output).
   */
  public void getDonorParameters(int index, IntByReference d1, IntByReference d2, IntByReference d3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_getDonorParameters(pointer, index, d1, d2, d3, parameters);
  }

  /**
   * Get the parameters for a donor.
   *
   * @param index      The index of the donor.
   * @param d1         The index of the first atom that defines the donor (output).
   * @param d2         The index of the second atom that defines the donor (output).
   * @param d3         The index of the third atom that defines the donor (output).
   * @param parameters The parameters for the donor (output).
   */
  public void getDonorParameters(int index, IntBuffer d1, IntBuffer d2, IntBuffer d3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_getDonorParameters(pointer, index, d1, d2, d3, parameters);
  }

  /**
   * Get the energy expression for the force.
   *
   * @return The energy expression for the force.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomHbondForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index    The index of the exclusion.
   * @param donor    The index of the donor (output).
   * @param acceptor The index of the acceptor (output).
   */
  public void getExclusionParticles(int index, IntByReference donor, IntByReference acceptor) {
    OpenMM_CustomHbondForce_getExclusionParticles(pointer, index, donor, acceptor);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index    The index of the exclusion.
   * @param donor    The index of the donor (output).
   * @param acceptor The index of the acceptor (output).
   */
  public void getExclusionParticles(int index, IntBuffer donor, IntBuffer acceptor) {
    OpenMM_CustomHbondForce_getExclusionParticles(pointer, index, donor, acceptor);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index    The index of the function.
   * @param name     The name of the function as it appears in expressions (output).
   * @param function A TabulatedFunction object defining the function (output).
   * @param min      The minimum value of the independent variable for which the function is defined (output).
   * @param max      The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference function, DoubleByReference min, DoubleByReference max) {
    OpenMM_CustomHbondForce_getFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index    The index of the function.
   * @param name     The name of the function as it appears in expressions (output).
   * @param function A TabulatedFunction object defining the function (output).
   * @param min      The minimum value of the independent variable for which the function is defined (output).
   * @param max      The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference function, DoubleBuffer min, DoubleBuffer max) {
    OpenMM_CustomHbondForce_getFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomHbondForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomHbondForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_CustomHbondForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of acceptors.
   *
   * @return The number of acceptors.
   */
  public int getNumAcceptors() {
    return OpenMM_CustomHbondForce_getNumAcceptors(pointer);
  }

  /**
   * Get the number of donors.
   *
   * @return The number of donors.
   */
  public int getNumDonors() {
    return OpenMM_CustomHbondForce_getNumDonors(pointer);
  }

  /**
   * Get the number of exclusions.
   *
   * @return The number of exclusions.
   */
  public int getNumExclusions() {
    return OpenMM_CustomHbondForce_getNumExclusions(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   * @deprecated This method exists only for backward compatibility. Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomHbondForce_getNumFunctions(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomHbondForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-acceptor parameters.
   *
   * @return The number of per-acceptor parameters.
   */
  public int getNumPerAcceptorParameters() {
    return OpenMM_CustomHbondForce_getNumPerAcceptorParameters(pointer);
  }

  /**
   * Get the number of per-donor parameters.
   *
   * @return The number of per-donor parameters.
   */
  public int getNumPerDonorParameters() {
    return OpenMM_CustomHbondForce_getNumPerDonorParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomHbondForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the name of a per-acceptor parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerAcceptorParameterName(int index) {
    Pointer p = OpenMM_CustomHbondForce_getPerAcceptorParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of a per-donor parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerDonorParameterName(int index) {
    Pointer p = OpenMM_CustomHbondForce_getPerDonorParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get a reference to a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function.
   * @return The TabulatedFunction object defining the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomHbondForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function.
   * @return The name of the function as it appears in expressions.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomHbondForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for an acceptor.
   *
   * @param index      The index of the acceptor.
   * @param a1         The index of the first atom that defines the acceptor.
   * @param a2         The index of the second atom that defines the acceptor.
   * @param a3         The index of the third atom that defines the acceptor.
   * @param parameters The parameters for the acceptor.
   */
  public void setAcceptorParameters(int index, int a1, int a2, int a3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_setAcceptorParameters(pointer, index, a1, a2, a3, parameters);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_CustomHbondForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the parameters for a donor.
   *
   * @param index      The index of the donor.
   * @param d1         The index of the first atom that defines the donor.
   * @param d2         The index of the second atom that defines the donor.
   * @param d3         The index of the third atom that defines the donor.
   * @param parameters The parameters for the donor.
   */
  public void setDonorParameters(int index, int d1, int d2, int d3, PointerByReference parameters) {
    OpenMM_CustomHbondForce_setDonorParameters(pointer, index, d1, d2, d3, parameters);
  }

  /**
   * Set the energy expression for the force.
   *
   * @param energy The energy expression for the force.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomHbondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the energy expression for the force.
   *
   * @param energy The energy expression for the force.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomHbondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the particles in an exclusion.
   *
   * @param index    The index of the exclusion.
   * @param donor    The index of the donor.
   * @param acceptor The index of the acceptor.
   */
  public void setExclusionParticles(int index, int donor, int acceptor) {
    OpenMM_CustomHbondForce_setExclusionParticles(pointer, index, donor, acceptor);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index    The index of the function.
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @param min      The minimum value of the independent variable for which the function is defined.
   * @param max      The maximum value of the independent variable for which the function is defined.
   */
  public void setFunctionParameters(int index, String name, PointerByReference function, double min, double max) {
    OpenMM_CustomHbondForce_setFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index    The index of the function.
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @param min      The minimum value of the independent variable for which the function is defined.
   * @param max      The maximum value of the independent variable for which the function is defined.
   */
  public void setFunctionParameters(int index, Pointer name, PointerByReference function, double min, double max) {
    OpenMM_CustomHbondForce_setFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomHbondForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomHbondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomHbondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_CustomHbondForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the name of a per-acceptor parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerAcceptorParameterName(int index, String name) {
    OpenMM_CustomHbondForce_setPerAcceptorParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-acceptor parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerAcceptorParameterName(int index, Pointer name) {
    OpenMM_CustomHbondForce_setPerAcceptorParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-donor parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerDonorParameterName(int index, String name) {
    OpenMM_CustomHbondForce_setPerDonorParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-donor parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerDonorParameterName(int index, Pointer name) {
    OpenMM_CustomHbondForce_setPerDonorParameterName(pointer, index, name);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomHbondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}