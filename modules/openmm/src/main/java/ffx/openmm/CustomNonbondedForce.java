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
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addInteractionGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_createExclusionsFromBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getInteractionGroupParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumComputedValues;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumInteractionGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumPerParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getUseLongRangeCorrection;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_getUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setInteractionGroupParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setUseLongRangeCorrection;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions;

/**
 * This class implements nonbonded interactions between particles.  Unlike NonbondedForce, the functional form
 * of the interaction is completely customizable, and may involve arbitrary algebraic expressions and tabulated
 * functions.  It may depend on the distance between particles, as well as on arbitrary global and
 * per-particle parameters.  It also optionally supports periodic boundary conditions and cutoffs for long range interactions.
 *
 * <p>To use this class, create a CustomNonbondedForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each pair of particles.  The expression may depend on r, the distance
 * between the particles, as well as on any parameters you choose.  Then call addPerParticleParameter() to define per-particle
 * parameters, and addGlobalParameter() to define global parameters.  The values of per-particle parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 *
 * <p>Next, call addParticle() once for each particle in the System to set the values of its per-particle parameters.
 * The number of particles for which you set parameters must be exactly equal to the number of particles in the
 * System, or else an exception will be thrown when you try to create a Context.  After a particle has been added,
 * you can modify its parameters by calling setParticleParameters().  This will have no effect on Contexts that already exist
 * unless you call updateParametersInContext().
 *
 * <p>CustomNonbondedForce also lets you specify "exclusions", particular pairs of particles whose interactions should be
 * omitted from force and energy calculations.  This is most often used for particles that are bonded to each other.
 *
 * <p>As an example, the following code creates a CustomNonbondedForce that implements a 12-6 Lennard-Jones potential:
 *
 * <pre>{@code
 * CustomNonbondedForce force = new CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)");
 * }</pre>
 *
 * <p>This force depends on two parameters: sigma and epsilon.  The following code defines these as per-particle parameters:
 *
 * <pre>{@code
 * force.addPerParticleParameter("sigma");
 * force.addPerParticleParameter("epsilon");
 * }</pre>
 *
 * <p>The expression <i>must</i> be symmetric with respect to the two particles.  It typically will only be evaluated once
 * for each pair of particles, and no guarantee is made about which particle will be identified as "particle 1".  In the
 * above example, the energy only depends on the products sigma1*sigma2 and epsilon1*epsilon2, both of which are unchanged
 * if the labels 1 and 2 are reversed.  In contrast, if it depended on the difference sigma1-sigma2, the results would
 * be undefined, because reversing the labels 1 and 2 would change the energy.
 *
 * <p>The energy also may depend on "computed values".  These are similar to per-particle parameters, but instead of being
 * specified in advance, their values are computed based on global and per-particle parameters.  For example, the following
 * code uses a global parameter (lambda) to interpolate between two different sigma values for each particle (sigmaA and sigmaB).
 *
 * <pre>{@code
 * CustomNonbondedForce force = new CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)");
 * force.addComputedValue("sigma", "(1-lambda)*sigmaA + lambda*sigmaB");
 * force.addGlobalParameter("lambda", 0);
 * force.addPerParticleParameter("sigmaA");
 * force.addPerParticleParameter("sigmaB");
 * force.addPerParticleParameter("epsilon");
 * }</pre>
 *
 * <p>You could, of course, embed the computation of sigma directly into the energy expression, but then it would need to be
 * repeated for every interaction.  By separating it out as a computed value, it only needs to be computed once for each
 * particle instead of once for each interaction, thus saving computation time.
 *
 * <p>CustomNonbondedForce can operate in two modes.  By default, it computes the interaction of every particle in the System
 * with every other particle.  Alternatively, you can restrict it to only a subset of particle pairs.  To do this, specify
 * one or more "interaction groups".  An interaction group consists of two sets of particles that should interact with
 * each other.  Every particle in the first set interacts with every particle in the second set.  For example, you might use
 * this feature to compute a solute-solvent interaction energy, while omitting all interactions between two solute atoms
 * or two solvent atoms.
 */
public class CustomNonbondedForce extends Force {

  /**
   * Constructor.
   *
   * @param energy The energy expression for the force.
   */
  public CustomNonbondedForce(String energy) {
    super(OpenMM_CustomNonbondedForce_create(energy));
  }

  /**
   * Add a computed value to the force.
   *
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   * @return The index of the computed value that was added.
   */
  public int addComputedValue(String name, String expression) {
    return OpenMM_CustomNonbondedForce_addComputedValue(pointer, name, expression);
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add an exclusion to the force.
   *
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   * @return The index of the exclusion that was added.
   */
  public int addExclusion(int particle1, int particle2) {
    return OpenMM_CustomNonbondedForce_addExclusion(pointer, particle1, particle2);
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
    return OpenMM_CustomNonbondedForce_addFunction(pointer, name, function, min, max);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name  The name of the parameter.
   * @param value The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double value) {
    return OpenMM_CustomNonbondedForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add an interaction group to the force.
   *
   * @param group1 The set of particles in the first group.
   * @param group2 The set of particles in the second group.
   * @return The index of the interaction group that was added.
   */
  public int addInteractionGroup(IntSet group1, IntSet group2) {
    return OpenMM_CustomNonbondedForce_addInteractionGroup(pointer, group1.getPointer(), group2.getPointer());
  }

  /**
   * Add a particle to the force.
   *
   * @param parameters The parameters for the new particle.
   * @return The index of the particle that was added.
   */
  public int addParticle(DoubleArray parameters) {
    return OpenMM_CustomNonbondedForce_addParticle(pointer, parameters.getPointer());
  }

  /**
   * Add a per-particle parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(String name) {
    return OpenMM_CustomNonbondedForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomNonbondedForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Create exclusions based on the molecular topology.
   *
   * @param bonds      The bonds to use for determining exclusions.
   * @param bondCutoff The number of bonds within which particles should be excluded.
   */
  public void createExclusionsFromBonds(BondArray bonds, int bondCutoff) {
    OpenMM_CustomNonbondedForce_createExclusionsFromBonds(pointer, bonds.getPointer(), bondCutoff);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomNonbondedForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a computed value.
   *
   * @param index      The index of the computed value to get.
   * @param name       The name of the computed value (output).
   * @param expression The expression for computing the value (output).
   */
  public void getComputedValueParameters(int index, PointerByReference name, PointerByReference expression) {
    OpenMM_CustomNonbondedForce_getComputedValueParameters(pointer, index, name, expression);
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_CustomNonbondedForce_getCutoffDistance(pointer);
  }

  /**
   * Get the energy expression for the force.
   *
   * @return The energy expression for the force.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomNonbondedForce_getEnergyFunction(pointer);
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
    Pointer p = OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion to get.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntByReference particle1, IntByReference particle2) {
    OpenMM_CustomNonbondedForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion to get.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntBuffer particle1, IntBuffer particle2) {
    OpenMM_CustomNonbondedForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index    The index of the function to get.
   * @param name     The name of the function as it appears in expressions (output).
   * @param function A TabulatedFunction object defining the function (output).
   * @param min      The minimum value of the independent variable for which the function is defined (output).
   * @param max      The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference function, DoubleByReference min, DoubleByReference max) {
    OpenMM_CustomNonbondedForce_getFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index    The index of the function to get.
   * @param name     The name of the function as it appears in expressions (output).
   * @param function A TabulatedFunction object defining the function (output).
   * @param min      The minimum value of the independent variable for which the function is defined (output).
   * @param max      The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference function, DoubleBuffer min, DoubleBuffer max) {
    OpenMM_CustomNonbondedForce_getFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter to get.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter to get.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomNonbondedForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the parameters for an interaction group.
   *
   * @param index  The index of the interaction group to get.
   * @param group1 The set of particles in the first group (output).
   * @param group2 The set of particles in the second group (output).
   */
  public void getInteractionGroupParameters(int index, PointerByReference group1, PointerByReference group2) {
    OpenMM_CustomNonbondedForce_getInteractionGroupParameters(pointer, index, group1, group2);
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_CustomNonbondedForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of computed values.
   *
   * @return The number of computed values.
   */
  public int getNumComputedValues() {
    return OpenMM_CustomNonbondedForce_getNumComputedValues(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of exclusions.
   *
   * @return The number of exclusions.
   */
  public int getNumExclusions() {
    return OpenMM_CustomNonbondedForce_getNumExclusions(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   * @deprecated This method exists only for backward compatibility. Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomNonbondedForce_getNumFunctions(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomNonbondedForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of interaction groups.
   *
   * @return The number of interaction groups.
   */
  public int getNumInteractionGroups() {
    return OpenMM_CustomNonbondedForce_getNumInteractionGroups(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_CustomNonbondedForce_getNumParticles(pointer);
  }

  /**
   * Get the number of per-particle parameters.
   *
   * @return The number of per-particle parameters.
   */
  public int getNumPerParticleParameters() {
    return OpenMM_CustomNonbondedForce_getNumPerParticleParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomNonbondedForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index      The index of the particle to get.
   * @param parameters The parameters for the particle (output).
   */
  public void getParticleParameters(int index, PointerByReference parameters) {
    OpenMM_CustomNonbondedForce_getParticleParameters(pointer, index, parameters);
  }

  /**
   * Get the name of a per-particle parameter.
   *
   * @param index The index of the parameter to get.
   * @return The name of the parameter.
   */
  public String getPerParticleParameterName(int index) {
    Pointer p = OpenMM_CustomNonbondedForce_getPerParticleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the switching distance.
   *
   * @return The switching distance, measured in nm.
   */
  public double getSwitchingDistance() {
    return OpenMM_CustomNonbondedForce_getSwitchingDistance(pointer);
  }

  /**
   * Get a reference to a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function to get.
   * @return The TabulatedFunction object defining the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomNonbondedForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function to get.
   * @return The name of the function as it appears in expressions.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomNonbondedForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get whether to use the long range correction.
   *
   * @return 1 if the long range correction is used, 0 otherwise.
   */
  public int getUseLongRangeCorrection() {
    return OpenMM_CustomNonbondedForce_getUseLongRangeCorrection(pointer);
  }

  /**
   * Get whether to use a switching function.
   *
   * @return 1 if a switching function is used, 0 otherwise.
   */
  public int getUseSwitchingFunction() {
    return OpenMM_CustomNonbondedForce_getUseSwitchingFunction(pointer);
  }

  /**
   * Set the parameters for a computed value.
   *
   * @param index      The index of the computed value to set.
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   */
  public void setComputedValueParameters(int index, String name, String expression) {
    OpenMM_CustomNonbondedForce_setComputedValueParameters(pointer, index, name, expression);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_CustomNonbondedForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the energy expression for the force.
   *
   * @param expression The energy expression for the force.
   */
  public void setEnergyFunction(String expression) {
    OpenMM_CustomNonbondedForce_setEnergyFunction(pointer, expression);
  }

  /**
   * Set the particles in an exclusion.
   *
   * @param index     The index of the exclusion to set.
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   */
  public void setExclusionParticles(int index, int particle1, int particle2) {
    OpenMM_CustomNonbondedForce_setExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index    The index of the function to set.
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @param min      The minimum value of the independent variable for which the function is defined.
   * @param max      The maximum value of the independent variable for which the function is defined.
   */
  public void setFunctionParameters(int index, String name, PointerByReference function, double min, double max) {
    OpenMM_CustomNonbondedForce_setFunctionParameters(pointer, index, name, function, min, max);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index The index of the parameter to set.
   * @param value The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue(pointer, index, value);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter to set.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomNonbondedForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the parameters for an interaction group.
   *
   * @param index  The index of the interaction group to set.
   * @param group1 The set of particles in the first group.
   * @param group2 The set of particles in the second group.
   */
  public void setInteractionGroupParameters(int index, PointerByReference group1, PointerByReference group2) {
    OpenMM_CustomNonbondedForce_setInteractionGroupParameters(pointer, index, group1, group2);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_CustomNonbondedForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index      The index of the particle to set.
   * @param parameters The parameters for the particle.
   */
  public void setParticleParameters(int index, PointerByReference parameters) {
    OpenMM_CustomNonbondedForce_setParticleParameters(pointer, index, parameters);
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter to set.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, String name) {
    OpenMM_CustomNonbondedForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Set the switching distance.
   *
   * @param distance The switching distance, measured in nm.
   */
  public void setSwitchingDistance(double distance) {
    OpenMM_CustomNonbondedForce_setSwitchingDistance(pointer, distance);
  }

  /**
   * Set whether to use the long range correction.
   *
   * @param use 1 to use the long range correction, 0 otherwise.
   */
  public void setUseLongRangeCorrection(int use) {
    OpenMM_CustomNonbondedForce_setUseLongRangeCorrection(pointer, use);
  }

  /**
   * Set whether to use a switching function.
   *
   * @param use 1 to use a switching function, 0 otherwise.
   */
  public void setUseSwitchingFunction(int use) {
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(pointer, use);
  }

  /**
   * Update the per-particle parameters and tabulated functions in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomNonbondedForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}