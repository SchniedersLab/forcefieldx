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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_createExclusionsFromBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumParticlesPerSet;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumPerParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getPermutationMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_getTypeFilter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setPermutationMode;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_setTypeFilter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions;

/**
 * This class supports a wide variety of nonbonded N-particle interactions, where N is user specified.  The
 * interaction energy is determined by an arbitrary, user specified algebraic expression that is evaluated for
 * every possible set of N particles in the system.  It may depend on the positions of the individual particles,
 * the distances between pairs of particles, the angles formed by sets of three particles, and the dihedral
 * angles formed by sets of four particles.
 *
 * <p>Be aware that the cost of evaluating an N-particle interaction increases very rapidly with N.  Values larger
 * than N=3 are rarely used.
 *
 * <p>We refer to a set of particles for which the energy is being evaluated  as p1, p2, p3, etc.  The energy expression
 * may depend on the following variables and functions:
 *
 * <ul>
 * <li>x1, y1, z1, x2, y2, z2, etc.: The x, y, and z coordinates of the particle positions.  For example, x1
 * is the x coordinate of particle p1, and y3 is the y coordinate of particle p3.</li>
 * <li>distance(p1, p2): the distance between particles p1 and p2 (where "p1" and "p2" may be replaced by the names
 * of whichever particles you want to calculate the distance between).</li>
 * <li>angle(p1, p2, p3): the angle formed by the three specified particles.</li>
 * <li>dihedral(p1, p2, p3, p4): the dihedral angle formed by the four specified particles.</li>
 * <li>arbitrary global and per-particle parameters that you define.</li>
 * </ul>
 *
 * <p>To use this class, create a CustomManyParticleForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy of each set of particles.  Then call addPerParticleParameter() to define per-particle
 * parameters, and addGlobalParameter() to define global parameters.  The values of per-particle parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 *
 * <p>Next, call addParticle() once for each particle in the System to set the values of its per-particle parameters.
 * The number of particles for which you set parameters must be exactly equal to the number of particles in the
 * System, or else an exception will be thrown when you try to create a Context.  After a particle has been added,
 * you can modify its parameters by calling setParticleParameters().  This will have no effect on Contexts that already exist
 * unless you call updateParametersInContext().
 *
 * <p>Multi-particle interactions can be very expensive to evaluate, so they are usually used with a cutoff distance.  The exact
 * interpretation of the cutoff depends on the permutation mode, as discussed below.
 *
 * <p>CustomManyParticleForce also lets you specify "exclusions", particular pairs of particles whose interactions should be
 * omitted from force and energy calculations.  This is most often used for particles that are bonded to each other.
 * If you specify a pair of particles as an exclusion, <i>all</i> sets that include those two particles will be omitted.
 *
 * <p>As an example, the following code creates a CustomManyParticleForce that implements an Axilrod-Teller potential.  This
 * is an interaction between three particles that depends on all three distances and angles formed by the particles.
 *
 * <pre>{@code
 * CustomManyParticleForce force = new CustomManyParticleForce(3,
 *     "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;" +
 *     "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);" +
 *     "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
 * force.setPermutationMode(CustomManyParticleForce.SinglePermutation);
 * }</pre>
 *
 * <p>This force depends on one parameter, C.  The following code defines it as a global parameter:
 *
 * <pre>{@code
 * force.addGlobalParameter("C", 1.0);
 * }</pre>
 *
 * <p>Notice that the expression is symmetric with respect to the particles.  It only depends on the products
 * cos(theta1)*cos(theta2)*cos(theta3) and r12*r13*r23, both of which are unchanged if the labels p1, p2, and p3 are permuted.
 * This is required because we specified SinglePermutation as the permutation mode.  (This is the default, so we did not
 * really need to set it, but doing so makes the example clearer.)  In this mode, the expression is only evaluated once for
 * each set of particles.  No guarantee is made about which particle will be identified as p1, p2, etc.  Therefore, the
 * energy <i>must</i> be symmetric with respect to exchange of particles.  Otherwise, the results would be undefined because
 * permuting the labels would change the energy.
 */
public class CustomManyParticleForce extends Force {

  /**
   * Create a CustomManyParticleForce.
   *
   * @param particlesPerSet The number of particles involved in each interaction.
   * @param energy          The algebraic expression that gives the interaction energy for each set of particles.
   */
  public CustomManyParticleForce(int particlesPerSet, String energy) {
    super(OpenMM_CustomManyParticleForce_create(particlesPerSet, energy));
  }

  /**
   * Add an exclusion for a pair of particles.
   *
   * @param particle1 The index of the first particle.
   * @param particle2 The index of the second particle.
   * @return The index of the exclusion that was added.
   */
  public int addExclusion(int particle1, int particle2) {
    return OpenMM_CustomManyParticleForce_addExclusion(pointer, particle1, particle2);
  }

  /**
   * Add a new global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomManyParticleForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a new global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomManyParticleForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the Force.
   *
   * @param parameters The list of parameters for the new particle.
   * @param type       The particle type.
   * @return The index of the particle that was added.
   */
  public int addParticle(PointerByReference parameters, int type) {
    return OpenMM_CustomManyParticleForce_addParticle(pointer, parameters, type);
  }

  /**
   * Add a new per-particle parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(String name) {
    return OpenMM_CustomManyParticleForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a new per-particle parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(Pointer name) {
    return OpenMM_CustomManyParticleForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomManyParticleForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(Pointer name, PointerByReference function) {
    return OpenMM_CustomManyParticleForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Identify exclusions based on the molecular topology.
   *
   * @param bonds      The set of bonds based on which to construct exclusions.
   * @param bondCutoff Pairs of particles that are separated by this many bonds or fewer are added as exclusions.
   */
  public void createExclusionsFromBonds(PointerByReference bonds, int bondCutoff) {
    OpenMM_CustomManyParticleForce_createExclusionsFromBonds(pointer, bonds, bondCutoff);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomManyParticleForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the cutoff distance (in nm) being used for interactions.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_CustomManyParticleForce_getCutoffDistance(pointer);
  }

  /**
   * Get the algebraic expression that gives the interaction energy for each set of particles.
   *
   * @return The energy expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomManyParticleForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion for which to get particles.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntByReference particle1, IntByReference particle2) {
    OpenMM_CustomManyParticleForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion for which to get particles.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntBuffer particle1, IntBuffer particle2) {
    OpenMM_CustomManyParticleForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter for which to get the default value.
   * @return The parameter default value.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter for which to get the name.
   * @return The parameter name.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomManyParticleForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the method used for handling long range nonbonded interactions.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_CustomManyParticleForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of exclusions.
   *
   * @return The number of exclusions.
   */
  public int getNumExclusions() {
    return OpenMM_CustomManyParticleForce_getNumExclusions(pointer);
  }

  /**
   * Get the number of global parameters that the interaction depends on.
   *
   * @return The number of parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomManyParticleForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles for which force field parameters have been defined.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_CustomManyParticleForce_getNumParticles(pointer);
  }

  /**
   * Get the number of particles involved in each interaction.
   *
   * @return The number of particles per set.
   */
  public int getNumParticlesPerSet() {
    return OpenMM_CustomManyParticleForce_getNumParticlesPerSet(pointer);
  }

  /**
   * Get the number of per-particle parameters that the interaction depends on.
   *
   * @return The number of parameters.
   */
  public int getNumPerParticleParameters() {
    return OpenMM_CustomManyParticleForce_getNumPerParticleParameters(pointer);
  }

  /**
   * Get the number of tabulated functions that have been defined.
   *
   * @return The number of functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomManyParticleForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index      The index of the particle for which to get parameters.
   * @param parameters The list of parameters (output).
   * @param type       The particle type (output).
   */
  public void getParticleParameters(int index, PointerByReference parameters, IntByReference type) {
    OpenMM_CustomManyParticleForce_getParticleParameters(pointer, index, parameters, type);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index      The index of the particle for which to get parameters.
   * @param parameters The list of parameters (output).
   * @param type       The particle type (output).
   */
  public void getParticleParameters(int index, PointerByReference parameters, IntBuffer type) {
    OpenMM_CustomManyParticleForce_getParticleParameters(pointer, index, parameters, type);
  }

  /**
   * Get the name of a per-particle parameter.
   *
   * @param index The index of the parameter for which to get the name.
   * @return The parameter name.
   */
  public String getPerParticleParameterName(int index) {
    Pointer p = OpenMM_CustomManyParticleForce_getPerParticleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the permutation mode.
   *
   * @return The permutation mode.
   */
  public int getPermutationMode() {
    return OpenMM_CustomManyParticleForce_getPermutationMode(pointer);
  }

  /**
   * Get a reference to a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function to get.
   * @return The TabulatedFunction object defining the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomManyParticleForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function that may appear in the energy expression.
   *
   * @param index The index of the function to get.
   * @return The name of the function as it appears in expressions.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomManyParticleForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the type filter for the specified type.
   *
   * @param index The particle type index.
   * @param types The allowed types for interactions (output).
   */
  public void getTypeFilter(int index, PointerByReference types) {
    OpenMM_CustomManyParticleForce_getTypeFilter(pointer, index, types);
  }

  /**
   * Set the cutoff distance (in nm) being used for interactions.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_CustomManyParticleForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the algebraic expression that gives the interaction energy for each set of particles.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomManyParticleForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the algebraic expression that gives the interaction energy for each set of particles.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomManyParticleForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the particles in an exclusion.
   *
   * @param index     The index of the exclusion for which to set particles.
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   */
  public void setExclusionParticles(int index, int particle1, int particle2) {
    OpenMM_CustomManyParticleForce_setExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter for which to set the default value.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomManyParticleForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomManyParticleForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the method used for handling long range nonbonded interactions.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_CustomManyParticleForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index      The index of the particle for which to set parameters.
   * @param parameters The list of parameters for the particle.
   * @param type       The particle type.
   */
  public void setParticleParameters(int index, PointerByReference parameters, int type) {
    OpenMM_CustomManyParticleForce_setParticleParameters(pointer, index, parameters, type);
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, String name) {
    OpenMM_CustomManyParticleForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, Pointer name) {
    OpenMM_CustomManyParticleForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Set the permutation mode.
   *
   * @param mode The permutation mode.
   */
  public void setPermutationMode(int mode) {
    OpenMM_CustomManyParticleForce_setPermutationMode(pointer, mode);
  }

  /**
   * Set the type filter for the specified type.
   *
   * @param index The particle type index.
   * @param types The allowed types for interactions.
   */
  public void setTypeFilter(int index, PointerByReference types) {
    OpenMM_CustomManyParticleForce_setTypeFilter(pointer, index, types);
  }

  /**
   * Update the per-particle parameters and tabulated functions in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomManyParticleForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}