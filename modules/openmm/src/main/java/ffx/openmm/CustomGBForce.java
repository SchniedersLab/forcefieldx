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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyTerm;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getEnergyTermParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumComputedValues;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumEnergyTerms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumPerParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setEnergyTermParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_usesPeriodicBoundaryConditions;

/**
 * This class implements complex, multiple stage nonbonded interactions between particles.  It is designed primarily
 * for implementing Generalized Born implicit solvation models, although it is not strictly limited to that purpose.
 * The interaction is specified as a series of computations, each defined by an arbitrary algebraic expression.
 * It also allows tabulated functions to be defined and used with the computations.  It optionally supports periodic boundary
 * conditions and cutoffs for long range interactions.
 * <p>
 * The computation consists of calculating some number of per-particle <b>computed values</b>, followed by one or more
 * <b>energy terms</b>.  A computed value is a scalar value that is computed for each particle in the system.  It may
 * depend on an arbitrary set of global and per-particle parameters, and well as on other computed values that have
 * been calculated before it.  Once all computed values have been calculated, the energy terms and their derivatives
 * are evaluated to determine the system energy and particle forces.  The energy terms may depend on global parameters,
 * per-particle parameters, and per-particle computed values.
 * <p>
 * When specifying a computed value or energy term, you provide an algebraic expression to evaluate and a <b>computation type</b>
 * describing how the expression is to be evaluated.  There are two main types of computations:
 * <ul>
 * <li><b>Single Particle</b>: The expression is evaluated once for each particle in the System.  In the case of a computed
 * value, this means the value for a particle depends only on other properties of that particle (its position, parameters, and other
 * computed values).  In the case of an energy term, it means each particle makes an independent contribution to the System
 * energy.</li>
 * <li><b>Particle Pairs</b>: The expression is evaluated for every pair of particles in the system.  In the case of a computed
 * value, the value for a particular particle is calculated by pairing it with every other particle in the system, evaluating
 * the expression for each pair, and summing them.  For an energy term, each particle pair makes an independent contribution to
 * the System energy.  (Note that energy terms are assumed to be symmetric with respect to the two interacting particles, and
 * therefore are evaluated only once per pair.  In contrast, expressions for computed values need not be symmetric and therefore are calculated
 * twice for each pair: once when calculating the value for the first particle, and again when calculating the value for the
 * second particle.)</li>
 * </ul>
 * <p>
 * Be aware that, although this class is extremely general in the computations it can define, particular Platforms may only support
 * more restricted types of computations.  In particular, all currently existing Platforms require that the first computed value
 * <i>must</i> be a particle pair computation, and all computed values after the first <i>must</i> be single particle computations.
 * This is sufficient for most Generalized Born models, but might not permit some other types of calculations to be implemented.
 * <p>
 * This is a complicated class to use, and an example may help to clarify it.  The following code implements the OBC variant
 * of the GB/SA solvation model, using the ACE approximation to estimate surface area:
 * <pre>
 *   {@code
 *    CustomGBForce* custom = new CustomGBForce();
 *    custom->addPerParticleParameter("q");
 *    custom->addPerParticleParameter("radius");
 *    custom->addPerParticleParameter("scale");
 *    custom->addGlobalParameter("solventDielectric", obc->getSolventDielectric());
 *    custom->addGlobalParameter("soluteDielectric", obc->getSoluteDielectric());
 *    custom->addComputedValue("I", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
 *                                  "U=r+sr2;"
 *                                  "C=2*(1/or1-1/L)*step(sr2-r-or1);"
 *                                  "L=max(or1, D);"
 *                                  "D=abs(r-sr2);"
 *                                  "sr2 = scale2*or2;"
 *                                  "or1 = radius1-0.009; or2 = radius2-0.009", CustomGBForce::ParticlePairNoExclusions);
 *    custom->addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
 *                                  "psi=I*or; or=radius-0.009", CustomGBForce::SingleParticle);
 *    custom->addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
 *                          CustomGBForce::SingleParticle);
 *    custom->addEnergyTerm("-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
 *                          "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePair);
 *   }
 * </pre>
 * <p>
 * It begins by defining three per-particle parameters (charge, atomic radius, and scale factor) and two global parameters
 * (the dielectric constants for the solute and solvent).  It then defines a computed value "I" of type ParticlePair.  The
 * expression for evaluating it is a complicated function of the distance between each pair of particles (r), their atomic
 * radii (radius1 and radius2), and their scale factors (scale1 and scale2).  Very roughly speaking, it is a measure of the
 * distance between each particle and other nearby particles.
 * <p>
 * Next a computation is defined for the Born Radius (B).  It is computed independently for each particle, and is a function of
 * that particle's atomic radius and the intermediate value I defined above.
 * <p>
 * Finally, two energy terms are defined.  The first one is computed for each particle and represents the surface area term,
 * as well as the self interaction part of the polarization energy.  The second term is calculated for each pair of particles,
 * and represents the screening of electrostatic interactions by the solvent.
 * <p>
 * After defining the force as shown above, you should then call addParticle() once for each particle in the System to set the
 * values of its per-particle parameters (q, radius, and scale).  The number of particles for which you set parameters must be
 * exactly equal to the number of particles in the System, or else an exception will be thrown when you try to create a Context.
 * After a particle has been added, you can modify its parameters by calling setParticleParameters().  This will have no effect
 * on Contexts that already exist unless you call updateParametersInContext().
 * <p>
 * CustomGBForce also lets you specify "exclusions", particular pairs of particles whose interactions should be
 * omitted from calculations.  This is most often used for particles that are bonded to each other.  Even if you specify exclusions,
 * however, you can use the computation type ParticlePairNoExclusions to indicate that exclusions should not be applied to a
 * particular piece of the computation.
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and &circ; (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.  In expressions for particle pair calculations, the names of per-particle parameters and computed values
 * have the suffix "1" or "2" appended to them to indicate the values for the two interacting particles.  As seen in the above example,
 * an expression may also involve intermediate quantities that are defined following the main expression, using ";" as a separator.
 * <p>
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in expressions.
 */
public class CustomGBForce extends Force {

  /**
   * Create a CustomGBForce.
   */
  public CustomGBForce() {
    super(OpenMM_CustomGBForce_create());
  }

  /**
   * Add a computed value to the force.
   *
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   * @param type       The type of computation to perform.
   * @return The index of the computed value that was added.
   */
  public int addComputedValue(String name, String expression, int type) {
    return OpenMM_CustomGBForce_addComputedValue(pointer, name, expression, type);
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomGBForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add an energy term to the force.
   *
   * @param expression The expression for computing the energy term.
   * @param type       The type of computation to perform.
   * @return The index of the energy term that was added.
   */
  public int addEnergyTerm(String expression, int type) {
    return OpenMM_CustomGBForce_addEnergyTerm(pointer, expression, type);
  }

  /**
   * Add an exclusion to the force.
   *
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   * @return The index of the exclusion that was added.
   */
  public int addExclusion(int particle1, int particle2) {
    return OpenMM_CustomGBForce_addExclusion(pointer, particle1, particle2);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name   The name of the function as it appears in expressions.
   * @param values The tabulated values of the function.
   * @param min    The minimum value of the independent variable for which the function is defined.
   * @param max    The maximum value of the independent variable for which the function is defined.
   * @return The index of the function that was added.
   * @deprecated This method exists only for backward compatibility. Use addTabulatedFunction() instead.
   */
  @Deprecated
  public int addFunction(String name, PointerByReference values, double min, double max) {
    return OpenMM_CustomGBForce_addFunction(pointer, name, values, min, max);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomGBForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the force.
   *
   * @param particleParameters The parameters for the particle.
   * @return The index of the particle that was added.
   */
  public int addParticle(DoubleArray particleParameters) {
    return OpenMM_CustomGBForce_addParticle(pointer, particleParameters.getPointer());
  }

  /**
   * Add a per-particle parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(String name) {
    return OpenMM_CustomGBForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomGBForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomGBForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value (output).
   * @param expression The expression for computing the value (output).
   * @param type       The type of computation to perform (output).
   */
  public void getComputedValueParameters(int index, PointerByReference name, PointerByReference expression, IntBuffer type) {
    OpenMM_CustomGBForce_getComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Get the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value (output).
   * @param expression The expression for computing the value (output).
   * @param type       The type of computation to perform (output).
   */
  public void getComputedValueParameters(int index, PointerByReference name, PointerByReference expression, IntByReference type) {
    OpenMM_CustomGBForce_getComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance.
   */
  public double getCutoffDistance() {
    return OpenMM_CustomGBForce_getCutoffDistance(pointer);
  }

  /**
   * Get the name of a parameter with respect to which the derivative of the energy should be computed.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomGBForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term (output).
   * @param type       The type of computation to perform (output).
   */
  public void getEnergyTermParameters(int index, PointerByReference expression, IntBuffer type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Get the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term (output).
   * @param type       The type of computation to perform (output).
   */
  public void getEnergyTermParameters(int index, PointerByReference expression, IntByReference type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntBuffer particle1, IntBuffer particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntByReference particle1, IntByReference particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function (output).
   * @param values The tabulated values of the function (output).
   * @param min    The minimum value of the independent variable for which the function is defined (output).
   * @param max    The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleBuffer min, DoubleBuffer max) {
    OpenMM_CustomGBForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function (output).
   * @param values The tabulated values of the function (output).
   * @param min    The minimum value of the independent variable for which the function is defined (output).
   * @param max    The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleByReference min, DoubleByReference max) {
    OpenMM_CustomGBForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomGBForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomGBForce_getGlobalParameterName(pointer, index);
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
    return OpenMM_CustomGBForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of computed values.
   *
   * @return The number of computed values.
   */
  public int getNumComputedValues() {
    return OpenMM_CustomGBForce_getNumComputedValues(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of energy terms.
   *
   * @return The number of energy terms.
   */
  public int getNumEnergyTerms() {
    return OpenMM_CustomGBForce_getNumEnergyTerms(pointer);
  }

  /**
   * Get the number of exclusions.
   *
   * @return The number of exclusions.
   */
  public int getNumExclusions() {
    return OpenMM_CustomGBForce_getNumExclusions(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   * @deprecated This method exists only for backward compatibility. Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomGBForce_getNumFunctions(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomGBForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_CustomGBForce_getNumParticles(pointer);
  }

  /**
   * Get the number of per-particle parameters.
   *
   * @return The number of per-particle parameters.
   */
  public int getNumPerParticleParameters() {
    return OpenMM_CustomGBForce_getNumPerParticleParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomGBForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particleParameters The parameters for the particle (output).
   */
  public void getParticleParameters(int index, DoubleArray particleParameters) {
    OpenMM_CustomGBForce_getParticleParameters(pointer, index, particleParameters.getPointer());
  }

  /**
   * Get the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerParticleParameterName(int index) {
    Pointer p = OpenMM_CustomGBForce_getPerParticleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get a reference to a tabulated function.
   *
   * @param index The index of the function.
   * @return A reference to the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomGBForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function.
   *
   * @param index The index of the function.
   * @return The name of the function.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomGBForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   * @param type       The type of computation to perform.
   */
  public void setComputedValueParameters(int index, String name, String expression, int type) {
    OpenMM_CustomGBForce_setComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_CustomGBForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term.
   * @param type       The type of computation to perform.
   */
  public void setEnergyTermParameters(int index, String expression, int type) {
    OpenMM_CustomGBForce_setEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Set the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   */
  public void setExclusionParticles(int index, int particle1, int particle2) {
    OpenMM_CustomGBForce_setExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function.
   * @param values The tabulated values of the function.
   * @param min    The minimum value of the independent variable for which the function is defined.
   * @param max    The maximum value of the independent variable for which the function is defined.
   */
  public void setFunctionParameters(int index, String name, PointerByReference values, double min, double max) {
    OpenMM_CustomGBForce_setFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomGBForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomGBForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_CustomGBForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particleParameters The parameters for the particle.
   */
  public void setParticleParameters(int index, DoubleArray particleParameters) {
    OpenMM_CustomGBForce_setParticleParameters(pointer, index, particleParameters.getPointer());
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, String name) {
    OpenMM_CustomGBForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomGBForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}