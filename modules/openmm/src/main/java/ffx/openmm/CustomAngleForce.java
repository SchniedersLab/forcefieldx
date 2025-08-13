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

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addAngle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addPerAngleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumAngles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumPerAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getPerAngleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setPerAngleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions;

/**
 * This class implements interactions between sets of three particles that depend on the angle between them.
 * Unlike HarmonicAngleForce, the functional form of the interaction is completely customizable, and may
 * involve arbitrary algebraic expressions.  In addition to the angle formed by the particles, it may depend
 * on arbitrary global and per-angle parameters.
 * <p>
 * To use this class, create a CustomAngleForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each set of particles.  The expression may depend on theta, the angle
 * formed by the particles, as well as on any parameters you choose.  Then call addPerAngleParameter() to define per-angle
 * parameters, and addGlobalParameter() to define global parameters.  The values of per-angle parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context.setParameter().
 * Finally, call addAngle() once for each angle.  After an angle has been added, you can modify its parameters by calling setAngleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 * <p>
 * As an example, the following code creates a CustomAngleForce that implements a harmonic potential:
 *
 * <pre>
 *    CustomAngleForce force = new CustomAngleForce("0.5*k*(theta-theta0)^2");
 * </pre>
 * <p>
 * This force depends on two parameters: the spring constant k and equilibrium angle theta0.  The following code defines these parameters:
 *
 * <pre>
 *    force.addPerAngleParameter("k");
 *    force.addPerAngleParameter("theta0");
 * </pre>
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */
public class CustomAngleForce extends Force {

  /**
   * Create a CustomAngleForce.
   *
   * @param energy an algebraic expression giving the interaction energy between three particles as a function
   *               of theta, the angle between them
   */
  public CustomAngleForce(String energy) {
    super(OpenMM_CustomAngleForce_create(energy));
  }

  /**
   * Add an angle term to the force field.
   *
   * @param i1         the index of the first particle connected by the angle
   * @param i2         the index of the second particle connected by the angle
   * @param i3         the index of the third particle connected by the angle
   * @param parameters the list of parameters for the new angle
   * @return the index of the angle that was added
   */
  public int addAngle(int i1, int i2, int i3, DoubleArray parameters) {
    return OpenMM_CustomAngleForce_addAngle(pointer, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   * The parameter must have already been added with addGlobalParameter().
   *
   * @param name the name of the parameter
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomAngleForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a new global parameter that the interaction may depend on.  The default value provided to
   * this method is the initial value of the parameter in newly created Contexts.  You can change
   * the value at any time by calling setParameter() on the Context.
   *
   * @param name  the name of the parameter
   * @param value the default value of the parameter
   * @return the index of the parameter that was added
   */
  public int addGlobalParameter(String name, double value) {
    return OpenMM_CustomAngleForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add a new per-angle parameter that the interaction may depend on.
   *
   * @param name the name of the parameter
   * @return the index of the parameter that was added
   */
  public int addPerAngleParameter(String name) {
    return OpenMM_CustomAngleForce_addPerAngleParameter(pointer, name);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomAngleForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the force field parameters for an angle term.
   *
   * @param index      the index of the angle for which to get parameters
   * @param i1         the index of the first particle connected by the angle
   * @param i2         the index of the second particle connected by the angle
   * @param i3         the index of the third particle connected by the angle
   * @param parameters the list of parameters for the angle
   */
  public void getAngleParameters(int index, IntBuffer i1, IntBuffer i2, IntBuffer i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Get the force field parameters for an angle term.
   *
   * @param index      the index of the angle for which to get parameters
   * @param i1         the index of the first particle connected by the angle
   * @param i2         the index of the second particle connected by the angle
   * @param i3         the index of the third particle connected by the angle
   * @param parameters the list of parameters for the angle
   */
  public void getAngleParameters(int index, IntByReference i1, IntByReference i2, IntByReference i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Get the algebraic expression that gives the interaction energy for each angle
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomAngleForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of a global parameter with respect to which this Force should compute the
   * derivative of the energy.
   *
   * @param index the index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives()
   * @return the parameter name
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(pointer, index);
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
    return OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomAngleForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of angles.
   *
   * @return The number of angles.
   */
  public int getNumAngles() {
    return OpenMM_CustomAngleForce_getNumAngles(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomAngleForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-angle parameters.
   *
   * @return The number of per-angle parameters.
   */
  public int getNumPerAngleParameters() {
    return OpenMM_CustomAngleForce_getNumPerAngleParameters(pointer);
  }

  /**
   * Get the name of a per-angle parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerAngleParameterName(int index) {
    Pointer p = OpenMM_CustomAngleForce_getPerAngleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
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
   * Set the energy function expression.
   *
   * @param energy The energy function expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomAngleForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @param value The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(pointer, index, value);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomAngleForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-angle parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerAngleParameterName(int index, String name) {
    OpenMM_CustomAngleForce_setPerAngleParameterName(pointer, index, name);
  }

  /**
   * Set whether this force uses periodic boundary conditions.
   *
   * @param periodic If true, periodic boundary conditions will be used.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
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

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}