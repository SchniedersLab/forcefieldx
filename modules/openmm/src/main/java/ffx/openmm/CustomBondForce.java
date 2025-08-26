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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumPerBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_usesPeriodicBoundaryConditions;

/**
 * This class implements bonded interactions between pairs of particles. Unlike HarmonicBondForce, the functional form
 * of the interaction is completely customizable, and may involve arbitrary algebraic expressions.
 * It may depend on the distance between particles, as well as on arbitrary global and
 * per-bond parameters.
 * <p>
 * To use this class, create a CustomBondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each pair of bonded particles. The expression may depend on r, the distance
 * between the particles, as well as on any parameters you choose. Then call addPerBondParameter() to define per-bond
 * parameters, and addGlobalParameter() to define global parameters. The values of per-bond parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 * Finally, call addBond() once for each bond. After a bond has been added, you can modify its parameters by calling setBondParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 * <p>
 * As an example, the following code creates a CustomBondForce that implements a harmonic potential:
 * <p>
 * &lt;pre&gt;
 * {@code
 * CustomBondForce force = new CustomBondForce("0.5*k*(r-r0)^2");
 * }
 * &lt;/pre&gt;
 * <p>
 * This force depends on two parameters: the spring constant k and equilibrium distance r0. The following code defines these parameters:
 * <p>
 * &lt;pre&gt;
 * {@code
 * force.addPerBondParameter("k");
 * force.addPerBondParameter("r0");
 * }
 * &lt;/pre&gt;
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed. You can then query its value in a Context by calling getState() on it.
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select. All trigonometric functions
 * are defined in radians, and log is the natural logarithm. step(x) = 0 if x is less than 0, 1 otherwise. delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */
public class CustomBondForce extends Force {

  /**
   * Create a CustomBondForce.
   *
   * @param energy an algebraic expression giving the interaction energy between two bonded particles as a function
   *               of r, the distance between them
   */
  public CustomBondForce(String energy) {
    super(OpenMM_CustomBondForce_create(energy));
  }

  /**
   * Add a bond term to the force field.
   *
   * @param i1         the index of the first particle connected by the bond
   * @param i2         the index of the second particle connected by the bond
   * @param parameters the list of parameters for the new bond
   * @return the index of the bond that was added
   */
  public int addBond(int i1, int i2, DoubleArray parameters) {
    return OpenMM_CustomBondForce_addBond(pointer, i1, i2, parameters.getPointer());
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   * The parameter must have already been added with addGlobalParameter().
   *
   * @param name the name of the parameter
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomBondForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a new global parameter that the interaction may depend on. The default value provided to
   * this method is the initial value of the parameter in newly created Contexts. You can change
   * the value at any time by calling setParameter() on the Context.
   *
   * @param name  the name of the parameter
   * @param value the default value of the parameter
   * @return the index of the parameter that was added
   */
  public int addGlobalParameter(String name, double value) {
    return OpenMM_CustomBondForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add a new per-bond parameter that the interaction may depend on.
   *
   * @param name the name of the parameter
   * @return the index of the parameter that was added
   */
  public int addPerBondParameter(String name) {
    return OpenMM_CustomBondForce_addPerBondParameter(pointer, name);
  }

  /**
   * Destroy the OpenMM CustomBondForce.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomBondForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the force field parameters for a bond term.
   *
   * @param index      the index of the bond for which to get parameters
   * @param i1         the index of the first particle connected by the bond
   * @param i2         the index of the second particle connected by the bond
   * @param parameters the list of parameters for the bond
   */
  public void getBondParameters(int index, IntBuffer i1, IntBuffer i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_getBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Get the force field parameters for a bond term.
   *
   * @param index      the index of the bond for which to get parameters
   * @param i1         the index of the first particle connected by the bond
   * @param i2         the index of the second particle connected by the bond
   * @param parameters the list of parameters for the bond
   */
  public void getBondParameters(int index, IntByReference i1, IntByReference i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_getBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Get the algebraic expression that gives the interaction energy for each bond
   *
   * @return the energy function expression
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomBondForce_getEnergyFunction(pointer);
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
    Pointer p = OpenMM_CustomBondForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomBondForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomBondForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of bonds for which force field parameters have been defined.
   *
   * @return the number of bonds
   */
  public int getNumBonds() {
    return OpenMM_CustomBondForce_getNumBonds(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomBondForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-bond parameters.
   *
   * @return The number of per-bond parameters.
   */
  public int getNumPerBondParameters() {
    return OpenMM_CustomBondForce_getNumPerBondParameters(pointer);
  }

  /**
   * Get the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerBondParameterName(int index) {
    Pointer p = OpenMM_CustomBondForce_getPerBondParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for one bond in the OpenMM System.
   *
   * @param index      The index of the bond.
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param parameters The bond parameters.
   */
  public void setBondParameters(int index, int i1, int i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_setBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Set the energy function expression.
   *
   * @param energy The energy function expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomBondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @param value The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMM_CustomBondForce_setGlobalParameterDefaultValue(pointer, index, value);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomBondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerBondParameterName(int index, String name) {
    OpenMM_CustomBondForce_setPerBondParameterName(pointer, index, name);
  }

  /**
   * Set whether this force uses periodic boundary conditions.
   *
   * @param periodic If true, periodic boundary conditions will be used.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}