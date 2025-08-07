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
import com.sun.jna.ptr.PointerByReference;

import java.nio.DoubleBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumParticlesPerBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumPerBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions;

/**
 * This class supports a wide variety of bonded interactions.  It defines a "bond" as a single energy term
 * that depends on the positions of a fixed set of particles.  The number of particles involved in a bond, and how
 * the energy depends on their positions, is configurable.  It may depend on the positions of individual particles,
 * the distances between pairs of particles, the angles formed by sets of three particles, and the dihedral
 * angles formed by sets of four particles.
 * <p>
 * We refer to the particles in a bond as p1, p2, p3, etc.  For each bond, CustomCompoundBondForce evaluates a
 * user supplied algebraic expression to determine the interaction energy.  The expression may depend on the
 * following variables and functions:
 * <ul>
 * <li>x1, y1, z1, x2, y2, z2, etc.: The x, y, and z coordinates of the particle positions.  For example, x1
 * is the x coordinate of particle p1, and y3 is the y coordinate of particle p3.</li>
 * <li>distance(p1, p2): the distance between particles p1 and p2 (where "p1" and "p2" may be replaced by the names
 * of whichever particles you want to calculate the distance between).</li>
 * <li>angle(p1, p2, p3): the angle formed by the three specified particles.</li>
 * <li>dihedral(p1, p2, p3, p4): the dihedral angle formed by the four specified particles, guaranteed to be in the range [-pi,+pi].</li>
 * </ul>
 * <p>
 * The expression also may involve tabulated functions, and may depend on arbitrary
 * global and per-bond parameters.
 * <p>
 * To use this class, create a CustomCompoundBondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy of each bond.  Then call addPerBondParameter() to define per-bond
 * parameters and addGlobalParameter() to define global parameters.  The values of per-bond parameters are specified
 * as part of the system definition, while values of global parameters may be modified during a simulation by calling
 * Context::setParameter().
 * <p>
 * Next, call addBond() to define bonds and specify their parameter values.  After a bond has been added, you can
 * modify its parameters by calling setBondParameters().  This will have no effect on Contexts that already exist unless
 * you call updateParametersInContext().
 * <p>
 * As an example, the following code creates a CustomCompoundBondForce that implements a Urey-Bradley potential.  This
 * is an interaction between three particles that depends on the angle formed by p1-p2-p3, and on the distance between
 * p1 and p3.
 * <pre>
 *   {@code
 *    CustomCompoundBondForce* force = new CustomCompoundBondForce(3, "0.5*(kangle*(angle(p1,p2,p3)-theta0)^2+kbond*(distance(p1,p3)-r0)^2)");
 *   }
 * </pre>
 * <p>
 * This force depends on four parameters: kangle, kbond, theta0, and r0.  The following code defines these as per-bond parameters:
 * <pre>
 *   {@code
 *    force->addPerBondParameter("kangle");
 *    force->addPerBondParameter("kbond");
 *    force->addPerBondParameter("theta0");
 *    force->addPerBondParameter("r0");
 *   }
 * </pre>
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and &circ; (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 * <p>
 * This class also supports the functions pointdistance(x1, y1, z1, x2, y2, z2),
 * pointangle(x1, y1, z1, x2, y2, z2, x3, y3, z3), and pointdihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4).
 * These functions are similar to distance(), angle(), and dihedral(), but the arguments are the
 * coordinates of points to perform the calculation based on rather than the names of particles.
 * This enables more flexible geometric calculations.  For example, the following computes the distance
 * from particle p1 to the midpoint between particles p2 and p3.
 * <pre>
 *   {@code
 *    CustomCompoundBondForce* force = new CustomCompoundBondForce(3, "pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)");
 *   }
 * </pre>
 * <p>
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in the expression.
 */
public class CustomCompoundBondForce extends Force {

  /**
   * Create a CustomCompoundBondForce.
   *
   * @param numParticles The number of particles per bond.
   * @param energy       The energy expression for the force.
   */
  public CustomCompoundBondForce(int numParticles, String energy) {
    super(OpenMM_CustomCompoundBondForce_create(numParticles, energy));
  }

  /**
   * Add a bond to the force.
   *
   * @param particles  The indices of the particles in the bond.
   * @param parameters The bond parameters.
   * @return The index of the bond that was added.
   */
  public int addBond(IntArray particles, DoubleArray parameters) {
    return OpenMM_CustomCompoundBondForce_addBond(pointer, particles.getPointer(), parameters.getPointer());
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative(pointer, name);
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
    return OpenMM_CustomCompoundBondForce_addFunction(pointer, name, values, min, max);
  }

  /**
   * Add a global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomCompoundBondForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a per-bond parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerBondParameter(String name) {
    return OpenMM_CustomCompoundBondForce_addPerBondParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomCompoundBondForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomCompoundBondForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a bond.
   *
   * @param index      The index of the bond.
   * @param particles  The indices of the particles in the bond.
   * @param parameters The bond parameters.
   */
  public void getBondParameters(int index, IntArray particles, DoubleArray parameters) {
    OpenMM_CustomCompoundBondForce_getBondParameters(pointer, index, particles.getPointer(), parameters.getPointer());
  }

  /**
   * Get the energy expression for the force.
   *
   * @return The energy expression for the force.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomCompoundBondForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of a parameter with respect to which the derivative of the energy should be computed.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function.
   * @param values The tabulated values.
   * @param min    The minimum value of the independent variable.
   * @param max    The maximum value of the independent variable.
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleBuffer min, DoubleBuffer max) {
    OpenMM_CustomCompoundBondForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function.
   * @param values The tabulated values.
   * @param min    The minimum value of the independent variable.
   * @param max    The maximum value of the independent variable.
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleByReference min, DoubleByReference max) {
    OpenMM_CustomCompoundBondForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomCompoundBondForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of bonds.
   *
   * @return The number of bonds.
   */
  public int getNumBonds() {
    return OpenMM_CustomCompoundBondForce_getNumBonds(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   * @deprecated This method exists only for backward compatibility. Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomCompoundBondForce_getNumFunctions(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomCompoundBondForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles per bond.
   *
   * @return The number of particles per bond.
   */
  public int getNumParticlesPerBond() {
    return OpenMM_CustomCompoundBondForce_getNumParticlesPerBond(pointer);
  }

  /**
   * Get the number of per-bond parameters.
   *
   * @return The number of per-bond parameters.
   */
  public int getNumPerBondParameters() {
    return OpenMM_CustomCompoundBondForce_getNumPerBondParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerBondParameterName(int index) {
    Pointer p = OpenMM_CustomCompoundBondForce_getPerBondParameterName(pointer, index);
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
    return OpenMM_CustomCompoundBondForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function.
   *
   * @param index The index of the function.
   * @return The name of the function.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomCompoundBondForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for a bond.
   *
   * @param index      The index of the bond.
   * @param particles  The indices of the particles in the bond.
   * @param parameters The bond parameters.
   */
  public void setBondParameters(int index, IntArray particles, DoubleArray parameters) {
    OpenMM_CustomCompoundBondForce_setBondParameters(pointer, index, particles.getPointer(), parameters.getPointer());
  }

  /**
   * Set the energy expression for the force.
   *
   * @param energy The energy expression for the force.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomCompoundBondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function.
   * @param values The tabulated values.
   * @param min    The minimum value of the independent variable.
   * @param max    The maximum value of the independent variable.
   */
  public void setFunctionParameters(int index, String name, PointerByReference values, double min, double max) {
    OpenMM_CustomCompoundBondForce_setFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerBondParameterName(int index, String name) {
    OpenMM_CustomCompoundBondForce_setPerBondParameterName(pointer, index, name);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be used.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context to update.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}