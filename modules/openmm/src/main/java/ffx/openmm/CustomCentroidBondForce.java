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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getGroupParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumGroupsPerBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumPerBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_getPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setGroupParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions;

/**
 * This class is similar to CustomCompoundBondForce, but instead of applying forces between individual particles,
 * it applies them between the centers of groups of particles.  This is useful for a variety of purposes, such as
 * restraints to keep two molecules from moving too far apart.
 * <p>
 * When using this class, you define groups of particles, and the center of each group is calculated as a weighted
 * average of the particle positions.  By default, the particle masses are used as weights, so the center position
 * is the center of mass.  You can optionally specify different weights to use.  You then add bonds just as with
 * CustomCompoundBondForce, but instead of specifying the particles that make up a bond, you specify the groups.
 * <p>
 * When creating a CustomCentroidBondForce, you specify the number of groups involved in a bond, and an expression
 * for the energy of each bond.  It may depend on the center positions of individual groups, the distances between
 * the centers of pairs of groups, the angles formed by sets of three groups, and the dihedral angles formed by
 * sets of four groups.
 * <p>
 * We refer to the groups in a bond as g1, g2, g3, etc.  For each bond, CustomCentroidBondForce evaluates a
 * user supplied algebraic expression to determine the interaction energy.  The expression may depend on the
 * following variables and functions:
 *
 * <ul>
 * <li>x1, y1, z1, x2, y2, z2, etc.: The x, y, and z coordinates of the centers of the groups.  For example, x1
 * is the x coordinate of the center of group g1, and y3 is the y coordinate of the center of group g3.</li>
 * <li>distance(g1, g2): the distance between the centers of groups g1 and g2 (where "g1" and "g2" may be replaced
 * by the names of whichever groups you want to calculate the distance between).</li>
 * <li>angle(g1, g2, g3): the angle formed by the centers of the three specified groups.</li>
 * <li>dihedral(g1, g2, g3, g4): the dihedral angle formed by the centers of the four specified groups.</li>
 * </ul>
 * <p>
 * The expression also may involve tabulated functions, and may depend on arbitrary global and per-bond parameters.
 * <p>
 * To use this class, create a CustomCentroidBondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy of each bond.  Then call addPerBondParameter() to define per-bond
 * parameters and addGlobalParameter() to define global parameters.  The values of per-bond parameters are specified
 * as part of the system definition, while values of global parameters may be modified during a simulation by calling
 * Context::setParameter().
 * <p>
 * Next call addGroup() to define the particle groups.  Each group is specified by the particles it contains, and
 * the weights to use when computing the center position.
 * <p>
 * Then call addBond() to define bonds and specify their parameter values.  After a bond has been added, you can
 * modify its parameters by calling setBondParameters().  This will have no effect on Contexts that already exist unless
 * you call updateParametersInContext().
 * <p>
 * As an example, the following code creates a CustomCentroidBondForce that implements a harmonic force between the
 * centers of mass of two groups of particles.
 *
 * <pre>
 *   {@code
 *    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "0.5*k*distance(g1,g2)^2");
 *    force->addPerBondParameter("k");
 *    force->addGroup(particles1);
 *    force->addGroup(particles2);
 *    vector<int> bondGroups;
 *    bondGroups.push_back(0);
 *    bondGroups.push_back(1);
 *    vector<double> bondParameters;
 *    bondParameters.push_back(k);
 *    force->addBond(bondGroups, bondParameters);
 *    }
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
 * <p>
 * This class also supports the functions pointdistance(x1, y1, z1, x2, y2, z2),
 * pointangle(x1, y1, z1, x2, y2, z2, x3, y3, z3), and pointdihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4).
 * These functions are similar to distance(), angle(), and dihedral(), but the arguments are the
 * coordinates of points to perform the calculation based on rather than the names of groups.
 * This enables more flexible geometric calculations.  For example, the following computes the distance
 * from group g1 to the midpoint between groups g2 and g3.
 *
 * <pre>
 *   {@code
 *    CustomCentroidBondForce* force = new CustomCentroidBondForce(3, "pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)");
 *   }
 * </pre>
 * <p>
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in the expression.
 */
public class CustomCentroidBondForce extends Force {

  /**
   * Create a CustomCentroidBondForce.
   *
   * @param numGroups the number of groups used to define each bond
   * @param energy    an algebraic expression giving the interaction energy of each bond as a function
   *                  of particle positions, inter-particle distances, angles, and dihedrals, and any global
   *                  and per-bond parameters
   */
  public CustomCentroidBondForce(int numGroups, String energy) {
    pointer = OpenMM_CustomCentroidBondForce_create(numGroups, energy);
  }

  /**
   * Get the number of groups used to define each bond.
   */
  public int getNumGroupsPerBond() {
    return OpenMM_CustomCentroidBondForce_getNumGroupsPerBond(pointer);
  }

  /**
   * Get the number of particle groups that have been defined.
   */
  public int getNumGroups() {
    return OpenMM_CustomCentroidBondForce_getNumGroups(pointer);
  }

  /**
   * Get the number of bonds for which force field parameters have been defined.
   */
  public int getNumBonds() {
    return OpenMM_CustomCentroidBondForce_getNumBonds(pointer);
  }

  /**
   * Get the number of per-bond parameters that the interaction depends on.
   */
  public int getNumPerBondParameters() {
    return OpenMM_CustomCentroidBondForce_getNumPerBondParameters(pointer);
  }

  /**
   * Get the number of global parameters that the interaction depends on.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomCentroidBondForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of global parameters with respect to which the derivative of the energy
   * should be computed.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of tabulated functions that have been defined.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the number of tabulated functions that have been defined.
   *
   * @deprecated This method exists only for backward compatibility.  Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomCentroidBondForce_getNumFunctions(pointer);
  }

  /**
   * Get the algebraic expression that gives the interaction energy of each bond
   */
  public String getEnergyFunction() {
    return OpenMM_CustomCentroidBondForce_getEnergyFunction(pointer).getString(0);
  }

  /**
   * Set the algebraic expression that gives the interaction energy of each bond
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomCentroidBondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Add a new per-bond parameter that the interaction may depend on.
   *
   * @param name the name of the parameter
   * @return the index of the parameter that was added
   */
  public int addPerBondParameter(String name) {
    return OpenMM_CustomCentroidBondForce_addPerBondParameter(pointer, name);
  }

  /**
   * Get the name of a per-bond parameter.
   *
   * @param index the index of the parameter for which to get the name
   * @return the parameter name
   */
  public String getPerBondParameterName(int index) {
    return OpenMM_CustomCentroidBondForce_getPerBondParameterName(pointer, index).getString(0);
  }

  /**
   * Set the name of a per-bond parameter.
   *
   * @param index the index of the parameter for which to set the name
   * @param name  the name of the parameter
   */
  public void setPerBondParameterName(int index, String name) {
    OpenMM_CustomCentroidBondForce_setPerBondParameterName(pointer, index, name);
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
    return OpenMM_CustomCentroidBondForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index the index of the parameter for which to get the name
   * @return the parameter name
   */
  public String getGlobalParameterName(int index) {
    return OpenMM_CustomCentroidBondForce_getGlobalParameterName(pointer, index).getString(0);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index the index of the parameter for which to set the name
   * @param name  the name of the parameter
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index the index of the parameter for which to get the default value
   * @return the parameter default value
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        the index of the parameter for which to set the default value
   * @param defaultValue the default value of the parameter
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   * The parameter must have already been added with addGlobalParameter().
   *
   * @param name the name of the parameter
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Get the name of a global parameter with respect to which this Force should compute the
   * derivative of the energy.
   *
   * @param index the index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives()
   * @return the parameter name
   */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName(pointer, index).getString(0);
  }

  /**
   * Add a particle group.
   *
   * @param particles the indices of the particles to include in the group
   * @param weights   the weight to use for each particle when computing the center position.
   *                  If this is omitted, then particle masses will be used as weights.
   * @return the index of the group that was added
   */
  public int addGroup(IntArray particles, DoubleArray weights) {
    return OpenMM_CustomCentroidBondForce_addGroup(pointer, particles.getPointer(), weights.getPointer());
  }

  /**
   * Get the properties of a group.
   *
   * @param index     the index of the group to get
   * @param particles the indices of the particles in the group
   * @param weights   the weight used for each particle when computing the center position.
   *                  If no weights were specified, this vector will be empty indicating that particle
   *                  masses should be used as weights.
   */
  public void getGroupParameters(int index, IntArray particles, DoubleArray weights) {
    OpenMM_CustomCentroidBondForce_getGroupParameters(pointer, index, particles.getPointer(), weights.getPointer());
  }

  /**
   * Set the properties of a group.
   *
   * @param index     the index of the group to set
   * @param particles the indices of the particles in the group
   * @param weights   the weight to use for each particle when computing the center position.
   *                  If this is omitted, then particle masses will be used as weights.
   */
  public void setGroupParameters(int index, IntArray particles, DoubleArray weights) {
    OpenMM_CustomCentroidBondForce_setGroupParameters(pointer, index, particles.getPointer(), weights.getPointer());
  }

  /**
   * Add a bond to the force
   *
   * @param groups     the indices of the groups the bond depends on
   * @param parameters the list of per-bond parameter values for the new bond
   * @return the index of the bond that was added
   */
  public int addBond(IntArray groups, DoubleArray parameters) {
    return OpenMM_CustomCentroidBondForce_addBond(pointer, groups.getPointer(), parameters.getPointer());
  }

  /**
   * Get the properties of a bond.
   *
   * @param index      the index of the bond to get
   * @param groups     the indices of the groups in the bond
   * @param parameters the list of per-bond parameter values for the bond
   */
  public void getBondParameters(int index, IntArray groups, DoubleArray parameters) {
    OpenMM_CustomCentroidBondForce_getBondParameters(pointer, index, groups.getPointer(), parameters.getPointer());
  }

  /**
   * Set the properties of a bond.
   *
   * @param index      the index of the bond to set
   * @param groups     the indices of the groups in the bond
   * @param parameters the list of per-bond parameter values for the bond
   */
  public void setBondParameters(int index, IntArray groups, DoubleArray parameters) {
    OpenMM_CustomCentroidBondForce_setBondParameters(pointer, index, groups.getPointer(), parameters.getPointer());
  }


  /**
   * Update the per-bond parameters and tabulated functions in a Context to match those stored in this Force object.  This method provides
   * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
   * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
   * to copy them over to the Context.
   * <p>
   * This method has several limitations.  The only information it updates is the values of per-bond parameters and tabulated
   * functions.  All other aspects of the Force (such as the energy function) are unaffected and can only be changed by reinitializing
   * the Context.  Neither the definitions of groups nor the set of groups involved in a bond can be changed, nor can new
   * bonds be added.  Also, while the tabulated values of a function can change, everything else about it (its dimensions,
   * the data range) must not be changed.
   *
   * @param context the OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomCentroidBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
   *
   * @param periodic 1 if periodic boundary conditions should be used, 0 if not.
   */
  public void setUsesPeriodicBoundaryConditions(int periodic) {
    OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic);
  }

  /**
   * Returns whether this force makes use of periodic boundary
   * conditions.
   *
   * @return true if force uses PBC and false otherwise
   */
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }

  /**
   * Destroy the OpenMM CustomCentroidBondForce.
   */
  public void destroy() {
    OpenMM_CustomCentroidBondForce_destroy(pointer);
  }

}
