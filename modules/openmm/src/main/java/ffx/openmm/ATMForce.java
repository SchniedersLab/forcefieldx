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
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import java.nio.DoubleBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_addForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getNumForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_getPerturbationEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_ATMForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * The ATMForce class implements the Alchemical Transfer Method (ATM) for OpenMM.
 * ATM is used to compute the binding free energies of molecular complexes and of other equilibrium processes.
 * ATM and its implementation are described in the open access article:
 * <p>
 * Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio.
 * Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.
 * J. Chem. Inf. Model.  62, 309 (2022)
 * https://doi.org/10.1021/acs.jcim.1c01129
 * <p>
 * Refer to the publication above for a detailed description of the ATM method and the parameters used in this API
 * and please cite it to support our work if you use this software in your research.
 * <p>
 * The ATMForce implements an arbitrary potential energy function that depends on the potential
 * energies (u0 and u1) of the system before and after a set of atoms are displaced by a specified amount.
 * For example, you might displace a molecule from the solvent bulk to a receptor binding site to simulate
 * a binding process.  The potential energy function typically also depends on one or more parameters that
 * are dialed to implement alchemical transformations.
 * <p>
 * To use this class, create an ATMForce object, passing an algebraic expression to the
 * constructor that defines the potential energy. This expression can be any combination
 * of the variables u0 and u1. Then call addGlobalParameter() to define the parameters on which the potential energy expression depends.
 * The values of global parameters may be modified during a simulation by calling Context::setParameter().
 * Next, call addForce() to add Force objects that define the terms of the potential energy function
 * that change upon displacement. Finally, call addParticle() to specify the displacement applied to
 * each particle. Displacements can be changed by calling setParticleParameters(). As any per-particle parameters,
 * changes in displacements take effect only after calling updateParametersInContext().
 * <p>
 * As an example, the following code creates an ATMForce based on the change in energy of
 * two particles when the second particle is displaced by 1 nm in the x direction.
 * The energy change is dialed using an alchemical parameter Lambda, which in this case is set to 1/2:
 *
 * <pre>
 *    ATMForce atmforce = new ATMForce("u0 + Lambda*(u1 - u0)");
 *    atm.addGlobalParameter("Lambda", 0.5);
 *    atm.addParticle(new OpenMM_Vec3(0, 0, 0));
 *    atm.addParticle(new OpenMM_Vec3(1, 0, 0));
 *    CustomBondForce force = new CustomBondForce("0.5*r^2");
 *    atm.addForce(force);
 * </pre>
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta,
 * select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 * <p>
 * If instead of the energy expression the ATMForce constructor specifies the values of a series of parameters,
 * the default energy expression is used:
 *
 * <pre>
 *    select(step(Direction), u0, u1) + ((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0;
 *    usc = select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);
 *    fsc = (z^Acore-1)/(z^Acore+1);
 *    z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;
 *    y = (u-Ubcore)/(Umax-Ubcore);
 *    u = select(step(Direction), 1, -1)*(u1-u0)
 * </pre>
 * <p>
 * which is the same as the soft-core softplus alchemical potential energy function in the Azimi et al. paper above.
 * <p>
 * The ATMForce is then added to the System as any other Force
 *
 * <pre>
 *  system.addForce(atmforce);
 * </pre>
 * <p>
 * after which it will be used for energy/force evaluations for molecular dynamics and energy optimization.
 * You can call getPerturbationEnergy() to query the values of u0 and u1, which are needed for computing
 * free energies.
 * <p>
 * In most cases, particles are only displaced in one of the two states evaluated by this force.  It computes the
 * change in energy between the current particle coordinates (as stored in the Context) and the displaced coordinates.
 * In some cases, it is useful to apply displacements to both states.  You can do this by providing two displacement
 * vectors to addParticle():
 *
 * <pre>
 *    atm.addParticle(new OpenMM_Vec3(1, 0, 0), new OpenMM_Vec3(-1, 0, 0));
 * </pre>
 * <p>
 * In this case, u1 will be computed after displacing the particle in the positive x direction, and
 * u0 will be computed after displacing it in the negative x direction.
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 */
public class ATMForce extends Force {

  /**
   * Create an ATMForce object.
   *
   * @param energy an algebraic expression giving the energy of the system as a function
   *               of u0 and u1, the energies before and after displacement
   */
  public ATMForce(String energy) {
    super(OpenMM_ATMForce_create(energy));
  }

  /**
   * Create an ATMForce object with the default softplus energy expression.  The values passed to
   * this constructor are the default values of the global parameters for newly created Contexts.
   * Their values can be changed by calling setParameter() on the Context using the parameter
   * names defined by the Lambda1(), Lambda2(), etc. methods below.
   *
   * @param lambda1   the default value of the Lambda1 parameter (dimensionless).  This should be
   *                  a number between 0 and 1.
   * @param lambda2   the default value of the Lambda2 parameter (dimensionless).  This should be
   *                  a number between 0 and 1.
   * @param alpha     the default value of the Alpha parameter (kJ/mol)^-1
   * @param uh        the default value of the Uh parameter (kJ/mol)
   * @param w0        the default value of the W0 parameter (kJ/mol)
   * @param umax      the default value of the Umax parameter (kJ/mol)
   * @param ubcore    the default value of the Ubcore parameter (kJ/mol)
   * @param acore     the default value of the Acore parameter dimensionless)
   * @param direction the default value of the Direction parameter (dimensionless).  This should be
   *                  either 1 for the forward transfer, or -1 for the backward transfer.
   */
  public ATMForce(double lambda1, double lambda2, double alpha, double uh, double w0,
                  double umax, double ubcore, double acore, double direction) {
    super(OpenMM_ATMForce_create_2(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore, direction));
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   * The parameter must have already been added with addGlobalParameter().
   *
   * @param name the name of the parameter
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_ATMForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Request that this Force compute the derivative of its energy with respect to a global parameter.
   * The parameter must have already been added with addGlobalParameter().
   *
   * @param name the name of the parameter
   */
  public void addEnergyParameterDerivative(Pointer name) {
    OpenMM_ATMForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a Force whose energy will be computed by the ATMForce.
   *
   * @param force the Force to the be added, which should have been created on the heap with the
   *              "new" operator.  The ATMForce takes over ownership of it, and deletes the Force when the
   *              ATMForce itself is deleted.
   * @return The index within ATMForce of the force that was added
   */
  public int addForce(Force force) {
    return OpenMM_ATMForce_addForce(pointer, force.getPointer());
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
    return OpenMM_ATMForce_addGlobalParameter(pointer, name, defaultValue);
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
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_ATMForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the force.
   * <p>
   * All of the particles in the System must be added to the ATMForce in the same order
   * as they appear in the System.
   *
   * @param displacement1 the displacement of the particle for the target state in nm
   * @param displacement0 the displacement of the particle for the initial state in nm
   * @return the index of the particle that was added
   */
  public int addParticle(OpenMM_Vec3 displacement1, OpenMM_Vec3 displacement0) {
    return OpenMM_ATMForce_addParticle(pointer, displacement1, displacement0);
  }


  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_ATMForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the energy function.
   *
   * @return The energy function.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_ATMForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of an energy parameter derivative.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter derivative.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_ATMForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get a force by index.
   *
   * @param index The index of the force.
   * @return The force at the specified index.
   */
  public PointerByReference getForce(int index) {
    return OpenMM_ATMForce_getForce(pointer, index);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_ATMForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_ATMForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of energy parameter derivatives.
   *
   * @return The number of energy parameter derivatives.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_ATMForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of forces.
   *
   * @return The number of forces.
   */
  public int getNumForces() {
    return OpenMM_ATMForce_getNumForces(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_ATMForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_ATMForce_getNumParticles(pointer);
  }

  /**
   * Get the parameters for a particle
   *
   * @param index         the index in the force for the particle for which to get parameters
   * @param displacement1 the displacement of the particle for the target state in nm
   * @param displacement0 the displacement of the particle for the initial state in nm
   */
  public void getParticleParameters(int index, OpenMM_Vec3 displacement1, OpenMM_Vec3 displacement0) {
    OpenMM_ATMForce_getParticleParameters(pointer, index, displacement1, displacement0);
  }

  /**
   * Returns the current perturbation energy.
   *
   * @param context the Context for which to return the energy
   * @param u1      on exit, the energy of the displaced state
   * @param u0      on exit, the energy of the non-displaced state
   * @param energy  on exit, the value of this force's energy function
   */
  public void getPerturbationEnergy(Context context, DoubleByReference u0, DoubleByReference u1, DoubleByReference energy) {
    OpenMM_ATMForce_getPerturbationEnergy(pointer, context.getPointer(), u0, u1, energy);
  }

  /**
   * Returns the current perturbation energy.
   *
   * @param context the Context for which to return the energy
   * @param u1      on exit, the energy of the displaced state
   * @param u0      on exit, the energy of the non-displaced state
   * @param energy  on exit, the value of this force's energy function
   */
  public void getPerturbationEnergy(Context context, DoubleBuffer u0, DoubleBuffer u1, DoubleBuffer energy) {
    OpenMM_ATMForce_getPerturbationEnergy(pointer, context.getPointer(), u0, u1, energy);
  }


  /**
   * Set the energy function.
   *
   * @param energy The energy function.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_ATMForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the energy function.
   *
   * @param energy The energy function.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_ATMForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_ATMForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_ATMForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_ATMForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the parameters for a particle
   *
   * @param index         the index in the force of the particle for which to set parameters
   * @param displacement1 the displacement of the particle for the target state in nm
   * @param displacement0 the displacement of the particle for the initial state in nm
   */
  public void setParticleParameters(int index, OpenMM_Vec3 displacement1, OpenMM_Vec3 displacement0) {
    OpenMM_ATMForce_setParticleParameters(pointer, index, displacement1, displacement0);
  }

  /**
   * Update the per-particle parameters in a Context to match those stored in this Force object.  This method
   * should be called after updating parameters with setParticleParameters() to copy them over to the Context.
   * The only information this method updates is the values of per-particle parameters.  The number of particles
   * cannot be changed.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_ATMForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Returns whether or not this force makes use of periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_ATMForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}