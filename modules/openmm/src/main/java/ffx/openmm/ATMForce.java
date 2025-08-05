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
 * This class implements the Alchemical Transfer Method (ATM) force for free energy calculations.
 * ATM is a method for computing free energy differences between different states of a system
 * by introducing alchemical parameters that smoothly transform one state into another.
 * <p>
 * The ATM force allows for the calculation of binding free energies and other thermodynamic
 * properties by defining a reaction coordinate that connects different chemical states.
 */
public class ATMForce extends Force {

  /**
   * Create a new ATMForce with an energy function.
   *
   * @param energy The energy function for the ATM force.
   */
  public ATMForce(String energy) {
    super(OpenMM_ATMForce_create(energy));
  }

  /**
   * Create a new ATMForce with specific parameters.
   *
   * @param lambda1   The first lambda parameter.
   * @param lambda2   The second lambda parameter.
   * @param alpha     The alpha parameter.
   * @param uh        The uh parameter.
   * @param w0        The w0 parameter.
   * @param umax      The umax parameter.
   * @param ubcore    The ubcore parameter.
   * @param acore     The acore parameter.
   * @param direction The direction parameter.
   */
  public ATMForce(double lambda1, double lambda2, double alpha, double uh, double w0,
                  double umax, double ubcore, double acore, double direction) {
    super(OpenMM_ATMForce_create_2(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore, direction));
  }

  /**
   * Add an energy parameter derivative.
   *
   * @param name The name of the parameter.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_ATMForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add an energy parameter derivative.
   *
   * @param name The name of the parameter.
   */
  public void addEnergyParameterDerivative(Pointer name) {
    OpenMM_ATMForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a force to the ATM force.
   *
   * @param force The force to add.
   * @return The index of the force that was added.
   */
  public int addForce(Force force) {
    return OpenMM_ATMForce_addForce(pointer, force.getPointer());
  }

  /**
   * Add a global parameter.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_ATMForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a global parameter.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_ATMForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the force.
   *
   * @param coordinates1 The coordinates for state 1.
   * @param coordinates2 The coordinates for state 2.
   * @return The index of the particle that was added.
   */
  public int addParticle(OpenMM_Vec3 coordinates1, OpenMM_Vec3 coordinates2) {
    return OpenMM_ATMForce_addParticle(pointer, coordinates1, coordinates2);
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
   * Get the parameters for a particle.
   *
   * @param index        The index of the particle.
   * @param coordinates1 The coordinates for state 1 (output).
   * @param coordinates2 The coordinates for state 2 (output).
   */
  public void getParticleParameters(int index, OpenMM_Vec3 coordinates1, OpenMM_Vec3 coordinates2) {
    OpenMM_ATMForce_getParticleParameters(pointer, index, coordinates1, coordinates2);
  }

  /**
   * Get the perturbation energy.
   *
   * @param context The context.
   * @param u0      The u0 energy (output).
   * @param u1      The u1 energy (output).
   * @param u2      The u2 energy (output).
   */
  public void getPerturbationEnergy(Context context, DoubleByReference u0, DoubleByReference u1, DoubleByReference u2) {
    OpenMM_ATMForce_getPerturbationEnergy(pointer, context.getPointer(), u0, u1, u2);
  }

  /**
   * Get the perturbation energy.
   *
   * @param context The context.
   * @param u0      The u0 energy (output).
   * @param u1      The u1 energy (output).
   * @param u2      The u2 energy (output).
   */
  public void getPerturbationEnergy(Context context, DoubleBuffer u0, DoubleBuffer u1, DoubleBuffer u2) {
    OpenMM_ATMForce_getPerturbationEnergy(pointer, context.getPointer(), u0, u1, u2);
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
   * Set the parameters for a particle.
   *
   * @param index        The index of the particle.
   * @param coordinates1 The coordinates for state 1.
   * @param coordinates2 The coordinates for state 2.
   */
  public void setParticleParameters(int index, OpenMM_Vec3 coordinates1, OpenMM_Vec3 coordinates2) {
    OpenMM_ATMForce_setParticleParameters(pointer, index, coordinates1, coordinates2);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_ATMForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_ATMForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}