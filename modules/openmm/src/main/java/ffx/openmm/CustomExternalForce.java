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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getNumPerParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_getPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_setPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions;

/**
 * Custom External Force.
 */
public class CustomExternalForce extends Force {

  /**
   * Create a CustomExternalForce.
   *
   * @param energy The energy expression for the force.
   */
  public CustomExternalForce(String energy) {
    super(OpenMM_CustomExternalForce_create(energy));
  }

  /**
   * Add a global parameter that the interaction may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomExternalForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the force.
   *
   * @param index              The particle index.
   * @param particleParameters The particle parameters.
   * @return The index of the particle that was added.
   */
  public int addParticle(int index, DoubleArray particleParameters) {
    return OpenMM_CustomExternalForce_addParticle(pointer, index, particleParameters.getPointer());
  }

  /**
   * Add a per-particle parameter that the interaction may depend on.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(String name) {
    return OpenMM_CustomExternalForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomExternalForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the energy expression for the force.
   *
   * @return The energy expression for the force.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomExternalForce_getEnergyFunction(pointer);
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
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomExternalForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomExternalForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_CustomExternalForce_getNumParticles(pointer);
  }

  /**
   * Get the number of per-particle parameters.
   *
   * @return The number of per-particle parameters.
   */
  public int getNumPerParticleParameters() {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particle           The index of the particle (output).
   * @param particleParameters The particle parameters (output).
   */
  public void getParticleParameters(int index, IntBuffer particle, DoubleArray particleParameters) {
    OpenMM_CustomExternalForce_getParticleParameters(pointer, index, particle, particleParameters.getPointer());
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particle           The index of the particle (output).
   * @param particleParameters The particle parameters (output).
   */
  public void getParticleParameters(int index, IntByReference particle, DoubleArray particleParameters) {
    OpenMM_CustomExternalForce_getParticleParameters(pointer, index, particle, particleParameters.getPointer());
  }

  /**
   * Get the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerParticleParameterName(int index) {
    Pointer p = OpenMM_CustomExternalForce_getPerParticleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the energy expression for the force.
   *
   * @param energy The energy expression for the force.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomExternalForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomExternalForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particle           The index of the particle.
   * @param particleParameters The particle parameters.
   */
  public void setParticleParameters(int index, int particle, DoubleArray particleParameters) {
    OpenMM_CustomExternalForce_setParticleParameters(pointer, index, particle, particleParameters.getPointer());
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, String name) {
    OpenMM_CustomExternalForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context to update.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomExternalForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}