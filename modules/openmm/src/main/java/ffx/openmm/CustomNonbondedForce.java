// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addInteractionGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setUseSwitchingFunction;

/**
 * Custom Non-bonded Force.
 */
public class CustomNonbondedForce extends Force {

  public CustomNonbondedForce(String energy) {
    pointer = OpenMM_CustomNonbondedForce_create(energy);
  }

  /**
   * Add a global parameter.
   *
   * @param name  The parameter name.
   * @param value The parameter value.
   */
  public void addGlobalParameter(String name, double value) {
    OpenMM_CustomNonbondedForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add a per particle parameter.
   *
   * @param name The parameter name.
   */
  public void addPerParticleParameter(String name) {
    OpenMM_CustomNonbondedForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a particle to the force.
   *
   * @param parameters The particle parameters.
   */
  public void addParticle(DoubleArray parameters) {
    OpenMM_CustomNonbondedForce_addParticle(pointer, parameters.getPointer());
  }

  /**
   * Add an interaction group.
   *
   * @param group1 The first group.
   * @param group2 The second group.
   */
  public void addInteractionGroup(IntSet group1, IntSet group2) {
    OpenMM_CustomNonbondedForce_addInteractionGroup(pointer, group1.getPointer(), group2.getPointer());
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
   * Set the cutoff distance.
   *
   * @param off The cutoff distance.
   */
  public void setCutoffDistance(double off) {
    OpenMM_CustomNonbondedForce_setCutoffDistance(pointer, off);
  }

  /**
   * Flag to contol use of a switching function.
   *
   * @param useSwitchingFunction If 1, the switching function is used.
   */
  public void setUseSwitchingFunction(int useSwitchingFunction) {
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(pointer, useSwitchingFunction);
  }

  /**
   * Set the switching distance.
   *
   * @param switchingDistance The switching distance.
   */
  public void setSwitchingDistance(double switchingDistance) {
    OpenMM_CustomNonbondedForce_setSwitchingDistance(pointer, switchingDistance);
  }

  /**
   * Add an exclusion.
   *
   * @param particle1 The first particle.
   * @param particle2 The second particle.
   */
  public void addExclusion(int particle1, int particle2) {
    OpenMM_CustomNonbondedForce_addExclusion(pointer, particle1, particle2);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomNonbondedForce_destroy(pointer);
      pointer = null;
    }
  }
}
