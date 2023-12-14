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
package ffx.potential.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions;

public class OpenMMCustomCentroidBondForce extends OpenMMForce {

  /**
   * OpenMM CustomCentroidBondForce constructor.
   *
   * @param nGroups The number of particles in the bond.
   * @param energy  The energy expression.
   */
  public OpenMMCustomCentroidBondForce(int nGroups, String energy) {
    pointer = OpenMM_CustomCentroidBondForce_create(nGroups, energy);
  }

  /**
   * Add a per bond parameters
   *
   * @param name The parameter name.
   */
  public void addPerBondParameter(String name) {
    OpenMM_CustomCentroidBondForce_addPerBondParameter(pointer, name);
  }

  /**
   * Add a atoms of atoms to the force.
   *
   * @param atoms  The group of atoms.
   * @param weight The weight of each atom.
   */
  public void addGroup(OpenMMIntArray atoms, OpenMMDoubleArray weight) {
    OpenMM_CustomCentroidBondForce_addGroup(pointer, atoms.getPointer(), weight.getPointer());
  }

  /**
   * Add a bond between two groups to the force.
   *
   * @param groups     The two groups.
   * @param parameters The parameters of each groups.
   */
  public void addBond(OpenMMIntArray groups, OpenMMDoubleArray parameters) {
    OpenMM_CustomCentroidBondForce_addBond(pointer, groups.getPointer(), parameters.getPointer());
  }

  /**
   * Set whether to use periodic boundary conditions.
   *
   * @param periodic 1 if periodic boundary conditions should be used, 0 if not.
   */
  public void setUsesPeriodicBoundaryConditions(int periodic) {
    OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic);
  }

  /**
   * Destroy the OpenMM CustomCentroidBondForce.
   */
  public void destroy() {
    OpenMM_CustomCentroidBondForce_destroy(pointer);
  }

}
