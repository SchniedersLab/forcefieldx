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

import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getDataTypes;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPeriodicBoxVolume;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getStepCount;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getTime;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getVelocities;

/**
 * A State object records a snapshot of the current state of a simulation at a point in time.
 * You create it by calling getState() on a Context.
 * <p>
 * When a State is created, you specify what information should be stored in it.  This saves
 * time and memory by only copying in the information that you actually want.  This is especially
 * important for forces and energies, since they may need to be calculated.  If you query a
 * State object for a piece of information which is not available (because it was not requested
 * when the State was created), it will throw an exception.
 */
public class State {

  /**
   * State pointer.
   */
  PointerByReference pointer;

  /**
   * Constructor.
   *
   * @param pointer Pointer to the state returned by a Context.
   */
  public State(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Destroy the state.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_State_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the data types.
   *
   * @return The data types.
   */
  public int getDataTypes() {
    return OpenMM_State_getDataTypes(pointer);
  }

  /**
   * Get the energy parameter derivatives.
   *
   * @return The energy parameter derivatives.
   */
  public PointerByReference getEnergyParameterDerivatives() {
    return OpenMM_State_getEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the forces.
   *
   * @return The forces.
   */
  public double[] getForces() {
    Vec3Array forces = new Vec3Array(OpenMM_State_getForces(pointer));
    return forces.getArray();
  }

  /**
   * Get the kinetic energy.
   *
   * @return The kinetic energy.
   */
  public double getKineticEnergy() {
    return OpenMM_State_getKineticEnergy(pointer);
  }

  /**
   * Get the parameters.
   *
   * @return The parameters.
   */
  public PointerByReference getParameters() {
    return OpenMM_State_getParameters(pointer);
  }

  /**
   * Get the periodic box vectors.
   *
   * @return The periodic box vectors.
   */
  public double[][] getPeriodicBoxVectors() {
    OpenMM_Vec3 a = new OpenMM_Vec3();
    OpenMM_Vec3 b = new OpenMM_Vec3();
    OpenMM_Vec3 c = new OpenMM_Vec3();
    OpenMM_State_getPeriodicBoxVectors(pointer, a, b, c);
    double[][] latticeVectors = new double[3][3];
    latticeVectors[0][0] = a.x;
    latticeVectors[0][1] = a.y;
    latticeVectors[0][2] = a.z;
    latticeVectors[1][0] = b.x;
    latticeVectors[1][1] = b.y;
    latticeVectors[1][2] = b.z;
    latticeVectors[2][0] = c.x;
    latticeVectors[2][1] = c.y;
    latticeVectors[2][2] = c.z;
    return latticeVectors;
  }

  /**
   * Get the periodic box volume.
   *
   * @return The periodic box volume.
   */
  public double getPeriodicBoxVolume() {
    return OpenMM_State_getPeriodicBoxVolume(pointer);
  }

  /**
   * Get the pointer to the state.
   *
   * @return The pointer to the state.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get the positions.
   *
   * @return The positions.
   */
  public double[] getPositions() {
    Vec3Array positions = new Vec3Array(OpenMM_State_getPositions(pointer));
    return positions.getArray();
  }

  /**
   * Get the potential energy.
   *
   * @return The potential energy.
   */
  public double getPotentialEnergy() {
    return OpenMM_State_getPotentialEnergy(pointer);
  }

  /**
   * Get the step count.
   *
   * @return The step count.
   */
  public long getStepCount() {
    return OpenMM_State_getStepCount(pointer);
  }

  /**
   * Get the time.
   *
   * @return The time.
   */
  public double getTime() {
    return OpenMM_State_getTime(pointer);
  }

  /**
   * Get the velocities.
   *
   * @return The velocities.
   */
  public double[] getVelocities() {
    Vec3Array velocities = new Vec3Array(OpenMM_State_getVelocities(pointer));
    return velocities.getArray();
  }

}