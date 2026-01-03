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
package ffx.potential.openmm;

import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.State;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.util.Arrays;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

/**
 * Retrieve state information from an OpenMM Simulation.
 */
public class OpenMMState extends State {

  /**
   * Potential energy (kcal/mol).
   */
  public final double potentialEnergy;
  /**
   * Kinetic energy (kcal/mol).
   */
  public final double kineticEnergy;
  /**
   * Total energy (kcal/mol).
   */
  public final double totalEnergy;
  /**
   * Mask of information to retrieve.
   */
  private final int dataTypes;

  /**
   * Construct an OpenMM State with the requested information.
   *
   * @param pointer Pointer to an OpenMM state.
   */
  protected OpenMMState(PointerByReference pointer) {
    super(pointer);

    // Set the data types mask using the super class method.
    this.dataTypes = super.getDataTypes();
    if (stateContains(OpenMM_State_Energy)) {
      // Set the energy fields using the super class method and convert units.
      potentialEnergy = super.getPotentialEnergy() * OpenMM_KcalPerKJ;
      kineticEnergy = super.getKineticEnergy() * OpenMM_KcalPerKJ;
      totalEnergy = potentialEnergy + kineticEnergy;
    } else {
      potentialEnergy = 0.0;
      kineticEnergy = 0.0;
      totalEnergy = 0.0;
    }
  }

  /**
   * The acceleration array will contain the acceleration information for all atoms. This
   * method will convert the OpenMM forces to accelerations using the atom masses.
   *
   * @param a     Acceleration components for all atoms.
   * @param atoms The array of atoms, which is needed to convert from force to acceleration.
   * @return The acceleration for each atom in units of Angstroms per picosecond squared.
   */
  public double[] getAccelerations(@Nullable double[] a, Atom[] atoms) {
    // Check if the state contains forces.
    if (!stateContains(OpenMM_State_Forces)) {
      return a;
    }
    double[] forces = getForces();
    int n = forces.length;

    // Validate the atoms array exists and is not empty.
    if (atoms == null || atoms.length == 0) {
      throw new IllegalArgumentException("Atoms array must not be null or empty.");
    }
    // Validate the number of degrees of freedom.
    if (atoms.length * 3 != n) {
      String message = format(" The number of atoms (%d) does not match the number of degrees of freedom (%d).", atoms.length, n);
      throw new IllegalArgumentException(message);
    }
    if (a == null || a.length != n) {
      a = new double[n];
    }

    int index = 0;
    for (Atom atom : atoms) {
      double mass = atom.getMass();
      double xx = forces[index] * OpenMM_AngstromsPerNm / mass;
      double yy = forces[index + 1] * OpenMM_AngstromsPerNm / mass;
      double zz = forces[index + 2] * OpenMM_AngstromsPerNm / mass;
      a[index] = xx;
      a[index + 1] = yy;
      a[index + 2] = zz;
      index += 3;
    }
    return a;
  }

  /**
   * The acceleration array will contain the acceleration information for all atoms. This
   * method will convert the OpenMM forces to accelerations using the atom masses.
   *
   * @param a     Acceleration components for all atoms.
   * @param atoms The array of atoms, which is needed to convert from force to acceleration.
   * @return The acceleration for each atom in units of Angstroms per picosecond squared.
   */
  public double[] getActiveAccelerations(@Nullable double[] a, Atom[] atoms) {
    if (!stateContains(OpenMM_State_Forces)) {
      return a;
    }
    return filterToActive(getAccelerations(null, atoms), a, atoms);
  }

  /**
   * The force array contains the OpenMM force information for all atoms.
   *
   * @param g Gradient array to use.
   * @return The gradient for all atoms in units of kcal/mol/Angstrom (in a new array if g is null or the wrong size).
   */
  public double[] getGradient(@Nullable double[] g) {
    // Check if the state contains forces.
    if (!stateContains(OpenMM_State_Forces)) {
      return g;
    }
    double[] forces = getForces();
    int n = forces.length;

    // Validate the gradient array exists and is the correct size.
    if (g == null || g.length != n) {
      g = new double[n];
    }

    for (int i = 0; i < n; i++) {
      double xx = -forces[i] * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
      if (isNaN(xx) || isInfinite(xx)) {
        throw new EnergyException(
            format(" The gradient of degree of freedom %d is %8.3f.", i, xx));
      }
      g[i] = xx;
    }
    return g;
  }

  /**
   * The force array contains the OpenMM force information for active atoms.
   *
   * @param g Gradient array to use.
   * @return The gradient for all atoms in units of kcal/mol/Angstrom (in a new array if g is null or the wrong size).
   */
  public double[] getActiveGradient(@Nullable double[] g, Atom[] atoms) {
    if (!stateContains(OpenMM_State_Forces)) {
      return g;
    }
    return filterToActive(getGradient(null), g, atoms);
  }

  /**
   * Read the periodic lattice vectors from a state.
   *
   * <p>The crystal instance will be updated, and passed to the ForceFieldEnergy instance.
   */
  public double[][] getPeriodicBoxVectors() {
    if (!stateContains(OpenMM_State_Positions)) {
      return null;
    }

    double[][] latticeVectors = super.getPeriodicBoxVectors();
    latticeVectors[0][0] *= OpenMM_AngstromsPerNm;
    latticeVectors[0][1] *= OpenMM_AngstromsPerNm;
    latticeVectors[0][2] *= OpenMM_AngstromsPerNm;
    latticeVectors[1][0] *= OpenMM_AngstromsPerNm;
    latticeVectors[1][1] *= OpenMM_AngstromsPerNm;
    latticeVectors[1][2] *= OpenMM_AngstromsPerNm;
    latticeVectors[2][0] *= OpenMM_AngstromsPerNm;
    latticeVectors[2][1] *= OpenMM_AngstromsPerNm;
    latticeVectors[2][2] *= OpenMM_AngstromsPerNm;
    return latticeVectors;
  }

  /**
   * The position array contains the OpenMM atomic position information for all atoms. The
   * returned array x is in units of Angstroms.
   *
   * @param x Atomic coordinates array to use.
   * @return The atomic coordinates (in a new array if x is null or the wrong size).
   */
  public double[] getPositions(@Nullable double[] x) {
    // Check if the state contains positions.
    if (!stateContains(OpenMM_State_Positions)) {
      return x;
    }

    double[] pos = getPositions();
    int n = pos.length;

    // Allocate x if null or the wrong size.
    if (x == null || x.length != n) {
      x = new double[n];
    }

    for (int i = 0; i < n; i++) {
      x[i] = pos[i] * OpenMM_AngstromsPerNm;
    }

    return x;
  }

  /**
   * The position array contains the OpenMM atomic position information for active atoms.
   * The returned array x is in units of Angstroms.
   *
   * @param x Atomic coordinates array to use.
   * @return The atomic coordinates (in a new array if x is null or the wrong size).
   */
  public double[] getActivePositions(@Nullable double[] x, Atom[] atoms) {
    if (!stateContains(OpenMM_State_Positions)) {
      return x;
    }
    return filterToActive(getPositions(null), x, atoms);
  }

  /**
   * The velocity array contains the OpenMM atomic position information for all atoms. The
   * returned array v is in units of Angstroms per picosecond.
   *
   * @param v Atomic velocity array to use.
   * @return The atomic velocities (in a new array if v is null or the wrong size).
   */
  public double[] getVelocities(@Nullable double[] v) {
    if (!stateContains(OpenMM_State_Velocities)) {
      return v;
    }

    double[] vel = getVelocities();
    int n = vel.length;

    // Validate the velocity array exists and is the correct size.
    if (v == null || v.length != n) {
      v = new double[n];
    }

    for (int i = 0; i < n; i++) {
      v[i] = vel[i] * OpenMM_AngstromsPerNm;
    }

    return v;
  }

  /**
   * The velocity array contains the OpenMM atomic position information for active atoms.
   * The returned array v is in units of Angstroms per picosecond.
   *
   * @param v Atomic velocity array to use.
   * @return The atomic velocities (in a new array if v is null or the wrong size).
   */
  public double[] getActiveVelocities(@Nullable double[] v, Atom[] atoms) {
    if (!stateContains(OpenMM_State_Velocities)) {
      return v;
    }
    return filterToActive(getVelocities(null), v, atoms);
  }

  /**
   * Get the periodic box volume.
   *
   * @return The periodic box volume.
   */
  @Override
  public double getPeriodicBoxVolume() {
    return super.getPeriodicBoxVolume()
        * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm;
  }

  /**
   * Get the potential energy. This field will be zero if the dataTypes mask did not include the energy.
   *
   * @return The potential energy.
   */
  @Override
  public double getPotentialEnergy() {
    return potentialEnergy;
  }

  /**
   * Get the kinetic energy. This field will be zero if the dataTypes mask did not include the energy.
   *
   * @return The kinetic energy.
   */
  @Override
  public double getKineticEnergy() {
    return kineticEnergy;
  }

  /**
   * Get the total energy. This field will be zero if the dataTypes mask did not include the energy.
   *
   * @return The total energy.
   */
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * Get the mask of information contained in the state.
   *
   * @return The mask of information contained in the state.
   */
  @Override
  public int getDataTypes() {
    return dataTypes;
  }

  /**
   * Check to see if the state contains the requested information.
   *
   * @param dataType Information to check for.
   * @return boolean indicating whether the state contains the requested information.
   */
  private boolean stateContains(int dataType) {
    return (dataTypes & dataType) == dataType;
  }

  /**
   * Filter an array to only include elements where the corresponding Atom is active.
   *
   * @param source The source array to filter.
   * @param target The target array to fill with filtered values (or null to create a new one).
   * @param atoms  The array of Atoms, which should have a corresponding active mask.
   */
  private static double[] filterToActive(double[] source, @Nullable double[] target, Atom[] atoms) {
    if (source == null || atoms == null) {
      throw new IllegalArgumentException("The arrays must be non-null.");
    }

    // Validate that the source array length is three times the number of atoms.
    if (source.length != atoms.length * 3) {
      throw new IllegalArgumentException("Source array length must be three times the number of atoms.");
    }

    // Count the number of active atoms.
    int count = (int) Arrays.stream(atoms).filter(Atom::isActive).count();

    // Ensure target is large enough to hold the filtered values.
    if (target == null || target.length < count * 3) {
      target = new double[count * 3];
    }

    // Fill the target array with values from the source array for active atoms.
    int sourceIndedx = 0;
    int targetIndex = 0;
    for (Atom atom : atoms) {
      if (atom.isActive()) {
        target[targetIndex] = source[sourceIndedx];
        target[targetIndex + 1] = source[sourceIndedx + 1];
        target[targetIndex + 2] = source[sourceIndedx + 2];
        targetIndex += 3;
      }
      // Always increment the source index by 3 for each atom.
      sourceIndedx += 3;
    }
    return target;
  }
}
