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
package ffx.potential.openmm;

import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.State;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;

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
   * Array of atoms.
   */
  private final Atom[] atoms;
  /**
   * Degrees of freedom.
   */
  private final int n;
  /**
   * Number of atoms.
   */
  private final int nAtoms;

  /**
   * Construct an OpenMM State with the requested information.
   *
   * @param pointer Pointer to an OpenMM state.
   * @param atoms   Array of atoms.
   * @param dof     Degrees of freedom.
   */
  protected OpenMMState(PointerByReference pointer, Atom[] atoms, int dof) {
    super(pointer);
    // Set the data types mask using the super class method.
    this.dataTypes = super.getDataTypes();
    this.atoms = atoms;
    this.n = dof;
    nAtoms = atoms.length;
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
   * The force array contains the OpenMM force information for all atoms. The returned array a
   * contains accelerations for active atoms only.
   *
   * @param a Acceleration components for only active atomic coordinates.
   * @return The acceleration for each active atomic coordinate.
   */
  public double[] getAccelerations(@Nullable double[] a) {
    if (!stateContains(OpenMM_State_Forces)) {
      return a;
    }

    if (a == null || a.length != n) {
      a = new double[n];
    }

    double[] forces = getForces();
    for (int i = 0, index = 0; i < nAtoms; i++, index += 3) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        double mass = atom.getMass();
        double xx = forces[index] * OpenMM_AngstromsPerNm / mass;
        double yy = forces[index + 1] * OpenMM_AngstromsPerNm / mass;
        double zz = forces[index + 2] * OpenMM_AngstromsPerNm / mass;
        a[index] = xx;
        a[index + 1] = yy;
        a[index + 2] = zz;
        atom.setAcceleration(xx, yy, zz);
      }
    }
    return a;
  }

  /**
   * The force array contains the OpenMM force information for all atoms. The returned array g
   * contains components for active atoms only.
   *
   * @param g Gradient components for only active atomic coordinates.
   * @return g The gradient includes only active atoms
   */
  public double[] getGradient(@Nullable double[] g) {
    if (!stateContains(OpenMM_State_Forces)) {
      return g;
    }

    if (g == null || g.length != n) {
      g = new double[n];
    }

    double[] forces = getForces();
    for (int i = 0, index = 0; i < nAtoms; i++, index += 3) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        double xx = -forces[index] * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        double yy = -forces[index + 1] * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        double zz = -forces[index + 2] * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        if (isNaN(xx) || isInfinite(xx) || isNaN(yy) || isInfinite(yy) || isNaN(zz) || isInfinite(zz)) {
          StringBuilder sb = new StringBuilder(
              format(" The gradient of atom %s is (%8.3f,%8.3f,%8.3f).", atom, xx, yy, zz));
          double[] vals = new double[3];
          atom.getVelocity(vals);
          sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          atom.getAcceleration(vals);
          sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          atom.getPreviousAcceleration(vals);
          sb.append(
              format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          throw new EnergyException(sb.toString());
        }
        g[index] = xx;
        g[index + 1] = yy;
        g[index + 2] = zz;
        atom.setXYZGradient(xx, yy, zz);
      }
    }
    return g;
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
   * The positions array contains the OpenMM atomic position information for all atoms. The
   * returned array x contains coordinates only for active atoms.
   *
   * @param x Atomic coordinates only for active atoms.
   * @return x The atomic coordinates for only active atoms.
   */
  public double[] getPositions(@Nullable double[] x) {
    if (!stateContains(OpenMM_State_Positions)) {
      return x;
    }

    if (x == null || x.length != n) {
      x = new double[n];
    }

    double[] pos = getPositions();
    for (int i = 0, index = 0; i < nAtoms; i++, index += 3) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        double xx = pos[index] * OpenMM_AngstromsPerNm;
        double yy = pos[index + 1] * OpenMM_AngstromsPerNm;
        double zz = pos[index + 2] * OpenMM_AngstromsPerNm;
        x[index] = xx;
        x[index + 1] = yy;
        x[index + 2] = zz;
        atom.moveTo(xx, yy, zz);
      }
    }
    return x;
  }

  /**
   * The positions array contains the OpenMM atomic position information for all atoms. The
   * returned array x contains coordinates for active atoms only.
   *
   * @param v Velocity only for active atomic coordinates.
   * @return v The velocity for each active atomic coordinate.
   */
  public double[] getVelocities(@Nullable double[] v) {
    if (!stateContains(OpenMM_State_Velocities)) {
      return v;
    }

    if (v == null || v.length != n) {
      v = new double[n];
    }

    double[] vel = getVelocities();
    for (int i = 0, index = 0; i < nAtoms; i++, index += 3) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        double xx = vel[index] * OpenMM_AngstromsPerNm;
        double yy = vel[index + 1] * OpenMM_AngstromsPerNm;
        double zz = vel[index + 2] * OpenMM_AngstromsPerNm;
        v[index] = xx;
        v[index + 1] = yy;
        v[index + 2] = zz;
        atom.setVelocity(xx, yy, zz);
      }
    }
    return v;
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
}
