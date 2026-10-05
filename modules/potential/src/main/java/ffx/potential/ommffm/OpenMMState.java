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
package ffx.potential.ommffm;

import ffx.openmm.ffm.State;
import ffx.openmm.ffm.System;
import ffx.openmm.ffm.bindings.OpenMMNative;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.lang.foreign.MemorySegment;

import static ffx.openmm.ffm.OpenMMUnits.ANGSTROMS_PER_NM;
import static ffx.openmm.ffm.OpenMMUnits.KCAL_PER_KJ;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

/**
 * Retrieve state information from an OpenMM Simulation using the FFM backend.
 */
public class OpenMMState extends State {

  public final double potentialEnergy;
  public final double kineticEnergy;
  public final double totalEnergy;
  private final int dataTypes;

  public OpenMMState(MemorySegment pointer) {
    super(pointer);

    this.dataTypes = super.getDataTypes();
    if (stateContains(OpenMMNative.OpenMM_State_Energy())) {
      potentialEnergy = super.getPotentialEnergy() * KCAL_PER_KJ;
      kineticEnergy = super.getKineticEnergy() * KCAL_PER_KJ;
      totalEnergy = potentialEnergy + kineticEnergy;
    } else {
      potentialEnergy = 0.0;
      kineticEnergy = 0.0;
      totalEnergy = 0.0;
    }
  }

  public double[] getAccelerations(@Nullable double[] a, Atom[] atoms) {
    if (!stateContains(OpenMMNative.OpenMM_State_Forces())) {
      return a;
    }
    double[] forces = getForces();
    int n = forces.length;

    if (atoms == null || atoms.length == 0) {
      throw new IllegalArgumentException("Atoms array must not be null or empty.");
    }
    if (atoms.length * 3 != n) {
      throw new IllegalArgumentException(
          format(" The number of atoms (%d) does not match the number of degrees of freedom (%d).", atoms.length, n));
    }
    if (a == null || a.length != n) {
      a = new double[n];
    }

    int index = 0;
    for (Atom atom : atoms) {
      double mass = atom.getMass();
      double xx = forces[index] * ANGSTROMS_PER_NM / mass;
      double yy = forces[index + 1] * ANGSTROMS_PER_NM / mass;
      double zz = forces[index + 2] * ANGSTROMS_PER_NM / mass;
      a[index] = xx;
      a[index + 1] = yy;
      a[index + 2] = zz;
      index += 3;
    }
    return a;
  }

  public double[] getActiveAccelerations(@Nullable double[] a, Atom[] atoms) {
    if (!stateContains(OpenMMNative.OpenMM_State_Forces())) {
      return a;
    }
    return filterToActive(getAccelerations(null, atoms), a, atoms);
  }

  public double[] getGradient(@Nullable double[] g) {
    if (!stateContains(OpenMMNative.OpenMM_State_Forces())) {
      return g;
    }
    double[] forces = getForces();
    int n = forces.length;

    if (g == null || g.length != n) {
      g = new double[n];
    }

    for (int i = 0; i < n; i++) {
      double xx = -forces[i] * NM_PER_ANGSTROM * KCAL_PER_KJ;
      if (isNaN(xx) || isInfinite(xx)) {
        throw new EnergyException(
            format(" The gradient of degree of freedom %d is %8.3f.", i, xx));
      }
      g[i] = xx;
    }
    return g;
  }

  public double[] getActiveGradient(@Nullable double[] g, Atom[] atoms) {
    if (!stateContains(OpenMMNative.OpenMM_State_Forces())) {
      return g;
    }
    return filterToActive(getGradient(null), g, atoms);
  }

  public double[][] getPeriodicBoxVectors(@Nullable double[][] latticeVectors) {
    if (!stateContains(OpenMMNative.OpenMM_State_Positions())) {
      return latticeVectors;
    }

    if (latticeVectors == null || latticeVectors.length != 3 || latticeVectors[0].length != 3) {
      latticeVectors = new double[3][3];
    }

    System.PeriodicBoxVectors vectors = super.getPeriodicBoxVectors();
    latticeVectors[0][0] = vectors.a().x() * ANGSTROMS_PER_NM;
    latticeVectors[0][1] = vectors.a().y() * ANGSTROMS_PER_NM;
    latticeVectors[0][2] = vectors.a().z() * ANGSTROMS_PER_NM;
    latticeVectors[1][0] = vectors.b().x() * ANGSTROMS_PER_NM;
    latticeVectors[1][1] = vectors.b().y() * ANGSTROMS_PER_NM;
    latticeVectors[1][2] = vectors.b().z() * ANGSTROMS_PER_NM;
    latticeVectors[2][0] = vectors.c().x() * ANGSTROMS_PER_NM;
    latticeVectors[2][1] = vectors.c().y() * ANGSTROMS_PER_NM;
    latticeVectors[2][2] = vectors.c().z() * ANGSTROMS_PER_NM;
    return latticeVectors;
  }

  public double[] getCoordinates(@Nullable double[] x) {
    if (!stateContains(OpenMMNative.OpenMM_State_Positions())) {
      return x;
    }
    double[] positions = getPositions();
    int n = positions.length;
    if (x == null || x.length != n) {
      x = new double[n];
    }
    for (int i = 0; i < n; i++) {
      x[i] = positions[i] * ANGSTROMS_PER_NM;
    }
    return x;
  }

  public double[] getActiveCoordinates(@Nullable double[] x, Atom[] atoms) {
    if (!stateContains(OpenMMNative.OpenMM_State_Positions())) {
      return x;
    }
    return filterToActive(getCoordinates(null), x, atoms);
  }

  public double[] getVelocities(@Nullable double[] v) {
    if (!stateContains(OpenMMNative.OpenMM_State_Velocities())) {
      return v;
    }
    double[] velocities = super.getVelocities();
    int n = velocities.length;
    if (v == null || v.length != n) {
      v = new double[n];
    }
    for (int i = 0; i < n; i++) {
      v[i] = velocities[i] * ANGSTROMS_PER_NM;
    }
    return v;
  }

  public double[] getActiveVelocities(@Nullable double[] v, Atom[] atoms) {
    if (!stateContains(OpenMMNative.OpenMM_State_Velocities())) {
      return v;
    }
    return filterToActive(getVelocities(null), v, atoms);
  }

  public double getPotentialEnergy() {
    return potentialEnergy;
  }

  public double getKineticEnergy() {
    return kineticEnergy;
  }

  public double getTotalEnergy() {
    return totalEnergy;
  }

  @Override
  public double getPeriodicBoxVolume() {
    return super.getPeriodicBoxVolume() * ANGSTROMS_PER_NM * ANGSTROMS_PER_NM * ANGSTROMS_PER_NM;
  }

  public int getDataTypes() {
    return dataTypes;
  }

  public boolean stateContains(int mask) {
    return (dataTypes & mask) == mask;
  }

  private double[] filterToActive(double[] values, @Nullable double[] activeValues, Atom[] atoms) {
    if (values == null) {
      return null;
    }
    int nActive = 0;
    for (Atom atom : atoms) {
      if (atom.isActive()) {
        nActive++;
      }
    }
    int n = nActive * 3;
    if (activeValues == null || activeValues.length != n) {
      activeValues = new double[n];
    }
    int index = 0;
    int activeIndex = 0;
    for (Atom atom : atoms) {
      if (atom.isActive()) {
        activeValues[activeIndex++] = values[index];
        activeValues[activeIndex++] = values[index + 1];
        activeValues[activeIndex++] = values[index + 2];
      }
      index += 3;
    }
    return activeValues;
  }
}
