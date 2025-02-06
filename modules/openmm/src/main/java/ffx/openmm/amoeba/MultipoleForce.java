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
package ffx.openmm.amoeba;

import ffx.openmm.Context;
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_addMultipole;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setAEwald;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCovalentMap;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMultipoleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPolarizationType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_updateParametersInContext;

/**
 * Amoeba Polarizable Multipole Force.
 */
public class MultipoleForce extends Force {

  public MultipoleForce() {
    pointer = OpenMM_AmoebaMultipoleForce_create();
  }

  /**
   * Set the polarization method.
   *
   * @param method The polarization method.
   */
  public void setPolarizationType(int method) {
    OpenMM_AmoebaMultipoleForce_setPolarizationType(pointer, method);
  }

  /**
   * Set extrapolation coefficients.
   *
   * @param exptCoefficients The extrapolation coefficients.
   */
  public void setExtrapolationCoefficients(DoubleArray exptCoefficients) {
    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(pointer, exptCoefficients.getPointer());
  }

  /**
   * Add a multipole.
   *
   * @param charge         The charge.
   * @param dipole         The dipole.
   * @param quadrupole     The quadrupole.
   * @param axisType       The axis type.
   * @param zaxis          The z-axis.
   * @param xaxis          The x-axis.
   * @param yaxis          The y-axis.
   * @param thole          The Thole parameter.
   * @param pdamp          The damping factor.
   * @param polarizability The polarizability.
   */
  public void addMultipole(double charge, DoubleArray dipole, DoubleArray quadrupole, int axisType,
                           int zaxis, int xaxis, int yaxis, double thole, double pdamp, double polarizability) {
    OpenMM_AmoebaMultipoleForce_addMultipole(pointer, charge, dipole.getPointer(), quadrupole.getPointer(),
        axisType, zaxis, xaxis, yaxis, thole, pdamp, polarizability);
  }

  /**
   * Set the multipole parameters.
   *
   * @param index          The atom index.
   * @param charge         The charge.
   * @param dipoles        The dipole.
   * @param quadrupoles    The quadrupole.
   * @param axisType       The axis type.
   * @param zaxis          The z-axis.
   * @param xaxis          The x-axis.
   * @param yaxis          The y-axis.
   * @param thole          The Thole parameter.
   * @param pdamp          The damping factor.
   * @param polarizability The polarizability.
   */
  public void setMultipoleParameters(int index, double charge, DoubleArray dipoles, DoubleArray quadrupoles,
                                     int axisType, int zaxis, int xaxis, int yaxis,
                                     double thole, double pdamp, double polarizability) {
    OpenMM_AmoebaMultipoleForce_setMultipoleParameters(pointer, index, charge,
        dipoles.getPointer(), quadrupoles.getPointer(), axisType, zaxis, xaxis, yaxis, thole, pdamp, polarizability);
  }

  /**
   * Set the nonbonded method for the multipole force.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_AmoebaMultipoleForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the cutoff distance for the multipole force.
   *
   * @param cutoff The cutoff distance.
   */
  public void setCutoffDistance(double cutoff) {
    OpenMM_AmoebaMultipoleForce_setCutoffDistance(pointer, cutoff);
  }

  /**
   * Set the Ewald coefficient for the multipole force.
   *
   * @param aewald The Ewald coefficient.
   */
  public void setAEwald(double aewald) {
    OpenMM_AmoebaMultipoleForce_setAEwald(pointer, aewald);
  }

  /**
   * Set the Ewald error tolerance for the multipole force.
   *
   * @param ewaldTolerance The Ewald error tolerance.
   */
  public void setEwaldErrorTolerance(double ewaldTolerance) {
    OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(pointer, ewaldTolerance);
  }

  /**
   * Set the PME grid dimensions for the multipole force.
   *
   * @param gridDimensions The PME grid dimensions.
   */
  public void setPmeGridDimensions(IntArray gridDimensions) {
    OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(pointer, gridDimensions.getPointer());
  }

  //    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce, pme.getPolarEps());

  /**
   * Set the mutual induced target maximum number of iterations.
   *
   * @param iterations The mutual induced max iterations.
   */
  public void setMutualInducedMaxIterations(int iterations) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(pointer, iterations);
  }

  /**
   * Set the mutual induced target epsilon.
   *
   * @param epsilon The mutual induced target epsilon.
   */
  public void setMutualInducedTargetEpsilon(double epsilon) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(pointer, epsilon);
  }

  /**
   * Set the covalent map.
   *
   * @param i            The atom index.
   * @param covalentType The covalent type.
   * @param covalentMap  The covalent map.
   */
  public void setCovalentMap(int i, int covalentType, IntArray covalentMap) {
    OpenMM_AmoebaMultipoleForce_setCovalentMap(pointer, i, covalentType, covalentMap.getPointer());
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaMultipoleForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaMultipoleForce_destroy(pointer);
      pointer = null;
    }
  }

}
