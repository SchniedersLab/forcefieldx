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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.Context;
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_addMultipole;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getAEwald;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getCovalentMap;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getCovalentMaps;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getElectrostaticPotential;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getInducedDipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getMultipoleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getNumMultipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getPMEParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getPmeGridDimensions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getPolarizationType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_getTotalDipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setAEwald;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCovalentMap;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMultipoleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPolarizationType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * This class implements the Amoeba multipole interaction.
 * <p>
 * To use it, create an AmoebaMultipoleForce object then call addMultipole() once for each atom.  After
 * an entry has been added, you can modify its force field parameters by calling setMultipoleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */
public class MultipoleForce extends Force {

  public MultipoleForce() {
    super(OpenMM_AmoebaMultipoleForce_create());
  }

  /**
   * Add multipole-related info for a particle
   *
   * @param charge              the particle's charge
   * @param molecularDipole     the particle's molecular dipole (vector of size 3)
   * @param molecularQuadrupole the particle's molecular quadrupole (vector of size 9)
   * @param axisType            the particle's axis type
   * @param multipoleAtomZ      index of first atom used in constructing lab<->molecular frames
   * @param multipoleAtomX      index of second atom used in constructing lab<->molecular frames
   * @param multipoleAtomY      index of second atom used in constructing lab<->molecular frames
   * @param thole               Thole parameter
   * @param dampingFactor       dampingFactor parameter
   * @param polarity            polarity parameter
   * @return the index of the particle that was added
   */
  public int addMultipole(double charge, DoubleArray molecularDipole, DoubleArray molecularQuadrupole, int axisType,
                          int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) {
    return OpenMM_AmoebaMultipoleForce_addMultipole(pointer, charge, molecularDipole.getPointer(), molecularQuadrupole.getPointer(),
        axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaMultipoleForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
   * based on the Ewald error tolerance.
   *
   * @return the Ewald alpha parameter
   * @deprecated This method exists only for backward compatibility.  Use getPMEParameters() instead.
   */
  public double getAEwald() {
    return OpenMM_AmoebaMultipoleForce_getAEwald(pointer);
  }

  /**
   * Get the covalent map for a given atom index and covalent type.
   *
   * @param i            The atom index.
   * @param covalentType The covalent type.
   * @return An IntArray representing the covalent map for the specified atom index and covalent type.
   */
  public IntArray getCovalentMap(int i, int covalentType) {
    IntArray covalentMap = new IntArray(0);
    if (pointer != null) {
      OpenMM_AmoebaMultipoleForce_getCovalentMap(pointer, i, covalentType, covalentMap.getPointer());
    }
    return covalentMap;
  }

  /**
   * Get all covalent maps for a given atom index.
   *
   * @param i The atom index.
   * @return An IntArray containing all covalent maps for the specified atom.
   */
  public IntArray getCovalentMaps(int i) {
    IntArray covalentMaps = new IntArray(0);
    if (pointer != null) {
      OpenMM_AmoebaMultipoleForce_getCovalentMaps(pointer, i, covalentMaps.getPointer());
    }
    return covalentMaps;
  }

  /**
   * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
   * is NoCutoff, this value will have no effect.
   *
   * @return the cutoff distance, measured in nm
   */
  public double getCutoffDistance() {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance(pointer);
  }

  /**
   * Get the electrostatic potential at specified points.
   *
   * @param context The OpenMM context.
   * @param points  The points at which to calculate the potential.
   * @return A DoubleArray containing the electrostatic potential at each point.
   */
  public DoubleArray getElectrostaticPotential(Context context, DoubleArray points) {
    DoubleArray potential = new DoubleArray(0);
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(pointer, context.getPointer(),
          points.getPointer(), potential.getPointer());
    }
    return potential;
  }

  /**
   * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
   * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
   * parameter so that the average error level will be less than the tolerance.  There is not a
   * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
   * <p>
   * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
   */
  public double getEwaldErrorTolerance() {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(pointer);
  }

  /**
   * Get the extrapolation coefficients.
   *
   * @return A PointerByReference to the extrapolation coefficients.
   */
  public PointerByReference getExtrapolationCoefficients() {
    return OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(pointer);
  }

  /**
   * Get the induced dipoles.
   *
   * @param context        The OpenMM context.
   * @param inducedDipoles The induced dipoles (output).
   */
  public void getInducedDipoles(Context context, DoubleArray inducedDipoles) {
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getInducedDipoles(pointer, context.getPointer(), inducedDipoles.getPointer());
    }
  }

  /**
   * Get the lab frame permanent dipoles.
   *
   * @param context The OpenMM context.
   * @param dipoles The lab frame permanent dipoles (output).
   */
  public void getLabFramePermanentDipoles(Context context, DoubleArray dipoles) {
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(pointer, context.getPointer(), dipoles.getPointer());
    }
  }

  /**
   * Get the multipole parameters for a particle.
   *
   * @param index               the index of the atom for which to get parameters
   * @param charge              the particle's charge
   * @param molecularDipole     the particle's molecular dipole (vector of size 3)
   * @param molecularQuadrupole the particle's molecular quadrupole (vector of size 9)
   * @param axisType            the particle's axis type
   * @param multipoleAtomZ      index of first atom used in constructing lab<->molecular frames
   * @param multipoleAtomX      index of second atom used in constructing lab<->molecular frames
   * @param multipoleAtomY      index of second atom used in constructing lab<->molecular frames
   * @param thole               Thole parameter
   * @param dampingFactor       dampingFactor parameter
   * @param polarity            polarity parameter
   */
  public void getMultipoleParameters(int index, DoubleByReference charge, PointerByReference molecularDipole,
                                     PointerByReference molecularQuadrupole, IntByReference axisType,
                                     IntByReference multipoleAtomZ, IntByReference multipoleAtomX, IntByReference multipoleAtomY,
                                     DoubleByReference thole, DoubleByReference dampingFactor,
                                     DoubleByReference polarity) {
    OpenMM_AmoebaMultipoleForce_getMultipoleParameters(pointer, index, charge, molecularDipole, molecularQuadrupole,
        axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
  }

  /**
   * Get the mutual induced max iterations.
   *
   * @return The mutual induced max iterations.
   */
  public int getMutualInducedMaxIterations() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(pointer);
  }

  /**
   * Get the mutual induced target epsilon.
   *
   * @return The mutual induced target epsilon.
   */
  public double getMutualInducedTargetEpsilon() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(pointer);
  }

  /**
   * Get the nonbonded method for the multipole force.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_AmoebaMultipoleForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of multipoles.
   *
   * @return The number of multipoles.
   */
  public int getNumMultipoles() {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles(pointer);
  }

  /**
   * Get the parameters to use for PME calculations.  If alpha is 0 (the default), these parameters are
   * ignored and instead their values are chosen based on the Ewald error tolerance.
   *
   * @param alpha the separation parameter
   * @param nx    the number of grid points along the X axis
   * @param ny    the number of grid points along the Y axis
   * @param nz    the number of grid points along the Z axis
   */
  public void getPMEParameters(DoubleByReference alpha, IntByReference nx, IntByReference ny, IntByReference nz) {
    OpenMM_AmoebaMultipoleForce_getPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Get the PME parameters in the context.
   *
   * @param context The OpenMM context.
   * @param alpha   The Ewald alpha parameter (output).
   * @param nx      The PME grid size in x (output).
   * @param ny      The PME grid size in y (output).
   * @param nz      The PME grid size in z (output).
   */
  public void getPMEParametersInContext(Context context, DoubleByReference alpha,
                                        IntByReference nx, IntByReference ny, IntByReference nz) {
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(pointer, context.getPointer(), alpha, nx, ny, nz);
    }
  }

  /**
   * Get the PME grid dimensions.
   *
   * @return An IntArray containing the PME grid dimensions.
   */
  public IntArray getPmeGridDimensions() {
    IntArray gridDimensions = new IntArray(0);
    if (pointer != null) {
      OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(pointer, gridDimensions.getPointer());
    }
    return gridDimensions;
  }

  /**
   * Get the PME B-spline order.
   *
   * @return The PME B-spline order.
   */
  public int getPmeBSplineOrder() {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(pointer);
  }

  /**
   * Get the polarization type.
   *
   * @return The polarization type.
   */
  public int getPolarizationType() {
    return OpenMM_AmoebaMultipoleForce_getPolarizationType(pointer);
  }

  /**
   * Get the system multipole moments.
   *
   * @param context The OpenMM context.
   * @param moments The system multipole moments (output).
   */
  public void getSystemMultipoleMoments(Context context, DoubleArray moments) {
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(pointer, context.getPointer(), moments.getPointer());
    }
  }

  /**
   * Get the total dipoles.
   *
   * @param context The OpenMM context.
   * @param dipoles The total dipoles (output).
   */
  public void getTotalDipoles(Context context, DoubleArray dipoles) {
    if (context.hasContextPointer() && pointer != null) {
      OpenMM_AmoebaMultipoleForce_getTotalDipoles(pointer, context.getPointer(), dipoles.getPointer());
    }
  }

  /**
   * Set the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
   * based on the Ewald error tolerance.
   *
   * @param aewald alpha parameter
   * @deprecated This method exists only for backward compatibility.  Use setPMEParameters() instead.
   */
  public void setAEwald(double aewald) {
    OpenMM_AmoebaMultipoleForce_setAEwald(pointer, aewald);
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
   * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
   * is NoCutoff, this value will have no effect.
   *
   * @param cutoff the cutoff distance, measured in nm
   */
  public void setCutoffDistance(double cutoff) {
    OpenMM_AmoebaMultipoleForce_setCutoffDistance(pointer, cutoff);
  }

  /**
   * Set the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
   * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
   * parameter so that the average error level will be less than the tolerance.  There is not a
   * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
   * <p>
   * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
   *
   * @param ewaldTolerance the error tolerance
   */
  public void setEwaldErrorTolerance(double ewaldTolerance) {
    OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(pointer, ewaldTolerance);
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
   * Set the multipole parameters for a particle.
   *
   * @param index               the index of the atom for which to set parameters
   * @param charge              the particle's charge
   * @param molecularDipole     the particle's molecular dipole (vector of size 3)
   * @param molecularQuadrupole the particle's molecular quadrupole (vector of size 9)
   * @param axisType            the particle's axis type
   * @param multipoleAtomZ      index of first atom used in constructing lab<->molecular frames
   * @param multipoleAtomX      index of second atom used in constructing lab<->molecular frames
   * @param multipoleAtomY      index of second atom used in constructing lab<->molecular frames
   * @param thole               thole parameter
   * @param dampingFactor       damping factor parameter
   * @param polarity            polarity parameter
   */
  public void setMultipoleParameters(int index, double charge, DoubleArray molecularDipole, DoubleArray molecularQuadrupole,
                                     int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY,
                                     double thole, double dampingFactor, double polarity) {
    OpenMM_AmoebaMultipoleForce_setMultipoleParameters(pointer, index, charge,
        molecularDipole.getPointer(), molecularQuadrupole.getPointer(), axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
  }

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
   * Set the nonbonded method for the multipole force.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_AmoebaMultipoleForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters to use for PME calculations.  If alpha is 0 (the default), these parameters are
   * ignored and instead their values are chosen based on the Ewald error tolerance.
   *
   * @param alpha the separation parameter
   * @param nx    the number of grid points along the X axis
   * @param ny    the number of grid points along the Y axis
   * @param nz    the number of grid points along the Z axis
   */
  public void setPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMM_AmoebaMultipoleForce_setPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Set the PME grid dimensions for the multipole force.
   *
   * @param gridDimensions The PME grid dimensions.
   */
  public void setPmeGridDimensions(IntArray gridDimensions) {
    OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(pointer, gridDimensions.getPointer());
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
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}