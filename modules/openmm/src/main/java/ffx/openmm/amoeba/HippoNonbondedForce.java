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
import ffx.openmm.Force;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_addException;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getDPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getDPMEParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getInducedDipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getNumExceptions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getPMEParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_getSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setDPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setPMEParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * This class implements the HIPPO (Hydrogen-like Intermolecular Polarizable POtential)
 * nonbonded force field. HIPPO is an advanced polarizable force field that includes
 * charge penetration, charge transfer, and many-body polarization effects.
 */
public class HippoNonbondedForce extends Force {

  /**
   * Create a new HippoNonbondedForce.
   */
  public HippoNonbondedForce() {
    super(OpenMM_HippoNonbondedForce_create());
  }

  /**
   * Add an exception to the force.
   *
   * @param particle1  The index of the first particle.
   * @param particle2  The index of the second particle.
   * @param chargeProd The charge product for the exception.
   * @param sigmaEps   The sigma epsilon parameter for the exception.
   * @param epsilon    The epsilon parameter for the exception.
   * @param damping    The damping parameter for the exception.
   * @param c6         The C6 dispersion parameter for the exception.
   * @param c8         The C8 dispersion parameter for the exception.
   * @param replace    Whether to replace an existing exception.
   * @return The index of the exception that was added.
   */
  public int addException(int particle1, int particle2, double chargeProd, double sigmaEps,
                          double epsilon, double damping, double c6, double c8, int replace) {
    return OpenMM_HippoNonbondedForce_addException(pointer, particle1, particle2, chargeProd,
        sigmaEps, epsilon, damping, c6, c8, replace);
  }

  /**
   * Add a particle to the force.
   *
   * @param charge             The charge of the particle.
   * @param dipole             The dipole moment of the particle.
   * @param quadrupole         The quadrupole moment of the particle.
   * @param c6                 The C6 dispersion parameter.
   * @param c8                 The C8 dispersion parameter.
   * @param c10                The C10 dispersion parameter.
   * @param c12                The C12 dispersion parameter.
   * @param sigma              The sigma parameter.
   * @param epsilon            The epsilon parameter.
   * @param damping            The damping parameter.
   * @param polarizability     The polarizability.
   * @param thole              The Thole damping parameter.
   * @param covalentMap        The covalent map.
   * @param polarizationGroup  The polarization group.
   * @param axisType           The axis type.
   * @param multipoleFrameType The multipole frame type.
   * @return The index of the particle that was added.
   */
  public int addParticle(double charge, PointerByReference dipole, PointerByReference quadrupole,
                         double c6, double c8, double c10, double c12, double sigma, double epsilon,
                         double damping, double polarizability, double thole, int covalentMap,
                         int polarizationGroup, int axisType, int multipoleFrameType) {
    return OpenMM_HippoNonbondedForce_addParticle(pointer, charge, dipole, quadrupole, c6, c8, c10,
        c12, sigma, epsilon, damping, polarizability, thole,
        covalentMap, polarizationGroup, axisType,
        multipoleFrameType);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_HippoNonbondedForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_HippoNonbondedForce_getCutoffDistance(pointer);
  }

  /**
   * Get the DPME parameters.
   *
   * @param alpha The Ewald alpha parameter (output).
   * @param nx    The number of grid points along the x axis (output).
   * @param ny    The number of grid points along the y axis (output).
   * @param nz    The number of grid points along the z axis (output).
   */
  public void getDPMEParameters(DoubleByReference alpha, IntByReference nx, IntByReference ny, IntByReference nz) {
    OpenMM_HippoNonbondedForce_getDPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Get the DPME parameters.
   *
   * @param alpha The Ewald alpha parameter (output).
   * @param nx    The number of grid points along the x axis (output).
   * @param ny    The number of grid points along the y axis (output).
   * @param nz    The number of grid points along the z axis (output).
   */
  public void getDPMEParameters(DoubleBuffer alpha, IntBuffer nx, IntBuffer ny, IntBuffer nz) {
    OpenMM_HippoNonbondedForce_getDPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Get the DPME parameters in context.
   *
   * @param context The context.
   * @param alpha   The Ewald alpha parameter (output).
   * @param nx      The number of grid points along the x axis (output).
   * @param ny      The number of grid points along the y axis (output).
   * @param nz      The number of grid points along the z axis (output).
   */
  public void getDPMEParametersInContext(Context context, DoubleByReference alpha, IntByReference nx,
                                         IntByReference ny, IntByReference nz) {
    OpenMM_HippoNonbondedForce_getDPMEParametersInContext(pointer, context.getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Get the DPME parameters in context.
   *
   * @param context The context.
   * @param alpha   The Ewald alpha parameter (output).
   * @param nx      The number of grid points along the x axis (output).
   * @param ny      The number of grid points along the y axis (output).
   * @param nz      The number of grid points along the z axis (output).
   */
  public void getDPMEParametersInContext(Context context, DoubleBuffer alpha, IntBuffer nx,
                                         IntBuffer ny, IntBuffer nz) {
    OpenMM_HippoNonbondedForce_getDPMEParametersInContext(pointer, context.getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Get the Ewald error tolerance.
   *
   * @return The Ewald error tolerance.
   */
  public double getEwaldErrorTolerance() {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(pointer);
  }

  /**
   * Get the parameters for an exception.
   *
   * @param index      The index of the exception.
   * @param particle1  The index of the first particle (output).
   * @param particle2  The index of the second particle (output).
   * @param chargeProd The charge product for the exception (output).
   * @param sigmaEps   The sigma epsilon parameter for the exception (output).
   * @param epsilon    The epsilon parameter for the exception (output).
   * @param damping    The damping parameter for the exception (output).
   * @param c6         The C6 dispersion parameter for the exception (output).
   * @param c8         The C8 dispersion parameter for the exception (output).
   */
  public void getExceptionParameters(int index, IntByReference particle1, IntByReference particle2,
                                     DoubleByReference chargeProd, DoubleByReference sigmaEps,
                                     DoubleByReference epsilon, DoubleByReference damping,
                                     DoubleByReference c6, DoubleByReference c8) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(pointer, index, particle1, particle2,
        chargeProd, sigmaEps, epsilon, damping, c6, c8);
  }

  /**
   * Get the parameters for an exception.
   *
   * @param index      The index of the exception.
   * @param particle1  The index of the first particle (output).
   * @param particle2  The index of the second particle (output).
   * @param chargeProd The charge product for the exception (output).
   * @param sigmaEps   The sigma epsilon parameter for the exception (output).
   * @param epsilon    The epsilon parameter for the exception (output).
   * @param damping    The damping parameter for the exception (output).
   * @param c6         The C6 dispersion parameter for the exception (output).
   * @param c8         The C8 dispersion parameter for the exception (output).
   */
  public void getExceptionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                     DoubleBuffer chargeProd, DoubleBuffer sigmaEps,
                                     DoubleBuffer epsilon, DoubleBuffer damping,
                                     DoubleBuffer c6, DoubleBuffer c8) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(pointer, index, particle1, particle2,
        chargeProd, sigmaEps, epsilon, damping, c6, c8);
  }

  /**
   * Get the extrapolation coefficients.
   *
   * @return The extrapolation coefficients.
   */
  public PointerByReference getExtrapolationCoefficients() {
    return OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(pointer);
  }

  /**
   * Get the induced dipoles.
   *
   * @param context The context.
   * @param dipoles The induced dipoles (output).
   */
  public void getInducedDipoles(Context context, PointerByReference dipoles) {
    OpenMM_HippoNonbondedForce_getInducedDipoles(pointer, context.getPointer(), dipoles);
  }

  /**
   * Get the lab frame permanent dipoles.
   *
   * @param context The context.
   * @param dipoles The lab frame permanent dipoles (output).
   */
  public void getLabFramePermanentDipoles(Context context, PointerByReference dipoles) {
    OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(pointer, context.getPointer(), dipoles);
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_HippoNonbondedForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of exceptions.
   *
   * @return The number of exceptions.
   */
  public int getNumExceptions() {
    return OpenMM_HippoNonbondedForce_getNumExceptions(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_HippoNonbondedForce_getNumParticles(pointer);
  }

  /**
   * Get the PME parameters.
   *
   * @param alpha The Ewald alpha parameter (output).
   * @param nx    The number of grid points along the x axis (output).
   * @param ny    The number of grid points along the y axis (output).
   * @param nz    The number of grid points along the z axis (output).
   */
  public void getPMEParameters(DoubleByReference alpha, IntByReference nx, IntByReference ny, IntByReference nz) {
    OpenMM_HippoNonbondedForce_getPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Get the PME parameters.
   *
   * @param alpha The Ewald alpha parameter (output).
   * @param nx    The number of grid points along the x axis (output).
   * @param ny    The number of grid points along the y axis (output).
   * @param nz    The number of grid points along the z axis (output).
   */
  public void getPMEParameters(DoubleBuffer alpha, IntBuffer nx, IntBuffer ny, IntBuffer nz) {
    OpenMM_HippoNonbondedForce_getPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Get the PME parameters in context.
   *
   * @param context The context.
   * @param alpha   The Ewald alpha parameter (output).
   * @param nx      The number of grid points along the x axis (output).
   * @param ny      The number of grid points along the y axis (output).
   * @param nz      The number of grid points along the z axis (output).
   */
  public void getPMEParametersInContext(Context context, DoubleByReference alpha, IntByReference nx,
                                        IntByReference ny, IntByReference nz) {
    OpenMM_HippoNonbondedForce_getPMEParametersInContext(pointer, context.getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Get the PME parameters in context.
   *
   * @param context The context.
   * @param alpha   The Ewald alpha parameter (output).
   * @param nx      The number of grid points along the x axis (output).
   * @param ny      The number of grid points along the y axis (output).
   * @param nz      The number of grid points along the z axis (output).
   */
  public void getPMEParametersInContext(Context context, DoubleBuffer alpha, IntBuffer nx,
                                        IntBuffer ny, IntBuffer nz) {
    OpenMM_HippoNonbondedForce_getPMEParametersInContext(pointer, context.getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param charge             The charge of the particle (output).
   * @param dipole             The dipole moment of the particle (output).
   * @param quadrupole         The quadrupole moment of the particle (output).
   * @param c6                 The C6 dispersion parameter (output).
   * @param c8                 The C8 dispersion parameter (output).
   * @param c10                The C10 dispersion parameter (output).
   * @param c12                The C12 dispersion parameter (output).
   * @param sigma              The sigma parameter (output).
   * @param epsilon            The epsilon parameter (output).
   * @param damping            The damping parameter (output).
   * @param polarizability     The polarizability (output).
   * @param thole              The Thole damping parameter (output).
   * @param covalentMap        The covalent map (output).
   * @param polarizationGroup  The polarization group (output).
   * @param axisType           The axis type (output).
   * @param multipoleFrameType The multipole frame type (output).
   */
  public void getParticleParameters(int index, DoubleByReference charge, PointerByReference dipole,
                                    PointerByReference quadrupole, DoubleByReference c6,
                                    DoubleByReference c8, DoubleByReference c10, DoubleByReference c12,
                                    DoubleByReference sigma, DoubleByReference epsilon,
                                    DoubleByReference damping, DoubleByReference polarizability,
                                    DoubleByReference thole, IntByReference covalentMap,
                                    IntByReference polarizationGroup, IntByReference axisType,
                                    IntByReference multipoleFrameType) {
    OpenMM_HippoNonbondedForce_getParticleParameters(pointer, index, charge, dipole, quadrupole,
        c6, c8, c10, c12, sigma, epsilon, damping,
        polarizability, thole, covalentMap,
        polarizationGroup, axisType, multipoleFrameType);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param charge             The charge of the particle (output).
   * @param dipole             The dipole moment of the particle (output).
   * @param quadrupole         The quadrupole moment of the particle (output).
   * @param c6                 The C6 dispersion parameter (output).
   * @param c8                 The C8 dispersion parameter (output).
   * @param c10                The C10 dispersion parameter (output).
   * @param c12                The C12 dispersion parameter (output).
   * @param sigma              The sigma parameter (output).
   * @param epsilon            The epsilon parameter (output).
   * @param damping            The damping parameter (output).
   * @param polarizability     The polarizability (output).
   * @param thole              The Thole damping parameter (output).
   * @param covalentMap        The covalent map (output).
   * @param polarizationGroup  The polarization group (output).
   * @param axisType           The axis type (output).
   * @param multipoleFrameType The multipole frame type (output).
   */
  public void getParticleParameters(int index, DoubleBuffer charge, PointerByReference dipole,
                                    PointerByReference quadrupole, DoubleBuffer c6,
                                    DoubleBuffer c8, DoubleBuffer c10, DoubleBuffer c12,
                                    DoubleBuffer sigma, DoubleBuffer epsilon,
                                    DoubleBuffer damping, DoubleBuffer polarizability,
                                    DoubleBuffer thole, IntBuffer covalentMap,
                                    IntBuffer polarizationGroup, IntBuffer axisType,
                                    IntBuffer multipoleFrameType) {
    OpenMM_HippoNonbondedForce_getParticleParameters(pointer, index, charge, dipole, quadrupole,
        c6, c8, c10, c12, sigma, epsilon, damping,
        polarizability, thole, covalentMap,
        polarizationGroup, axisType, multipoleFrameType);
  }

  /**
   * Get the switching distance.
   *
   * @return The switching distance, measured in nm.
   */
  public double getSwitchingDistance() {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance(pointer);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_HippoNonbondedForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the DPME parameters.
   *
   * @param alpha The Ewald alpha parameter.
   * @param nx    The number of grid points along the x axis.
   * @param ny    The number of grid points along the y axis.
   * @param nz    The number of grid points along the z axis.
   */
  public void setDPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMM_HippoNonbondedForce_setDPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Set the Ewald error tolerance.
   *
   * @param tolerance The Ewald error tolerance.
   */
  public void setEwaldErrorTolerance(double tolerance) {
    OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(pointer, tolerance);
  }

  /**
   * Set the parameters for an exception.
   *
   * @param index      The index of the exception.
   * @param particle1  The index of the first particle.
   * @param particle2  The index of the second particle.
   * @param chargeProd The charge product for the exception.
   * @param sigmaEps   The sigma epsilon parameter for the exception.
   * @param epsilon    The epsilon parameter for the exception.
   * @param damping    The damping parameter for the exception.
   * @param c6         The C6 dispersion parameter for the exception.
   * @param c8         The C8 dispersion parameter for the exception.
   */
  public void setExceptionParameters(int index, int particle1, int particle2, double chargeProd,
                                     double sigmaEps, double epsilon, double damping, double c6, double c8) {
    OpenMM_HippoNonbondedForce_setExceptionParameters(pointer, index, particle1, particle2,
        chargeProd, sigmaEps, epsilon, damping, c6, c8);
  }

  /**
   * Set the extrapolation coefficients.
   *
   * @param coefficients The extrapolation coefficients.
   */
  public void setExtrapolationCoefficients(PointerByReference coefficients) {
    OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(pointer, coefficients);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_HippoNonbondedForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the PME parameters.
   *
   * @param alpha The Ewald alpha parameter.
   * @param nx    The number of grid points along the x axis.
   * @param ny    The number of grid points along the y axis.
   * @param nz    The number of grid points along the z axis.
   */
  public void setPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMM_HippoNonbondedForce_setPMEParameters(pointer, alpha, nx, ny, nz);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param charge             The charge of the particle.
   * @param dipole             The dipole moment of the particle.
   * @param quadrupole         The quadrupole moment of the particle.
   * @param c6                 The C6 dispersion parameter.
   * @param c8                 The C8 dispersion parameter.
   * @param c10                The C10 dispersion parameter.
   * @param c12                The C12 dispersion parameter.
   * @param sigma              The sigma parameter.
   * @param epsilon            The epsilon parameter.
   * @param damping            The damping parameter.
   * @param polarizability     The polarizability.
   * @param thole              The Thole damping parameter.
   * @param covalentMap        The covalent map.
   * @param polarizationGroup  The polarization group.
   * @param axisType           The axis type.
   * @param multipoleFrameType The multipole frame type.
   */
  public void setParticleParameters(int index, double charge, PointerByReference dipole,
                                    PointerByReference quadrupole, double c6, double c8, double c10,
                                    double c12, double sigma, double epsilon, double damping,
                                    double polarizability, double thole, int covalentMap,
                                    int polarizationGroup, int axisType, int multipoleFrameType) {
    OpenMM_HippoNonbondedForce_setParticleParameters(pointer, index, charge, dipole, quadrupole,
        c6, c8, c10, c12, sigma, epsilon, damping,
        polarizability, thole, covalentMap,
        polarizationGroup, axisType, multipoleFrameType);
  }

  /**
   * Set the switching distance.
   *
   * @param distance The switching distance, measured in nm.
   */
  public void setSwitchingDistance(double distance) {
    OpenMM_HippoNonbondedForce_setSwitchingDistance(pointer, distance);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_HippoNonbondedForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}