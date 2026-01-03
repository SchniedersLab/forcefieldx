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
 * This class implements all nonbonded interactions in the HIPPO force field: electrostatics,
 * induction, charge transfer, dispersion, and repulsion. Although some of these are
 * conceptually distinct, they share parameters in common and are most efficiently computed
 * together. For example, the same multipole definitions are used for both electrostatics
 * and Pauli repulsion. Therefore, all of them are computed by a single Force object.
 * <p>
 * To use it, create a HippoNonbondedForce object, then call addParticle() once for each particle.
 * After an entry has been added, you can modify its force field parameters by calling setParticleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 * <p>
 * You also can specify "exceptions", particular pairs of particles whose interactions should be
 * reduced or completely omitted. Call addException() to define exceptions.
 */
public class HippoNonbondedForce extends Force {

  /**
   * Create a new HippoNonbondedForce.
   */
  public HippoNonbondedForce() {
    super(OpenMM_HippoNonbondedForce_create());
  }

  /**
   * Add an interaction to the list of exceptions that should be calculated differently from other interactions.
   * If all scale factors are set to 0, this will cause the interaction to be completely omitted from
   * force and energy calculations.
   *
   * @param particle1               the index of the first particle involved in the interaction
   * @param particle2               the index of the second particle involved in the interaction
   * @param multipoleMultipoleScale the factor by which to scale the Coulomb interaction between fixed multipoles
   * @param dipoleMultipoleScale    the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
   * @param dipoleDipoleScale       the factor by which to scale the Coulomb interaction between induced dipoles
   * @param dispersionScale         the factor by which to scale the dispersion interaction
   * @param repulsionScale          the factor by which to scale the Pauli repulsion
   * @param chargeTransferScale     the factor by which to scale the charge transfer interaction
   * @param replace                 determines the behavior if there is already an exception for the same two particles. If true, the existing one is replaced. If false, an exception is thrown.
   * @return the index of the exception that was added
   */
  public int addException(int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale,
                          double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, int replace) {
    return OpenMM_HippoNonbondedForce_addException(pointer, particle1, particle2, multipoleMultipoleScale,
        dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
  }

  /**
   * Add the nonbonded force parameters for a particle. This should be called once for each particle
   * in the System. When it is called for the i'th time, it specifies the parameters for the i'th particle.
   *
   * @param charge         the particle's charge
   * @param dipole         the particle's molecular dipole (vector of size 3)
   * @param quadrupole     the particle's molecular quadrupole (vector of size 9)
   * @param coreCharge     the charge of the atomic core
   * @param alpha          controls the width of the particle's electron density
   * @param epsilon        sets the magnitude of charge transfer
   * @param damping        sets the length scale for charge transfer
   * @param c6             the coefficient of the dispersion interaction
   * @param pauliK         the coefficient of the Pauli repulsion interaction
   * @param pauliQ         the charge used in computing the Pauli repulsion interaction
   * @param pauliAlpha     the width of the particle's electron density for computing the Pauli repulsion interaction
   * @param polarizability atomic polarizability
   * @param axisType       the particle's axis type
   * @param multipoleAtomZ index of first atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomX index of second atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomY index of third atom used in defining the local coordinate system for multipoles
   * @return the index of the particle that was added
   */
  public int addParticle(double charge, PointerByReference dipole, PointerByReference quadrupole, double coreCharge,
                         double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                         double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    return OpenMM_HippoNonbondedForce_addParticle(pointer, charge, dipole, quadrupole, coreCharge,
        alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
        polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
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
   * Get the scale factors for an interaction that should be calculated differently from others.
   *
   * @param index                   the index of the interaction for which to get parameters
   * @param particle1               the index of the first particle involved in the interaction
   * @param particle2               the index of the second particle involved in the interaction
   * @param multipoleMultipoleScale the factor by which to scale the Coulomb interaction between fixed multipoles
   * @param dipoleMultipoleScale    the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
   * @param dipoleDipoleScale       the factor by which to scale the Coulomb interaction between induced dipoles
   * @param dispersionScale         the factor by which to scale the dispersion interaction
   * @param repulsionScale          the factor by which to scale the Pauli repulsion
   * @param chargeTransferScale     the factor by which to scale the charge transfer interaction
   */
  public void getExceptionParameters(int index, IntByReference particle1, IntByReference particle2,
                                     DoubleByReference multipoleMultipoleScale, DoubleByReference dipoleMultipoleScale,
                                     DoubleByReference dipoleDipoleScale, DoubleByReference dispersionScale,
                                     DoubleByReference repulsionScale, DoubleByReference chargeTransferScale) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(pointer, index, particle1, particle2,
        multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
  }

  /**
   * Get the scale factors for an interaction that should be calculated differently from others.
   *
   * @param index                   the index of the interaction for which to get parameters
   * @param particle1               the index of the first particle involved in the interaction
   * @param particle2               the index of the second particle involved in the interaction
   * @param multipoleMultipoleScale the factor by which to scale the Coulomb interaction between fixed multipoles
   * @param dipoleMultipoleScale    the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
   * @param dipoleDipoleScale       the factor by which to scale the Coulomb interaction between induced dipoles
   * @param dispersionScale         the factor by which to scale the dispersion interaction
   * @param repulsionScale          the factor by which to scale the Pauli repulsion
   * @param chargeTransferScale     the factor by which to scale the charge transfer interaction
   */
  public void getExceptionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                     DoubleBuffer multipoleMultipoleScale, DoubleBuffer dipoleMultipoleScale,
                                     DoubleBuffer dipoleDipoleScale, DoubleBuffer dispersionScale,
                                     DoubleBuffer repulsionScale, DoubleBuffer chargeTransferScale) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(pointer, index, particle1, particle2,
        multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
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
   * Get the nonbonded force parameters for a particle.
   *
   * @param index          the index of the particle for which to get parameters
   * @param charge         the particle's charge
   * @param dipole         the particle's molecular dipole (vector of size 3)
   * @param quadrupole     the particle's molecular quadrupole (vector of size 9)
   * @param coreCharge     the charge of the atomic core
   * @param alpha          controls the width of the particle's electron density
   * @param epsilon        sets the magnitude of charge transfer
   * @param damping        sets the length scale for charge transfer
   * @param c6             the coefficient of the dispersion interaction
   * @param pauliK         the coefficient of the Pauli repulsion interaction
   * @param pauliQ         the charge used in computing the Pauli repulsion interaction
   * @param pauliAlpha     the width of the particle's electron density for computing the Pauli repulsion interaction
   * @param polarizability atomic polarizability
   * @param axisType       the particle's axis type
   * @param multipoleAtomZ index of first atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomX index of second atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomY index of third atom used in defining the local coordinate system for multipoles
   */
  public void getParticleParameters(int index, DoubleByReference charge, PointerByReference dipole,
                                    PointerByReference quadrupole, DoubleByReference coreCharge,
                                    DoubleByReference alpha, DoubleByReference epsilon,
                                    DoubleByReference damping, DoubleByReference c6,
                                    DoubleByReference pauliK, DoubleByReference pauliQ, DoubleByReference pauliAlpha,
                                    DoubleByReference polarizability, IntByReference axisType,
                                    IntByReference multipoleAtomZ, IntByReference multipoleAtomX, IntByReference multipoleAtomY) {
    OpenMM_HippoNonbondedForce_getParticleParameters(pointer, index, charge, dipole, quadrupole,
        coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
        polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
  }

  /**
   * Get the nonbonded force parameters for a particle.
   *
   * @param index          the index of the particle for which to get parameters
   * @param charge         the particle's charge
   * @param dipole         the particle's molecular dipole (vector of size 3)
   * @param quadrupole     the particle's molecular quadrupole (vector of size 9)
   * @param coreCharge     the charge of the atomic core
   * @param alpha          controls the width of the particle's electron density
   * @param epsilon        sets the magnitude of charge transfer
   * @param damping        sets the length scale for charge transfer
   * @param c6             the coefficient of the dispersion interaction
   * @param pauliK         the coefficient of the Pauli repulsion interaction
   * @param pauliQ         the charge used in computing the Pauli repulsion interaction
   * @param pauliAlpha     the width of the particle's electron density for computing the Pauli repulsion interaction
   * @param polarizability atomic polarizability
   * @param axisType       the particle's axis type
   * @param multipoleAtomZ index of first atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomX index of second atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomY index of third atom used in defining the local coordinate system for multipoles
   */
  public void getParticleParameters(int index, DoubleBuffer charge, PointerByReference dipole,
                                    PointerByReference quadrupole, DoubleBuffer coreCharge,
                                    DoubleBuffer alpha, DoubleBuffer epsilon,
                                    DoubleBuffer damping, DoubleBuffer c6,
                                    DoubleBuffer pauliK, DoubleBuffer pauliQ, DoubleBuffer pauliAlpha,
                                    DoubleBuffer polarizability, IntBuffer axisType,
                                    IntBuffer multipoleAtomZ, IntBuffer multipoleAtomX, IntBuffer multipoleAtomY) {
    OpenMM_HippoNonbondedForce_getParticleParameters(pointer, index, charge, dipole, quadrupole,
        coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
        polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
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
   * Set the scale factors for an interaction that should be calculated differently from others.
   *
   * @param index                   the index of the interaction for which to set parameters
   * @param particle1               the index of the first particle involved in the interaction
   * @param particle2               the index of the second particle involved in the interaction
   * @param multipoleMultipoleScale the factor by which to scale the Coulomb interaction between fixed multipoles
   * @param dipoleMultipoleScale    the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
   * @param dipoleDipoleScale       the factor by which to scale the Coulomb interaction between induced dipoles
   * @param dispersionScale         the factor by which to scale the dispersion interaction
   * @param repulsionScale          the factor by which to scale the Pauli repulsion
   * @param chargeTransferScale     the factor by which to scale the charge transfer interaction
   */
  public void setExceptionParameters(int index, int particle1, int particle2, double multipoleMultipoleScale,
                                     double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale,
                                     double repulsionScale, double chargeTransferScale) {
    OpenMM_HippoNonbondedForce_setExceptionParameters(pointer, index, particle1, particle2,
        multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
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
   * Set the nonbonded force parameters for a particle.
   *
   * @param index          the index of the particle for which to set parameters
   * @param charge         the particle's charge
   * @param dipole         the particle's molecular dipole (vector of size 3)
   * @param quadrupole     the particle's molecular quadrupole (vector of size 9)
   * @param coreCharge     the charge of the atomic core
   * @param alpha          controls the width of the particle's electron density
   * @param epsilon        sets the magnitude of charge transfer
   * @param damping        sets the length scale for charge transfer
   * @param c6             the coefficient of the dispersion interaction
   * @param pauliK         the coefficient of the Pauli repulsion interaction
   * @param pauliQ         the charge used in computing the Pauli repulsion interaction
   * @param pauliAlpha     the width of the particle's electron density for computing the Pauli repulsion interaction
   * @param polarizability atomic polarizability
   * @param axisType       the particle's axis type
   * @param multipoleAtomZ index of first atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomX index of second atom used in defining the local coordinate system for multipoles
   * @param multipoleAtomY index of third atom used in defining the local coordinate system for multipoles
   */
  public void setParticleParameters(int index, double charge, PointerByReference dipole,
                                    PointerByReference quadrupole, double coreCharge, double alpha, double epsilon,
                                    double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                                    double polarizability, int axisType, int multipoleAtomZ,
                                    int multipoleAtomX, int multipoleAtomY) {
    OpenMM_HippoNonbondedForce_setParticleParameters(pointer, index, charge, dipole, quadrupole,
        coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
        polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
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