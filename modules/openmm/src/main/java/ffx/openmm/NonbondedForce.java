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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_createExceptionsFromBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getNumExceptions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setPMEParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_updateParametersInContext;

/**
 * Nonbonded Force.
 */
public class NonbondedForce extends Force {

  public NonbondedForce() {
    pointer = OpenMM_NonbondedForce_create();
  }

  /**
   * Add a particle.
   *
   * @param charge The atomic charge.
   * @param sigma  The vdW sigma.
   * @param eps    The vdW eps.
   */
  public void addParticle(double charge, double sigma, double eps) {
    OpenMM_NonbondedForce_addParticle(pointer, charge, sigma, eps);
  }

  /**
   * Set the particle parameters.
   *
   * @param index  The particle index.
   * @param charge The atomic charge.
   * @param sigma  The vdW sigma.
   * @param eps    The vdW eps.
   */
  public void setParticleParameters(int index, double charge, double sigma, double eps) {
    OpenMM_NonbondedForce_setParticleParameters(pointer, index, charge, sigma, eps);
  }

  /**
   * Get the particle parameters.
   *
   * @param index  The particle index.
   * @param charge The atomic charge.
   * @param sigma  The vdW sigma.
   * @param eps    The vdW eps.
   */
  public void getParticleParameters(int index, DoubleByReference charge, DoubleByReference sigma, DoubleByReference eps) {
    OpenMM_NonbondedForce_getParticleParameters(pointer, index, charge, sigma, eps);
  }

  /**
   * Create exceptions from bonds.
   *
   * @param bondArray      The bond array.
   * @param coulomb14Scale The coulomb 1-4 scale.
   * @param lj14Scale      The LJ 1-4 scale.
   */
  public void createExceptionsFromBonds(BondArray bondArray, double coulomb14Scale, double lj14Scale) {
    OpenMM_NonbondedForce_createExceptionsFromBonds(pointer, bondArray.getPointer(), coulomb14Scale, lj14Scale);
  }

  /**
   * Get the exception parameters.
   *
   * @param index      The exception index.
   * @param particle1  The first particle.
   * @param particle2  The second particle.
   * @param chargeProd The charge product.
   * @param sigma      The sigma vdW parameter.
   * @param eps        The eps vdW parameter.
   */
  public void getExceptionParameters(int index, IntByReference particle1, IntByReference particle2,
                                     DoubleByReference chargeProd, DoubleByReference sigma, DoubleByReference eps) {
    OpenMM_NonbondedForce_getExceptionParameters(pointer, index, particle1, particle2, chargeProd, sigma, eps);
  }

  /**
   * Set the exception parameters.
   *
   * @param index      The exception index.
   * @param particle1  The first particle.
   * @param particle2  The second particle.
   * @param chargeProd The charge product.
   * @param sigma      The sigma vdW parameter.
   * @param eps        The eps vdW parameter.
   */
  public void setExceptionParameters(int index, int particle1, int particle2, double chargeProd, double sigma, double eps) {
    OpenMM_NonbondedForce_setExceptionParameters(pointer, index, particle1, particle2, chargeProd, sigma, eps);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_NonbondedForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the PME parameters.
   *
   * @param aEwald The Ewald alpha.
   * @param nx     The PME grid size in x.
   * @param ny     The PME grid size in y.
   * @param nz     The PME grid size in z.
   */
  public void setPMEParameters(double aEwald, int nx, int ny, int nz) {
    OpenMM_NonbondedForce_setPMEParameters(pointer, aEwald, nx, ny, nz);
  }

  /**
   * Set the cutoff distance.
   *
   * @param cutoffDistance The cutoff distance.
   */
  public void setCutoffDistance(double cutoffDistance) {
    OpenMM_NonbondedForce_setCutoffDistance(pointer, cutoffDistance);
  }

  /**
   * Set if a switching function will be used.
   *
   * @param useSwitchingFunction The switching distance flag.
   */
  public void setUseSwitchingFunction(int useSwitchingFunction) {
    OpenMM_NonbondedForce_setUseSwitchingFunction(pointer, useSwitchingFunction);
  }

  /**
   * Set the switching distance.
   *
   * @param switchingDistance The switching distance.
   */
  public void setSwitchingDistance(double switchingDistance) {
    OpenMM_NonbondedForce_setSwitchingDistance(pointer, switchingDistance);
  }

  /**
   * Set if a dispersion correction will be used.
   *
   * @param useDispersionCorrection The dispersion correction flag.
   */
  public void setUseDispersionCorrection(int useDispersionCorrection) {
    OpenMM_NonbondedForce_setUseDispersionCorrection(pointer, useDispersionCorrection);
  }

  /**
   * Get the number of exceptions.
   *
   * @return The number of exceptions.
   */
  public int getNumExceptions() {
    return OpenMM_NonbondedForce_getNumExceptions(pointer);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_NonbondedForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  public void destroy() {
    if (pointer != null) {
      OpenMM_NonbondedForce_destroy(pointer);
      pointer = null;
    }
  }

}
