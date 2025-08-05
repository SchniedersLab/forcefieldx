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
import ffx.openmm.Context;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * Amoeba Generalized Kirkwood Force.
 */
public class GeneralizedKirkwoodForce extends Force {

  public GeneralizedKirkwoodForce() {
    super(OpenMM_AmoebaGeneralizedKirkwoodForce_create());
  }

  /**
   * Add a particle to the force with 3 parameters.
   *
   * @param charge   The charge of the particle.
   * @param radius   The radius of the particle.
   * @param hctScale The hctScale of the particle.
   * @return The index of the added particle.
   */
  public int addParticle(double charge, double radius, double hctScale) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(pointer, charge, radius, hctScale);
  }

  /**
   * Add a particle to the force.
   *
   * @param charge   The charge of the particle.
   * @param radius   The radius of the particle.
   * @param hctScale The hctScale of the particle.
   * @param descreen The descreen of the particle.
   * @param neck     The neck of the particle.
   */
  public void addParticle_1(double charge, double radius, double hctScale, double descreen, double neck) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(pointer, charge, radius, hctScale, descreen, neck);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the descreen offset.
   *
   * @return The descreen offset.
   */
  public double getDescreenOffset() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(pointer);
  }

  /**
   * Get the dielectric offset.
   *
   * @return The dielectric offset.
   */
  public double getDielectricOffset() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(pointer);
  }

  /**
   * Get the include cavity term.
   *
   * @return The include cavity term.
   */
  public int getIncludeCavityTerm() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(pointer);
  }

  /**
   * Get the number of particles in the force.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(pointer);
  }

  /**
   * Get the particle parameters.
   *
   * @param index    The index of the particle.
   * @param charge   The charge of the particle (output).
   * @param radius   The radius of the particle (output).
   * @param hctScale The hctScale of the particle (output).
   * @param descreen The descreen of the particle (output).
   * @param neck     The neck of the particle (output).
   */
  public void getParticleParameters(int index, DoubleByReference charge, DoubleByReference radius,
                                    DoubleByReference hctScale, DoubleByReference descreen,
                                    DoubleByReference neck) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(pointer, index, charge, radius,
        hctScale, descreen, neck);
  }

  /**
   * Get the probe radius.
   *
   * @return The probe radius.
   */
  public double getProbeRadius() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(pointer);
  }

  /**
   * Get the solute dielectric constant.
   *
   * @return The solute dielectric constant.
   */
  public double getSoluteDielectric() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(pointer);
  }

  /**
   * Get the solvent dielectric constant.
   *
   * @return The solvent dielectric constant.
   */
  public double getSolventDielectric() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(pointer);
  }

  /**
   * Get the surface area factor.
   *
   * @return The surface area factor.
   */
  public double getSurfaceAreaFactor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(pointer);
  }

  /**
   * Get the tanh parameters.
   *
   * @param beta0 The tanh parameter beta0 (output).
   * @param beta1 The tanh parameter beta1 (output).
   * @param beta2 The tanh parameter beta2 (output).
   */
  public void getTanhParameters(DoubleByReference beta0, DoubleByReference beta1, DoubleByReference beta2) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(pointer, beta0, beta1, beta2);
  }

  /**
   * Get the tanh rescaling.
   *
   * @return The tanh rescaling.
   */
  public int getTanhRescaling() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(pointer);
  }

  /**
   * Set the descreen offset.
   *
   * @param offset The descreen offset.
   */
  public void setDescreenOffset(double offset) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(pointer, offset);
  }

  /**
   * Set the dielectric offset.
   *
   * @param offset The dielectric offset.
   */
  public void setDielectricOffset(double offset) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(pointer, offset);
  }

  /**
   * Set the include cavity term.
   *
   * @param includeCavityTerm The include cavity term.
   */
  public void setIncludeCavityTerm(int includeCavityTerm) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(pointer, includeCavityTerm);
  }

  /**
   * Set the particle parameters.
   *
   * @param index    The index of the particle.
   * @param charge   The charge of the particle.
   * @param radius   The radius of the particle.
   * @param hctScale The hctScale of the particle.
   * @param descreen The descreen of the particle.
   * @param neck     The neck of the particle.
   */
  public void setParticleParameters_1(int index, double charge, double radius, double hctScale, double descreen, double neck) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(pointer, index, charge, radius, hctScale, descreen, neck);
  }

  /**
   * Set the probe radius.
   *
   * @param radius The probe radius.
   */
  public void setProbeRadius(double radius) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(pointer, radius);
  }

  /**
   * Set the solute dielectric constant.
   *
   * @param dielectric The solute dielectric constant.
   */
  public void setSoluteDielectric(double dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(pointer, dielectric);
  }

  /**
   * Set the solvent dielectric constant.
   *
   * @param dielectric The solvent dielectric constant.
   */
  public void setSolventDielectric(double dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(pointer, dielectric);
  }

  /**
   * Set the surface area factor.
   *
   * @param surfaceAreaFactor The surface area factor.
   */
  public void setSurfaceAreaFactor(double surfaceAreaFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(pointer, surfaceAreaFactor);
  }

  /**
   * Set the tanh parameters.
   *
   * @param beta0 The tanh parameter beta0.
   * @param beta1 The tanh parameter beta1.
   * @param beta2 The tanh parameter beta2.
   */
  public void setTanhParameters(double beta0, double beta1, double beta2) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(pointer, beta0, beta1, beta2);
  }

  /**
   * Set the tanh rescaling.
   *
   * @param tanhRescale The tanh rescaling.
   */
  public void setTanhRescaling(int tanhRescale) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(pointer, tanhRescale);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }

}
