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
package ffx.openmm.amoeba;

import ffx.openmm.Context;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_destroy;
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

/**
 * Amoeba Generalized Kirkwood Force.
 */
public class GeneralizedKirkwoodForce extends Force {

  public GeneralizedKirkwoodForce() {
    pointer = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
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
   * Set the solute dielectric constant.
   *
   * @param dielectric The solute dielectric constant.
   */
  public void setSoluteDielectric(double dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(pointer, dielectric);
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
   * Set the tanh rescaling.
   *
   * @param tanhRescale The tanh rescaling.
   */
  public void setTanhRescaling(int tanhRescale) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(pointer, tanhRescale);
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
   * Set the probe radius.
   *
   * @param radius The probe radius.
   */
  public void setProbeRadius(double radius) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(pointer, radius);
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
   * Set the surface area factor.
   *
   * @param surfaceAreaFactor The surface area factor.
   */
  public void setSurfaceAreaFactor(double surfaceAreaFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(pointer, surfaceAreaFactor);
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
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(pointer);
      pointer = null;
    }
  }

}
