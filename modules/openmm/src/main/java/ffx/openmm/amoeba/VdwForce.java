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
import ffx.openmm.IntArray;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticleType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addTypePair;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcoreAlpha;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcorePower;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;

/**
 * Amoeba van der Waals Force.
 */
public class VdwForce extends Force {

  /**
   * The Amoeba vdW Force constructor.
   */
  public VdwForce() {
    pointer = OpenMM_AmoebaVdwForce_create();
  }

  /**
   * Add a particle type to the vdW Force.
   *
   * @param rad The radius.
   * @param eps The well depth.
   * @return The type.
   */
  public int addParticleType(double rad, double eps) {
    return OpenMM_AmoebaVdwForce_addParticleType(pointer, rad, eps);
  }

  /**
   * Add a type pair to the vdW Force.
   *
   * @param type1 The first type.
   * @param type2 The second type.
   * @param rad   The radius.
   * @param eps   The well depth.
   */
  public void addTypePair(int type1, int type2, double rad, double eps) {
    OpenMM_AmoebaVdwForce_addTypePair(pointer, type1, type2, rad, eps);
  }

  /**
   * Add a particle to the vdW Force.
   *
   * @param ired            The particle ired.
   * @param type            The particle type.
   * @param reductionFactor The reduction factor.
   * @param isAlchemical    The alchemical flag.
   * @param scaleFactor     The scale factor.
   */
  public void addParticle_1(int ired, int type, double reductionFactor, int isAlchemical, double scaleFactor) {
    OpenMM_AmoebaVdwForce_addParticle_1(pointer, ired, type, reductionFactor, isAlchemical, scaleFactor);
  }

  /**
   * Set the particle parameters.
   *
   * @param index           The particle index.
   * @param ired            The particle reduction index.
   * @param rad             The radius.
   * @param eps             The well depth.
   * @param reductionFactor The reduction factor.
   * @param isAlchemical    The alchemical flag.
   * @param type            The type.
   * @param scaleFactor     The scale factor.
   */
  public void setParticleParameters(int index, int ired, double rad, double eps, double reductionFactor,
                                    int isAlchemical, int type, double scaleFactor) {
    OpenMM_AmoebaVdwForce_setParticleParameters(pointer, index, ired, rad, eps, reductionFactor, isAlchemical, type, scaleFactor);
  }

  /**
   * Set the cutoff distance.
   *
   * @param cutoff The cutoff distance.
   */
  public void setCutoffDistance(double cutoff) {
    OpenMM_AmoebaVdwForce_setCutoffDistance(pointer, cutoff);
  }

  /**
   * Set the vdW force to use a long-range dispersion correction.
   *
   * @param value The flag.
   */
  public void setUseDispersionCorrection(int value) {
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection(pointer, value);
  }

  /**
   * Set the non-bonded method.
   *
   * @param method The non-bonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_AmoebaVdwForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the alchemical method.
   *
   * @param method The alchemical method.
   */
  public void setAlchemicalMethod(int method) {
    OpenMM_AmoebaVdwForce_setAlchemicalMethod(pointer, method);
  }

  /**
   * Set the softcore power.
   *
   * @param vdWSoftcoreAlpha The softcore power.
   */
  public void setSoftcoreAlpha(double vdWSoftcoreAlpha) {
    OpenMM_AmoebaVdwForce_setSoftcoreAlpha(pointer, vdWSoftcoreAlpha);
  }

  /**
   * Set the softcore power.
   *
   * @param vdwSoftcorePower The softcore power.
   */
  public void setSoftcorePower(int vdwSoftcorePower) {
    OpenMM_AmoebaVdwForce_setSoftcorePower(pointer, vdwSoftcorePower);
  }

  /**
   * Set the particle exclusions.
   *
   * @param i          The particle index.
   * @param exclusions The exclusions.
   */
  public void setParticleExclusions(int i, IntArray exclusions) {
    OpenMM_AmoebaVdwForce_setParticleExclusions(pointer, i, exclusions.getPointer());
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaVdwForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Destroy the vdW Force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaVdwForce_destroy(pointer);
      pointer = null;
    }
  }

}
