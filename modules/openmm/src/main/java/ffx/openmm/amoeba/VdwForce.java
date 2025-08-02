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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import ffx.openmm.Context;
import ffx.openmm.Force;
import ffx.openmm.IntArray;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_Lambda;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticleType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addTypePair;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getEpsilonCombiningRule;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getNumParticleTypes;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getNumTypePairs;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getParticleExclusions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getParticleTypeParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getPotentialFunction;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getSigmaCombiningRule;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getSoftcoreAlpha;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getSoftcorePower;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getTypePairParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getUseLambdaComplement;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getUseParticleTypes;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setEpsilonCombiningRule;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleTypeParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setPotentialFunction;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSigmaCombiningRule;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcoreAlpha;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcorePower;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setTypePairParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setUseLambdaComplement;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * Amoeba van der Waals Force.
 */
public class VdwForce extends Force {

  /**
   * Create an OpenMM VdwForce.
   */
  public VdwForce() {
    super(OpenMM_AmoebaVdwForce_create());
  }

  /**
   * Add a particle to the force.
   *
   * @param ired            The index of the particle that this particle is reduced to.
   * @param rad             The radius of the particle.
   * @param eps             The epsilon of the particle.
   * @param reductionFactor The reduction factor.
   * @param isAlchemical    Whether the particle is alchemical.
   * @param scaleFactor     The scale factor.
   * @return The index of the added particle.
   */
  public int addParticle(int ired, double rad, double eps, double reductionFactor, int isAlchemical, double scaleFactor) {
    return OpenMM_AmoebaVdwForce_addParticle(pointer, ired, rad, eps, reductionFactor, isAlchemical, scaleFactor);
  }

  /**
   * Add a particle to the force.
   *
   * @param ired            The index of the particle that this particle is reduced to.
   * @param type            The type of the particle.
   * @param reductionFactor The reduction factor.
   * @param isAlchemical    Whether the particle is alchemical.
   * @param scaleFactor     The scale factor.
   */
  public void addParticle_1(int ired, int type, double reductionFactor, int isAlchemical, double scaleFactor) {
    OpenMM_AmoebaVdwForce_addParticle_1(pointer, ired, type, reductionFactor, isAlchemical, scaleFactor);
  }

  /**
   * Add a particle type to the force.
   *
   * @param rad The radius of the particle type.
   * @param eps The epsilon of the particle type.
   * @return The index of the added particle type.
   */
  public int addParticleType(double rad, double eps) {
    return OpenMM_AmoebaVdwForce_addParticleType(pointer, rad, eps);
  }

  /**
   * Add a type pair to the force.
   *
   * @param type1 The first type.
   * @param type2 The second type.
   * @param rad   The radius.
   * @param eps   The epsilon.
   * @return The index of the added type pair.
   */
  public int addTypePair(int type1, int type2, double rad, double eps) {
    return OpenMM_AmoebaVdwForce_addTypePair(pointer, type1, type2, rad, eps);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaVdwForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the alchemical method.
   *
   * @return The alchemical method.
   */
  public int getAlchemicalMethod() {
    return OpenMM_AmoebaVdwForce_getAlchemicalMethod(pointer);
  }

  /**
   * Get the cutoff.
   *
   * @return The cutoff.
   */
  public double getCutoff() {
    return OpenMM_AmoebaVdwForce_getCutoff(pointer);
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance.
   */
  public double getCutoffDistance() {
    return OpenMM_AmoebaVdwForce_getCutoffDistance(pointer);
  }

  /**
   * Get the epsilon combining rule.
   *
   * @return The epsilon combining rule.
   */
  public String getEpsilonCombiningRule() {
    Pointer rule = OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(pointer);
    if (rule == null) {
      return null;
    }
    return rule.getString(0);
  }

  /**
   * Get the lambda parameter.
   *
   * @return The lambda parameter.
   */
  public Pointer getLambda() {
    return OpenMM_AmoebaVdwForce_Lambda();
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_AmoebaVdwForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_AmoebaVdwForce_getNumParticles(pointer);
  }

  /**
   * Get the number of particle types.
   *
   * @return The number of particle types.
   */
  public int getNumParticleTypes() {
    return OpenMM_AmoebaVdwForce_getNumParticleTypes(pointer);
  }

  /**
   * Get the number of type pairs.
   *
   * @return The number of type pairs.
   */
  public int getNumTypePairs() {
    return OpenMM_AmoebaVdwForce_getNumTypePairs(pointer);
  }

  /**
   * Get the particle exclusions.
   *
   * @param i The index of the particle.
   * @return An IntArray containing the exclusions.
   */
  public IntArray getParticleExclusions(int i) {
    IntArray exclusions = new IntArray(0);
    if (pointer != null) {
      OpenMM_AmoebaVdwForce_getParticleExclusions(pointer, i, exclusions.getPointer());
    }
    return exclusions;
  }

  /**
   * Get the particle parameters.
   *
   * @param index           The index of the particle.
   * @param ired            The index of the particle that this particle is reduced to (output).
   * @param rad             The radius of the particle (output).
   * @param eps             The epsilon of the particle (output).
   * @param reductionFactor The reduction factor (output).
   * @param isAlchemical    Whether the particle is alchemical (output).
   * @param type            The type of the particle (output).
   * @param scaleFactor     The scale factor (output).
   */
  public void getParticleParameters(int index, IntByReference ired, DoubleByReference rad,
                                    DoubleByReference eps, DoubleByReference reductionFactor,
                                    IntByReference isAlchemical, IntByReference type,
                                    DoubleByReference scaleFactor) {
    OpenMM_AmoebaVdwForce_getParticleParameters(pointer, index, ired, rad, eps, reductionFactor,
        isAlchemical, type, scaleFactor);
  }

  /**
   * Get the particle type parameters.
   *
   * @param index The index of the particle type.
   * @param rad   The radius of the particle type (output).
   * @param eps   The epsilon of the particle type (output).
   */
  public void getParticleTypeParameters(int index, DoubleByReference rad, DoubleByReference eps) {
    OpenMM_AmoebaVdwForce_getParticleTypeParameters(pointer, index, rad, eps);
  }

  /**
   * Get the potential function.
   *
   * @return The potential function.
   */
  public int getPotentialFunction() {
    return OpenMM_AmoebaVdwForce_getPotentialFunction(pointer);
  }

  /**
   * Get the sigma combining rule.
   *
   * @return The sigma combining rule.
   */
  public String getSigmaCombiningRule() {
    Pointer rule = OpenMM_AmoebaVdwForce_getSigmaCombiningRule(pointer);
    if (rule == null) {
      return null;
    }
    return rule.getString(0);
  }

  /**
   * Get the softcore alpha.
   *
   * @return The softcore alpha.
   */
  public double getSoftcoreAlpha() {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha(pointer);
  }

  /**
   * Get the softcore power.
   *
   * @return The softcore power.
   */
  public int getSoftcorePower() {
    return OpenMM_AmoebaVdwForce_getSoftcorePower(pointer);
  }

  /**
   * Get the type pair parameters.
   *
   * @param index The index of the type pair.
   * @param type1 The first type (output).
   * @param type2 The second type (output).
   * @param rad   The radius (output).
   * @param eps   The epsilon (output).
   */
  public void getTypePairParameters(int index, IntByReference type1, IntByReference type2,
                                    DoubleByReference rad, DoubleByReference eps) {
    OpenMM_AmoebaVdwForce_getTypePairParameters(pointer, index, type1, type2, rad, eps);
  }

  /**
   * Get whether to use dispersion correction.
   *
   * @return 1 if dispersion correction is used, 0 otherwise.
   */
  public int getUseDispersionCorrection() {
    return OpenMM_AmoebaVdwForce_getUseDispersionCorrection(pointer);
  }

  /**
   * Get whether to use lambda complement.
   *
   * @return 1 if lambda complement is used, 0 otherwise.
   */
  public int getUseLambdaComplement() {
    return OpenMM_AmoebaVdwForce_getUseLambdaComplement(pointer);
  }

  /**
   * Get whether to use particle types.
   *
   * @return 1 if particle types are used, 0 otherwise.
   */
  public int getUseParticleTypes() {
    return OpenMM_AmoebaVdwForce_getUseParticleTypes(pointer);
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
   * Set the cutoff.
   *
   * @param cutoff The cutoff.
   */
  public void setCutoff(double cutoff) {
    OpenMM_AmoebaVdwForce_setCutoff(pointer, cutoff);
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
   * Set the epsilon combining rule.
   *
   * @param rule The epsilon combining rule.
   */
  public void setEpsilonCombiningRule(String rule) {
    OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(pointer, rule);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_AmoebaVdwForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the particle exclusions.
   *
   * @param i          The index of the particle.
   * @param exclusions The exclusions.
   */
  public void setParticleExclusions(int i, IntArray exclusions) {
    OpenMM_AmoebaVdwForce_setParticleExclusions(pointer, i, exclusions.getPointer());
  }

  /**
   * Set the particle parameters.
   *
   * @param index           The index of the particle.
   * @param ired            The index of the particle that this particle is reduced to.
   * @param rad             The radius of the particle.
   * @param eps             The epsilon of the particle.
   * @param reductionFactor The reduction factor.
   * @param isAlchemical    Whether the particle is alchemical.
   * @param type            The type of the particle.
   * @param scaleFactor     The scale factor.
   */
  public void setParticleParameters(int index, int ired, double rad, double eps, double reductionFactor, int isAlchemical, int type, double scaleFactor) {
    OpenMM_AmoebaVdwForce_setParticleParameters(pointer, index, ired, rad, eps, reductionFactor, isAlchemical, type, scaleFactor);
  }

  /**
   * Set the particle type parameters.
   *
   * @param index The index of the particle type.
   * @param rad   The radius of the particle type.
   * @param eps   The epsilon of the particle type.
   */
  public void setParticleTypeParameters(int index, double rad, double eps) {
    OpenMM_AmoebaVdwForce_setParticleTypeParameters(pointer, index, rad, eps);
  }

  /**
   * Set the potential function.
   *
   * @param function The potential function.
   */
  public void setPotentialFunction(int function) {
    OpenMM_AmoebaVdwForce_setPotentialFunction(pointer, function);
  }

  /**
   * Set the sigma combining rule.
   *
   * @param rule The sigma combining rule.
   */
  public void setSigmaCombiningRule(String rule) {
    OpenMM_AmoebaVdwForce_setSigmaCombiningRule(pointer, rule);
  }

  /**
   * Set the softcore alpha.
   *
   * @param vdWSoftcoreAlpha The softcore alpha.
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
   * Set the type pair parameters.
   *
   * @param index The index of the type pair.
   * @param type1 The first type.
   * @param type2 The second type.
   * @param rad   The radius.
   * @param eps   The epsilon.
   */
  public void setTypePairParameters(int index, int type1, int type2, double rad, double eps) {
    OpenMM_AmoebaVdwForce_setTypePairParameters(pointer, index, type1, type2, rad, eps);
  }

  /**
   * Set whether to use dispersion correction.
   *
   * @param value 1 to use dispersion correction, 0 otherwise.
   */
  public void setUseDispersionCorrection(int value) {
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection(pointer, value);
  }

  /**
   * Set whether to use lambda complement.
   *
   * @param useLambdaComplement 1 to use lambda complement, 0 otherwise.
   */
  public void setUseLambdaComplement(int useLambdaComplement) {
    OpenMM_AmoebaVdwForce_setUseLambdaComplement(pointer, useLambdaComplement);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaVdwForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}