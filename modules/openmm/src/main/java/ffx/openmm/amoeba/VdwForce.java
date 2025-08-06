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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_getUseParticleTypes;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setEpsilonCombiningRule;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setLambdaName;
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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * This class models van der Waals forces in the AMOEBA force field.  It can use
 * either buffered 14-7 potential or a Lennard-Jones 12-6 potential.
 * <p>
 * This class can operate in two different modes.  In one mode, force field parameters
 * are defined for each particle.  When two particles interact, a combining rule is
 * used to calculate the interaction parameters based on the parameters for the two
 * particles.  To use the class in this mode, call the version of addParticle() that
 * takes sigma and epsilon values.  It should be called once for each particle in the
 * System.
 * <p>
 * In the other mode, each particle has a type index, and parameters are specified for
 * each type rather than each individual particle.  By default this mode also uses a
 * combining rule, but you can override it by defining alternate parameters to use for
 * specific pairs of particle types.  To use the class in this mode, call the version of
 * addParticle() that takes a type index.  It should be called once for each particle
 * in the System.  You also must call addParticleType() once for each type.  If you
 * wish to override the combining for particular pairs of types, do so by calling
 * addTypePair().
 * <p>
 * A unique feature of this class is that the interaction site for a particle does not need to be
 * exactly at the particle's location.  Instead, it can be placed a fraction of the distance from that
 * particle to another one.  This is typically done for hydrogens to place the interaction site slightly
 * closer to the parent atom.  The fraction is known as the "reduction factor", since it reduces the distance
 * from the parent atom to the interaction site.
 * <p>
 * Support is also available for softcore interactions based on setting a per particle alchemical flag and
 * setting the AmoebaVdwForce to use an "AlchemicalMethod" -- either Decouple or Annihilate.
 * For Decouple, two alchemical atoms interact normally. For Annihilate, all interactions involving an
 * alchemical atom are influenced. The softcore state is specified by setting a single
 * Context parameter "AmoebaVdwLambda" between 0.0 and 1.0.
 * <p>
 * The softcore functional form can be modified by setting the softcore power (default of 5) and the softcore
 * alpha (default of 0,7). For more information on the softcore functional form see Eq. 2 from:
 * Jiao, D.;  Golubkov, P. A.;  Darden, T. A.; Ren, P.,
 * Calculation of protein-ligand binding free energy by using a polarizable potential.
 * Proc. Natl. Acad. Sci. U.S.A. 2008, 105 (17), 6290-6295.
 * https://www.pnas.org/content/105/17/6290.
 */
public class VdwForce extends Force {

  /**
   * Create an Amoeba VdwForce.
   */
  public VdwForce() {
    super(OpenMM_AmoebaVdwForce_create());
  }

  /**
   * Add the force field parameters for a vdw particle.  This version is used when parameters
   * are defined for each particle.
   *
   * @param parentIndex     the index of the parent particle
   * @param sigma           vdw sigma
   * @param epsilon         vdw epsilon
   * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
   *                        at which the interaction site should be placed
   * @param isAlchemical    if true, this vdW particle is undergoing an alchemical change.
   * @param scaleFactor     a scale factor to apply to all interactions involving this particle (used for CpHMD).
   * @return index of added particle
   */
  public int addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor, int isAlchemical, double scaleFactor) {
    return OpenMM_AmoebaVdwForce_addParticle(pointer, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, scaleFactor);
  }

  /**
   * Add the force field parameters for a vdw particle. This version is used when parameters
   * are defined by particle type.
   *
   * @param parentIndex     the index of the parent particle
   * @param typeIndex       the index of the particle type for this particle
   * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
   *                        at which the interaction site should be placed
   * @param isAlchemical    if true, this vdW particle is undergoing an alchemical change.
   * @param scaleFactor     a scale factor to apply to all interactions involving this particle (used for CpHMD).
   * @return index of added particle
   */
  public int addParticle(int parentIndex, int typeIndex, double reductionFactor, int isAlchemical, double scaleFactor) {
    return OpenMM_AmoebaVdwForce_addParticle_1(pointer, parentIndex, typeIndex, reductionFactor, isAlchemical, scaleFactor);
  }

  /**
   * Add a particle type.
   *
   * @param sigma   the sigma value for particles of this type
   * @param epsilon the epsilon value for particles of this type
   * @return the index of the particle type that was just added.
   */
  public int addParticleType(double sigma, double epsilon) {
    return OpenMM_AmoebaVdwForce_addParticleType(pointer, sigma, epsilon);
  }

  /**
   * Add a type pair.  This overrides the standard combining rule for interactions
   * between particles of two particular types.
   *
   * @param type1   the index of the first particle type
   * @param type2   the index of the second particle type
   * @param sigma   the sigma value for interactions between particles of these two types
   * @param epsilon the epsilon  value for interactions between particles of these two types
   * @return the index of the type pair that was just added.
   */
  public int addTypePair(int type1, int type2, double sigma, double epsilon) {
    return OpenMM_AmoebaVdwForce_addTypePair(pointer, type1, type2, sigma, epsilon);
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
   * Get the cutoff distance.
   *
   * @deprecated This method exists only for backward compatibility.  Use getCutoffDistance() instead.
   */
  public double getCutoff() {
    return OpenMM_AmoebaVdwForce_getCutoff(pointer);
  }

  /**
   * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
   * is NoCutoff, this value will have no effect.
   *
   * @return the cutoff distance, measured in nm
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
    return OpenMM_AmoebaVdwForce_Lambda(pointer);
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
   * Get whether to add a contribution to the energy that approximately represents the effect of VdW
   * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
   * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
   * this contribution can improve the quality of results.
   */
  public boolean getUseDispersionCorrection() {
    return OpenMM_AmoebaVdwForce_getUseDispersionCorrection(pointer) != 0;
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
   * Set the cutoff distance.
   *
   * @deprecated This method exists only for backward compatibility.  Use setCutoffDistance() instead.
   */
  public void setCutoff(double cutoff) {
    OpenMM_AmoebaVdwForce_setCutoff(pointer, cutoff);
  }

  /**
   * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
   * is NoCutoff, this value will have no effect.
   *
   * @param distance the cutoff distance, measured in nm
   */
  public void setCutoffDistance(double distance) {
    OpenMM_AmoebaVdwForce_setCutoffDistance(pointer, distance);
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
   * Set the lambda parameter name.
   *
   * @param name The name of the lambda parameter.
   */
  public void setLambdaName(String name) {
    OpenMM_AmoebaVdwForce_setLambdaName(pointer, name);
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
   * Set whether to add a contribution to the energy that approximately represents the effect of VdW
   * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
   * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
   * this contribution can improve the quality of results.
   */
  public void setUseDispersionCorrection(boolean useCorrection) {
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection(pointer, useCorrection ? 1 : 0);
  }

  /**
   * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
   * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
   * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
   * to copy them over to the Context.
   * <p>
   * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
   * (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be changed by reinitializing the Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaVdwForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Returns whether or not this force makes use of periodic boundary
   * conditions.
   *
   * @returns true if nonbondedMethod uses PBC and false otherwise
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}