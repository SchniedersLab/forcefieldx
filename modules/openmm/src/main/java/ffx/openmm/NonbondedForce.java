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
package ffx.openmm;

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_addException;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_createExceptionsFromBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getNumExceptions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getNumParticles;
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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_usesPeriodicBoundaryConditions;

/**
 * This class implements nonbonded interactions between particles, including a Coulomb force to represent
 * electrostatics and a Lennard-Jones force to represent van der Waals interactions. It optionally supports
 * periodic boundary conditions and cutoffs for long range interactions. Lennard-Jones interactions are
 * calculated with the Lorentz-Berthelot combining rule: it uses the arithmetic mean of the sigmas and the
 * geometric mean of the epsilons for the two interacting particles.
 * <p>
 * To use this class, create a NonbondedForce object, then call addParticle() once for each particle in the
 * System to define its parameters. The number of particles for which you define nonbonded parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context. After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters(). This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().
 * <p>
 * NonbondedForce also lets you specify "exceptions", particular pairs of particles whose interactions should be
 * computed based on different parameters than those defined for the individual particles. This can be used to
 * completely exclude certain interactions from the force calculation, or to alter how they interact with each other.
 * <p>
 * Many molecular force fields omit Coulomb and Lennard-Jones interactions between particles separated by one
 * or two bonds, while using modified parameters for those separated by three bonds (known as "1-4 interactions").
 * This class provides a convenience method for this case called createExceptionsFromBonds(). You pass to it
 * a list of bonds and the scale factors to use for 1-4 interactions. It identifies all pairs of particles which
 * are separated by 1, 2, or 3 bonds, then automatically creates exceptions for them.
 * <p>
 * When using a cutoff, by default Lennard-Jones interactions are sharply truncated at the cutoff distance.
 * Optionally you can instead use a switching function to make the interaction smoothly go to zero over a finite
 * distance range. To enable this, call setUseSwitchingFunction(). You must also call setSwitchingDistance()
 * to specify the distance at which the interaction should begin to decrease. The switching distance must be
 * less than the cutoff distance.
 * <p>
 * Another optional feature of this class (enabled by default) is to add a contribution to the energy which approximates
 * the effect of all Lennard-Jones interactions beyond the cutoff in a periodic system. When running a simulation
 * at constant pressure, this can improve the quality of the result. Call setUseDispersionCorrection() to set whether
 * this should be used.
 * <p>
 * In some applications, it is useful to be able to inexpensively change the parameters of small groups of particles.
 * Usually this is done to interpolate between two sets of parameters. For example, a titratable group might have
 * two states it can exist in, each described by a different set of parameters for the atoms that make up the
 * group. You might then want to smoothly interpolate between the two states. This is done by first calling
 * addGlobalParameter() to define a Context parameter, then addParticleParameterOffset() to create a "parameter offset"
 * that depends on the Context parameter. Each offset defines the following:
 * <ul>
 * <li>A Context parameter used to interpolate between the states.</li>
 * <li>A single particle whose parameters are influenced by the Context parameter.</li>
 * <li>Three scale factors (chargeScale, sigmaScale, and epsilonScale) that specify how the Context parameter
 * affects the particle.</li>
 * </ul>
 * <p>
 * The "effective" parameters for a particle (those used to compute forces) are given by
 * <pre>{@code
 * charge = baseCharge + param*chargeScale
 * sigma = baseSigma + param*sigmaScale
 * epsilon = baseEpsilon + param*epsilonScale
 * }</pre>
 * where the "base" values are the ones specified by addParticle() and "param" is the current value
 * of the Context parameter. A single Context parameter can apply offsets to multiple particles,
 * and multiple parameters can be used to apply offsets to the same particle. Parameters can also be used
 * to modify exceptions in exactly the same way by calling addExceptionParameterOffset().
 */
public class NonbondedForce extends Force {

  /**
   * Create a new NonbondedForce.
   */
  public NonbondedForce() {
    super(OpenMM_NonbondedForce_create());
  }

  /**
   * Add an exception to the force.
   *
   * @param particle1  The index of the first particle.
   * @param particle2  The index of the second particle.
   * @param chargeProd The charge product.
   * @param sigma      The sigma vdW parameter.
   * @param epsilon    The epsilon vdW parameter.
   * @param replace    If true, replace any existing exception.
   * @return The index of the exception that was added.
   */
  public int addException(int particle1, int particle2, double chargeProd, double sigma, double epsilon, boolean replace) {
    return OpenMM_NonbondedForce_addException(pointer, particle1, particle2, chargeProd, sigma, epsilon, replace ? 1 : 0);
  }

  /**
   * Add a particle.
   *
   * @param charge The atomic charge.
   * @param sigma  The vdW sigma.
   * @param eps    The vdW eps.
   * @return The index of the particle that was added.
   */
  public int addParticle(double charge, double sigma, double eps) {
    return OpenMM_NonbondedForce_addParticle(pointer, charge, sigma, eps);
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
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_NonbondedForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance.
   */
  public double getCutoffDistance() {
    return OpenMM_NonbondedForce_getCutoffDistance(pointer);
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
   * Get the number of exceptions.
   *
   * @return The number of exceptions.
   */
  public int getNumExceptions() {
    return OpenMM_NonbondedForce_getNumExceptions(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_NonbondedForce_getNumParticles(pointer);
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
   * Set the cutoff distance.
   *
   * @param cutoffDistance The cutoff distance.
   */
  public void setCutoffDistance(double cutoffDistance) {
    OpenMM_NonbondedForce_setCutoffDistance(pointer, cutoffDistance);
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
   * Set if a switching function will be used.
   *
   * @param useSwitchingFunction The switching distance flag.
   */
  public void setUseSwitchingFunction(int useSwitchingFunction) {
    OpenMM_NonbondedForce_setUseSwitchingFunction(pointer, useSwitchingFunction);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_NonbondedForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_NonbondedForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}