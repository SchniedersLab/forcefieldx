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
package ffx.openmm;

import com.sun.jna.ptr.DoubleByReference;

import java.nio.DoubleBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_getSurfaceAreaEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_setSurfaceAreaEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions;

/**
 * This class implements an implicit solvation force using the GBSA-OBC model.
 * <p>
 * To use this class, create a GBSAOBCForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define GBSA parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().  This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().
 * <p>
 * When using this Force, the System should also include a NonbondedForce, and both objects must specify
 * identical charges for all particles.  Otherwise, the results will not be correct.  Furthermore, if the
 * nonbonded method is set to CutoffNonPeriodic or CutoffPeriodic, you should call setReactionFieldDielectric(1.0)
 * on the NonbondedForce to turn off the reaction field approximation, which does not produce correct results
 * when combined with GBSA.
 */
public class GBSAOBCForce extends Force {

  /**
   * Create a new GBSAOBCForce.
   */
  public GBSAOBCForce() {
    super(OpenMM_GBSAOBCForce_create());
  }

  /**
   * Add a particle to the force field.
   *
   * @param charge The charge of the particle, measured in units of the proton charge.
   * @param radius The GBSA radius of the particle, measured in nm.
   * @param scale  The OBC scaling parameter for the particle.
   * @return The index of the particle that was added.
   */
  public int addParticle(double charge, double radius, double scale) {
    return OpenMM_GBSAOBCForce_addParticle(pointer, charge, radius, scale);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_GBSAOBCForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the cutoff distance (in nm) being used for nonbonded interactions.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_GBSAOBCForce_getCutoffDistance(pointer);
  }

  /**
   * Get the method used for handling long range nonbonded interactions.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_GBSAOBCForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of particles in the force field.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_GBSAOBCForce_getNumParticles(pointer);
  }

  /**
   * Get the force field parameters for a particle.
   *
   * @param index  The index of the particle for which to get parameters.
   * @param charge The charge of the particle, measured in units of the proton charge (output).
   * @param radius The GBSA radius of the particle, measured in nm (output).
   * @param scale  The OBC scaling parameter for the particle (output).
   */
  public void getParticleParameters(int index, DoubleByReference charge, DoubleByReference radius, DoubleByReference scale) {
    OpenMM_GBSAOBCForce_getParticleParameters(pointer, index, charge, radius, scale);
  }

  /**
   * Get the force field parameters for a particle.
   *
   * @param index  The index of the particle for which to get parameters.
   * @param charge The charge of the particle, measured in units of the proton charge (output).
   * @param radius The GBSA radius of the particle, measured in nm (output).
   * @param scale  The OBC scaling parameter for the particle (output).
   */
  public void getParticleParameters(int index, DoubleBuffer charge, DoubleBuffer radius, DoubleBuffer scale) {
    OpenMM_GBSAOBCForce_getParticleParameters(pointer, index, charge, radius, scale);
  }

  /**
   * Get the dielectric constant for the solute.
   *
   * @return The solute dielectric constant.
   */
  public double getSoluteDielectric() {
    return OpenMM_GBSAOBCForce_getSoluteDielectric(pointer);
  }

  /**
   * Get the dielectric constant for the solvent.
   *
   * @return The solvent dielectric constant.
   */
  public double getSolventDielectric() {
    return OpenMM_GBSAOBCForce_getSolventDielectric(pointer);
  }

  /**
   * Get the energy scale for the surface area term, measured in kJ/mol/nm&circ;2.
   *
   * @return The surface area energy scale.
   */
  public double getSurfaceAreaEnergy() {
    return OpenMM_GBSAOBCForce_getSurfaceAreaEnergy(pointer);
  }

  /**
   * Set the cutoff distance (in nm) being used for nonbonded interactions.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_GBSAOBCForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the method used for handling long range nonbonded interactions.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_GBSAOBCForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the force field parameters for a particle.
   *
   * @param index  The index of the particle for which to set parameters.
   * @param charge The charge of the particle, measured in units of the proton charge.
   * @param radius The GBSA radius of the particle, measured in nm.
   * @param scale  The OBC scaling parameter for the particle.
   */
  public void setParticleParameters(int index, double charge, double radius, double scale) {
    OpenMM_GBSAOBCForce_setParticleParameters(pointer, index, charge, radius, scale);
  }

  /**
   * Set the dielectric constant for the solute.
   *
   * @param dielectric The solute dielectric constant.
   */
  public void setSoluteDielectric(double dielectric) {
    OpenMM_GBSAOBCForce_setSoluteDielectric(pointer, dielectric);
  }

  /**
   * Set the dielectric constant for the solvent.
   *
   * @param dielectric The solvent dielectric constant.
   */
  public void setSolventDielectric(double dielectric) {
    OpenMM_GBSAOBCForce_setSolventDielectric(pointer, dielectric);
  }

  /**
   * Set the energy scale for the surface area term, measured in kJ/mol/nm&circ;2.
   *
   * @param energy The surface area energy scale.
   */
  public void setSurfaceAreaEnergy(double energy) {
    OpenMM_GBSAOBCForce_setSurfaceAreaEnergy(pointer, energy);
  }

  /**
   * Update the particle parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_GBSAOBCForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}