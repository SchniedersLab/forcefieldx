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
package ffx.openmm.drude;

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import ffx.openmm.Context;
import ffx.openmm.Force;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_addScreenedPair;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_getNumScreenedPairs;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_getScreenedPairParameters;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_setScreenedPairParameters;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * This class implements forces that are specific to Drude oscillators.  There are two distinct forces
 * it applies: an anisotropic harmonic force connecting each Drude particle to its parent particle; and
 * a screened Coulomb interaction between specific pairs of dipoles.  The latter is typically used between
 * closely bonded particles whose Coulomb interaction would otherwise be fully excluded.
 * <p>
 * To use this class, create a DrudeForce object, then call addParticle() once for each Drude particle in the
 * System to define its parameters.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().  This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().  Likewise, call addScreenedPair() for each pair of dipoles (each dipole
 * consisting of a Drude particle and its parent) that should be computed.
 */
public class DrudeForce extends Force {

  /**
   * Create a DrudeForce.
   */
  public DrudeForce() {
    super(OpenMM_DrudeForce_create());
  }

  /**
   * Add a Drude particle to which forces should be applied.
   *
   * @param particle       the index within the System of the Drude particle
   * @param particle1      the index within the System of the particle to which the Drude particle is attached
   * @param particle2      the index within the System of the second particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso12 will be ignored.
   * @param particle3      the index within the System of the third particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param particle4      the index within the System of the fourth particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param charge         The charge on the Drude particle
   * @param polarizability The isotropic polarizability
   * @param aniso12        The scale factor for the polarizability along the direction defined by particle1 and particle2
   * @param aniso34        The scale factor for the polarizability along the direction defined by particle3 and particle4
   * @return the index of the particle that was added
   */
  public int addParticle(int particle, int particle1, int particle2, int particle3, int particle4,
                         double charge, double polarizability, double aniso12, double aniso34) {
    return OpenMM_DrudeForce_addParticle(pointer, particle, particle1, particle2, particle3, particle4,
        charge, polarizability, aniso12, aniso34);
  }

  /**
   * Add an interaction to the list of screened pairs.
   *
   * @param particle1 the index within this Force of the first particle involved in the interaction
   * @param particle2 the index within this Force of the second particle involved in the interaction
   * @param thole     the Thole screening factor
   * @return the index of the screenedPair that was added
   */
  public int addScreenedPair(int particle1, int particle2, double thole) {
    return OpenMM_DrudeForce_addScreenedPair(pointer, particle1, particle2, thole);
  }

  /**
   * Destroy the force.
   * <p>
   * This method releases the memory associated with the DrudeForce object.
   * After calling this method, the force should not be used.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_DrudeForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the number of particles for which force field parameters have been defined.
   */
  public int getNumParticles() {
    return OpenMM_DrudeForce_getNumParticles(pointer);
  }

  /**
   * Get the number of special interactions that should be calculated differently from other interactions.
   */
  public int getNumScreenedPairs() {
    return OpenMM_DrudeForce_getNumScreenedPairs(pointer);
  }

  /**
   * Get the parameters for a Drude particle.
   *
   * @param index          the index of the Drude particle for which to get parameters
   * @param particle       the index within the System of the Drude particle
   * @param particle1      the index within the System of the particle to which the Drude particle is attached
   * @param particle2      the index within the System of the second particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso12 will be ignored.
   * @param particle3      the index within the System of the third particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param particle4      the index within the System of the fourth particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param charge         The charge on the Drude particle
   * @param polarizability The isotropic polarizability
   * @param aniso12        The scale factor for the polarizability along the direction defined by particle1 and particle2
   * @param aniso34        The scale factor for the polarizability along the direction defined by particle3 and particle4
   */
  public void getParticleParameters(int index, IntByReference particle, IntByReference particle1,
                                    IntByReference particle2, IntByReference particle3, IntByReference particle4,
                                    DoubleByReference charge, DoubleByReference polarizability,
                                    DoubleByReference aniso12, DoubleByReference aniso34) {
    OpenMM_DrudeForce_getParticleParameters(pointer, index, particle, particle1, particle2,
        particle3, particle4, charge, polarizability, aniso12, aniso34);
  }

  /**
   * Get the force field parameters for a screened pair.
   *
   * @param index     the index of the pair for which to get parameters
   * @param particle1 the index within this Force of the first particle involved in the interaction
   * @param particle2 the index within this Force of the second particle involved in the interaction
   * @param thole     the Thole screening factor
   */
  public void getScreenedPairParameters(int index, IntByReference particle1, IntByReference particle2,
                                        DoubleByReference thole) {
    OpenMM_DrudeForce_getScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Get the force field parameters for a screened pair.
   *
   * @param index     the index of the pair for which to get parameters
   * @param particle1 the index within this Force of the first particle involved in the interaction
   * @param particle2 the index within this Force of the second particle involved in the interaction
   * @param thole     the Thole screening factor
   */
  public void getScreenedPairParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                        DoubleBuffer thole) {
    OpenMM_DrudeForce_getScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Set the parameters for a Drude particle.
   *
   * @param index          the index of the Drude particle for which to set parameters
   * @param particle       the index within the System of the Drude particle
   * @param particle1      the index within the System of the particle to which the Drude particle is attached
   * @param particle2      the index within the System of the second particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso12 will be ignored.
   * @param particle3      the index within the System of the third particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param particle4      the index within the System of the fourth particle used for defining anisotropic polarizability.
   *                       This may be set to -1, in which case aniso34 will be ignored.
   * @param charge         The charge on the Drude particle
   * @param polarizability The isotropic polarizability
   * @param aniso12        The scale factor for the polarizability along the direction defined by particle1 and particle2
   * @param aniso34        The scale factor for the polarizability along the direction defined by particle3 and particle4
   */
  public void setParticleParameters(int index, int particle, int particle1, int particle2,
                                    int particle3, int particle4, double charge, double polarizability,
                                    double aniso12, double aniso34) {
    OpenMM_DrudeForce_setParticleParameters(pointer, index, particle, particle1, particle2,
        particle3, particle4, charge, polarizability, aniso12, aniso34);
  }

  /**
   * Set the force field parameters for a screened pair.
   *
   * @param index     the index of the pair for which to get parameters
   * @param particle1 the index within this Force of the first particle involved in the interaction
   * @param particle2 the index within this Force of the second particle involved in the interaction
   * @param thole     the Thole screening factor
   */
  public void setScreenedPairParameters(int index, int particle1, int particle2, double thole) {
    OpenMM_DrudeForce_setScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
   * <p>
   * Periodic boundary conditions are only applied to screened pairs.  They are never used for the
   * force between a Drude particle and its parent particle, regardless of this setting.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the particle and screened pair parameters in a Context to match those stored in this Force object.  This method
   * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
   * Simply call setParticleParameters() and setScreenedPairParameters() to modify this object's parameters, then call
   * updateParametersInContext() to copy them over to the Context.
   * <p>
   * This method has several limitations.  It can be used to modify the numeric parameters associated with a particle or
   * screened pair (polarizability, thole, etc.), but not the identities of the particles they involve.  It also cannot
   * be used to add new particles or screenedPairs, only to change the parameters of existing ones.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_DrudeForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Returns whether or not this force makes use of periodic boundary
   * conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_DrudeForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}