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
 * DrudeForce implements the Drude oscillator model for electronic polarization.
 * <p>
 * The Drude model represents electronic polarization by attaching a fictitious charged
 * particle (Drude particle) to each polarizable atom via a harmonic spring. The Drude
 * particle carries a small negative charge, while the core atom carries a compensating
 * positive charge. When an external electric field is applied, the Drude particle is
 * displaced from the core, creating an induced dipole moment.
 * <p>
 * This force implements the harmonic restraint between Drude particles and their
 * parent atoms, as well as the electrostatic interactions involving the Drude charges.
 * The total energy includes:
 * <ul>
 * <li>Harmonic restraint energy: E = 0.5 * k * r²</li>
 * <li>Electrostatic interactions between Drude particles and other charges</li>
 * <li>Screening interactions to prevent overpolarization at short distances</li>
 * </ul>
 * <p>
 * The Drude model provides a classical treatment of electronic polarization that
 * can significantly improve the accuracy of molecular simulations, particularly
 * for systems involving ions, polar molecules, or strong electric fields.
 */
public class DrudeForce extends Force {

  /**
   * Create a new DrudeForce.
   * <p>
   * This constructor initializes a new Drude force object that will handle
   * the interactions between Drude particles and their parent atoms, as well
   * as the electrostatic interactions involving Drude charges.
   */
  public DrudeForce() {
    super(OpenMM_DrudeForce_create());
  }

  /**
   * Add a Drude particle to the force.
   * <p>
   * This method adds a Drude oscillator consisting of a parent atom and its
   * associated Drude particle. The Drude particle is connected to the parent
   * atom by a harmonic spring with the specified force constant.
   *
   * @param particle1      The index of the parent atom to which the Drude particle is attached.
   * @param particle2      The index of the Drude particle.
   * @param particle3      The index of an additional particle (if applicable, -1 if not used).
   * @param particle4      The index of an additional particle (if applicable, -1 if not used).
   * @param particle5      The index of an additional particle (if applicable, -1 if not used).
   * @param charge         The charge of the Drude particle (typically negative).
   * @param polarizability The polarizability of the atom (in units of nm³).
   * @param aniso12        The anisotropy parameter for the 1-2 direction.
   * @param aniso34        The anisotropy parameter for the 3-4 direction.
   * @return The index of the added Drude particle pair.
   */
  public int addParticle(int particle1, int particle2, int particle3, int particle4, int particle5,
                         double charge, double polarizability, double aniso12, double aniso34) {
    return OpenMM_DrudeForce_addParticle(pointer, particle1, particle2, particle3, particle4, particle5,
        charge, polarizability, aniso12, aniso34);
  }

  /**
   * Add a screened pair to the force.
   * <p>
   * This method adds a screened pair interaction between two particles to prevent
   * overpolarization at short distances. The screening is typically applied between
   * atoms that are close in the molecular topology to avoid unphysical behavior.
   *
   * @param particle1 The index of the first particle in the screened pair.
   * @param particle2 The index of the second particle in the screened pair.
   * @param thole     The Thole screening parameter that controls the strength of screening.
   * @return The index of the added screened pair.
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
   * Get the number of Drude particle pairs in the force.
   * <p>
   * This method returns the total number of Drude oscillators (parent-Drude
   * particle pairs) that have been added to this force.
   *
   * @return The number of Drude particle pairs.
   */
  public int getNumParticles() {
    return OpenMM_DrudeForce_getNumParticles(pointer);
  }

  /**
   * Get the number of screened pairs in the force.
   * <p>
   * This method returns the total number of screened pair interactions that
   * have been added to this force to prevent overpolarization at short distances.
   *
   * @return The number of screened pairs.
   */
  public int getNumScreenedPairs() {
    return OpenMM_DrudeForce_getNumScreenedPairs(pointer);
  }

  /**
   * Get the parameters for a Drude particle pair.
   * <p>
   * This method retrieves the parameters that define a specific Drude oscillator,
   * including the particle indices, charge, polarizability, and anisotropy parameters.
   *
   * @param index          The index of the Drude particle pair to query.
   * @param particle1      The index of the parent atom (output).
   * @param particle2      The index of the Drude particle (output).
   * @param particle3      The index of an additional particle (output, -1 if not used).
   * @param particle4      The index of an additional particle (output, -1 if not used).
   * @param particle5      The index of an additional particle (output, -1 if not used).
   * @param charge         The charge of the Drude particle (output).
   * @param polarizability The polarizability of the atom (output).
   * @param aniso12        The anisotropy parameter for the 1-2 direction (output).
   * @param aniso34        The anisotropy parameter for the 3-4 direction (output).
   */
  public void getParticleParameters(int index, IntByReference particle1, IntByReference particle2,
                                    IntByReference particle3, IntByReference particle4, IntByReference particle5,
                                    DoubleByReference charge, DoubleByReference polarizability,
                                    DoubleByReference aniso12, DoubleByReference aniso34) {
    OpenMM_DrudeForce_getParticleParameters(pointer, index, particle1, particle2, particle3,
        particle4, particle5, charge, polarizability, aniso12, aniso34);
  }

  /**
   * Get the parameters for a screened pair.
   * <p>
   * This method retrieves the parameters that define a specific screened pair
   * interaction, including the particle indices and Thole screening parameter.
   *
   * @param index     The index of the screened pair to query.
   * @param particle1 The index of the first particle (output).
   * @param particle2 The index of the second particle (output).
   * @param thole     The Thole screening parameter (output).
   */
  public void getScreenedPairParameters(int index, IntByReference particle1, IntByReference particle2,
                                        DoubleByReference thole) {
    OpenMM_DrudeForce_getScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Get the parameters for a screened pair.
   * <p>
   * This method retrieves the parameters that define a specific screened pair
   * interaction, including the particle indices and Thole screening parameter.
   * This overloaded version uses IntBuffer and DoubleBuffer for output parameters.
   *
   * @param index     The index of the screened pair to query.
   * @param particle1 The index of the first particle (output).
   * @param particle2 The index of the second particle (output).
   * @param thole     The Thole screening parameter (output).
   */
  public void getScreenedPairParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                        DoubleBuffer thole) {
    OpenMM_DrudeForce_getScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Set the parameters for a Drude particle pair.
   * <p>
   * This method modifies the parameters of an existing Drude oscillator.
   * The particle indices typically should not be changed after creation,
   * but the charge, polarizability, and anisotropy parameters can be adjusted.
   *
   * @param index          The index of the Drude particle pair to modify.
   * @param particle1      The index of the parent atom.
   * @param particle2      The index of the Drude particle.
   * @param particle3      The index of an additional particle (-1 if not used).
   * @param particle4      The index of an additional particle (-1 if not used).
   * @param particle5      The index of an additional particle (-1 if not used).
   * @param charge         The charge of the Drude particle.
   * @param polarizability The polarizability of the atom.
   * @param aniso12        The anisotropy parameter for the 1-2 direction.
   * @param aniso34        The anisotropy parameter for the 3-4 direction.
   */
  public void setParticleParameters(int index, int particle1, int particle2, int particle3,
                                    int particle4, int particle5, double charge, double polarizability,
                                    double aniso12, double aniso34) {
    OpenMM_DrudeForce_setParticleParameters(pointer, index, particle1, particle2, particle3,
        particle4, particle5, charge, polarizability, aniso12, aniso34);
  }

  /**
   * Set the parameters for a screened pair.
   * <p>
   * This method modifies the parameters of an existing screened pair interaction.
   * The particle indices typically should not be changed after creation,
   * but the Thole screening parameter can be adjusted.
   *
   * @param index     The index of the screened pair to modify.
   * @param particle1 The index of the first particle.
   * @param particle2 The index of the second particle.
   * @param thole     The Thole screening parameter.
   */
  public void setScreenedPairParameters(int index, int particle1, int particle2, double thole) {
    OpenMM_DrudeForce_setScreenedPairParameters(pointer, index, particle1, particle2, thole);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   * <p>
   * This method controls whether the DrudeForce takes into account periodic boundary
   * conditions when calculating interactions. This is typically important for systems
   * in periodic boxes to ensure proper treatment of long-range electrostatic interactions.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   * <p>
   * This method updates the Drude force parameters in the specified OpenMM context.
   * This is necessary when parameters have been modified after the context was created.
   * The context must have been created with this force included in the system.
   *
   * @param context The OpenMM Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_DrudeForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   * <p>
   * This method returns whether the DrudeForce takes into account periodic
   * boundary conditions when calculating interactions. For Drude forces,
   * this typically depends on the underlying electrostatic treatment.
   *
   * @return True if the force uses periodic boundary conditions, false otherwise.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_DrudeForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}