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
package ffx.potential.openmm;

import ffx.potential.bonded.Atom;

import javax.annotation.Nullable;

/**
 * An interface for classes that provide an OpenMM potential energy implementation.
 */
public interface OpenMMPotential {

  /**
   * Returns the Context instance.
   *
   * @return context
   */
  OpenMMContext getContext();

  /**
   * Update the OpenMM Context.
   *
   * @param integratorName Integrator to use.
   * @param timeStep       Time step.
   * @param temperature    Temperature (K).
   * @param forceCreation  Force a new Context to be created, even if the existing one matches the
   *                       request.
   */
  void updateContext(String integratorName, double timeStep, double temperature, boolean forceCreation);

  /**
   * Create an immutable OpenMM State.
   *
   * <p>State.free() must be called to free OpenMM memory.
   *
   * @param mask The State mask.
   * @return Returns the State.
   */
  OpenMMState getOpenMMState(int mask);

  /**
   * Get a reference to the System instance.
   *
   * @return a reference to the OpenMMSystem.
   */
  OpenMMSystem getSystem();

  /**
   * Update active atoms.
   *
   * @return True if there are inactive atoms.
   */
  boolean setActiveAtoms();

  /**
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  void updateParameters(@Nullable Atom[] atoms);
}
