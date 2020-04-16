// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;
import java.io.File;

/**
 * AlgorithmFunctions, on top of the core functionality of PotentialsFunctions, describes additional
 * functionality such as molecular dynamics and L-BFGS local optimization.
 *
 * <p>This is implemented in two locations: UIUtils in the User Interfaces package, and in
 * AlgorithmsUtils in the Algorithm package.
 *
 * <p>The UIUtils implementation is the default for Force Field X; on top of the core functionality,
 * it also updates the FFX graphical user interface and tree structure.
 *
 * <p>The AlgorithmUtils implementation lacks the extra functionality of the UIUtils implementation,
 * and simply accomplishes the required task. This is used by our tests, and is also potentially
 * useful for third parties who would like to use FFX without its GUI.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface AlgorithmFunctions extends PotentialsFunctions {

  /**
   * Returns a default Listener if available (null by default).
   *
   * @return An AlgorithmListener or null.
   */
  default AlgorithmListener getDefaultListener() {
    return null;
  }

  /**
   * Runs molecular dynamics on an assembly.
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param nStep Timesteps
   * @param timeStep Time per step
   * @param printInterval a double.
   * @param saveInterval a double.
   * @param temperature a double.
   * @param initVelocities Initialize velocities from Maxwell-Boltzmann distribution
   * @param dyn Dynamics file
   */
  void md(
      MolecularAssembly assembly,
      int nStep,
      double timeStep,
      double printInterval,
      double saveInterval,
      double temperature,
      boolean initVelocities,
      File dyn);

  /**
   * Relax the coordinates of a MolecularAssembly and minimize its potential energy
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param eps RMS gradient convergence criteria
   * @return A <code>Potential</code>
   */
  Potential minimize(MolecularAssembly assembly, double eps);
}
