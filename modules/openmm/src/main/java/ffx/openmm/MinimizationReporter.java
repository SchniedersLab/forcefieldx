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
package ffx.openmm;

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MinimizationReporter_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MinimizationReporter_destroy;

/**
 * A MinimizationReporter can be passed to LocalEnergyMinimizer::minimize() to provide
 * periodic information on the progress of minimization, and to give you the chance to
 * stop minimization early.  Define a subclass that overrides report() and implement it
 * to take whatever action you want.
 * <p>
 * To correctly interpret the information passed to the reporter, you need to know a bit
 * about how the minimizer works.  The L-BFGS algorithm used by the minimizer does not
 * support constraints.  The minimizer therefore replaces all constraints with harmonic
 * restraints, then performs unconstrained minimization of a combined objective function
 * that is the sum of the system's potential energy and the restraint energy.  Once
 * minimization completes, it checks whether all constraints are satisfied to an acceptable
 * tolerance.  It not, it increases the strength of the harmonic restraints and performs
 * additional minimization.  If the error in constrained distances is especially large,
 * it may choose to throw out all work that has been done so far and start over with
 * stronger restraints.  This has several important consequences.
 *
 * <ul>
 * <li>The objective function being minimized not actually the same as the potential energy.</li>
 * <li>The objective function and the potential energy can both increase between iterations.</li>
 * <li>The total number of iterations performed could be larger than the number specified
 *     by the maxIterations argument, if that many iterations leaves unacceptable constraint errors.</li>
 * <li>All work is provisional.  It is possible for the minimizer to throw it out and start over.</li>
 * </ul>
 */
public class MinimizationReporter {

  /**
   * Reporter pointer.
   */
  private PointerByReference pointer;

  /**
   * Constructor.
   */
  public MinimizationReporter() {
    pointer = OpenMM_MinimizationReporter_create();
  }

  /**
   * Pointer to the reporter.
   *
   * @return the pointer.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Destroy the reporter.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_MinimizationReporter_destroy(pointer);
      pointer = null;
    }
  }

}
