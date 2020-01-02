//******************************************************************************
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
//******************************************************************************
package ffx.numerics.optimization;

import ffx.numerics.optimization.LineSearch.LineSearchResult;

/**
 * This interface allows the optimizer to notify registered instances of
 * successful steps. Currently only one listener is supported.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface OptimizationListener {

    /**
     * This method is called by the optimizer after each step.
     * <p>
     * It can be used to log status messages, update the user interface, or
     * gracefully terminate the optimizer.
     *
     * @param iter  Number of iterations.
     * @param nfun  Number of function evaluations.
     * @param grms  RMS gradient at current solution.
     * @param xrms  RMS coordinate change at current solution.
     * @param f     Function value at current solution.
     * @param df    Change in the function value compared to the previous solution.
     * @param angle Current angle between gradient and search direction.
     * @param info  Result of the line search (null at iteraction == 0).
     * @return A return value of false will terminate the optimization.
     * @since 1.0
     */
    boolean optimizationUpdate(int iter, int nfun, double grms,
                               double xrms, double f, double df, double angle, LineSearchResult info);
}
