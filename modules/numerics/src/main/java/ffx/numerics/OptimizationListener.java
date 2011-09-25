/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.numerics;

import ffx.numerics.LineSearch.LineSearchResult;

/**
 * This interface allows the optimizer to notify registered instances of
 * successful steps. Currently only one listener is supported.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public interface OptimizationListener {

    /**
     * This method is called by the optimizer after each step.
     *
     * It can be used to log status messages, update the user interface, or
     * gracefully terminate the optimizer.
     *
     * @param iter Number of iterations.
     * @param nfun Number of function evaluations.
     * @param grms RMS gradient at current solution.
     * @param xrms RMS coordinate change at current solution.
     * @param f Function value at current solution.
     * @param df Change in the function value compared to the previous solution.
     * @param angle Current angle between gradient and search direction.
     * @param info Result of the line search (null at iteraction == 0).
     * @return A return value of false will terminate the optimization.
     * @since 1.0
     */
    public abstract boolean optimizationUpdate(int iter, int nfun, double grms,
            double xrms, double f, double df, double angle, LineSearchResult info);

}
