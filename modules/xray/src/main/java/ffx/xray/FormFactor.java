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

package ffx.xray;

import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>FormFactor interface.</p>
 *
 * @author Tim Fenn
 * @version $Id: $
 */
public interface FormFactor {
    /**
     * Compute the real space density rho
     *
     * @param f the current density to modify
     * @param lambda the state variable
     * @param xyz the requested point for evaluating density
     * @return the real space density value at xyz
     */
    double rho(double f, double lambda, double xyz[]);

    /**
     * Compute the real space gradient
     *
     * @param xyz the requested point for evaluating gradient
     * @param dfc the multiplier to apply to the gradient
     * @param refinementmode {@link RefinementMinimize.RefinementMode}
     * determines which gradients will be computed
     */
    void rho_grad(double xyz[], double dfc, RefinementMode refinementmode);

    /**
     * update the coordinates to the current position
     *
     * @param xyz an array of double.
     */
    void update(double xyz[]);

    /**
     * update the coordinates to the current position and Badd
     *
     * @param xyz an array of double.
     * @param badd a double.
     */
    void update(double xyz[], double badd);
}
