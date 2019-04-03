//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.xray;

import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>
 * FormFactor interface.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public interface FormFactor {

    /**
     * Compute the real space density rho
     *
     * @param f      the current density to modify
     * @param lambda the state variable
     * @param xyz    the requested point for evaluating density
     * @return the real space density value at xyz
     */
    double rho(double f, double lambda, double[] xyz);

    /**
     * Compute the real space gradient
     *
     * @param xyz            the requested point for evaluating gradient
     * @param dfc            the multiplier to apply to the gradient
     * @param refinementmode {@link ffx.xray.RefinementMinimize.RefinementMode}
     *                       determines which gradients will be computed
     */
    void rhoGrad(double[] xyz, double dfc, RefinementMode refinementmode);

    /**
     * update the coordinates to the current position
     *
     * @param xyz an array of double.
     */
    void update(double[] xyz);

    /**
     * update the coordinates to the current position and Badd
     *
     * @param xyz  an array of double.
     * @param badd a double.
     */
    void update(double[] xyz, double badd);
}
