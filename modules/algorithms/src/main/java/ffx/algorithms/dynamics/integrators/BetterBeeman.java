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
package ffx.algorithms.dynamics.integrators;

import ffx.numerics.Potential;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;

/**
 * Integrate Newton's equations of motion using a Beeman multistep recursion
 * formula; the actual coefficients are Brooks' "Better Beeman" values.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BetterBeeman extends Integrator {

    private double dt2_8;
    private double dt_8;

    /**
     * Constructor for BetterBeeman.
     *
     * @param nVariables number of Variables.
     * @param x          Cartesian coordinates (Angstroms).
     * @param v          Velocities.
     * @param a          Accelerations.
     * @param aPrevious  Previous Accelerations.
     * @param mass       Mass.
     */
    public BetterBeeman(int nVariables, double[] x, double[] v, double[] a,
                        double[] aPrevious, double[] mass) {
        super(nVariables, x, v, a, aPrevious, mass);
        dt_8 = 0.125 * dt;
        dt2_8 = dt * dt_8;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Store the current atom positions, then find new atom positions and
     * half-step velocities via Beeman recursion.
     */
    @Override
    public void preForce(Potential potential) {
        for (int i = 0; i < nVariables; i++) {
            double temp = 5.0 * a[i] - aPrevious[i];
            x[i] += v[i] * dt + temp * dt2_8;
            v[i] += temp * dt_8;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Use Newton's second law to get the next acceleration and find the
     * full-step velocities using the Beeman recusion.
     */
    @Override
    public void postForce(double[] gradient) {
        copyAccelerationToPrevious();
        for (int i = 0; i < nVariables; i++) {
            a[i] = -KCAL_TO_GRAM_ANG2_PER_PS2 * gradient[i] / mass[i];
            v[i] += (3.0 * a[i] + aPrevious[i]) * dt_8;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setTimeStep(double dt) {
        this.dt = dt;
        dt_8 = 0.125 * dt;
        dt2_8 = dt * dt_8;
    }
}
