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
package ffx.numerics.switching;

/**
 * The ConstantSwitch returns 1 for all input values x. This is useful for
 * having a single code path that accomodates both "real" switching behavior and no
 * switching behavior.
 *
 *  @author Jacob M. Litman
 *  @author Michael J. Schnieders
 */
public class ConstantSwitch implements UnivariateSwitchingFunction {

    /**
     * {@inheritDoc}
     */
    @Override
    public double getZeroBound() {
        return Double.NaN;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getOneBound() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean constantOutsideBounds() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean validOutsideBounds() {
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getHighestOrderZeroDerivative() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean symmetricToUnity() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        return 1.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double firstDerivative(double x) {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double secondDerivative(double x) {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        if (order < 0) {
            throw new IllegalArgumentException(String.format(" Order must be > 0, was %d", order));
        }
        return 0.0;
    }

    @Override
    public String toString() {
        return "Constant-value f(x) = 1, with no switching behavior (i.e. a dummy switch)";
    }
}
