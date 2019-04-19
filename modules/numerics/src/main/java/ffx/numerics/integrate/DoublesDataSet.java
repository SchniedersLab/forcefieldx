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
package ffx.numerics.integrate;

import static java.lang.String.format;
import static java.lang.System.arraycopy;

import static ffx.numerics.integrate.FunctionDataCurve.approxEquals;

/**
 * Descibes a set of x, f(x) obtained by some mechanism; intended for numerical
 * integration.
 *
 * @author Jacob M. Litman
 */
public class DoublesDataSet implements DataSet {

    /**
     * An array of input points.
     */
    private final double[] x;
    /**
     * An array of results f(x).
     */
    private final double[] fX;
    /**
     * Lower bound.
     */
    private final double lb;
    /**
     * Upper bound.
     */
    private final double ub;
    /**
     * Number of data points.
     */
    private final int nX;
    /**
     * Bin size.
     */
    private final double sep;
    /**
     * If ends should have 1/2 regular separation.
     */
    private final boolean halfWidthEnd;

    /**
     * Constructs a DataSet from actual data, with no known underlying function (or at least none with an analytically
     * solved integral). Assumes no half-width end bins (such as found in OSRW).
     *
     * @param x  Points where f(x) is known
     * @param fX Values/estimates of f(x)
     */
    public DoublesDataSet(double[] x, double[] fX) {
        this(x, fX, false);
    }

    /**
     * Constructs a DataSet from actual data, with no known underlying function (or at least none with an analytically
     * solved integral).
     *
     * @param x          Points where f(x) is known
     * @param fX         Values/estimates of f(x)
     * @param halvedEnds Whether the first and last bins are half-width (such as OSRW)
     */
    public DoublesDataSet(double[] x, double[] fX, boolean halvedEnds) {
        nX = x.length;
        assert nX == fX.length;

        this.x = new double[nX];
        arraycopy(x, 0, this.x, 0, nX);

        this.fX = new double[nX];
        arraycopy(fX, 0, this.fX, 0, nX);

        lb = x[0];
        ub = x[nX - 1];
        //sep = ((ub - lb) / ((double) nX));
        halfWidthEnd = halvedEnds;
        double sepDist = ub - lb;
        sep = halfWidthEnd ? (sepDist / ((double) nX - 2)) : (sepDist / ((double) nX - 1));
        assertXIntegrity(this.x);
    }

    /**
     * Constructs a DataSet from another DataSet, effectively masquerading a test set such as a sine wave as data from
     * an "unknown" function. Used primarily for testing purposes.
     *
     * @param set DataSet to cast
     */
    public DoublesDataSet(DataSet set) {
        nX = set.numPoints();

        this.x = new double[nX];
        arraycopy(set.getX(), 0, this.x, 0, nX);

        this.fX = new double[nX];
        arraycopy(set.getAllFxPoints(), 0, this.fX, 0, nX);

        lb = x[0];
        ub = x[nX - 1];
        this.halfWidthEnd = set.halfWidthEnds();
        sep = set.binWidth();
        // Possibly add assertion for separation distance
        assertXIntegrity(x);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double lowerBound() {
        return lb;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double upperBound() {
        return ub;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numPoints() {
        return x.length;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double binWidth() {
        return sep;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFxPoint(int index) {
        return fX[index];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getAllFxPoints() {
        double[] pts = new double[nX];
        arraycopy(fX, 0, pts, 0, nX);
        return pts;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean halfWidthEnds() {
        return halfWidthEnd;
    }

    /**
     * Used to check that the passed-in x array is composed of equally-spaced
     * points from lb to ub.
     *
     * @param x
     */
    private void assertXIntegrity(double[] x) {
        assert ub > lb;
        if (halfWidthEnd) {
            assert x.length == nX;
            assert lb == x[0];
            assert ub == x[nX - 1];

            assert approxEquals(x[1], lb + 0.5 * sep);
            assert approxEquals(x[nX - 2], (ub - 0.5 * sep));

            for (int i = 2; i < (nX - 2); i++) {
                double target = lb + 0.5 * sep;
                target += ((i - 1) * sep);
                assert approxEquals(x[i], target);
            }
        } else {
            for (int i = 0; i < x.length; i++) {
                assert approxEquals(x[i], x[0] + i * sep);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getX() {
        double[] copyX = new double[x.length];
        arraycopy(x, 0, copyX, 0, x.length);
        return copyX;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(format("Data set with %d points from lower bound %9.3g and upper bound %9.3g", nX, lb, ub));
        if (halfWidthEnd) {
            sb.append(" and half-width start/end bins");
        }
        sb.append(".");
        return sb.toString();
    }
}
