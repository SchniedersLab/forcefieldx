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
package ffx.numerics.integrate;

/**
 * A DataSet represents a set of points along a single dimension, and is able
 * to be numerically integrated.
 *
 * @author Jacob M. Litman
 */
public interface DataSet {

    /**
     * Separation between points along x; should be uniform.
     *
     * @return a double.
     */
    double binWidth();

    /**
     * Returns copy of the array of points f(x) to integrate (y-axis).
     *
     * @return an array of {@link double} objects.
     */
    double[] getAllFxPoints();

    /**
     * Point f(x) at index.
     *
     * @param index a int.
     * @return a double.
     */
    double getFxPoint(int index);

    /**
     * Returns copy of the array of points x (x-axis).
     *
     * @return an array of {@link double} objects.
     */
    double[] getX();

    /**
     * Does this data set have half-width start/end bins. Intended for OST,
     * where the first and last bins are half the regular width.
     *
     * @return a boolean.
     */
    boolean halfWidthEnds();

    /**
     * Lower bound of the points along x.
     *
     * @return a double.
     */
    double lowerBound();

    /**
     * Number of points along x.
     *
     * @return a int.
     */
    int numPoints();

    /**
     * Upper bound of the points along x.
     *
     * @return a double.
     */
    double upperBound();
}
