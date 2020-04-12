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
package ffx.algorithms.optimize.anneal;

/**
 * Temperature schedule for simulated annealing
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface AnnealingSchedule {
    /**
     * Gets the starting temperature.
     *
     * @return Starting temperature in Kelvin.
     */
    double getHighTemp();

    /**
     * Gets the final temperature.
     *
     * @return Final temperature in Kelvin.
     */
    double getLowTemp();

    /**
     * Gets the number of annealing windows (including repeat windows).
     *
     * @return Number of annealing windows.
     */
    int getNumWindows();

    /**
     * Get the temperature for annealing step i.
     *
     * @param i An annealing step [0-nWindows)
     * @return Associated temperature.
     */
    double getTemperature(int i);

    /**
     * Get all temperatures this schedule specifies.
     *
     * @return An array of temperatures specified.
     */
    double[] getTemperatures();

    /**
     * Returns the longest window to be used (normalized to the number of MD steps in a "regular" window).
     *
     * @return Maximum normalized window length.
     */
    double maxWindowLength();

    /**
     * Returns the shortest window to be used (normalized to the number of MD steps in a "regular" window).
     *
     * @return Minimum normalized window length.
     */
    double minWindowLength();

    /**
     * Returns the sum of window lengths to be used (normalized to the number of MD steps in a "regular" window).
     *
     * @return Total normalized window length.
     */
    double totalWindowLength();

    /**
     * Get the relative size of a window (normalized to the number of MD steps in a "regular" window).
     *
     * @param window Window to check.
     * @return Normalized length of the window.
     */
    double windowLength(int window);
}
