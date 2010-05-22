/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.potential.nonbonded;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;
import java.util.logging.Logger;

/**
 * A fixed schedule that load balances work chunks across threads.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class SpatialDensitySchedule extends IntegerSchedule {

    private static Logger logger = Logger.getLogger(SpatialDensitySchedule.class.getName());
    private final int nThreads;
    private final int nAtoms;
    private final int atomsPerChunk[];
    private final double loadBalancePercentage;
    private final boolean threadDone[];
    private final Range ranges[];

    public SpatialDensitySchedule(int nThreads, int nAtoms,
            int atomsPerChunk[], double loadBalancePercentage) {
        this.nAtoms = nAtoms;
        this.atomsPerChunk = atomsPerChunk;
        this.nThreads = nThreads;
        ranges = new Range[nThreads];
        threadDone = new boolean[nThreads];

        if (loadBalancePercentage > 0.01 && loadBalancePercentage <= 1.0) {
            this.loadBalancePercentage = loadBalancePercentage;
        } else  {
            this.loadBalancePercentage = 1.0;
        }
    }

    /**
     * This is a fixed schedule.
     * @return true
     */
    @Override
    public boolean isFixedSchedule() {
        return true;
    }

    @Override
    public void start(int nThreads, Range range) {

        assert(nThreads == this.nThreads);

        int lb = range.lb();
        int ub = range.ub();

        // Null out the thread ranges.
        for (int i=0; i < nThreads; i++) {
            ranges[i] = null;
            threadDone[i] = false;
        }

        int threadID = 0;
        int start = 0;
        int total = 0;

        int goal = (int) (nAtoms * loadBalancePercentage / nThreads);
        for (int i=lb; i<=ub; i++) {
            // Count the number of atoms in each work chunk.
            total += atomsPerChunk[i];
            // Check if the load balancing goal has been reached.
            if (total > goal) {
                int stop = i;
                // Define the range for this thread.
                ranges[threadID++] = new Range(start, stop);
                // Initialization for all, but the last thread.
                if (threadID < nThreads - 1) {
                    start = i+1;
                    total = 0;
                } else {
                    // The last thread gets the rest of the work chunks.
                    start = i + 1;
                    stop = ub;
                    ranges[threadID] = new Range(start, stop);
                    break;
                }
            }
        }
    }

    @Override
    public Range next(int threadID) {
        if (!threadDone[threadID]) {
            threadDone[threadID] = true;
            return ranges[threadID];
        }
        return null;
    }
}
