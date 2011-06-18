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
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;

/**
 * A fixed schedule that balances pairwise work across threads.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class PairwiseSchedule extends IntegerSchedule {

    private static final Logger logger = Logger.getLogger(PairwiseSchedule.class.getName());
    private final int nAtoms;
    private int nThreads;
    private Range ranges[];
    private boolean threadDone[];

    public PairwiseSchedule(int nThreads, int nAtoms) {
        this.nAtoms = nAtoms;
        this.nThreads = nThreads;
        ranges = new Range[nThreads];
        threadDone = new boolean[nThreads];
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
    public void start(int nThreads, Range chunkRange) {
        assert (nThreads == this.nThreads);
        assert (chunkRange.lb() == 0);
        assert (chunkRange.ub() == nAtoms - 1);

        for (int i = 0; i < nThreads; i++) {
            threadDone[i] = false;
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

    public void setRanges(Range ranges[]) {
        this.ranges = ranges;
    }
}
