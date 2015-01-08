/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;

/**
 * A fixed schedule that balances pairwise work across threads.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class PairwiseSchedule extends IntegerSchedule {

    private static final Logger logger = Logger.getLogger(PairwiseSchedule.class.getName());
    private int nAtoms;
    private final int nThreads;
    private final Range ranges[];
    private final boolean threadDone[];

    /**
     * <p>
     * Constructor for PairwiseSchedule.</p>
     *
     * @param nThreads a int.
     * @param nAtoms a int.
     * @param ranges an array of {@link edu.rit.util.Range} objects.
     */
    public PairwiseSchedule(int nThreads, int nAtoms, Range ranges[]) {
        this.nAtoms = nAtoms;
        this.nThreads = nThreads;
        this.ranges = ranges;
        threadDone = new boolean[nThreads];
    }

    public void setAtoms(int nAtoms) {
        this.nAtoms = nAtoms;
    }

    /**
     * {@inheritDoc}
     *
     * This is a fixed schedule.
     */
    @Override
    public boolean isFixedSchedule() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void start(int nThreads, Range chunkRange) {
        assert (nThreads == this.nThreads);
        assert (chunkRange.lb() == 0);
        assert (chunkRange.ub() == nAtoms - 1);

        for (int i = 0; i < nThreads; i++) {
            threadDone[i] = false;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Range next(int threadID) {
        if (!threadDone[threadID]) {
            threadDone[threadID] = true;
            return ranges[threadID];
        }
        return null;
    }

    /**
     * <p>
     * updateRanges</p>
     *
     * @param totalInteractions a int.
     * @param listCount an array of int.
     */
    public void updateRanges(int totalInteractions, int listCount[]) {
        int id = 0;
        int goal = totalInteractions / nThreads;
        int num = 0;
        int start = 0;
        for (int i = 0; i < nAtoms; i++) {
            num += listCount[i];
            if (num >= goal) {
                /**
                 * Last thread gets the remaining atoms.
                 */
                if (id == nThreads - 1) {
                    ranges[id] = new Range(start, nAtoms - 1);
                    break;
                }

                ranges[id] = new Range(start, i);

                // Zero out the interaction counter.
                num = 0;
                // Next thread.
                id++;
                // Next range starts at i+1.
                start = i + 1;

                /**
                 * Out of atoms. Threads remaining get a null range.
                 */
                if (start == nAtoms) {
                    for (int j = id; j < nThreads; j++) {
                        ranges[j] = null;
                    }
                    break;
                }
            } else if (i == nAtoms - 1) {
                /**
                 * Last atom without reaching goal for current thread.
                 */
                ranges[id] = new Range(start, nAtoms - 1);
                for (int j = id + 1; j < nThreads; j++) {
                    ranges[j] = null;
                }
            }
        }
    }
}
