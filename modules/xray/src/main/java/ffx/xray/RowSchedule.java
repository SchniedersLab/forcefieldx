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
package ffx.xray;

import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;

/**
 * <p>RowSchedule class.</p>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RowSchedule extends IntegerSchedule {

    private final int[] lowerBounds;
    private final int fftZ;
    private final int fftY;
    private int nThreads;
    private boolean[] threadDone;
    private Range[] ranges;
    private int[] weights;

    /**
     * <p>Constructor for RowSchedule.</p>
     *
     * @param nThreads a int.
     * @param fftZ     a int.
     * @param fftY     a int.
     */
    protected RowSchedule(int nThreads, int fftZ, int fftY) {
        this.nThreads = nThreads;
        threadDone = new boolean[nThreads];
        ranges = new Range[nThreads];
        lowerBounds = new int[nThreads + 1];
        this.fftY = fftY;
        this.fftZ = fftZ;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isFixedSchedule() {
        return true;
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
     * {@inheritDoc}
     */
    @Override
    public void start(int nThreads, Range chunkRange) {
        this.nThreads = nThreads;

        if (nThreads != threadDone.length) {
            threadDone = new boolean[nThreads];
        }
        fill(threadDone, false);

        if (nThreads != ranges.length) {
            ranges = new Range[nThreads];
        }
        fill(lowerBounds, 0);
        defineRanges();
    }

    /**
     * <p>updateWeights.</p>
     *
     * @param weights an array of {@link int} objects.
     */
    void updateWeights(int[] weights) {
        this.weights = weights;
    }

    private int totalWeight() {
        int totalWeight = 0;
        for (int i = 0; i < fftZ * fftY; i++) {
            totalWeight += weights[i];
        }
        return totalWeight;
    }

    private void defineRanges() {
        double totalWeight = totalWeight();

        // Infrequent edge case where the total weight is less than or equal to the number of threads.
        if (totalWeight <= nThreads) {
            Range temp = new Range(0, fftZ * fftY - 1);
            ranges = temp.subranges(nThreads);
            return;
        }

        // Handle the case where we only have a single thread, which will receive all the rows.
        if (nThreads == 1) {
            ranges[0] = new Range(0, fftZ * fftY - 1);
            return;
        }

        double targetWeight = (totalWeight / nThreads);
        int lastRow = fftZ * fftY - 1;

        int currentRow = 0;
        lowerBounds[0] = 0;
        int currentThread = 0;
        while (currentThread < nThreads) {
            int threadWeight = 0;
            while (threadWeight < targetWeight && currentRow < lastRow) {
                threadWeight += weights[currentRow];
                currentRow++;
            }
            currentThread++;
            if (currentRow < lastRow) {
                lowerBounds[currentThread] = currentRow;
            } else {
                lowerBounds[currentThread] = lastRow;
                break;
            }
        }


        int lastThread = currentThread;

        // Loop over all threads that will receive work except the final one.
        for (currentThread = 0; currentThread < lastThread - 1; currentThread++) {
            ranges[currentThread] = new Range(lowerBounds[currentThread], lowerBounds[currentThread + 1] - 1);
            //logger.info(String.format("Range for thread %d %s.", currentThread, ranges[currentThread]));

        }

        // Final range for the last thread that will receive work.
        ranges[lastThread - 1] = new Range(lowerBounds[lastThread - 1], lastRow);
        //logger.info(String.format("Range for thread %d %s.", lastThread - 1, ranges[lastThread - 1]));

        // Left-over threads with null ranges.
        for (int it = lastThread; it < nThreads; it++) {
            ranges[it] = null;
        }
    }

    /**
     * <p>getThreadWeights.</p>
     *
     * @return an array of {@link int} objects.
     */
    int[] getThreadWeights() {
        int[] weightsToReturn = new int[nThreads];
        arraycopy(weights, 0, weightsToReturn, 0, nThreads);
        return weightsToReturn;
    }

    /**
     * <p>Getter for the field <code>lowerBounds</code>.</p>
     *
     * @return an array of {@link int} objects.
     */
    int[] getLowerBounds() {
        int[] boundsToReturn = new int[nThreads];
        arraycopy(lowerBounds, 1, boundsToReturn, 0, nThreads);
        return boundsToReturn;
    }

}
