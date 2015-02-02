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
package ffx.xray;

import java.util.logging.Logger;

import static java.util.Arrays.fill;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;

/**
 * @author Armin Avdic
 */
public class SliceSchedule extends IntegerSchedule {

    private static final Logger logger = Logger.getLogger(SliceSchedule.class.getName());

    private int nThreads;
    private boolean threadDone[];
    private Range ranges[];
    private final int lowerBounds[];
    private final int fftZ;
    private int weights[];

    protected SliceSchedule(int nThreads, int fftZ) {
        this.nThreads = nThreads;
        threadDone = new boolean[nThreads];
        ranges = new Range[nThreads];
        lowerBounds = new int[nThreads + 1];
        this.fftZ = fftZ;
    }

    public void updateWeights(int weights[]) {
        this.weights = weights;
    }

    @Override
    public boolean isFixedSchedule() {
        return true;
    }

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

    @Override
    public Range next(int threadID) {
        if (!threadDone[threadID]) {
            threadDone[threadID] = true;
            return ranges[threadID];
        }
        return null;
    }

    private int totalWeight() {
        int totalWeight = 0;
        for (int i = 0; i < fftZ; i++) {
            totalWeight += weights[i];
        }
        return totalWeight;
    }

    private void defineRanges() {
        double totalWeight = totalWeight();

        /**
         * Infrequent edge case where the total weight is less than or equal to
         * the number of threads.
         */
        if (totalWeight <= nThreads) {
            Range temp = new Range(0, fftZ - 1);
            ranges = temp.subranges(nThreads);
            return;
        }

        /**
         * Handle the case where we only have a single thread, which will
         * receive all the slices.
         */
        if (nThreads == 1) {
            ranges[0] = new Range(0, fftZ - 1);
            return;
        }

        double targetWeight = (totalWeight / nThreads) * .95;
        int lastSlice = fftZ - 1;

        int currentSlice = 0;
        lowerBounds[0] = 0;
        int currentThread = 0;
        while (currentThread < nThreads) {
            int threadWeight = 0;
            while (threadWeight < targetWeight && currentSlice < lastSlice) {
                threadWeight += weights[currentSlice];
                currentSlice++;
            }
            currentThread++;
            if (currentSlice < lastSlice) {
                lowerBounds[currentThread] = currentSlice;
            } else {
                lowerBounds[currentThread] = lastSlice;
                break;
            }
        }

        int lastThread = currentThread;

        /**
         * Loop over all threads that will receive work except the final one.
         */
        for (currentThread = 0; currentThread < lastThread - 1; currentThread++) {
            ranges[currentThread] = new Range(lowerBounds[currentThread], lowerBounds[currentThread + 1] - 1);
            //logger.info(String.format("Range for thread %d %s.", currentThread, ranges[currentThread]));
        }
        /**
         * Final range for the last thread that will receive work.
         */
        ranges[lastThread - 1] = new Range(lowerBounds[lastThread - 1], lastSlice);
        //logger.info(String.format("Range for thread %d %s.", lastThread - 1, ranges[lastThread - 1]));

        /**
         * Left-over threads with null ranges.
         */
        for (int it = lastThread; it < nThreads; it++) {
            ranges[it] = null;
        }
    }
    public int[] getThreadWeights(){
        if(lowerBounds != null){
            int[] weightsToReturn = new int[nThreads];
            for (int i = 0; i<nThreads; i++){
                weightsToReturn[i] = weights[i];
            }
            return weightsToReturn;
        } else {
            return null;
        }
    }

    public int[] getLowerBounds() {
        if (lowerBounds != null) {
            int[] boundsToReturn = new int[nThreads];
            for (int i = 0; i < nThreads; i++) {
                boundsToReturn[i] = lowerBounds[i + 1];
            }
            return boundsToReturn;
        } else {
            return null;
        }
    }
}
