/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.xray;

import java.util.logging.Logger;

import static java.util.Arrays.fill;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;

public class RowSchedule extends IntegerSchedule {

    private static final Logger logger = Logger.getLogger(RowSchedule.class.getName());

    private int nThreads;
    private boolean threadDone[];
    private Range ranges[];
    private final int lowerBounds[];
    private final int fftZ;
    private final int fftY;
    private int weights[];

    protected RowSchedule(int nThreads, int fftZ, int fftY) {
        this.nThreads = nThreads;
        threadDone = new boolean[nThreads];
        ranges = new Range[nThreads];
        lowerBounds = new int[nThreads + 1];
        this.fftY = fftY;
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
        for (int i = 0; i < fftZ * fftY; i++) {
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
            Range temp = new Range(0, fftZ*fftY - 1);
            ranges = temp.subranges(nThreads);
            return;
        }

        /**
         * Handle the case where we only have a single thread, which will
         * receive all the rows.
         */
        if (nThreads == 1) {
            ranges[0] = new Range(0, fftZ*fftY - 1);
            return;
        }

        double targetWeight = (totalWeight / nThreads) * .95;
        int lastRow = fftZ*fftY - 1;

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
        ranges[lastThread - 1] = new Range(lowerBounds[lastThread - 1], lastRow);
        //logger.info(String.format("Range for thread %d %s.", lastThread - 1, ranges[lastThread - 1]));

        /**
         * Left-over threads with null ranges.
         */
        for (int it = lastThread; it < nThreads; it++) {
            ranges[it] = null;
        }
    }

    public int[] getWeightPerThread() {
        if (lowerBounds != null) {
            int[] weightToReturn = new int[nThreads];
            for (int i = 0; i < nThreads; i++) {
                weightToReturn[i] = lowerBounds[i + 1];
            }
            return weightToReturn;
        } else {
            return null;
        }
    }
}
