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

/**
 *
 * @author avdic
 */
public class SliceSchedule extends IntegerSchedule {

    private static final Logger logger = Logger.getLogger(SliceSchedule.class.getName());

    private int nThreads;
    private boolean threadDone[];
    private Range ranges[];
    private final int intervals[];
    private final int ub[];
    private final int lb[];
    private final int fftZ;
    private int weights[];

    protected SliceSchedule(int nThreads, int fftZ) {
        this.nThreads = nThreads;
        threadDone = new boolean[nThreads];
        ranges = new Range[nThreads];
        intervals = new int[nThreads + 1];
        ub = new int[nThreads];
        lb = new int[nThreads];
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
        fill(intervals, 0);
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
        if (totalWeight > nThreads) {
            double targetWeight = (totalWeight / (nThreads)) * .95;
            int j = 0;
            boolean quit = false;
            intervals[0] = 0;
            int i = 0;

            while (i < (nThreads) && !quit) {
                int threadWeight = 0;
                while (threadWeight < targetWeight && j < (fftZ - 1)) {
                    threadWeight += weights[j];
                    j++;
                }
                if (j < (fftZ - 1)) {
                    intervals[i + 1] = j;
                } else {
                    quit = true;
                }
                i++;
            }

            /**
             * Check if final slices remain to be assigned.
             */
            boolean terminatorFound = false;
            int terminator = 1;

            while (!terminatorFound) {
                if (intervals[terminator] != 0) {
                    terminator++;
                } else {
                    terminatorFound = true;
                }

                if (terminator == nThreads) {
                    terminatorFound = true;
                }
            }

            intervals[terminator] = fftZ - 1;

            int iThreads = 0;
            while (iThreads < (terminator - 1)) {
                ranges[iThreads] = new Range(intervals[iThreads], intervals[iThreads + 1] - 1);
                //      logger.info(String.format("Range for thread %d %s %d.", iThreads, ranges[iThreads], fftZ));
                iThreads++;
            }
            ranges[terminator - 1] = new Range(intervals[terminator - 1], intervals[terminator]);
            //   logger.info(String.format("Range for thread %d %s %d.", terminator - 1, ranges[terminator - 1], fftZ));

            for (int it = terminator; it < nThreads; it++) {
                ranges[it] = null;
            }

        } else {
            Range temp = new Range(0, fftZ - 1);
            ranges = temp.subranges(nThreads);
        }
    }

    public int[] getWeightPerThread() {
        if (intervals != null) {
            int[] weightToReturn = new int[nThreads];
            for (int i = 0; i < nThreads; i++) {
                weightToReturn[i] = intervals[i + 1];
            }
            return weightToReturn;
        } else {
            return null;
        }
    }
}
