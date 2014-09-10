/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.xray;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;
import java.util.Arrays;
import java.util.logging.Logger;

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
    private final int fftZ;
    private int weights[];

    protected SliceSchedule(int nThreads, int fftZ) {
        this.nThreads = nThreads;
        threadDone = new boolean[nThreads];
        ranges = new Range[nThreads];
        intervals = new int[nThreads + 1];
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
        Arrays.fill(threadDone, false);

        if (nThreads != ranges.length) {
            ranges = new Range[nThreads];
        }

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
            double targetWeight = totalWeight / nThreads;
            int j = 0;
            boolean risk;
            int diff;
            int prevDiff;
            int prevWeight;
            double prop;
            double tolerance = (fftZ / nThreads) / (fftZ * fftZ);

            intervals[0] = 0;
            for (int i = 0; i < nThreads - 1; i++) {
                int threadWeight = 0;
                prop = 0;
                diff = 0;
                prevDiff = 0;
                prevWeight = 0;
                risk = false;

                while (threadWeight < targetWeight && j < fftZ - 1 && !risk) {
                    threadWeight += weights[j];

                    prevDiff = diff;
                    diff = weights[j] - prevWeight;
                    prevWeight = weights[j];
                    prop = (targetWeight - threadWeight) / targetWeight;
                    if (prop > (1.0 - tolerance)) {
                        if (diff * prevDiff > 0) {
                            risk = true;
                        }
                    }
                    j++;
                }
                intervals[i + 1] = j;
            }
            /**
             * Check if final slices remain to be assigned.
             */
            if (j < fftZ - 1) {
                intervals[nThreads] = fftZ - 1;
            }

            for (int i = 0; i < nThreads - 1; i++) {
                ranges[i] = new Range(intervals[i], intervals[i + 1] - 1);
                //     logger.info(String.format("Range for thread %d %s %d.", i, ranges[i], fftZ));
            }
            ranges[nThreads - 1] = new Range(intervals[nThreads - 1], intervals[nThreads]);
            //   logger.info(String.format("Range for thread %d %s %d.", nThreads-1, ranges[nThreads - 1], fftZ));
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
