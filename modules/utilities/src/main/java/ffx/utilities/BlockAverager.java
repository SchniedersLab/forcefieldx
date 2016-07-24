/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.utilities;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Computes the true uncertainty for metadynamics histogram data.
 * Takes a log file containing histogram over time; calculates the block-averaged 
 * uncertainty of each lambda bin. First-order error propagation is used to combine
 * bin uncertainties into a total uncertainty.
 * 
 * @author S. LuCore
 */
public class BlockAverager {
    private static final Logger logger = Logger.getLogger(BlockAverager.class.getName());
    private static int histoIndexer = 0;
    
    private final int numObs;
    private final int numBins;
    private final int maxBlockSize;
    private final int blockSizeStep;
    private final double psPerHisto;
    private final List<Histogram> histoList = new ArrayList<>();
    private double[] stdError;
    
    /** Parallel Stuff                **/
    private final ParallelTeam parallelTeam;
    private final int numThreads;
    
    /** Debugging and Testing Options **/
    private final MODE mode = (System.getProperty("ba-mode") == null) ? 
            MODE.dG : MODE.valueOf(System.getProperty("ba-mode"));
    private final String preGrep = System.getProperty("ba-preGrep");
    private final FITTER fitter = (System.getProperty("ba-fitter") == null) ? 
            FITTER.LOG : FITTER.valueOf(System.getProperty("ba-fitter"));
    private final int polyDegree = 2;
    private boolean TEST = (System.getProperty("ba-test") != null);
    
    private boolean blockByBin = (System.getProperty("ba-byBin") != null);
    
    public enum MODE {
        avgFL, dG, G;
    }
    public enum FITTER {
        POLYNOMIAL, POWER, LOG, LINEAR;
    }
    
    /**
     * Constructor grabs all histograms from the file and loads them into data structures.
     *  TODO: figure out how to disregard histogram-bin combos that aren't (currently) changing per time.
     */
    public BlockAverager(String filename,    boolean testMode, 
            Optional<String> grepCmd,        Optional<Double> psPerHisto, 
            Optional<Integer> blockSizeStep, Optional<Integer> maxBlockSize) throws IOException {
        this.TEST = testMode;
        this.psPerHisto = (psPerHisto.isPresent()) ? psPerHisto.get() : 1.0;
        this.blockSizeStep = (blockSizeStep.isPresent()) ? blockSizeStep.get() : 100;
        int linesPerHistogram = (System.getProperty("ba-lph") == null) ? 
            201 : Integer.parseInt(System.getProperty("ba-lph"));
        
        if (TEST) {
            logger.info(" Testing Mode ");
            linesPerHistogram = 1;
        }
        
        File parallelInFile = new File(filename);
        int nThreads = ParallelTeam.getDefaultThreadCount();
        parallelTeam = new ParallelTeam(nThreads);
        numThreads = parallelTeam.getThreadCount();
        BlockRegion parallelBlock = new BlockRegion(parallelInFile);
        try {
            parallelTeam.execute(parallelBlock);
        } catch (Exception ex) {
            Logger.getLogger(BlockAverager.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        // Step 1: Find histograms and create a stream.
        Scanner scan = null;
        File outFile = null;
        if (preGrep != null) {
            File file = new File(preGrep);
            BufferedReader br = new BufferedReader(new FileReader(file));
            scan = new Scanner(br);
        } else {
            outFile = new File(filename + "-ba.tmp");
            if (outFile.exists()) {
                logger.info(format(" Previous temp file exists: %s", outFile.getName()));
                if (!outFile.canWrite()) {
                    logger.severe(format("Lacked write permissions to temp file."));
                }
                System.out.print(format("   Delete it? (Y/N) "));
                Scanner kb = new Scanner(System.in);
                if (kb.nextLine().toUpperCase().startsWith("Y")) {
                    outFile.delete();
                    logger.info("");
                } else {
                    logger.severe("Aborted by user.");
                }
            }

            // Manually accomplish a 'grep -A 201 Bins filename'.
            File inFile = new File(filename);
            BufferedReader br = new BufferedReader(new FileReader(inFile));
            scan = new Scanner(br);
            BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
            
            logger.info(" Parsing logfile... ");
            int numFound = 0;
            while(scan.hasNextLine()) {
                String line = scan.nextLine();
                if (TEST) {     // No headers in test data.
                    if (++numFound % 100 == 0) {
                        logger.info(format("    Parsed %d histograms.", numFound));
                    }
                    bw.write(line);
                    bw.newLine();
                    continue;
                }
                if (line.contains("Lambda Bins")) {
                    if (++numFound % 100 == 0) {
                        logger.info(format("    Parsed %d histograms.", numFound));
                    }
                    bw.write(line);
                    bw.newLine();
                    for (int i = 0; i < linesPerHistogram; i++) {
                        if (!scan.hasNextLine() && i < linesPerHistogram) {
                            logger.warning(format("Found incomplete histogram: %d, %s", 
                                    numFound, line));
                        }
                        bw.write(scan.nextLine());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            scan = new Scanner(outFile);
        }
        
        // Parse stream into data structures.
        List<Bin> binList = new ArrayList<>();
        Histogram histo = null;
        while (scan.hasNextLine()) {
            String line = scan.nextLine();
            String[] tokens = line.split("\\s+");
            // Catch grep flotsam.
            if (tokens[0].startsWith("--")) {
                continue;
            }
            // Header line signals time for a new histogram.
            if (line.contains("Lambda Bins") || TEST) {
                if (histo != null) {
                    histoList.add(histo);
                }
                histo = new Histogram(++histoIndexer);
                if (histoIndexer % 100 == 0) {
                    if (psPerHisto.isPresent()) {
                        logger.info(format(" BlockAverager loaded %d histograms (simTime %.2f ps).", 
                                histoIndexer, histoIndexer * this.psPerHisto));
                    } else {
                        logger.info(format(" BlockAverager loaded %d histograms.", histoIndexer));
                    }
                }
                if (TEST) {     // No headers in test data.
                    histo.bins.add(new Bin(tokens));
                }
                continue;
            }
            histo.bins.add(new Bin(tokens));
        }
        histoList.add(histo);
        Collections.sort(histoList);
        logger.info(format(""));
        
        numObs = histoList.size();
        this.maxBlockSize = (maxBlockSize.isPresent()) ? maxBlockSize.get() : numObs;
        
        // Validate
        for (int i = 1; i < histoList.size(); i++) {
            if (histoList.get(i).index != histoList.get(i-1).index + 1
                || histoList.get(i).bins.size() != histoList.get(i-1).bins.size()) {
                logger.warning(format("Improper indexing or bin size mismatch. i,i-1,binsi,binsi-1: %d %d %d %d", 
                        histoList.get(i).index,         histoList.get(i-1).index, 
                        histoList.get(i).bins.size(),   histoList.get(i-1).bins.size()));
                throw new ArithmeticException();
            }
        }
        
        if (outFile != null && outFile.exists()) {
            outFile.delete();
        }
        numBins = histoList.get(0).bins.size();
        this.describe();
    }
    
    /**
     * Use first-order error propagation to combine bin uncertainties into a total std error.
     */
    public double computeTotalUncertainty() {
        logger.info(format(" Total Combined StdError of %s:", mode.toString()));
        double totalStdError;
        double sumSq = 0.0;
        for (int bin = 0; bin < numBins; bin++) {
            sumSq += stdError[bin] * stdError[bin];
        }
        totalStdError = Math.sqrt(sumSq);
        logger.info(format("    Log stdErr: %12.10g ", totalStdError));
        return totalStdError;
    }
    
    public final void describe() {
        StringBuilder sb = new StringBuilder();
        sb.append(format(" BlockAverager over %s: \n", mode.toString()));
        sb.append(format("    histograms:     %4d \n", numObs));
        sb.append(format("    blockSizeStep:  %4d \n", blockSizeStep));
        sb.append(format("    maxBlockSize:   %4d \n", maxBlockSize));
//        if (TEST) {
//            sb.append(format("\n HistoList: \n"));
//            for (int i = 0; i < histoList.size(); i++) {
//                Histogram histo = histoList.get(i);
//                sb.append(format("    %4d (%d)    %6.4f    %6.4f\n", 
//                        histo.index, histo.bins.size(), 
//                        histo.bins.get(0).count, histo.bins.get(0).avgFL));
//            }
//        }
        logger.info(sb.toString());
    }
    
    /**
     * Compute the statistical uncertainty of G in each histogram bin and overall.
     *  Loop over increasing values of block size.
     *  For each, calculate the block means and their standard deviation.
     *  Then limit(blockStdErr, blockSize->entireTraj) == trajStdErr.
     * 
     * @return aggregate standard error of the total free energy change
     */
    public double[] computeBinUncertainties() {
        double[][] sems = new double[numBins][maxBlockSize+1];
        BinDev[][] binStDevs = new BinDev[numBins][maxBlockSize+1];
        
        List<WeightedObservedPoint>[] obsDev = new ArrayList[numBins];
        List<WeightedObservedPoint>[] obsErr = new ArrayList[numBins];
        
        for (int binIndex = 0; binIndex < numBins; binIndex++) {
            logger.info(format(" Computing stdError for bin %d...", binIndex));
            obsDev[binIndex] = new ArrayList<>();
            obsErr[binIndex] = new ArrayList<>();
            for (int blockSize = 1; blockSize <= maxBlockSize; blockSize += blockSizeStep) {
                int numBlocks = (int) Math.floor(numObs / blockSize);
                binStDevs[binIndex][blockSize] = new BinDev(binIndex, blockSize);
                sems[binIndex][blockSize] = binStDevs[binIndex][blockSize].stdev / Math.sqrt(numBlocks - 1);
                obsDev[binIndex].add(new WeightedObservedPoint(1.0, blockSize, binStDevs[binIndex][blockSize].stdev));
                obsErr[binIndex].add(new WeightedObservedPoint(1.0, blockSize, sems[binIndex][blockSize]));
                if (TEST) {
                    logger.info(format("  bin,blockSize,stdev,sem: %d %d %.6g %.6g", 
                            binIndex, blockSize, 
                            binStDevs[binIndex][blockSize].stdev, 
                            sems[binIndex][blockSize]));
                }
            }
        }
        
        // Fit a function to (blockSize v. stdError) and extrapolate to blockSize == entire trajectory.
        // This is our correlation-corrected estimate of the std error for this lambda bin.
        stdError = new double[numBins];
        for (int binIndex = 0; binIndex < numBins; binIndex++) {
            logger.info(format("\n Bin %d : fitting & extrapolating blockSize v. stdError", binIndex));
            if (fitter == FITTER.POLYNOMIAL) {
                // Fit a polynomial (shitty).
                double[] coeffsPoly = PolynomialCurveFitter.create(polyDegree).fit(obsErr[binIndex]);
                PolynomialFunction poly = new PolynomialFunction(coeffsPoly);
                logger.info(format("    Poly %d:   %12.10g     %s", 
                        polyDegree, poly.value(numObs), Arrays.toString(poly.getCoefficients())));
            } else if (fitter == FITTER.POWER) {
                // Fit a power function (better).
                double[] coeffsPow = powerFit(obsErr[binIndex]);
                double powerExtrapolated = coeffsPow[0]*pow(numObs,coeffsPow[1]);
                logger.info(format("    Power:     %12.10g     %s", powerExtrapolated, Arrays.toString(coeffsPow)));
            }
            // Fit a log function (best).
            double[] logCoeffs = logFit(obsErr[binIndex]);
            double logExtrap = logCoeffs[0] + logCoeffs[1]*log(numObs);
            logger.info(format("    Log sem:   %12.10g     Residual: %12.10g     Coeffs: %6.4f, %6.4f \n", 
                    logExtrap, logCoeffs[0], logCoeffs[1]));
            
            // Also try fitting a linear function for the case of uncorrelated or extremely well-converged data.
            double [] linearCoef = linearFit(obsErr[binIndex]);
            double linearExtrap = linearCoef[0] + linearCoef[1]*numObs;
            logger.info(format("    Lin. sem:  %12.10g     Residual: %12.10g     Coeffs: %6.4f, %6.4f \n", 
                    linearExtrap, linearCoef[0], linearCoef[1]));
            
            stdError[binIndex] = logExtrap;
        }
        return stdError;
    }
    
    /**
     * Returns [A,B] such that f(x) = A*(x^B) is minimized for the input points.
     * As described at http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
     */
    private double[] powerFit(List<WeightedObservedPoint> obs) {
        int n = obs.size();
        double sumlnxlny = 0.0;
        double sumlnx = 0.0;
        double sumlny = 0.0;
        double sumlnxsq = 0.0;
        for (int i = 0; i < n; i++) {
            final double x = obs.get(i).getX() * obs.get(i).getWeight();
            final double y = obs.get(i).getY() * obs.get(i).getWeight();
            final double lnx = log(x);
            final double lny = log(y);
            sumlnxlny += lnx*lny;
            sumlnx += lnx;
            sumlny += lny;
            sumlnxsq += lnx*lnx;
        }
        final double b = (n*sumlnxlny - sumlnx*sumlny) / (n*sumlnxsq - sumlnx*sumlnx);
        final double a = (sumlny - b*sumlnx) / n;
        final double B = b;
        final double A = Math.exp(a);
        double[] ret = {A,B};
        return ret;
    }
    
    /**
     * Returns [A,B] such that f(x)= a + b*ln(x) is minimized for the input points.
     * As described at http://mathworld.wolfram.com/LeastSquaresFittingLogarithmic.html
     */
    private double[] logFit(List<WeightedObservedPoint> obs) {
        int n = obs.size();
        double sumylnx = 0.0;
        double sumy = 0.0;
        double sumlnx = 0.0;
        double sumlnxsq = 0.0;
        for (int i = 0; i < n; i++) {
            final double x = obs.get(i).getX() * obs.get(i).getWeight();
            final double y = obs.get(i).getY() * obs.get(i).getWeight();
            final double lnx = log(x);
            final double lny = log(y);
            sumylnx += y*lnx;
            sumy += y;
            sumlnx += lnx;
            sumlnxsq += lnx*lnx;
        }
        final double b = (n*sumylnx - sumy*sumlnx) / (n*sumlnxsq - sumlnx*sumlnx);
        final double a = (sumy - b*sumlnx) / n;
        double[] ret = {a,b};
        return ret;
    }
    
    /**
     * Returns [A,B,R^2] such that f(x)= a + b*x is minimized for the input points.
     * As described at http://mathworld.wolfram.com/LeastSquaresFitting.html
     */
    private double[] linearFit(List<WeightedObservedPoint> obs) {
        int n = obs.size();
        double sumx = 0.0;
        double sumy = 0.0;
        double sumxsq = 0.0;
        double sumysq = 0.0;
        double sumxy = 0.0;        
        for (int i = 0; i < n; i++) {
            final double x = obs.get(i).getX() * obs.get(i).getWeight();
            final double y = obs.get(i).getY() * obs.get(i).getWeight();
            sumx += x;
            sumy += y;
            sumxsq += x*x;
            sumysq += y*y;
            sumxy += x*y;
        }
        final double xbar = sumx / n;
        final double ybar = sumy / n;
        
        final double ssxx = sumxsq - n*xbar*xbar;
        final double ssyy = sumysq - n*ybar*ybar;
        final double ssxy = sumxy - n*xbar*ybar;
        final double rsq = (ssxy*ssxy) / (ssxx*ssyy);
        
        final double a = (ybar*sumxsq - xbar*sumxy) / (sumxsq - n*xbar*xbar);
        final double b = (sumxy - n*xbar*ybar) / (sumxsq - n*xbar*xbar);
        double[] ret = {a,b,rsq};
        return ret;
    }
    
    /**
     * Computes a residual to the given points for the provided fit type and coefficients.
     */
    private double residual(List<WeightedObservedPoint> obs, FITTER fitter, double[] coeffs) {
        int n = obs.size();
        double sumydiffsq = 0.0;
        for (int i = 0; i < n; i++) {
            final double x = obs.get(i).getX() * obs.get(i).getWeight();
            final double y = obs.get(i).getY() * obs.get(i).getWeight();
            double value;
            switch (fitter) {
                case LINEAR:
                    value = coeffs[0] + coeffs[1]*x;
                    break;
                case LOG:
                    value = coeffs[0] + coeffs[1]*log(x);
                    break;
                case POWER:
                    value = coeffs[0] * pow(x, coeffs[1]);
                    break;
                default:
                    throw new UnsupportedOperationException();
            }
            sumydiffsq += (y - value)*(y - value);
        }
        double residual = sqrt(sumydiffsq);
        return residual;
    }
    
    /**
     * Computes stdev of one lambda bin at one block size.
     */
    private class BinDev {
        public final int binIndex;
        public final int blockSize;
        public final double[] mean;     // Mean of each block.
        public final double stdev;      // Stdev of the BLOCK MEANS.
                
        private BinDev(int binIndex, int blockSize) {
            this.binIndex = binIndex;
            this.blockSize = blockSize;
            
            int numBlocks;
            double meanSum = 0.0;
            if (!blockByBin) {
                // This blocks out every bin using histoList.size(), i.e. total simulation time.
                // Assumes an evenly-sampled histogram; error from undersampled bins is underestimated.
                numBlocks = (int) Math.floor(histoList.size() / blockSize);
                mean = new double[numBlocks];

                // Compute the mean of "avgFL" in each block; find the average mean.
                for (int block = 0; block < numBlocks; block++) {
                    double sum = 0.0;
                    for (int histo = block*blockSize; histo < block*blockSize + blockSize; histo++) {
                        switch (mode) {
                            case avgFL:
                                sum += histoList.get(histo).bins.get(binIndex).avgFL;
                                break;
                            case dG:
                                sum += histoList.get(histo).bins.get(binIndex).dG;
                                break;
                            case G:
                                sum += histoList.get(histo).bins.get(binIndex).G;
                                break;
                        }
                    }
                    mean[block] = sum / blockSize;
                    meanSum += mean[block];
                    logger.info(format("    Block mean: %8.4f          Mean sum: %8.4f", mean[block], meanSum));
                }
            } else {
                // Give EACH BIN its own blocking.
                logger.info(format(" Blocking for bin %d...", binIndex));
                int totalBinCounts = (int) floor(histoList.get(histoList.size() - 1).bins.get(binIndex).count);
                numBlocks = (int) Math.floor(totalBinCounts / blockSize);
                mean = new double[numBlocks];
                logger.info(format("    totalBinCounts,numBlocks,countsPerBlock: %d, %d, %d", 
                        totalBinCounts, numBlocks, blockSize));
                
                // Compute the mean of requested property in each block; find the average mean.
                for (int block = 0; block < numBlocks; block++) {
                    
                    // Find which histograms contribute to this bin's block.
                    int blockCountsLow = blockSize*block;
                    int blockCountsHigh = blockSize*block + blockSize;
                    logger.info(format("    Summing for block {%d , %d}: ", blockCountsLow, blockCountsHigh));
                    
                    double sum = 0.0;
                    int previousCount = -1;
                    double previousValue = 0.0;
                    for (int histo = 0; histo < histoList.size(); histo++) {
                        int count = (int) floor(histoList.get(histo).bins.get(binIndex).count);
                        if (count > blockCountsHigh) {
                            break;
                        }
                        if (count >= blockCountsLow) {
                            // So now this value is part of the block bin, but this is the tricky part...
                            // There may not have been NEW counts in this bin since our previous histogram.
                            if (previousCount == -1) {
                                previousCount = count;
                                continue;
                            }
                            if (count > previousCount) {
                                // There may also have been SEVERAL counts in this bin since our previous histogram.
                                // We'll have to average the property in question over the change in counts.
                                int deltaCount = count - previousCount;
                                double value = 0.0;
                                switch (mode) {
                                    case avgFL:
                                        value = histoList.get(histo).bins.get(binIndex).avgFL;
                                        break;
                                    case dG:
                                        value = histoList.get(histo).bins.get(binIndex).dG;
                                        break;
                                    case G:
                                        value = histoList.get(histo).bins.get(binIndex).G;
                                        break;
                                }
                                sum += value * deltaCount;
                                previousCount = count;
//                                logger.info(format("       Count changed! histo,count,deltaCount,deltaValue,addToSum: %d, %8.4f, %8.4f", 
//                                        histo, count, deltaCount, deltaValue, avgValue * deltaCount));
                            }
                        }
                    }
                    // Record the mean and add to the average mean.
                    mean[block] = sum / blockSize;
                    meanSum += mean[block];
                    logger.info(format("    Block mean: %8.4f          Mean sum: %8.4f", mean[block], meanSum));
                }
            }
            double meanOfMeans = meanSum / numBlocks;

            double sumSqDiff = 0.0;       // sum of the squared difference of each mean to the mean of means
            for (int block = 0; block < numBlocks; block++) {
                sumSqDiff += Math.pow(mean[block] - meanOfMeans, 2);
            }
            stdev = Math.sqrt(sumSqDiff / numBlocks);
            logger.info(format(" StDev of block means for bin %d: %8.4f", binIndex, stdev));
        }
    }
    
    /**
     * Uncorrelated process: x = 5 + 2 * rand(N,1) - 1
     * Correlated process:   x(1) = 0; x(t+1) = 0.95 * x(t) + 2 * rand(N,1) - 1; shift all up by 5
     */
    public static void generateTestData(String filename, int size) throws IOException {
        logger.info(format(" Generating test data set of size: %d.", size));
        File outFile = new File(filename);
        if (outFile.exists()) {
            logger.severe("Outfile already exists.");
        }
        Random rng = new Random();
        double[] uncor = new double[size];
        double[] corr = new double[size];   corr[0] = 0.0;
        for (int i = 0; i < size; i++) {
            uncor[i] = 5 + (2*rng.nextDouble() - 1);
            if (i > 0) {
                corr[i] = 0.95*corr[i-1] + (2*rng.nextDouble() - 1);
            }
        }
        for (int i = 0; i < size; i++) {
            corr[i] += 5.0;
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
        for (int i = 0; i < size; i++) {
            bw.write(format(" %4d  %6.4f  %6.4f", i, uncor[i], corr[i]));
            bw.newLine();
        }
        bw.close();
        logger.info(format("    Data saved to: %s", filename));
    }
    
    private class Histogram implements Comparable {
        public final int index;                             // int
        public final List<Bin> bins = new ArrayList<>();    // entries
        
        private Histogram(int index) {
            this.index = histoIndexer;
        }
        
        @Override
        public int compareTo(Object other) {
            if (!(other instanceof Histogram)) {
                throw new UnsupportedOperationException();
            }
            return Integer.compare(index, ((Histogram) other).index);
        }
    }
    
    private class Bin implements Comparable {
        public final double count;
        public final double binStart, binEnd;
        public final double FLbinStart, FLbinEnd;
        public final double avgFL;
        public final double dG, G;
        
        private String[] shift(String[] tokens) {
            String[] newTok = new String[tokens.length - 1];
            for (int i = 1; i < tokens.length; i++) {
                newTok[i-1] = tokens[i];
            }
            return newTok;
        }
        
        private Bin(String[] tokens) {
            // Remove empty starting tokens and process identifiers.
            if (tokens[0].equals("")) tokens = shift(tokens);
            if (tokens[0].startsWith("[")) tokens = shift(tokens);
            if (tokens[0].equals("")) tokens = shift(tokens);
            
            if (TEST) {
                count = (tokens.length > 1) ? Integer.parseInt(tokens[0]) : 0;
                avgFL = (tokens.length > 1) ? Double.parseDouble(tokens[1]) : Double.parseDouble(tokens[0]);
                dG = (tokens.length > 2) ? Double.parseDouble(tokens[2]) : 0.0;
                G = (tokens.length > 3) ? Double.parseDouble(tokens[3]) : 0.0;
                binStart = 0.0; binEnd = 0.0; FLbinStart = 0.0; FLbinEnd = 0.0;
                return;
            }
            
            if (tokens.length != 8) {
                logger.warning(format("Incorrect number of tokens on histogram line: %s", Arrays.toString(tokens)));
            }
            
            try {
                count = (tokens[0].contains(".")) ? Double.parseDouble(tokens[0]) : Integer.parseInt(tokens[0]);
                binStart = Double.parseDouble(tokens[1]);    // ^^ number of walker visits to this bin
                binEnd = Double.parseDouble(tokens[2]);      // defines range of the lambda bin
                FLbinStart = Double.parseDouble(tokens[3]);  // defines range of the lambda force bin
                FLbinEnd = Double.parseDouble(tokens[4]);
                avgFL = Double.parseDouble(tokens[5]);       // average force along lambda from this bin
                dG = Double.parseDouble(tokens[6]);          // free energy from this bin
                G = Double.parseDouble(tokens[7]);           // cumulative free energy sum
            } catch (NumberFormatException ex) {
                logger.warning(format("Bin creation failed for tokens: %s", 
                        Arrays.toString(tokens)));
                throw ex;
            }
        }
        
        @Override
        public int compareTo(Object o) {
            if (!(o instanceof Bin)) {
                throw new UnsupportedOperationException();
            }
            Bin ob = (Bin) o;
            if (this.binStart != ob.binStart) {
                return Double.compare(this.binStart, ob.binStart);
            } else {
                return Double.compare(this.count, ob.count);
            }
        }
        
        @Override
        public boolean equals(Object o) {
            if (o == null || !(o instanceof Bin)) {
                return false;
            }
            Bin ob = (Bin) o;
            if (this.binStart == ob.binStart && this.count == ob.count) {
                if (this.dG != ob.dG) {
                    logger.severe(format("Inconsistent bin information: binA %s, binB %s", this, ob));
                }
                return true;
            }
            return false;
        }
        
        @Override
        public String toString() {
            return format(" %5.3f %5.3f %5.3f %10.3f %10.3f", 
                    count, binStart, binEnd, avgFL, dG);
        }
    }
    
    private final static HashMap<Integer,Double> binLookup = new HashMap<>();
    static {
        binLookup.put(0,0.000);
        binLookup.put(1,0.003);
        binLookup.put(2,0.008);
        binLookup.put(3,0.012);
        binLookup.put(4,0.018);
        binLookup.put(5,0.023);
        binLookup.put(6,0.028);
        binLookup.put(7,0.033);
        binLookup.put(8,0.038);
        binLookup.put(9,0.042);
        binLookup.put(10,0.048);
        binLookup.put(11,0.053);
        binLookup.put(12,0.057);
        binLookup.put(13,0.063);
        binLookup.put(14,0.068);
        binLookup.put(15,0.073);
        binLookup.put(16,0.078);
        binLookup.put(17,0.083);
        binLookup.put(18,0.088);
        binLookup.put(19,0.093);
        binLookup.put(20,0.098);
        binLookup.put(21,0.103);
        binLookup.put(22,0.108);
        binLookup.put(23,0.113);
        binLookup.put(24,0.118);
        binLookup.put(25,0.123);
        binLookup.put(26,0.128);
        binLookup.put(27,0.133);
        binLookup.put(28,0.138);
        binLookup.put(29,0.143);
        binLookup.put(30,0.148);
        binLookup.put(31,0.153);
        binLookup.put(32,0.158);
        binLookup.put(33,0.163);
        binLookup.put(34,0.168);
        binLookup.put(35,0.173);
        binLookup.put(36,0.178);
        binLookup.put(37,0.183);
        binLookup.put(38,0.188);
        binLookup.put(39,0.193);
        binLookup.put(40,0.198);
        binLookup.put(41,0.203);
        binLookup.put(42,0.208);
        binLookup.put(43,0.213);
        binLookup.put(44,0.218);
        binLookup.put(45,0.223);
        binLookup.put(46,0.228);
        binLookup.put(47,0.233);
        binLookup.put(48,0.238);
        binLookup.put(49,0.243);
        binLookup.put(50,0.248);
        binLookup.put(51,0.253);
        binLookup.put(52,0.258);
        binLookup.put(53,0.263);
        binLookup.put(54,0.268);
        binLookup.put(55,0.273);
        binLookup.put(56,0.278);
        binLookup.put(57,0.283);
        binLookup.put(58,0.288);
        binLookup.put(59,0.293);
        binLookup.put(60,0.298);
        binLookup.put(61,0.303);
        binLookup.put(62,0.308);
        binLookup.put(63,0.313);
        binLookup.put(64,0.318);
        binLookup.put(65,0.323);
        binLookup.put(66,0.328);
        binLookup.put(67,0.333);
        binLookup.put(68,0.338);
        binLookup.put(69,0.343);
        binLookup.put(70,0.348);
        binLookup.put(71,0.353);
        binLookup.put(72,0.358);
        binLookup.put(73,0.363);
        binLookup.put(74,0.368);
        binLookup.put(75,0.373);
        binLookup.put(76,0.378);
        binLookup.put(77,0.383);
        binLookup.put(78,0.388);
        binLookup.put(79,0.393);
        binLookup.put(80,0.398);
        binLookup.put(81,0.403);
        binLookup.put(82,0.408);
        binLookup.put(83,0.413);
        binLookup.put(84,0.418);
        binLookup.put(85,0.423);
        binLookup.put(86,0.428);
        binLookup.put(87,0.433);
        binLookup.put(88,0.438);
        binLookup.put(89,0.443);
        binLookup.put(90,0.448);
        binLookup.put(91,0.453);
        binLookup.put(92,0.458);
        binLookup.put(93,0.463);
        binLookup.put(94,0.468);
        binLookup.put(95,0.473);
        binLookup.put(96,0.478);
        binLookup.put(97,0.483);
        binLookup.put(98,0.488);
        binLookup.put(99,0.493);
        binLookup.put(100,0.498);
        binLookup.put(101,0.503);
        binLookup.put(102,0.508);
        binLookup.put(103,0.513);
        binLookup.put(104,0.518);
        binLookup.put(105,0.523);
        binLookup.put(106,0.528);
        binLookup.put(107,0.533);
        binLookup.put(108,0.538);
        binLookup.put(109,0.543);
        binLookup.put(110,0.548);
        binLookup.put(111,0.553);
        binLookup.put(112,0.558);
        binLookup.put(113,0.563);
        binLookup.put(114,0.568);
        binLookup.put(115,0.573);
        binLookup.put(116,0.578);
        binLookup.put(117,0.583);
        binLookup.put(118,0.588);
        binLookup.put(119,0.593);
        binLookup.put(120,0.598);
        binLookup.put(121,0.603);
        binLookup.put(122,0.608);
        binLookup.put(123,0.613);
        binLookup.put(124,0.618);
        binLookup.put(125,0.623);
        binLookup.put(126,0.628);
        binLookup.put(127,0.633);
        binLookup.put(128,0.638);
        binLookup.put(129,0.643);
        binLookup.put(130,0.648);
        binLookup.put(131,0.653);
        binLookup.put(132,0.658);
        binLookup.put(133,0.663);
        binLookup.put(134,0.668);
        binLookup.put(135,0.673);
        binLookup.put(136,0.678);
        binLookup.put(137,0.683);
        binLookup.put(138,0.688);
        binLookup.put(139,0.693);
        binLookup.put(140,0.698);
        binLookup.put(141,0.703);
        binLookup.put(142,0.708);
        binLookup.put(143,0.713);
        binLookup.put(144,0.718);
        binLookup.put(145,0.723);
        binLookup.put(146,0.728);
        binLookup.put(147,0.733);
        binLookup.put(148,0.738);
        binLookup.put(149,0.743);
        binLookup.put(150,0.748);
        binLookup.put(151,0.753);
        binLookup.put(152,0.758);
        binLookup.put(153,0.763);
        binLookup.put(154,0.768);
        binLookup.put(155,0.773);
        binLookup.put(156,0.778);
        binLookup.put(157,0.783);
        binLookup.put(158,0.788);
        binLookup.put(159,0.793);
        binLookup.put(160,0.798);
        binLookup.put(161,0.803);
        binLookup.put(162,0.808);
        binLookup.put(163,0.813);
        binLookup.put(164,0.818);
        binLookup.put(165,0.823);
        binLookup.put(166,0.828);
        binLookup.put(167,0.833);
        binLookup.put(168,0.838);
        binLookup.put(169,0.843);
        binLookup.put(170,0.848);
        binLookup.put(171,0.853);
        binLookup.put(172,0.858);
        binLookup.put(173,0.863);
        binLookup.put(174,0.868);
        binLookup.put(175,0.873);
        binLookup.put(176,0.878);
        binLookup.put(177,0.883);
        binLookup.put(178,0.888);
        binLookup.put(179,0.893);
        binLookup.put(180,0.898);
        binLookup.put(181,0.903);
        binLookup.put(182,0.908);
        binLookup.put(183,0.913);
        binLookup.put(184,0.918);
        binLookup.put(185,0.923);
        binLookup.put(186,0.928);
        binLookup.put(187,0.933);
        binLookup.put(188,0.938);
        binLookup.put(189,0.943);
        binLookup.put(190,0.948);
        binLookup.put(191,0.953);
        binLookup.put(192,0.958);
        binLookup.put(193,0.963);
        binLookup.put(194,0.968);
        binLookup.put(195,0.973);
        binLookup.put(196,0.978);
        binLookup.put(197,0.983);
        binLookup.put(198,0.988);
        binLookup.put(199,0.993);
        binLookup.put(200,0.998);
    }
        
    private class BlockRegion extends ParallelRegion {

        private final File file;
        private final List<Bin>[] binLists;
        private final ParsingLoop parsingLoop;
        private final UncertaintyLoop binningLoop;
        private final int maxEntriesPerBin = (System.getProperty("ba-maxEntries") == null) ? 
                Integer.MAX_VALUE : Integer.parseInt(System.getProperty("ba-maxEntries"));

        public BlockRegion(File file) {
            // Make loops.
            this.file = file;
            parsingLoop = new ParsingLoop();
            binningLoop = new UncertaintyLoop();
            binLists = new ArrayList[binLookup.size()];
        }
        
        private void logIfMaster(String msg) {
            if (getThreadIndex() == 0) {
                logger.info(msg);
            }
        }

        @Override
        public void start() {
            // Single-threaded; initialize shared variables.
        }

        @Override
        public void finish() {
            // Single-threaded; cleanup.
        }

        @Override
        public void run() throws Exception {
            int threadIndex = getThreadIndex();
            logIfMaster(" Calling parse()...");
            execute(0, binLookup.size() - 1, parsingLoop);
            logIfMaster(" Thread zero finished parsing. Setting up barrier.");
            barrier();
            logIfMaster(" Barrier passed!");
            logIfMaster(" Launch binning loop...");
            execute(0, binLookup.size() - 1, binningLoop);
            logIfMaster(" Thread zero finsihed binning.");
        }
        
        /**
         * Search through a giant log file and create Bin objects.
         * Specialized for by-bin deviation calculation; when time happens on a 
         *  per-bin basis, the organization of the log file is of no consequence.
         */
        private class ParsingLoop extends IntegerForLoop {
            
            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }
            
            /**
             * lb,ub == binIndex. i.e. [0,201]
             */
            @Override
            public void run(int lb, int ub) {
                int thread = getThreadIndex();
                // Loop over bin indices assigned to this thread.
                for (int i = lb; i <= ub; i++) {
                    // Loop over all histograms and build BinDev objects for this binIndex.
                    binLists[i] = new ArrayList<>();
                    double targetBinStart = binLookup.get(i);
                    String target = format("%5.3f",targetBinStart);
                    String line = "";
                    int found = 0;
                    Bin previousBin = null;
                    try {
                        // Start reading the file again from the beginning.
                        BufferedReader br = new BufferedReader(new FileReader(file));
                        while ((line = br.readLine()) != null) {
                            String tokens[] = line.split("\\s+");
                            if (tokens != null && tokens.length >= 8) {
                                // Remove empty starting tokens and process identifiers.
                                if (tokens[0].equals("")) tokens = shift(tokens);
                                if (tokens[0].startsWith("[")) tokens = shift(tokens);
                                if (tokens[0].equals("")) tokens = shift(tokens);
                                if (tokens.length == 8 && tokens[1].equals(target)) {
                                    // We've found a histogram entry of our target bin.
                                    Bin bin = new Bin(tokens);
                                    // Only record entries that add NEW information about THIS bin.
                                    if (bin.equals(previousBin)) {
                                        binLists[i].add(bin);
                                        previousBin = bin;
                                    }
//                                    logger.info(format(" thread,i,count,dg: %d %d %d %.4f", 
//                                            getThreadIndex(), i, (int) bin.count, bin.dG));
                                    if (++found % 1000 == 0) {
                                        logger.info(format(" Thread %2d, bin %3d (%s) parsing %6d...", 
                                                thread, i, target, found));
                                        if (found >= maxEntriesPerBin) {
                                            logger.info(format(" Maximum bin entries reached by thread %d, bin %d (%s).", 
                                                    thread, i, target));
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        logger.info(format(" (Total) Thread %2d found %6d entries for bin %3d (%s).", 
                                thread, found, i, target));
                        Collections.sort(binLists[i]);
                        File outFile = new File(file.getName() + format(".%d.tmp", i));
                        BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
                        for (int bin = 0; bin < binLists[i].size(); bin++) {
                            bw.write(binLists[i].get(bin).toString());
                            bw.newLine();
                        }
                        bw.close();
                        logger.info(format(" Wrote evolution of bin %d to file: %s", i, outFile.getName()));
                    } catch (IOException ex) {
                        logger.severe(format(" IOException in ParsingLoop thread %d, line %s, exception %s", 
                                thread, line, ex.getMessage()));
                    }
                }
            }
            
            private String[] shift(String[] tokens) {
                String[] newTok = new String[tokens.length - 1];
                for (int i = 1; i < tokens.length; i++) {
                    newTok[i-1] = tokens[i];
                }
                return newTok;
            }
        }
        
        /**
         * Calculate bin deviations. Parallelized on a per-bin basis.
         */
        private class UncertaintyLoop extends IntegerForLoop {
            
            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }
            
            /**
             * lb,ub == binIndex. i.e. [0,201]
             */
            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    // Loop over all histograms and build BinDev objects for this binIndex.
                    
                }
            }
        }
    }
    
}
