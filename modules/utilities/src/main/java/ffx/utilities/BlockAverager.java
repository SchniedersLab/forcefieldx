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

import java.util.ArrayList;
import java.util.Arrays;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Collections;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import java.util.Optional;
import java.util.Random;
import java.util.Scanner;
import static java.lang.String.format;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.pow;

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
    
    /** Debugging and Testing Options **/
    private final MODE mode = (System.getProperty("ba-mode") == null) ? 
            MODE.dG : MODE.valueOf(System.getProperty("ba-mode"));
    private final String preGrep = System.getProperty("ba-preGrep");
    private final FITTER fitter = (System.getProperty("ba-fitter") == null) ? 
            FITTER.LOG : FITTER.valueOf(System.getProperty("ba-fitter"));
    private final int polyDegree = 2;
    private boolean TEST = (System.getProperty("ba-test") != null);
    
    public enum MODE {
        avgFL, dG, G;
    }
    public enum FITTER {
        POLYNOMIAL, POWER, LOG;
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
            // Should only be actual histogram entries left at this point.
            if (tokens.length != 8) {
                logger.warning(format("Incorrect number of tokens on histogram line: %s", Arrays.toString(tokens)));
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
            logger.info(format("    Log sem:   %12.10g     Coeffs: %6.4f, %6.4f \n", 
                    logExtrap, logCoeffs[0], logCoeffs[1]));
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
            int numBlocks = (int) Math.floor(histoList.size() / blockSize);
            mean = new double[numBlocks];
            
            // Compute the mean of "avgFL" in each block; find the average mean.
            double meanSum = 0.0;
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
//                // This stdev within one block is unnecessary. We need the stdev OF THE MEANS.
//                double sumSqDiffFromMean = 0.0;
//                for (int histo = block*blockSize; histo < block*blockSize + blockSize; histo++) {
//                    sumSqDiffFromMean += Math.pow(histoList.get(histo).bins.get(binIndex).avgFL - mean[block], 2);
//                }
//                stdev[block] = Math.sqrt(sumSqDiffFromMean / (blockSize - 1));
            }
            double meanOfMeans = meanSum / numBlocks;
            
            double sumSqDiff = 0.0;       // sum of the squared difference of each mean to the mean of means
            for (int block = 0; block < numBlocks; block++) {
                sumSqDiff += Math.pow(mean[block] - meanOfMeans, 2);
            }
            stdev = Math.sqrt(sumSqDiff / numBlocks);
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
            return Double.compare(this.binStart, ((Bin) o).binStart);
        }
    }
    
    /*
 Count   Lambda Bins    F_Lambda Bins   <   F_L  >       dG        G
      0  0.000 0.003      -5.0     5.0      0.002     0.000    0.000
      0  0.003 0.008       0.0     0.0      0.000     0.000    0.000
      0  0.008 0.013       0.0     0.0      0.000     0.000    0.000
      0  0.012 0.017       0.0     0.0      0.000     0.000    0.000
      0  0.018 0.023       0.0     0.0      0.000     0.000    0.000
      0  0.023 0.028       0.0     0.0      0.000     0.000    0.000
      0  0.028 0.033       0.0     0.0      0.000     0.000    0.000
      0  0.033 0.038       0.0     0.0      0.000     0.000    0.000
      0  0.038 0.042       0.0     0.0      0.000     0.000    0.000
      0  0.042 0.047       0.0     0.0      0.000     0.000    0.000
      0  0.048 0.053       0.0     0.0      0.000     0.000    0.000
      0  0.053 0.057       0.0     0.0      0.000     0.000    0.000
      0  0.057 0.062       0.0     0.0      0.000     0.000    0.000
      0  0.063 0.068       0.0     0.0      0.000     0.000    0.000
      0  0.068 0.073       0.0     0.0      0.000     0.000    0.000
      0  0.073 0.078       0.0     0.0      0.000     0.000    0.000
      0  0.078 0.083       0.0     0.0      0.000     0.000    0.000
      0  0.083 0.088       0.0     0.0      0.000     0.000    0.000
      0  0.088 0.093       0.0     0.0      0.000     0.000    0.000
      0  0.093 0.098       0.0     0.0      0.000     0.000    0.000
      0  0.098 0.103       0.0     0.0      0.000     0.000    0.000
      0  0.103 0.108       0.0     0.0      0.000     0.000    0.000
      0  0.108 0.113       0.0     0.0      0.000     0.000    0.000
      0  0.113 0.118       0.0     0.0      0.000     0.000    0.000
      0  0.118 0.123       0.0     0.0      0.000     0.000    0.000
      0  0.123 0.128       0.0     0.0      0.000     0.000    0.000
      0  0.128 0.133       0.0     0.0      0.000     0.000    0.000
      0  0.133 0.138       0.0     0.0      0.000     0.000    0.000
      0  0.138 0.143       0.0     0.0      0.000     0.000    0.000
      0  0.143 0.148       0.0     0.0      0.000     0.000    0.000
      0  0.148 0.153       0.0     0.0      0.000     0.000    0.000
      0  0.153 0.158       0.0     0.0      0.000     0.000    0.000
      0  0.158 0.163       0.0     0.0      0.000     0.000    0.000
      0  0.163 0.168       0.0     0.0      0.000     0.000    0.000
      0  0.168 0.173       0.0     0.0      0.000     0.000    0.000
      0  0.173 0.178       0.0     0.0      0.000     0.000    0.000
      0  0.178 0.183       0.0     0.0      0.000     0.000    0.000
      0  0.183 0.188       0.0     0.0      0.000     0.000    0.000
      0  0.188 0.193       0.0     0.0      0.000     0.000    0.000
      0  0.193 0.198       0.0     0.0      0.000     0.000    0.000
      0  0.198 0.203       0.0     0.0      0.000     0.000    0.000
      0  0.203 0.208       0.0     0.0      0.000     0.000    0.000
      0  0.208 0.213       0.0     0.0      0.000     0.000    0.000
      0  0.213 0.218       0.0     0.0      0.000     0.000    0.000
      0  0.218 0.223       0.0     0.0      0.000     0.000    0.000
      0  0.223 0.228       0.0     0.0      0.000     0.000    0.000
      0  0.228 0.233       0.0     0.0      0.000     0.000    0.000
      0  0.233 0.238       0.0     0.0      0.000     0.000    0.000
      0  0.238 0.243       0.0     0.0      0.000     0.000    0.000
      0  0.243 0.248       0.0     0.0      0.000     0.000    0.000
      0  0.248 0.253       0.0     0.0      0.000     0.000    0.000
      0  0.253 0.258       0.0     0.0      0.000     0.000    0.000
      0  0.258 0.263       0.0     0.0      0.000     0.000    0.000
      0  0.263 0.268       0.0     0.0      0.000     0.000    0.000
      0  0.268 0.273       0.0     0.0      0.000     0.000    0.000
      0  0.273 0.278       0.0     0.0      0.000     0.000    0.000
      0  0.278 0.283       0.0     0.0      0.000     0.000    0.000
      0  0.283 0.288       0.0     0.0      0.000     0.000    0.000
      0  0.288 0.293       0.0     0.0      0.000     0.000    0.000
      0  0.293 0.298       0.0     0.0      0.000     0.000    0.000
      0  0.298 0.303       0.0     0.0      0.000     0.000    0.000
      0  0.303 0.308       0.0     0.0      0.000     0.000    0.000
      0  0.308 0.313       0.0     0.0      0.000     0.000    0.000
      0  0.313 0.318       0.0     0.0      0.000     0.000    0.000
      0  0.318 0.323       0.0     0.0      0.000     0.000    0.000
      0  0.323 0.328       0.0     0.0      0.000     0.000    0.000
      0  0.328 0.333       0.0     0.0      0.000     0.000    0.000
      0  0.333 0.338       0.0     0.0      0.000     0.000    0.000
      0  0.338 0.343       0.0     0.0      0.000     0.000    0.000
      0  0.343 0.348       0.0     0.0      0.000     0.000    0.000
      0  0.348 0.353       0.0     0.0      0.000     0.000    0.000
      0  0.353 0.358       0.0     0.0      0.000     0.000    0.000
      0  0.358 0.363       0.0     0.0      0.000     0.000    0.000
      0  0.363 0.368       0.0     0.0      0.000     0.000    0.000
      0  0.368 0.373       0.0     0.0      0.000     0.000    0.000
      0  0.373 0.378       0.0     0.0      0.000     0.000    0.000
      0  0.378 0.383       0.0     0.0      0.000     0.000    0.000
      0  0.383 0.388       0.0     0.0      0.000     0.000    0.000
      0  0.388 0.393       0.0     0.0      0.000     0.000    0.000
      0  0.393 0.398       0.0     0.0      0.000     0.000    0.000
      0  0.398 0.403       0.0     0.0      0.000     0.000    0.000
      0  0.403 0.408       0.0     0.0      0.000     0.000    0.000
      0  0.408 0.413       0.0     0.0      0.000     0.000    0.000
      0  0.413 0.418       0.0     0.0      0.000     0.000    0.000
      0  0.418 0.423       0.0     0.0      0.000     0.000    0.000
      0  0.423 0.428       0.0     0.0      0.000     0.000    0.000
      0  0.428 0.433       0.0     0.0      0.000     0.000    0.000
      0  0.433 0.438       0.0     0.0      0.000     0.000    0.000
      0  0.438 0.443       0.0     0.0      0.000     0.000    0.000
      0  0.443 0.448       0.0     0.0      0.000     0.000    0.000
      0  0.448 0.453       0.0     0.0      0.000     0.000    0.000
      0  0.453 0.458       0.0     0.0      0.000     0.000    0.000
      0  0.458 0.463       0.0     0.0      0.000     0.000    0.000
      0  0.463 0.468       0.0     0.0      0.000     0.000    0.000
      0  0.468 0.473       0.0     0.0      0.000     0.000    0.000
      0  0.473 0.478       0.0     0.0      0.000     0.000    0.000
      0  0.478 0.483       0.0     0.0      0.000     0.000    0.000
      0  0.483 0.488       0.0     0.0      0.000     0.000    0.000
      0  0.488 0.493       0.0     0.0      0.000     0.000    0.000
      0  0.493 0.498       0.0     0.0      0.000     0.000    0.000
      0  0.498 0.503       0.0     0.0      0.000     0.000    0.000
      0  0.503 0.508       0.0     0.0      0.000     0.000    0.000
      0  0.508 0.513       0.0     0.0      0.000     0.000    0.000
      0  0.513 0.518       0.0     0.0      0.000     0.000    0.000
      0  0.518 0.523       0.0     0.0      0.000     0.000    0.000
      0  0.523 0.528       0.0     0.0      0.000     0.000    0.000
      0  0.528 0.533       0.0     0.0      0.000     0.000    0.000
      0  0.533 0.538       0.0     0.0      0.000     0.000    0.000
      0  0.538 0.543       0.0     0.0      0.000     0.000    0.000
      0  0.543 0.548       0.0     0.0      0.000     0.000    0.000
      0  0.548 0.553       0.0     0.0      0.000     0.000    0.000
      0  0.553 0.558       0.0     0.0      0.000     0.000    0.000
      0  0.558 0.563       0.0     0.0      0.000     0.000    0.000
      0  0.563 0.568       0.0     0.0      0.000     0.000    0.000
      0  0.568 0.573       0.0     0.0      0.000     0.000    0.000
      0  0.573 0.578       0.0     0.0      0.000     0.000    0.000
      0  0.578 0.583       0.0     0.0      0.000     0.000    0.000
      0  0.583 0.588       0.0     0.0      0.000     0.000    0.000
      0  0.588 0.593       0.0     0.0      0.000     0.000    0.000
      0  0.593 0.598       0.0     0.0      0.000     0.000    0.000
      0  0.598 0.603       0.0     0.0      0.000     0.000    0.000
      0  0.603 0.608       0.0     0.0      0.000     0.000    0.000
      0  0.608 0.613       0.0     0.0      0.000     0.000    0.000
      0  0.613 0.618       0.0     0.0      0.000     0.000    0.000
      0  0.618 0.623       0.0     0.0      0.000     0.000    0.000
      0  0.623 0.628       0.0     0.0      0.000     0.000    0.000
      0  0.628 0.633       0.0     0.0      0.000     0.000    0.000
      0  0.633 0.638       0.0     0.0      0.000     0.000    0.000
      0  0.638 0.643       0.0     0.0      0.000     0.000    0.000
      0  0.643 0.648       0.0     0.0      0.000     0.000    0.000
      0  0.648 0.653       0.0     0.0      0.000     0.000    0.000
      0  0.653 0.658       0.0     0.0      0.000     0.000    0.000
      0  0.658 0.663       0.0     0.0      0.000     0.000    0.000
      0  0.663 0.668       0.0     0.0      0.000     0.000    0.000
      0  0.668 0.673       0.0     0.0      0.000     0.000    0.000
      0  0.673 0.678       0.0     0.0      0.000     0.000    0.000
      0  0.678 0.683       0.0     0.0      0.000     0.000    0.000
      0  0.683 0.688       0.0     0.0      0.000     0.000    0.000
      0  0.688 0.693       0.0     0.0      0.000     0.000    0.000
      0  0.693 0.698       0.0     0.0      0.000     0.000    0.000
      0  0.698 0.703       0.0     0.0      0.000     0.000    0.000
      0  0.703 0.708       0.0     0.0      0.000     0.000    0.000
      0  0.708 0.713       0.0     0.0      0.000     0.000    0.000
      0  0.713 0.718       0.0     0.0      0.000     0.000    0.000
      0  0.718 0.723       0.0     0.0      0.000     0.000    0.000
      0  0.723 0.728       0.0     0.0      0.000     0.000    0.000
      0  0.728 0.733       0.0     0.0      0.000     0.000    0.000
      0  0.733 0.738       0.0     0.0      0.000     0.000    0.000
      0  0.738 0.743       0.0     0.0      0.000     0.000    0.000
      0  0.743 0.748       0.0     0.0      0.000     0.000    0.000
      0  0.748 0.753       0.0     0.0      0.000     0.000    0.000
      0  0.753 0.758       0.0     0.0      0.000     0.000    0.000
      0  0.758 0.763       0.0     0.0      0.000     0.000    0.000
      0  0.763 0.768       0.0     0.0      0.000     0.000    0.000
      0  0.768 0.773       0.0     0.0      0.000     0.000    0.000
      0  0.773 0.778       0.0     0.0      0.000     0.000    0.000
      0  0.778 0.783       0.0     0.0      0.000     0.000    0.000
      0  0.783 0.788       0.0     0.0      0.000     0.000    0.000
      0  0.788 0.793       0.0     0.0      0.000     0.000    0.000
      0  0.793 0.798       0.0     0.0      0.000     0.000    0.000
      0  0.798 0.803       0.0     0.0      0.000     0.000    0.000
      0  0.803 0.808       0.0     0.0      0.000     0.000    0.000
      0  0.808 0.813       0.0     0.0      0.000     0.000    0.000
      0  0.813 0.818       0.0     0.0      0.000     0.000    0.000
      0  0.818 0.823       0.0     0.0      0.000     0.000    0.000
      0  0.823 0.828       0.0     0.0      0.000     0.000    0.000
      0  0.828 0.833       0.0     0.0      0.000     0.000    0.000
      0  0.833 0.838       0.0     0.0      0.000     0.000    0.000
      0  0.838 0.843       0.0     0.0      0.000     0.000    0.000
      0  0.843 0.848       0.0     0.0      0.000     0.000    0.000
      0  0.848 0.853       0.0     0.0      0.000     0.000    0.000
      0  0.853 0.858       0.0     0.0      0.000     0.000    0.000
      0  0.858 0.863       0.0     0.0      0.000     0.000    0.000
      0  0.863 0.868       0.0     0.0      0.000     0.000    0.000
      0  0.868 0.873       0.0     0.0      0.000     0.000    0.000
      0  0.873 0.878       0.0     0.0      0.000     0.000    0.000
      0  0.878 0.883       0.0     0.0      0.000     0.000    0.000
      0  0.883 0.888       0.0     0.0      0.000     0.000    0.000
      0  0.888 0.893       0.0     0.0      0.000     0.000    0.000
      0  0.893 0.898       0.0     0.0      0.000     0.000    0.000
      0  0.898 0.903       0.0     0.0      0.000     0.000    0.000
      0  0.903 0.908       0.0     0.0      0.000     0.000    0.000
      0  0.908 0.913       0.0     0.0      0.000     0.000    0.000
      0  0.913 0.918       0.0     0.0      0.000     0.000    0.000
      0  0.918 0.923       0.0     0.0      0.000     0.000    0.000
      0  0.923 0.928       0.0     0.0      0.000     0.000    0.000
      0  0.928 0.933       0.0     0.0      0.000     0.000    0.000
      0  0.933 0.938       0.0     0.0      0.000     0.000    0.000
      0  0.938 0.943       0.0     0.0      0.000     0.000    0.000
      0  0.943 0.948       0.0     0.0      0.000     0.000    0.000
      0  0.948 0.953       0.0     0.0      0.000     0.000    0.000
      0  0.953 0.958       0.0     0.0      0.000     0.000    0.000
      0  0.958 0.963       0.0     0.0      0.000     0.000    0.000
      0  0.963 0.968       0.0     0.0      0.000     0.000    0.000
      0  0.968 0.973       0.0     0.0      0.000     0.000    0.000
      0  0.973 0.978       0.0     0.0      0.000     0.000    0.000
      0  0.978 0.983       0.0     0.0      0.000     0.000    0.000
      0  0.983 0.988       0.0     0.0      0.000     0.000    0.000
      0  0.988 0.993       0.0     0.0      0.000     0.000    0.000
      0  0.993 0.998       0.0     0.0      0.000     0.000    0.000
      0  0.998 1.000       0.0     0.0      0.000     0.000    0.000
    */
    
}
