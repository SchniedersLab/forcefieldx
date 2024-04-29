package ffx.numerics.estimator;

import ffx.numerics.estimator.MultistateBennettAcceptanceRatio.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import static org.apache.commons.lang3.math.NumberUtils.max;
import static org.apache.commons.lang3.math.NumberUtils.min;

/**
 * The MBARFilter class parses mbar (*.mbar or *.bar) files. Expected file format is a header
 * line including the number of snapshots contained, a name, and a temperature. Following the header
 * is a list of energies for each snapshot at each lambda value being considered with an index to start
 * the line. Then the energies go from least (0) to greatest (1) lambda value.
 *
 * Files < numLambda states are handled. Users should simply generate MBAR files with desired number of lambda
 * windows, this filter should handle the rest and warn about potential issues.
 *
 * @author Matthew J. Speranza
 * @since 1.0
 *
 */
public class MBARFilter {
    private static final Logger logger = Logger.getLogger(MBARFilter.class.getName());
    private File tempBarFile;
    private File[] barFiles;
    private File fileLocation;
    private ArrayList<ArrayList<Double>> tempFileEnergies;
    private double[][] fileEnergies;
    private double[][][] eAll;
    private double[] temperatures;
    private int[] snaps;
    private int[] numLambdas;
    private int windowsRead;
    private int windows;


    public MBARFilter(File fileLocation) {
        this.fileLocation = fileLocation;
        barFiles = fileLocation.listFiles((dir, name) -> name.matches("energy_\\d+.mbar") || name.matches("energy_\\d+.bar"));
        assert barFiles != null;
        if (barFiles.length == 0 || barFiles.length == 1) {
            String message = barFiles.length == 0 ?
                    " No files matching 'energy_\\d+.mbar' or 'energy_\\d+.bar' found in " +
                            fileLocation.getAbsolutePath() :
                    " Only one file matching 'energy_\\d+.mbar' or 'energy_\\d+.bar' found in " +
                            fileLocation.getAbsolutePath() + ". At least two are required.";
            logger.severe(message);
        }
        // Sort files by state number
        Arrays.sort(barFiles, (f1, f2) -> {
            int state1 = Integer.parseInt(f1.getName().split("\\.")[0].split("_")[1]);
            int state2 = Integer.parseInt(f2.getName().split("\\.")[0].split("_")[1]);
            return Integer.compare(state1, state2);
        });
        windows = barFiles.length;
        temperatures = new double[windows];
        snaps = new int[windows];
        numLambdas = new int[windows];
        this.parseFiles();
    }

    public MultistateBennettAcceptanceRatio getMBAR(SeedType seedType){
        return getMBAR(seedType, 1e-7);
    }

    public MultistateBennettAcceptanceRatio getMBAR(SeedType seedType, double tolerance) {
        double[] lambda = new double[windows];
        for (int i = 0; i < windows; i++) {
            lambda[i] = i / (windows - 1.0);
        }
        return new MultistateBennettAcceptanceRatio(lambda, eAll, temperatures, tolerance, seedType);
    }
    private void parseFiles(){
        eAll = new double[windows][][];
        for (int i = 0; i < windows; i++) {
            eAll[i] = readFile(barFiles[i].getName(), i);
        }
        if (windowsRead != windows) {
            logger.severe("Failed to read all files in " + fileLocation.getAbsolutePath());
        }
        int minSnaps = min(snaps);
        int maxSnaps = max(snaps);

        // Basically just make sure eAll isn't jagged
        boolean warn = minSnaps != maxSnaps;
        if (warn) {
            logger.warning("NOT ALL FILES CONTAINED THE SAME NUMBER OF SNAPSHOTS. ");
            logger.warning("SAMPLES PER WINDOW: " + Arrays.toString(snaps));
            double[][][] temp = new double[eAll.length][eAll[0].length][maxSnaps];
            for(int j = 0; j < windows; j++){
                for(int k = 0; k < windows; k++){
                    System.arraycopy(eAll[j][k], 0, temp[j][k], 0, snaps[j]);
                    for(int l = snaps[j]; l < maxSnaps; l++){ // Fill in the rest with NaNs
                        temp[j][k][l] = Double.NaN;
                    }
                }
            }
            eAll = temp;
        }
        // Fail if not all files have the same number of energy evaluations across lambda
        int maxLambdas = max(numLambdas);
        for (int i = 0; i < windows; i++) {
            if(numLambdas[i] != maxLambdas){
                logger.severe(" Number of lambda evaluations in file " + barFiles[i].getName() +
                        " does not match the number of lambda evaluations in the other files. This is unrecoverable.");
            }
        }
        // Handle files with more lambda windows than actual trajectories
        warn = maxLambdas != windows;
        if (warn) {
            logger.warning("NOT ALL FILES CONTAINED THE SAME NUMBER OF LAMBDA EVALUATIONS. FILLING IN MISSING " +
                    "LAMBDA EVALUATIONS WITH NaN (*end-states assumed to have sampling).");
            logger.warning("SEEDING WITH ZEROS");
            MultistateBennettAcceptanceRatio.FORCE_ZEROS_SEED = true;
            logger.warning("LAMBDA EVALUATIONS PER WINDOW: " + Arrays.toString(numLambdas));
            int diff = maxLambdas - windows;
            int gaps = windows-1;
            int numBetween = diff/gaps;
            if (gaps * numBetween + windows != maxLambdas) {
                logger.info("Gaps: " + gaps + " NumBetween: " + numBetween + " Windows: " + windows + " MaxLambdas: " + maxLambdas);
                logger.severe("Failed to fill in missing lambda evaluations evenly. Calculation: " + gaps + " * " + numBetween + " + " + windows + " != " + maxLambdas);
            }
            // Arraylist takes longer but is easier to write
            ArrayList<double[][]> e = new ArrayList<>(maxLambdas);
            ArrayList<Double> t = new ArrayList<>(maxLambdas);
            ArrayList<Integer> s = new ArrayList<>(maxLambdas);
            for (int i = 0; i < gaps; i++){
                e.add(eAll[i]);
                t.add(temperatures[i]);
                s.add(snaps[i]);
                for (int j = 0; j < numBetween; j++){
                     double[][] gapE = new double[maxLambdas][maxSnaps];
                     for (int k = 0; k < maxLambdas; k++){
                         for (int l = 0; l < maxSnaps; l++){
                             gapE[k][l] = Double.NaN;
                         }
                     }
                     e.add(gapE);
                    //TODO: Fix atrocious setting of temperatures
                     t.add(298.0);
                     s.add(0);
                }
                if (i == gaps-1){ // Add the last one that got missed by looping over gaps
                    e.add(eAll[i+1]);
                    t.add(temperatures[i+1]);
                    s.add(snaps[i+1]);
                }
            }
            temperatures = t.stream().mapToDouble(Double::doubleValue).toArray();
            snaps = s.stream().mapToInt(Integer::intValue).toArray();
            double[][][] temp = new double[maxLambdas][][];
            for (int i = 0; i < maxLambdas; i++) {
                temp[i] = e.get(i);
            }
            eAll = temp;
            windows = maxLambdas;
        }
    }

    public void writeFiles(File mbarFileLoc, double[][][] energies, double[] temperatures) {
        if (temperatures.length != windows) {
            double temp = temperatures[0];
            temperatures = new double[windows];
            for (int i = 0; i < windows; i++) {
                temperatures[i] = temp;
            }
        }
        for (int i = 0; i < windows; i++) {
            File file = new File(mbarFileLoc, "energy_" + i + ".mbar");
            writeFile(energies[i], file, temperatures[i]);
        }
    }

    /**
     * Parses the file matching the name given in the directory of 'fileLocation'.
     * @param fileName the name of the file to be parsed matching 'energy_\d+.mbar' or 'energy_\d+.bar'.
     * @return a double[][] of the energies for each snapshot at each lambda value
     */
    private double[][] readFile(String fileName, int state) {
        tempBarFile = new File(fileLocation, fileName);
        if (!tempBarFile.exists()) {
            logger.severe("File " + tempBarFile.getAbsolutePath() + " does not exist.");
        }
        tempFileEnergies = new ArrayList<>();
        for(int i = 0; i < windows; i++) {
            tempFileEnergies.add(new ArrayList<>());
        }
        try (FileReader fr1 = new FileReader(tempBarFile);
             BufferedReader br1 = new BufferedReader(fr1);) {
            // Read header
            String line = br1.readLine();
            String[] tokens = line.trim().split("\\t *| +");
            temperatures[state] = Double.parseDouble(tokens[1]);
            // Read energies (however many there are)
            int count = 0;
            int numLambda = 0;
            line = br1.readLine();
            while (line != null) {
                tokens = line.trim().split("\\t *| +");
                numLambda = tokens.length - 1;
                for (int i = 1; i < tokens.length; i++) {
                    if (tempFileEnergies.size() < i){
                        tempFileEnergies.add(new ArrayList<>());
                    }
                    tempFileEnergies.get(i-1).add(Double.parseDouble(tokens[i]));
                }
                count++;
                line = br1.readLine();
            }
            numLambdas[state] = numLambda;
            snaps[state] = count;
        } catch(IOException e){
            logger.info("Failed to read MBAR file: " + tempBarFile.getAbsolutePath());
            throw new RuntimeException(e);
        }
        // Convert to double[][]
        fileEnergies = new double[tempFileEnergies.size()][];
        for (int i = 0; i < tempFileEnergies.size(); i++) {
            fileEnergies[i] = new double[tempFileEnergies.get(i).size()];
            for (int j = 0; j < tempFileEnergies.get(i).size(); j++) {
                fileEnergies[i][j] = tempFileEnergies.get(i).get(j);
            }
        }
        windowsRead++;
        return fileEnergies;
    }

    public void writeFile(double[][] energies, File file, double temperature) {
        MultistateBennettAcceptanceRatio.writeFile(energies, file, temperature);
    }
}
