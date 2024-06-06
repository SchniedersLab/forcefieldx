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
    private MultistateBennettAcceptanceRatio mbar;


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
        this.mbar = new MultistateBennettAcceptanceRatio(lambda, eAll, temperatures, tolerance, seedType);
        return this.mbar;
    }

    /**
     * 10% of the total samples at different time points.
     * @param seedType
     * @param tol
     * @return an array of MBAR objects
     */
    public MultistateBennettAcceptanceRatio[] getPeriodComparisonMBAR(SeedType seedType, double tol){
        double[] lambda = new double[windows];
        for (int i = 0; i < windows; i++) {
            lambda[i] = i / (windows - 1.0);
        }
        MultistateBennettAcceptanceRatio[] mbar = new MultistateBennettAcceptanceRatio[10];
        for (int i = 0; i < 10; i++){
            double[][][] e = new double[windows][][];
            int maxSamples = max(snaps);
            int timePeriod = maxSamples / 10;
            for (int j = 0; j < windows; j++) {
                e[j] = new double[windows][];
                for (int k = 0; k < windows; k++) {
                    e[j][k] = new double[timePeriod];
                    if (timePeriod * (i + 1) > maxSamples) {
                        System.arraycopy(eAll[j][k], timePeriod * i, e[j][k], 0, maxSamples - timePeriod * i);
                    } else {
                        System.arraycopy(eAll[j][k], timePeriod * i, e[j][k], 0, timePeriod);
                    }
                }
            }
            logger.info(" Period: " + (timePeriod*i) + " - " + (timePeriod*(i+1)) + " samples calculation.");
            mbar[i] = new MultistateBennettAcceptanceRatio(lambda, e, temperatures, tol, seedType);
        }
        return mbar;
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
            logger.warning("FILES CONTAIN MORE LAMBDA WINDOWS THAN ACTUAL TRAJECTORIES. ");
            logger.severe("Create completely empty files (zero lines) to fill in the gaps.");
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
            if (line == null) { // Empty file
                for(int i = 0; i < windows; i++){
                    tempFileEnergies.get(i).add(Double.NaN);
                }
                snaps[state] = 0;
                temperatures[state] = 298; // Assumed default temp since 0 leads to division by zero
                MultistateBennettAcceptanceRatio.FORCE_ZEROS_SEED = true;
                if (state != 0) {
                    numLambdas[state] = numLambdas[state - 1];
                }
                windowsRead++;
                fileEnergies = new double[windows][];
                for (int i = 0; i < windows; i++) {
                    fileEnergies[i] = new double[tempFileEnergies.get(i).size()];
                    for (int j = 0; j < tempFileEnergies.get(i).size(); j++) {
                        fileEnergies[i][j] = tempFileEnergies.get(i).get(j);
                    }
                }
                return fileEnergies;
            }
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
            if(state != 0 && numLambdas[0] == 0){ // If the zeroth window is missing this wasn't set yet
                numLambdas[0] = numLambda;
            }
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

    public void setStartSnapshot(int startIndex) {
        for (int i = 0; i < windows; i++) {
            for (int j = 0; j < windows; j++) {
                try{
                    eAll[i][j] = Arrays.copyOfRange(eAll[i][j], startIndex, eAll[i][j].length);
                } catch (ArrayIndexOutOfBoundsException e){
                    logger.severe("Start index " + startIndex + " is out of bounds for file " + barFiles[i].getName());
                }
            }
        }
    }

    public void setEndSnapshot(int endIndex) {
        for (int i = 0; i < windows; i++) {
            for (int j = 0; j < windows; j++) {
                try{
                    eAll[i][j] = Arrays.copyOfRange(eAll[i][j], 0, endIndex);
                } catch (ArrayIndexOutOfBoundsException e){
                    logger.severe("End index " + endIndex + " is out of bounds for file " + barFiles[i].getName());
                }
            }
        }
    }

    /**
     * Read in observable data, try to leave as many fields in-tact as possible.
     * @param parentDirectory
     * @param multiDataObservable
     */
    public void readObservableData(File parentDirectory, boolean multiDataObservable) {
        barFiles = fileLocation.listFiles((dir, name) -> name.matches("derivative_\\d+.mbar") ||
                name.matches("derivative_\\d+.bar") ||
                name.matches("derivatives_\\d+.mbar") ||
                name.matches("derivatives_\\d+.bar"));
        assert barFiles != null && barFiles.length > 0;
        // Sort files by state number
        Arrays.sort(barFiles, (f1, f2) -> {
            int state1 = Integer.parseInt(f1.getName().split("\\.")[0].split("_")[1]);
            int state2 = Integer.parseInt(f2.getName().split("\\.")[0].split("_")[1]);
            return Integer.compare(state1, state2);
        });
        eAll = new double[windows][][];
        for (int i = 0; i < windows; i++) {
            eAll[i] = readFile(barFiles[i].getName(), i);
        }
        int minSnaps = min(snaps);
        int maxSnaps = max(snaps);

        // Basically just make sure eAll isn't jagged
        boolean warn = minSnaps != maxSnaps;
        if (warn) {
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
        // Set observable data and compute observable averages
        mbar.setObservableData(eAll, multiDataObservable);
    }
}
