package ffx.potential.parsers;

import ffx.numerics.estimator.MultistateBennettAcceptanceRatio;
import ffx.numerics.estimator.MultistateBennettAcceptanceRatio.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;
import static java.util.Arrays.stream;
import static org.apache.commons.lang3.math.NumberUtils.min;

/**
 * The MBARFilter class parses mbar (*.mbar or *.bar) files. Expected file format is a header
 * line including the number of snapshots contained, a name, and a temperature. Following the header
 * is a list of energies for each snapshot at each lambda value being considered with an index to start
 * the line. Then the energies go from least (0) to greatest (1) lambda value.
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
        windows = barFiles.length;
        temperatures = new double[windows];
        snaps = new int[windows];
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
            eAll[i] = readFile(barFiles[i].getName());
        }
        if (windowsRead != windows) {
            logger.severe("Failed to read all files in " + fileLocation.getAbsolutePath());
        }
    }

    /**
     * Parses the file matching the name given in the directory of 'fileLocation'.
     * @param fileName the name of the file to be parsed matching 'energy_\d+.mbar' or 'energy_\d+.bar'.
     * @return a double[][] of the energies for each snapshot at each lambda value
     */
    private double[][] readFile(String fileName) {
        tempBarFile = new File(fileLocation, fileName);
        if (!tempBarFile.exists()) {
            logger.severe("File " + tempBarFile.getAbsolutePath() + " does not exist.");
        }
        int state = Integer.parseInt(fileName.split("\\.")[0].split("_")[1]);
        tempFileEnergies = new ArrayList<>();
        for(int i = 0; i < windows; i++) {
            tempFileEnergies.add(new ArrayList<>());
        }
        try (FileReader fr1 = new FileReader(tempBarFile);
             BufferedReader br1 = new BufferedReader(fr1);) {
            String line = br1.readLine();
            String[] tokens = line.trim().split("\\t *| +");
            int numSnaps = Integer.parseInt(tokens[0]);
            temperatures[state] = Double.parseDouble(tokens[1]);
            int count = 0;
            while (line != null) {
                tokens = line.trim().split("\\t *| +");
                for (int i = 1; i < tokens.length; i++) {
                    tempFileEnergies.get(i).add(Double.parseDouble(tokens[i + 1]));
                }
                count++;
                line = br1.readLine();
            }
            snaps[state] = min(numSnaps, count);
        } catch(IOException e){
            logger.info("Failed to read MBAR file: " + tempBarFile.getAbsolutePath());
            throw new RuntimeException(e);
        }
        fileEnergies = new double[windows][];
        for (int i = 0; i < windows; i++) {
            fileEnergies[i] = new double[tempFileEnergies.get(i).size()];
            for (int j = 0; j < tempFileEnergies.get(i).size(); j++) {
                fileEnergies[i][j] = tempFileEnergies.get(i).get(j);
            }
        }
        windowsRead++;
        return fileEnergies;
    }
}
