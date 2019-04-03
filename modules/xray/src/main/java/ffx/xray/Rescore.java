//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import edu.rit.pj.ParallelTeam;

import ffx.algorithms.AlgorithmFunctions;
import ffx.potential.MolecularAssembly;
import ffx.realspace.RealSpaceData;
import ffx.realspace.parsers.RealSpaceFile;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.Keyword;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.parsers.DiffractionFile;
import static ffx.algorithms.misc.ClusterStructures.generatePath;
import static ffx.xray.Rescore.RescoreStrategy.NO_RESCORE;

/**
 * This class performs rescoring on a provided list of structure files.
 * Rescoring can be based on energy evaluations, minimization, x-ray
 * minimization, or real-space minimization.
 *
 * @author Jacob Litman
 * @since 1.0
 */
public class Rescore {

    private static final Logger logger = Logger.getLogger(Rescore.class.getName());
    private final RefinementMode refinementMode = RefinementMode.COORDINATES;
    private final AlgorithmFunctions utils;
    private RescoreStrategy rscType = NO_RESCORE;
    private List<DiffractionFile> diffractionFiles;
    private List<RealSpaceFile> mapFiles;

    private double eps = -1.0;
    private int maxiter = 1000;
    private double acceptThreshold = 0.0;
    private int numToAccept = 10;
    private boolean includeRejected = true;

    private String fileSuffix = "_rsc";
    private final Path pwdPath;
    private File resultsFile;
    private File resultDir;
    private Path resultPath;
    private boolean printModels = false;

    /**
     * <p>Constructor for Rescore.</p>
     *
     * @param utils a {@link ffx.algorithms.AlgorithmFunctions} object.
     */
    public Rescore(AlgorithmFunctions utils) {
        this.pwdPath = generatePath(new File(""));
        diffractionFiles = new ArrayList<>();
        mapFiles = new ArrayList<>();
        this.utils = utils;
    }

    /**
     * <p>setRefineEps.</p>
     *
     * @param eps a double.
     */
    public void setRefineEps(double eps) {
        this.eps = eps;
    }

    /**
     * <p>setMaxIter.</p>
     *
     * @param maxiter a int.
     */
    public void setMaxIter(int maxiter) {
        this.maxiter = maxiter;
    }

    /**
     * <p>setResultDirectory.</p>
     *
     * @param directory a {@link java.io.File} object.
     * @throws java.io.IOException if any.
     */
    public void setResultDirectory(File directory) throws IOException {
        if (directory == null) {
            resultDir = null;
        } else if (directory.exists() && !directory.isDirectory()) {
            throw new IOException(" Results directory could not be recognized as folder or directory.");
        } else {
            directory.mkdirs();
            resultDir = directory;
            resultPath = pwdPath.relativize(generatePath(resultDir));
        }
    }

    /**
     * <p>Setter for the field <code>acceptThreshold</code>.</p>
     *
     * @param acceptThreshold a double.
     */
    public void setAcceptThreshold(double acceptThreshold) {
        this.acceptThreshold = acceptThreshold;
    }

    /**
     * <p>Setter for the field <code>numToAccept</code>.</p>
     *
     * @param numToAccept a int.
     */
    public void setNumToAccept(int numToAccept) {
        this.numToAccept = numToAccept;
    }

    /**
     * <p>Setter for the field <code>includeRejected</code>.</p>
     *
     * @param includeRejected a boolean.
     */
    public void setIncludeRejected(boolean includeRejected) {
        this.includeRejected = includeRejected;
    }

    /**
     * <p>Setter for the field <code>fileSuffix</code>.</p>
     *
     * @param fileSuffix a {@link java.lang.String} object.
     */
    public void setFileSuffix(String fileSuffix) {
        this.fileSuffix = fileSuffix;
    }

    /**
     * <p>Setter for the field <code>resultsFile</code>.</p>
     *
     * @param resultsFile a {@link java.io.File} object.
     */
    public void setResultsFile(File resultsFile) {
        if (resultsFile != null) {
            this.resultsFile = resultsFile;
        }
    }

    /**
     * <p>Setter for the field <code>printModels</code>.</p>
     *
     * @param printModels a boolean.
     */
    public void setPrintModels(boolean printModels) {
        this.printModels = printModels;
    }

    /**
     * <p>setRescoreStrategy.</p>
     *
     * @param rscType a {@link ffx.xray.Rescore.RescoreStrategy} object.
     */
    public void setRescoreStrategy(RescoreStrategy rscType) {
        this.rscType = rscType;
    }

    /**
     * <p>Setter for the field <code>diffractionFiles</code>.</p>
     *
     * @param diffractionFiles a {@link java.util.List} object.
     */
    public void setDiffractionFiles(List<DiffractionFile> diffractionFiles) {
        this.diffractionFiles.addAll(diffractionFiles);
    }

    /**
     * <p>Setter for the field <code>mapFiles</code>.</p>
     *
     * @param mapFiles a {@link java.util.List} object.
     */
    public void setMapFiles(List<RealSpaceFile> mapFiles) {
        this.mapFiles.addAll(mapFiles);
    }

    /**
     * Launch the rescoring algorithm on provided files. Assumes it has been
     * given valid files to be run on; use CoordinateFileFilter.acceptDeep(File
     * file) before sending files to this method.
     *
     * @param modelFiles Files to xray.groovy.test.rescore.
     */
    public void runRsc(File[] modelFiles) {
        int numFiles = modelFiles.length;
        logger.info(String.format(" Rescoring %d files", numFiles));
        logger.info(String.format(" Rescore algorithm: %s", rscType.toString()));
        rescore(modelFiles);
    }

    private File rescoreSingle(File modelFile, RescoreStrategy rscType, DoubleIndexPair[] energies, int i) {
        Path filepath = generatePath(modelFile);
        if (filepath == null) {
            logger.warning(String.format(" Could not generate path to file %s", modelFile.toPath()));
            return null;
        }
        String filename = pwdPath.relativize(filepath).toString();
        File retFile = modelFile;
        try {
            MolecularAssembly assembly = utils.open(filename);
            switch (rscType) {
                case NO_RESCORE:
                    logger.warning(" Rescore is being called with rscType = NO_RESCORE");
                    break;
                case ENERGY_EVAL:
                    break;
                case MINIMIZE:
                    logger.info(String.format("\n Running minimize on %s", filename));
                    logger.info(String.format(" RMS gradient convergence criteria: %f", eps));
                    utils.energy(assembly);
                    utils.minimize(assembly, eps);

                    String ext = FilenameUtils.getExtension(filename);
                    ext = ".".concat(ext);

                    if (resultDir != null) {
                        filename = FilenameUtils.getBaseName(filename);
                        filename = FilenameUtils.concat(resultPath.toString(), filename);
                    } else {
                        filename = FilenameUtils.removeExtension(filename);
                    }
                    filename = filename.concat(fileSuffix).concat(ext);
                    retFile = new File(filename);
                    if (ext.toUpperCase().contains("XYZ")) {
                        utils.saveAsXYZ(assembly, retFile);
                    } else {
                        utils.saveAsPDB(assembly, retFile);
                    }
                    break;
                case XRAY_MIN:
                    logger.info(String.format("\n Running x-ray minimize on %s", filename));

                    DiffractionFile diffractionFile = null;
                    if (diffractionFiles.isEmpty()) {
                        diffractionFile = new DiffractionFile(assembly, 1.0, false);
                        diffractionFiles.add(diffractionFile);
                    }
                    CompositeConfiguration properties = Keyword.loadProperties(modelFile);
                    DiffractionData diffractionData = new DiffractionData(assembly, properties,
                            SolventModel.POLYNOMIAL, diffractionFiles.toArray(new DiffractionFile[0]));

                    diffractionData.scaleBulkFit();
                    diffractionData.printStats();
                    utils.energy(assembly);

                    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData, refinementMode);
                    if (eps < 0.0) {
                        eps = refinementMinimize.getEps();
                    }

                    if (maxiter < Integer.MAX_VALUE) {
                        logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps, maxiter));
                    } else {
                        logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", eps));
                    }

                    refinementMinimize.minimize(eps, maxiter);
                    diffractionData.scaleBulkFit();
                    diffractionData.printStats();

                    ext = FilenameUtils.getExtension(filename);
                    ext = ".".concat(ext);

                    if (resultDir != null) {
                        filename = FilenameUtils.getBaseName(filename);
                        filename = FilenameUtils.concat(resultPath.toString(), filename);
                    } else {
                        filename = FilenameUtils.removeExtension(filename);
                    }
                    filename = filename.concat(fileSuffix);
                    diffractionData.writeData(filename + ".mtz");
                    filename = filename.concat(ext);
                    diffractionData.writeModel(filename);

                    retFile = new File(filename);
                    if (diffractionFile != null) {
                        try {
                            diffractionFiles.remove(diffractionFile);
                        } catch (UnsupportedOperationException ex) {
                            // This should never occur, because diffractionFiles should be of a List type supporting remove(object).
                            diffractionFiles = new ArrayList<>();
                        }
                    }
                    break;
                case RS_MIN:
                    logger.info(String.format("\n Running real-space minimize on %s", filename));

                    RealSpaceFile realspaceFile = null;
                    if (mapFiles.isEmpty()) {
                        realspaceFile = new RealSpaceFile(assembly);
                        mapFiles.add(realspaceFile);
                    }
                    properties = Keyword.loadProperties(modelFile);
                    RealSpaceData realspaceData = new RealSpaceData(assembly,
                            properties, new ParallelTeam(),
                            mapFiles.toArray(new RealSpaceFile[0]));
                    utils.energy(assembly);

                    refinementMinimize = new RefinementMinimize(realspaceData, refinementMode);
                    if (eps < 0.0) {
                        eps = 1.0;
                    }

                    if (maxiter < Integer.MAX_VALUE) {
                        logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps, maxiter));
                    } else {
                        logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", eps));
                    }

                    ext = FilenameUtils.getExtension(filename);
                    ext = ".".concat(ext);
                    if (resultDir != null) {
                        filename = FilenameUtils.getBaseName(filename);
                        filename = FilenameUtils.concat(resultPath.toString(), filename);
                    } else {
                        filename = FilenameUtils.removeExtension(filename);
                    }
                    filename = filename.concat(fileSuffix).concat(ext);
                    retFile = new File(filename);
                    if (ext.toUpperCase().contains("XYZ")) {
                        utils.saveAsXYZ(assembly, retFile);
                    } else {
                        utils.saveAsPDB(assembly, retFile);
                    }

                    if (realspaceFile != null) {
                        try {
                            mapFiles.remove(realspaceFile);
                        } catch (UnsupportedOperationException ex) {
                            // This should never occur, because diffractionFiles should be of a List type supporting remove(object).
                            mapFiles = new ArrayList<>();
                        }
                    }
                    break;
                default:
                    logger.severe(" No valid xray.groovy.test.rescore type: FFX will not continue.");
            }
            double e = utils.returnEnergy(assembly);
            energies[i] = new DoubleIndexPair(i, e);
            utils.close(assembly);
        } catch (Exception ex) {
            logger.warning(String.format(" Exception rescoring on file %s", filename));
            logger.info(ex.toString());
        }
        return retFile;
    }

    /**
     * <p>rescore.</p>
     *
     * @param modelFiles an array of {@link java.io.File} objects.
     * @return an array of {@link java.io.File} objects.
     */
    public File[] rescore(File[] modelFiles) {
        int numFiles = modelFiles.length;
        DoubleIndexPair[] energies = new DoubleIndexPair[numFiles];
        File[] rescoredFiles = new File[numFiles];
        for (int i = 0; i < numFiles; i++) {
            rescoredFiles[i] = rescoreSingle(modelFiles[i], rscType, energies, i);
        }

        Arrays.sort(energies);
        ArrayList<Integer> acceptedFileIndices = new ArrayList<>();
        int numAccepted = 1;
        acceptedFileIndices.add(energies[0].getIndex());
        double minEnergy = energies[0].getDoubleValue();

        if (acceptThreshold > 0.0) {
            for (int i = 1; i < numFiles; i++) {
                DoubleIndexPair currentPair = energies[i];
                double e = currentPair.getDoubleValue();
                if (e - minEnergy <= acceptThreshold) {
                    acceptedFileIndices.add(currentPair.getIndex());
                    ++numAccepted;
                } else {
                    break;
                }
            }
        } else {
            numAccepted = numToAccept < numFiles ? numToAccept : numFiles;
            for (int i = 1; i < numAccepted; i++) {
                acceptedFileIndices.add(energies[i].getIndex());
            }
        }

        for (int i = 0; i < numAccepted; i++) {
            File filei = rescoredFiles[energies[i].getIndex()];
            Path pathi;
            pathi = generatePath(filei);
            String relPath = pwdPath.relativize(pathi).toString();
            logger.info(String.format(" Accepted: %s at %10.6f kcal/mol", relPath, energies[i].getDoubleValue()));
        }
        for (int i = numAccepted; i < numFiles; i++) {
            File filei = rescoredFiles[energies[i].getIndex()];
            Path pathi;
            pathi = generatePath(filei);
            String relPath = pwdPath.relativize(pathi).toString();
            logger.info(String.format(" Rejected: %s at %10.6f kcal/mol", relPath, energies[i].getDoubleValue()));
        }

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(resultsFile, true))) {
            Path rscFilePath = generatePath(resultsFile);
            String rscFileName = pwdPath.relativize(rscFilePath).toString();

            logger.info(String.format(" Printing accepted files to xray.groovy.test.rescore file %s", rscFileName));

            if (acceptThreshold > 0.0) {
                String message = String.format("Minimum potential energy: %f, threshold = %6.4f", minEnergy, acceptThreshold);
                bw.write(message);
                message = " " + message + "\n";
                logger.info(message);
            } else {
                double maxEnergy = energies[numAccepted - 1].getDoubleValue();
                String message = String.format("Minimum potential energy: %f, maximum accepted energy %f", minEnergy, maxEnergy);
                bw.write(message);
                message = " " + message + "\n";
                logger.info(message);
            }
            bw.newLine();
            bw.newLine();
            String message = String.format("Number of files accepted: %d", numAccepted);
            bw.write(message);
            bw.newLine();
            logger.info(String.format(" %s", message));

            for (int i = 0; i < numAccepted; i++) {
                int fileIndex = energies[i].getIndex();
                File pointedFile = rescoredFiles[fileIndex];
                Path pointedPath = generatePath(pointedFile);
                String relPath = pwdPath.relativize(pointedPath).toString();
                double thisEnergy = energies[i].getDoubleValue();

                logger.info(String.format(" Accepted file %d energy %9.3f < %9.3f kcal/mol", (i + 1), thisEnergy, minEnergy + acceptThreshold));
                logger.info(String.format(" %s", relPath));
                try {
                    bw.write(String.format("Accepted file: %s rank %d energy %f\n", relPath, (i + 1), thisEnergy));
                    if (printModels) {
                        bw.newLine();
                        try (BufferedReader br = new BufferedReader(new FileReader(pointedFile))) {
                            String line = br.readLine();
                            while (line != null) {
                                bw.write(line);
                                bw.newLine();
                                line = br.readLine();
                            }
                        }
                        bw.newLine();
                    }
                } catch (IOException ex) {
                    logger.warning(String.format(" File %s had exception printing to xray.groovy.test.rescore file %s", relPath, ex.toString()));
                }
            }
            message = String.format("\n Number of files not accepted: %d", numFiles - numAccepted);
            logger.info(String.format(" %s", message));
            bw.newLine();
            bw.write(message);
            bw.newLine();
            if (includeRejected) {
                for (int i = numAccepted; i < numFiles; i++) {
                    int fileIndex = energies[i].getIndex();
                    File pointedFile = rescoredFiles[fileIndex];
                    Path pointedPath = generatePath(pointedFile);
                    String relPath = pwdPath.relativize(pointedPath).toString();
                    double thisEnergy = energies[i].getDoubleValue();

                    message = String.format("Non-accepted file: %s rank %d energy %f", relPath, (i + 1), thisEnergy);
                    logger.info(String.format(" %s", message));
                    try {
                        bw.write(message);
                        bw.newLine();
                    } catch (IOException ex) {
                        logger.warning(String.format(" File %s had exception printing to xray.groovy.test.rescore file %s", relPath, ex.toString()));
                    }
                }
            }
        } catch (IOException ex) {
            logger.warning(String.format(" Exception in writing xray.groovy.test.rescore file: %s", ex.toString()));
        }
        return rescoredFiles;
    }

    public enum RescoreStrategy {

        NO_RESCORE {
            @Override
            public String toString() {
                return "none";
            }
        },
        ENERGY_EVAL {
            @Override
            public String toString() {
                return "potential energy";
            }
        },
        MINIMIZE {
            @Override
            public String toString() {
                return "force field minimization";
            }
        },
        XRAY_MIN {
            @Override
            public String toString() {
                return "x-ray hybrid target minimization";
            }
        },
        RS_MIN {
            @Override
            public String toString() {
                return "real-space hybrid target minimization";
            }
        };
    }
}
