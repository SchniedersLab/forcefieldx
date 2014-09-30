/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */

// Copyright license for hierarchical-clustering-java
/*******************************************************************************
 * Copyright 2013 Lars Behnke
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

// Copyright license for BioJava
/* 
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 7/8/2014
 *
 */

package ffx.xray;

import com.apporiented.algorithm.clustering.Cluster;
import ffx.algorithms.AlgorithmFunctions;
import ffx.potential.bonded.MolecularAssembly;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.Keyword;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;
import static ffx.xray.RescoreAndCluster.ClustAlg.*;
import static ffx.xray.RescoreAndCluster.RescoreStrategy.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.io.PDBFileReader;

/**
 * This class performs rescoring and clustering on a provided list of structure 
 * files. Rescoring can be based on energy evaluations, minimization, x-ray 
 * minimization, or real-space minimization. Clustering is based on all-atom 
 * RMSD (using BioJava libraries) and hierarchical clustering from Lars Behnke's
 * Java clustering package.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class RescoreAndCluster {
    private static final Logger logger = Logger.getLogger(RescoreAndCluster.class.getName());
    private final RefinementMode refinementMode = RefinementMode.COORDINATES;
    private final AlgorithmFunctions utils;
    private RescoreStrategy rscType = NO_RESCORE;
    private ClustAlg clusterAlg = NO_CLUSTERS;
    private List<DiffractionFile> diffractionFiles; 
    private List<RealSpaceFile> mapFiles;
    
    private boolean doRescore = false; //
    private boolean doCluster = false; //
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
    
    public RescoreAndCluster(AlgorithmFunctions utils) {
        this.pwdPath = generatePath(new File(""));
        diffractionFiles = new ArrayList<>();
        mapFiles = new ArrayList<>();
        this.utils = utils;
    }
    
    public void setRefineEps (double eps) {
        this.eps = eps;
    }
    
    public void setMaxIter (int maxiter) {
        this.maxiter = maxiter;
    }
    
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
    
    public void setAcceptThreshold (double acceptThreshold) {
        this.acceptThreshold = acceptThreshold;
    }
    
    public void setNumToAccept (int numToAccept) {
        this.numToAccept = numToAccept;
    }
    
    public void setIncludeRejected (boolean includeRejected) {
        this.includeRejected = includeRejected;
    }
    
    public void setFileSuffix (String fileSuffix) {
        this.fileSuffix = fileSuffix;
    }
    
    public void setResultsFile(File resultsFile) {
        if (resultsFile != null) {
            this.resultsFile = resultsFile;
        }
    }
    
    public void setPrintModels (boolean printModels) {
        this.printModels = printModels;
    }
    
    public void setRescoreStrategy(RescoreStrategy rscType) {
        this.rscType = rscType;
        if (rscType != NO_RESCORE) {
            //energies = new DoubleIndexPair[numFiles];
        } else {
            //rescoredFiles = modelFiles;
        }
    }
    
    public void setClusterAlg(ClustAlg clusterAlg) {
        this.clusterAlg = clusterAlg;
    }
    
    public void setDiffractionFiles(List<DiffractionFile> diffractionFiles) {
        this.diffractionFiles.addAll(diffractionFiles);
    }
    
    public void setMapFiles(List<RealSpaceFile> mapFiles) {
        this.mapFiles.addAll(mapFiles);
    }
    
    /**
     * Utility method which attempts to generate a file Path using the canonical
     * path string, else uses the absolute path.
     * @param file To find path of
     * @return Canonical or absolute path.
     */
    private Path generatePath (File file) {
        Path path;
        try {
            path = Paths.get(file.getCanonicalPath());
        } catch (IOException ex) {
            path = Paths.get(file.getAbsolutePath());
        }
        return path;
    }
    
    public Cluster cluster(File[] clustFiles) {
        int numFiles = clustFiles.length;
        double[][] distanceMatrix = new double[numFiles][numFiles];
        PDBFileReader reader = new PDBFileReader();
        StructurePairAligner aligner = new StructurePairAligner();
        try {
            for (int i = 0; i < numFiles; i++) {
                Structure strucI = reader.getStructure(clustFiles[i]);
                for (int j = i + 1; j < numFiles; j++) {
                    Structure strucJ = reader.getStructure(clustFiles[j]);
                    aligner.align(strucI, strucJ);
                    AlternativeAlignment[] alignments = aligner.getAlignments();
                    double bestRMSD = alignments[0].getRmsd();
                    for (int k = 1; k < alignments.length; k++) {
                        double alignRMSD = alignments[k].getRmsd();
                        bestRMSD = alignRMSD < bestRMSD ? alignRMSD : bestRMSD;
                    }
                    distanceMatrix[i][j] = bestRMSD;
                }
            }
        } catch (StructureException ex) {
            
        } catch (IOException ex) {
            Logger.getLogger(RescoreAndCluster.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }
    
    /**
     * Launch the rescoring/clustering algorithms on provided files. Assumes it 
     * has been given valid files to be run on; use 
     * CoordinateFileFilter.acceptDeep(File file) before sending files to this method.
     * @param modelFiles Files to rescore and/or cluster.
     */
    public void runRsc(File[] modelFiles) {
        int numFiles = modelFiles.length;
        if (rscType != NO_RESCORE) {
            if (clusterAlg != NO_CLUSTERS) {
                logger.info(String.format(" Rescoring and clustering %d files", numFiles));
                logger.info(String.format(" Rescore algorithm: %s, clustering algorithm: %s", rscType.toString(), clusterAlg.toString()));
                cluster(rescore(modelFiles)); // Cluster the files returned by rescore.
            } else {
                logger.info(String.format(" Rescoring %d files", numFiles));
                logger.info(String.format(" Rescore algorithm: %s", rscType.toString()));
                rescore(modelFiles);
            }
        } else if (clusterAlg != NO_CLUSTERS) {
            logger.info(String.format(" Clustering %d files", numFiles));
            logger.info(String.format(" Clustering algorithm: %s", clusterAlg.toString()));
            cluster(modelFiles);
        } else {
            logger.info(" No rescoring or clustering algorithm selected.");
        }
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
            MolecularAssembly[] openedAssemblies = utils.open(filename);
            MolecularAssembly assembly = openedAssemblies[0];
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
                    DiffractionData diffractionData = new DiffractionData(assembly, properties, SolventModel.POLYNOMIAL, diffractionFiles.toArray(new DiffractionFile[diffractionFiles.size()]));

                    diffractionData.scaleBulkFit();
                    diffractionData.printStats();
                    utils.energy(assembly);

                    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData, refinementMode);
                    if (eps < 0.0) {
                        eps = refinementMinimize.getEps();
                    }
                    logger.info(String.format("\n RMS gradient convergence criteria: %8.5f max number of iterations %d", eps, maxiter));
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
                    RealSpaceData realspaceData = new RealSpaceData(assembly, properties, mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
                    utils.energy(assembly);

                    refinementMinimize = new RefinementMinimize(realspaceData, refinementMode);
                    if (eps < 0.0) {
                        eps = 1.0;
                    }
                    logger.info(String.format("\n RMS gradient convergence criteria: %8.5f max number of iterations %d", eps, maxiter));
                    refinementMinimize.minimize(eps, maxiter);

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
                    logger.severe(" No valid rescore type: FFX will not continue.");
            }
            double e = utils.returnEnergy(assembly);
            energies[i] = new DoubleIndexPair(i, e);
            utils.close();
        } catch (Exception ex) {
            logger.warning(String.format(" Exception rescoring on file %s", filename));
            logger.info(ex.toString());
        }
        return retFile;
    }
    
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

            logger.info(String.format(" Printing accepted files to rescore file %s", rscFileName));

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
                
                logger.info(String.format(" Accepted file %d energy %9.3f < %9.3f kcal/mol", (i+1), thisEnergy, minEnergy + acceptThreshold));
                logger.info(String.format(" %s", relPath));
                try {
                    bw.write(String.format("Accepted file: %s rank %d energy %f\n", relPath, (i+1), thisEnergy));
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
                    logger.warning(String.format(" File %s had exception printing to rescore file %s", relPath, ex.toString()));
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
                    
                    message = String.format("Non-accepted file: %s rank %d energy %f", relPath, (i+1), thisEnergy);
                    logger.info(String.format(" %s", message));
                    try {
                        bw.write(message);
                        bw.newLine();
                    } catch (IOException ex) {
                        logger.warning(String.format(" File %s had exception printing to rescore file %s", relPath, ex.toString()));
                    }
                }
            }
        } catch (IOException ex) {
            logger.warning(String.format(" Exception in writing rescore file: %s", ex.toString()));
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
    public enum ClustAlg {
        /**
         * All algorithms start with each point a cluster, and then join the closest
         * clusters together until everything is one cluster. SLINK is Single 
         * Linkage; cluster-cluster distance is defined by the nearest two points.
         * This is vulnerable to chaining; two clusters might be joined by a handful
         * of intermediate points. CLINK is Complete Linkage; CLINK uses the greatest
         * distance between points in two clusters. AV_LINK (average link) is the
         * UPGMA (Unweighted Pair Group Method with Arithmetic Mean) function, which
         * takes the mean distance between points in a cluster.
         * 
         * Makes me wonder if there's a WPGMA algorithm which does weight one
         * way or the other, or perhaps a RPGMA RMSD-like algorithm.
         */
        NO_CLUSTERS {
            @Override
            public String toString() {
                return "none";
            }
        },
        SLINK {
            @Override
            public String toString() {
                return "single linkage";
            }
        },
        AV_LINK {
            @Override
            public String toString() {
                return "average linkage (UPGMA)";
            }
        },
        CLINK {
            @Override
            public String toString() {
                return "complete linkage";
            }
        };
    }
}