/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import com.apporiented.algorithm.clustering.AverageLinkageStrategy;
import com.apporiented.algorithm.clustering.Cluster;
import com.apporiented.algorithm.clustering.ClusteringAlgorithm;
import com.apporiented.algorithm.clustering.CompleteLinkageStrategy;
import com.apporiented.algorithm.clustering.DefaultClusteringAlgorithm;
import com.apporiented.algorithm.clustering.LinkageStrategy;
import com.apporiented.algorithm.clustering.SingleLinkageStrategy;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructurePairAligner;
import org.biojava.nbio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.nbio.structure.io.PDBFileReader;

/**
 * <p>ClusterStructures class.</p>
 *
 * @author Jacob M. Litman
 */
public class ClusterStructures {

    private static final Logger logger = Logger.getLogger(ClusterStructures.class.getName());
    private final AlgorithmFunctions utils;
    private File[] files;
    private Structure[] structureCache;
    private ClusterDistanceFunction distFunction = ClusterDistanceFunction.RMSD;
    private ClustAlg algorithm = ClustAlg.AV_LINK;
    private int numClusters = 0; // Over-rides rmsdCutoff if > 0.
    private int cacheSize = 1000;
    private int cacheStart = 0; // First structure to be cached.

    private int nFiles;
    private double rmsdCutoff = 1.0;
    private boolean copyFiles = true;
    private boolean parallel = true;
    private final Path pwdPath;
    private String outputPrefix = "ffx_cluster_";
    private File[] outputDirectories;
    private Path[] outputPaths;

    /**
     * Default constructor for ClusterStructures
     *
     * @param utils AlgorithmFunctions object to use
     */
    public ClusterStructures(AlgorithmFunctions utils) {
        if (utils != null) {
            this.utils = utils;
        } else {
            this.utils = new AlgorithmUtils();
        }
        pwdPath = generatePath(new File(""));
    }

    /**
     * Constructor, including files, for ClusterStructures
     *
     * @param utils AlgorithmFunctions object to use
     * @param files Files to cluster
     */
    public ClusterStructures(AlgorithmFunctions utils, File[] files) {
        this(utils);
        this.files = files;
        nFiles = files.length;
    }

    /**
     * <p>Setter for the field <code>files</code>.</p>
     *
     * @param files an array of {@link java.io.File} objects.
     */
    public void setFiles(File[] files) {
        this.files = files;
        nFiles = files.length;
    }

    /**
     * <p>Setter for the field <code>numClusters</code>.</p>
     *
     * @param numClusters a int.
     */
    public void setNumClusters(int numClusters) {
        this.numClusters = numClusters;
    }

    /**
     * <p>setDistanceFunction.</p>
     *
     * @param distFunction a {@link ffx.algorithms.ClusterStructures.ClusterDistanceFunction} object.
     */
    public void setDistanceFunction(ClusterDistanceFunction distFunction) {
        this.distFunction = distFunction;
    }

    /**
     * <p>Setter for the field <code>algorithm</code>.</p>
     *
     * @param algorithm a {@link ffx.algorithms.ClusterStructures.ClustAlg} object.
     */
    public void setAlgorithm(ClustAlg algorithm) {
        this.algorithm = algorithm;
    }

    /**
     * <p>setOutputDirectoryPrefix.</p>
     *
     * @param prefix a {@link java.lang.String} object.
     */
    public void setOutputDirectoryPrefix(String prefix) {
        this.outputPrefix = prefix;
    }

    /**
     * <p>Setter for the field <code>copyFiles</code>.</p>
     *
     * @param copyFiles a boolean.
     */
    public void setCopyFiles(boolean copyFiles) {
        this.copyFiles = copyFiles;
    }

    /**
     * <p>setClusterParallel.</p>
     *
     * @param parallel a boolean.
     */
    public void setClusterParallel(boolean parallel) {
        this.parallel = parallel;
    }

    /**
     * <p>Setter for the field <code>cacheSize</code>.</p>
     *
     * @param cacheSize a int.
     */
    public void setCacheSize(int cacheSize) {
        this.cacheSize = cacheSize;
    }

    /**
     * <p>Setter for the field <code>rmsdCutoff</code>.</p>
     *
     * @param rmsdCutoff a double.
     */
    public void setRmsdCutoff(double rmsdCutoff) {
        this.rmsdCutoff = rmsdCutoff;
    }

    /**
     * Generate directories to output cluster info
     *
     * @param nClusters Number of directories to generate.
     */
    private void generateOutputDirectories(int nClusters) {
        File outDir = new File(String.format("%s%d", outputPrefix, 1));
        Path relPath = pwdPath;
        outputDirectories = new File[nClusters];
        outputPaths = new Path[nClusters];
        if (outDir.exists()) {
            for (int i = 2; i < 1000; i++) {
                outDir = new File(String.format("ffx_clusters_%d", i));
                if (!outDir.exists()) {
                    outDir.mkdir();
                    Path outPath = generatePath(outDir);
                    relPath = pwdPath.relativize(outPath);
                    break;
                }
            }
            if (outDir.exists()) {
                logger.severe(" Could not make output directories for clustering: all directory names taken.");
            }
        }
        for (int i = 1; i <= nClusters; i++) {
            String namei = String.format("%s%s_%d", relPath.toString(), outputPrefix, i);
            File diri = new File(namei);
            diri.mkdirs();
            outputDirectories[i - 1] = diri;
            outputPaths[i - 1] = pwdPath.relativize(generatePath(diri));
        }
    }

    /**
     * Gets the specified structure, either from the array of stored Structures,
     * or by reading from disk (depending on array size).
     *
     * @param index  Structure index
     * @param reader PDB file reader to use
     * @return the corresponding Structure
     * @throws IOException If the PDBFileReader encounters an error
     */
    private Structure accessStructure(int index, PDBFileReader reader) throws IOException {
        if (index < cacheStart) {
            return reader.getStructure(files[index]);
        } else {
            return structureCache[index - cacheStart];
        }
    }

    /**
     * Main execution method for ClusterStructures.
     *
     * @return A list of the final clusters
     */
    public List<Cluster> cluster() {
        cacheStart = nFiles - cacheSize;
        List<Cluster> clusters;
        if (parallel) {
            clusters = clusterParallel();
        } else {
            clusters = clusterSequential();
        }

        int nClusters = clusters.size();
        generateOutputDirectories(nClusters);
        for (int i = 0; i < nClusters; i++) {
            Cluster cluster = clusters.get(i);
            String filename = String.format("%sffx_cluster_%d_summary", outputPaths[i].toString(), i);
            File summaryFile = new File(filename.concat(".txt"));

            for (int j = 1; j < 1000; j++) {
                if (summaryFile.exists()) {
                    summaryFile = new File(String.format("%s_%d.txt", filename, j));
                } else {
                    break;
                }
            }
            if (summaryFile.exists()) {
                logger.warning(String.format(" Could not make valid cluster summary "
                        + "file name for %s", filename));
                continue;
            }

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(summaryFile))) {
                bw.write(String.format("PDB files for cluster %d", i));
                bw.newLine();
                bw.newLine();

                List<Cluster> childClusters = getSubclusters(cluster, 0);
                for (Cluster child : childClusters) {
                    int index = Integer.parseInt(child.getName());
                    bw.write(files[index].getName());
                    bw.newLine();
                }

                if (copyFiles) {
                    for (Cluster child : childClusters) {
                        int index = Integer.parseInt(child.getName());
                        File childFile = files[index];
                        filename = outputPaths[i].toString().concat(childFile.getName());
                        File writeTo = new File(filename);
                        try {
                            FileUtils.copyFile(childFile, writeTo, false);
                        } catch (IOException ex) {
                            logger.warning(String.format(" Could not copy file %s", filename));
                        }
                    }
                }
            } catch (IOException ex) {
                logger.warning(String.format(" Failed to properly write summary "
                        + "file for cluster %d", i));
            }
        }
        return clusters;
    }

    /**
     * Performs clustering in parallel.
     *
     * @return Final clusters.
     */
    private List<Cluster> clusterParallel() {
        String[] names = new String[nFiles];
        double[][] rmsdDistances = new double[nFiles][nFiles];

        for (int i = 0; i < nFiles; i++) {
            rmsdDistances[i][i] = 0.0; // Ensure the diagonal is filled.
            names[i] = String.format("%d", i);
        }
        /* This stuff should go in the ParallelRegion start() method.
         nThreads = ParallelTeam.getDefaultThreadCount();
         fileReaders = new PDBFileReader[nThreads];
         for (int i = 0; i < nThreads; i++) {
         fileReaders[i] = new PDBFileReader();
         }*/
        return null;
    }

    /**
     * Performs clustering
     *
     * @return Final clusters.
     */
    private List<Cluster> clusterSequential() {
        String[] names = new String[nFiles];
        double[][] rmsdDistances = new double[nFiles][nFiles];
        PDBFileReader fileReader = new PDBFileReader();
        LinkageStrategy ls;
        switch (algorithm) {
            case CLINK:
                ls = new CompleteLinkageStrategy();
                break;
            case SLINK:
                ls = new SingleLinkageStrategy();
                break;
            case AV_LINK:
            default:
                ls = new AverageLinkageStrategy();
                break;
        }

        for (int i = 0; i < nFiles; i++) {
            rmsdDistances[i][i] = 0.0; // Ensure the diagonal is filled.
            names[i] = String.format("%d", i);
            if (i >= cacheStart) {
                try {
                    structureCache[i - cacheStart] = fileReader.getStructure(files[i]);
                } catch (IOException ex) {
                    logger.severe(String.format(" Error in reading file %s: %s",
                            files[i].getName(), ex.toString()));
                }
            }
        }

        StructurePairAligner aligner = new StructurePairAligner();
        for (int i = 0; i < nFiles; i++) {
            Structure structI = null;
            try {
                structI = accessStructure(i, fileReader);
            } catch (IOException ex) {
                logger.severe(String.format(" Error in reading file %s: %s",
                        files[i].getName(), ex.toString()));
            }
            for (int j = i; j < nFiles; j++) {
                Structure structJ = null;
                try {
                    structJ = accessStructure(j, fileReader);
                } catch (IOException ex) {
                    logger.severe(String.format(" Error in reading file %s: %s",
                            files[j].getName(), ex.toString()));
                }

                try {
                    aligner.align(structI, structJ);
                } catch (StructureException ex) {
                    logger.severe(String.format(" Exception aligning structures "
                            + "%d and %d: %s", i, j, ex.toString()));
                }
                AlternativeAlignment[] alignments = aligner.getAlignments();
                double minRMSD = alignments[0].getRmsd();
                for (int k = 1; k < alignments.length; k++) {
                    double rmsdK = alignments[k].getRmsd();
                    minRMSD = rmsdK < minRMSD ? rmsdK : minRMSD;
                }
                rmsdDistances[i][j] = minRMSD;
                rmsdDistances[j][i] = minRMSD;
            }
        }

        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(rmsdDistances, names, ls);
        List<Cluster> subClusters;
        int nClusters = 1;

        if (numClusters > 0) {
            subClusters = new ArrayList<>(Arrays.asList(cluster));

            while (nClusters < numClusters) {
                double maxDist = subClusters.get(0).getDistanceValue();
                Cluster maxCluster = subClusters.get(0);
                for (Cluster subcluster : subClusters) {
                    double dist = subcluster.getDistanceValue();
                    if (dist > maxDist) {
                        maxDist = dist;
                        maxCluster = subcluster;
                    }
                }

                List<Cluster> newClusters = maxCluster.getChildren();
                nClusters += (newClusters.size() - 1);
                subClusters.addAll(newClusters);
                subClusters.remove(maxCluster);
            }
            logger.severe(" Num clusters not implemented yet.");
        } else {
            subClusters = getSubclusters(cluster, rmsdCutoff);
            nClusters = subClusters.size();
        }

        assert nClusters == subClusters.size() : " nClusters != subClusters.size()";

        return subClusters;
    }

    /**
     * Recursively returns all subclusters in this Cluster with a distance
     * greater than RMSD cutoff.
     *
     * @param cluster
     * @param cutoff
     * @return A List of Cluster instances.
     */
    private List<Cluster> getSubclusters(Cluster cluster, double cutoff) {
        if (cluster.getDistanceValue() < cutoff || cluster.isLeaf()) {
            return Arrays.asList(cluster);
        } else {
            List<Cluster> clusters = new ArrayList<>();
            for (Cluster subcluster : cluster.getChildren()) {
                clusters.addAll(getSubclusters(subcluster, cutoff));
            }
            return clusters;
        }
    }

    /**
     * Recursively returns all leaf clusters under this cluster.
     *
     * @param cluster Cluster to search under
     * @return All child leaves
     */
    private List<Cluster> getLeafClusters(Cluster cluster) {
        if (cluster.isLeaf()) {
            return Arrays.asList(cluster);
        } else {
            List<Cluster> clusters = new ArrayList<>();
            for (Cluster subcluster : cluster.getChildren()) {
                clusters.addAll(getLeafClusters(subcluster));
            }
            return clusters;
        }
    }

    /**
     * Utility method which attempts to generate a file Path using the canonical
     * path string, else uses the absolute path.
     *
     * @param file To find path of
     * @return Canonical or absolute path.
     */
    public static Path generatePath(File file) {
        Path path;
        try {
            path = Paths.get(file.getCanonicalPath());
        } catch (IOException ex) {
            path = Paths.get(file.getAbsolutePath());
        }
        return path;
    }

    public enum ClusterDistanceFunction {

        RMSD, CA_RMSD, DIHEDRALS, BACKBONE_DIHEDRALS
    }

    public enum ClustAlg {

        /**
         * All algorithms start with each point a cluster, and then join the
         * closest clusters together until everything is one cluster. SLINK is
         * Single Linkage; cluster-cluster distance is defined by the nearest
         * two points. This is vulnerable to chaining; two clusters might be
         * joined by a handful of intermediate points. CLINK is Complete
         * Linkage; CLINK uses the greatest distance between points in two
         * clusters. AV_LINK (average link) is the UPGMA (Unweighted Pair Group
         * Method with Arithmetic Mean) function, which takes the mean distance
         * between points in a cluster.
         * <p>
         * Makes me wonder if there's a WPGMA algorithm which does weight one
         * way or the other, or perhaps a RPGMA RMSD-like algorithm.
         */
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

    // Copyright license for hierarchical-clustering-java
    /**
     * *****************************************************************************
     * Copyright 2013 Lars Behnke
     *
     * Licensed under the Apache License, Version 2.0 (the "License"); you may
     * not use this file except in compliance with the License. You may obtain a
     * copy of the License at
     *
     * http://www.apache.org/licenses/LICENSE-2.0
     *
     * Unless required by applicable law or agreed to in writing, software
     * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
     * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
     * License for the specific language governing permissions and limitations
     * under the License.
     * ****************************************************************************
     */
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
}
