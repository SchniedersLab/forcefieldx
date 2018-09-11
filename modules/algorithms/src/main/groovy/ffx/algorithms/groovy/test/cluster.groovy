/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

package ffx.algorithms.groovy.test;

import groovy.cli.picocli.CliBuilder

import ffx.potential.parsers.CoordinateFileFilter
import ffx.potential.parsers.PDBFileFilter
import static ffx.algorithms.ClusterStructures.ClustAlg.AV_LINK
import static ffx.algorithms.ClusterStructures.ClustAlg.CLINK
import static ffx.algorithms.ClusterStructures.ClustAlg.SLINK
import static ffx.algorithms.ClusterStructures.ClusterDistanceFunction.BACKBONE_DIHEDRALS
import static ffx.algorithms.ClusterStructures.ClusterDistanceFunction.CA_RMSD
import static ffx.algorithms.ClusterStructures.ClusterDistanceFunction.DIHEDRALS
import static ffx.algorithms.ClusterStructures.ClusterDistanceFunction.RMSD

boolean copyFiles = true;
boolean parallel = true;
String[] sourceFileNames;
File[] clusterFiles;
ffx.algorithms.AlgorithmFunctions utils;
String outputDirectoryName = "ffx_cluster_";
ffx.algorithms.ClusterStructures.ClustAlg algorithm = ffx.algorithms.ClusterStructures.ClustAlg.AV_LINK;
ffx.algorithms.ClusterStructures.ClusterDistanceFunction distFunction = RMSD;
int numClusters = 0;
int cacheSize = 1000;
double rmsdCutoff = 1.0;

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc test.cluster [options] <pdbfilename>');
cli.h(longOpt: 'help', 'Print this help message.');
cli.w(longOpt: 'write', args: 1, argName: 'true', 'Write copies of PDB files to cluster directories.');
cli.o(longOpt: 'outputDirectories', argName: 'ffx_cluster_', 'Prefix of cluster output directories (followed by number)');
cli.d(longOpt: 'distanceFunction', argName: '1', 'Cluster based on all-atom RMSD (1), CA RMSD (2), all-torsion RMSD (3), or backbone torsion RMSD (4)');
cli.a(longOpt: 'algorithm', argName: 'average', 'Make clusters using single, average, or complete linkage (SLINK, UPGMA, CLINK)');
cli.r(longOpt: 'rmsdCutoff', argName: '1.0', 'RMSD at which to separate clusters.');
cli.n(longOpt: 'numClusters', argName: '0', 'Number of clusters to generate; over-rides distance.');
cli.p(longOpt: 'parallel', argName: 'true', 'Clusters in parallel.');
cli.c(longOpt: 'cacheSize', argName: '1000', 'Number of structures to retain in memory.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    logger.info(" The filename can be a directory of files to cluster, a text file of files to cluster, or a file to cluster.");
    return cli.usage();
}

try {
    utils = getAlgorithmUtils();
} catch (MissingMethodException ex) {
    utils = new ffx.algorithms.AlgorithmUtils();
}

if (options.w) {
    copyFiles = Boolean.parseBoolean(options.w);
}

if (options.o) {
    outputDirectoryName = options.o;
}

if (options.c) {
    cacheSize = Integer.parseInt(options.c);
}

if (options.r) {
    rmsdCutoff = Double.parseDouble(options.r);
}

if (options.d) {
    int value = Integer.parseInt(options.d);
    switch (value) {
        case 1:
            distFunction = RMSD;
            break;
        case 2:
            distFunction = CA_RMSD;
            break;
        case 3:
            distFunction = DIHEDRALS;
            break;
        case 4:
            distFunction = BACKBONE_DIHEDRALS;
            break;
        default:
            logger.warning(String.format(" Invalid selection %d for distance function; must be 1, 2, 3, or 4. Defaulting to average linkage.", value));
            distFunction = RMSD;
            break;
    }
}

if (options.a) {
    String algo = options.a.toUpperCase();
    switch (algo) {
        case "SLINK":
        case "SINGLE":
        case "SINGLE_LINKAGE":
            algorithm = SLINK;
            break;
        case "AV_LINK":
        case "AVERAGE":
        case "UPGMA":
        case "AVERAGE_LINKAGE":
            algorithm = AV_LINK;
            break;
        case "CLINK":
        case "COMPLETE":
        case "COMPLETE_LINKAGE":
            algorithm = CLINK;
            break;
        default:
            logger.warning(String.format(" Invalid algorithm selection %s; must be SLINK, AV_LINK, or CLINK", algo));
            break;
    }
}

if (options.n) {
    numClusters = Integer.parseInt(options.n);
}

if (options.p) {
    parallel = Boolean.parseBoolean(options.p);
}

if (numClusters < 1 && rmsdCutoff <= 0) {
    logger.severe(" Invalid options: number of distances and RMSD cutoff are both 0 or less.");
}

PDBFileFilter pdbFilter = new PDBFileFilter();
CoordinateFileFilter coordinateFilter = new CoordinateFileFilter();
List<File> modelFiles = new ArrayList<>();
int numTempFiles = 0;
File tempDirectory = null;
String tempDirName = "";

// Temporary, simplified system for loading files: one directory of PDBs.
String dirname = arguments.get(i);
File dir = new File(dirname);
if (!dir.exists() || !dir.isDirectory()) {
    logger.severe(" Argument was not a directory (preliminary loading system only takes a directory)");
}
File[] allFiles = dir.listFiles();
for (File file : allFiles) {
    if (pdbFilter.acceptDeep(file)) {
        modelFiles.add(file);
    } else {
        logger.info(String.format(" File %s not accepted.", file.getName()));
    }
}

/**
 * Essentially: it loops through direct file arguments, the files listed by non-
 * coordinate files, and files in directories. If they are not PDB format, a temporary
 * file is written to a temporary directory. At the end of the try-catch 
 * block, the temporary directory (if it was ever created) is deleted.
 *
 * At each point where a non-PDB coordinate file is detected, it checks to see 
 * if the temporary directory has been created (creating it if necessary), opens 
 * the file, saves it as a PDB in the temporary directory, adds that temporary 
 * file to modelFiles, and closes the MolecularAssembly[].
 *
 * I'll test this out later.
 */
/*try {
    for (int i = 0; i < arguments.size(); i++) {
        String filename = arguments.get(i);
        File file = new File(filename);
        if (!file.exists) {
            logger.warning(String.format(" File %s does not exist", filename));
        } else if (file.isDirectory()) {
            File[] allFiles = file.listFiles();
            for (File modelFile : allFiles) {
                if (coordinateFilter.accept(modelFile)) {
                    if (pdbFilter.acceptDeep(modelFile)) {
                        modelFiles.add(modelFile);
                    } else {
                        if (tempDirectory == null) {
                            String tempDirectoryName = "ffx-tempfiles";
                            tempDirectory = new File(tempDirectoryName);
                            for (int j = 1; j < 1000; j++) {
                                if (tempDirectory.exists()) {
                                    tempDirectory = new File(String.format("%s-%d", tempDirectoryName, j));
                                } else {
                                    tempDirectory.mkdir();
                                    break;
                                }
                            }
                            if (tempDirectory.exists()) {
                                logger.severe(" Unable to make temporary directory to convert non-PDB files to PDB format");
                            }
                            tempDirName = tempDirectory.getName();
                        }
                        MolecularAssembly[] assemblies = utils.open(modelFile);
                        String newFileName = FilenameUtils.removeExtension(modelFile.getName()).concat("_tmp.pdb");
                        File tempFile = File.createTempFile(String.format("%s%c%s", tempDirName, File.separatorChar, newFileName));
                        utils.saveAsPDB(assemblies, tempFile);
                        modelFiles.add(tempFile);
                        utils.closeAll(assemblies);
                    }
                }
            }
        } else if (coordinateFilter.accept(file)) {
            if (pdbFilter.acceptDeep(file)) {
                modelFiles.add(file);
            } else {
                if (tempDirectory == null) {
                    String tempDirectoryName = "ffx-tempfiles";
                    tempDirectory = new File(tempDirectoryName);
                    for (int j = 1; j < 1000; j++) {
                        if (tempDirectory.exists()) {
                            tempDirectory = new File(String.format("%s-%d", tempDirectoryName, j));
                        } else {
                            tempDirectory.mkdir();
                            break;
                        }
                    }
                    if (tempDirectory.exists()) {
                        logger.severe(" Unable to make temporary directory to convert non-PDB files to PDB format");
                    }
                    tempDirName = tempDirectory.getName();
                }
                MolecularAssembly[] assemblies = utils.open(file);
                String newFileName = FilenameUtils.removeExtension(file.getName()).concat("_tmp.pdb");
                File tempFile = File.createTempFile(String.format("%s%c%s", tempDirName, File.separatorChar, newFileName));
                utils.saveAsPDB(assemblies, tempFile);
                modelFiles.add(tempFile);
                utils.closeAll(assemblies);
            }
        } else {
            BufferedReader sfReader;
            try {
                sfReader = new BufferedReader(new FileReader(file));
                String line = sfReader.readLine();
                while (line != null) {
                    File model = new File(line);
                    if (!model.exists()) {
                        logger.warning(String.format(" File %s does not exist", line));
                    } else if (coordinateFilter.accept(model)) {
                        if (pdbFilter.acceptDeep(model)) {
                            modelFiles.add(model);
                        } else {
                            if (tempDirectory == null) {
                                String tempDirectoryName = "ffx-tempfiles";
                                tempDirectory = new File(tempDirectoryName);
                                for (int j = 1; j < 1000; j++) {
                                    if (tempDirectory.exists()) {
                                        tempDirectory = new File(String.format("%s-%d", tempDirectoryName, j));
                                    } else {
                                        tempDirectory.mkdir();
                                        break;
                                    }
                                }
                                if (tempDirectory.exists()) {
                                    logger.severe(" Unable to make temporary directory to convert non-PDB files to PDB format");
                                }
                                tempDirName = tempDirectory.getName();
                            }
                            MolecularAssembly[] assemblies = utils.open(model);
                            String newFileName = FilenameUtils.removeExtension(model.getName()).concat("_tmp.pdb");
                            File tempFile = File.createTempFile(String.format("%s%c%s", tempDirName, File.separatorChar, newFileName));
                            utils.saveAsPDB(assemblies, tempFile);
                            modelFiles.add(tempFile);
                            utils.closeAll(assemblies);
                        }
                    } else {
                        logger.info(String.format(" File %s not recognized as valid PDB file.", line));
                    }
                    line = sfReader.readLine();
                }
            } catch (IOException ex) {
                sfReader.close();
                logger.warning(String.format(" File reader for input file %s shut down: exception %s", filename, ex.toString()));
            }
        }
    }
    for (int i = 0; i < arguments.size(); i++) {
        String filename = arguments.get(i);
        File file = new File(filename);
        if (!file.exists) {
            logger.warning(String.format(" File %s does not exist", filename));
        } else if (file.isDirectory()) {
            File[] allFiles = file.listFiles();
        } else if (pdbFilter.acceptDeep(file)) {
                
        } else {
            BufferedReader sfReader;
            try {
                sfReader = new BufferedReader(new FileReader(file));
                String line = sfReader.readLine();
                while (line != null) {
                    File model = new File(line);
                    if (!model.exists()) {
                        logger.warning(String.format(" File %s does not exist", line));
                    } else if (filter.acceptDeep(model)) {
                        clusterFiles.add(model);
                    } else {
                        logger.info(String.format(" File %s not recognized as valid PDB file.", line));
                    }
                    line = sfReader.readLine();
                }
            } catch (IOException ex) {
                sfReader.close();
                logger.warning(String.format(" File reader for input file %s shut down: exception %s", filename, ex.toString()));
            }
        }
    }
} catch (Exception ex) {
    if (tempDirectory != null) {
        //FileUtils.deleteDirectory(tempDirectory);
    }
    logger.severe(String.format(" Exception %s in reading files to cluster", ex.toString()));
}
if (tempDirectory != null) {
    //FileUtils.deleteDirectory(tempDirectory);
}*/

ffx.algorithms.ClusterStructures clusterer = new ffx.algorithms.ClusterStructures(utils, modelFiles);
clusterer.setAlgorithm(algorithm);
clusterer.setNumClusters(numClusters);
clusterer.setOutputDirectoryPrefix(outputDirectoryName);
clusterer.setDistanceFunction(distFunction);
clusterer.setCopyFiles(copyFiles);
clusterer.setClusterParallel(parallel);
clusterer.setCacheSize(cacheSize);
clusterer.setRmsdCutoff(rmsdCutoff);
clusterer.cluster();
