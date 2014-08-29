/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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

package ffx

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Java Imports
import java.nio.file.Paths;
import java.nio.file.Path;

// Force Field X X-ray Imports
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// Force Field X Minimize imports
import ffx.algorithms.Minimize;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.CoordinateFileFilter;

// Force Field X real-space Imports
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// Other Force Field X imports.
import ffx.utilities.DoubleIndexPair;

// RMS gradient per atom convergence criteria
double eps = -1.0;

// maximum number of refinement cycles
int maxiter = 1000;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;

String suffix = "_rsc";
boolean addListFile = true;
int refineMode = 0;
List diffractionfiles = new ArrayList();
List mapFiles = new ArrayList();
ArrayList<File> modelFiles = new ArrayList<>();
double acceptThreshold = 0.0;
int numToAccept = 10;
File resultsFile;
File[] rescoredFiles;
String rescoreFileName = "ffx_rescore.txt";
boolean versioned = true;
File directory;
boolean useDirectory = false;
Path directoryPath = null;
boolean printModels = false;
boolean includeRejected = true;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rescore [options] <modelfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.m(longOpt: 'minimize', args: 1, argName:'0', 'Do not minimize coordinates before rescore (0), minimize (1), minimize with X-ray target (2), or minimize with real space target (3).');
cli.e(longOpt:'eps', args:1, argName:'-1.0', 'RMS gradient convergence criteria (negative: automatically determine based on refinement type)');
cli.xd(longOpt:'xrayData', args:3, valueSeparator:',', argName:'data,1.0,false', 'specify input xray data filename or filenames (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment.');
cli.rd(longOpt:'realspaceData', args:2, valueSeparator:',', argName:'data,1.0', 'specify input map data filename or filenames (or simply provide the datafilename argument after the PDB file) and weight applied to the data.');
cli.i(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed refinement iterations');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');
cli.s(longOpt:'suffix', args:1, argName:'_rsc', 'output suffix');
cli.r(longOpt:'resultsFile', args:1, argName:'ffx_rescore.txt', 'File to output accepted models and scores to.');
cli.pm(longOpt:'printModels', args:1, argName:'false', 'Prints models to rescore file if true; else only prints file names');
cli.t(longOpt:'acceptThreshold', args:1, argName:'0.0', 'Accepts models within this distance (kcal/mol) of best structure; over-rides number to accept if >0.0');
cli.n(longOpt:'numberToAccept', args:1, argName:'10', 'Accepts this many of the best structures.');
cli.dr(longOpt:'directory', args:1, argName:'directory', 'If set, rescore files are written to this directory: is made if it does not exist.');
cli.ir(longOpt:'includeRejected', args:1, argName:'true', 'If set, includes names of files not accepted.');
//cli.dt(longOpt:'dataType', args:1, argName:'false', 'If true, refine modes 2 and 3 will parse data arguments as a line-separated list of {pdbfile,datafile} to match data files to PDBs. Preferred meth.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    logger.info(" modelfilename can either be a file with a list of files to rescore, or a directory/folder of files to rescore.");
    logger.info(" To minimize multiple structures with different data/map files, enter only a directory containing models and identically-named data/map files (ex: abc.pdb and abc.mtz).\n");
    return cli.usage();
}
/*if (options.dt) {
    matchDataFiles = Boolean.parseBoolean(options.dt);
}*/
if (options.m) {
    refineMode = Integer.parseInt(options.m);
    switch (refineMode) {
    case 0:
    case 1:
        if (arguments.size() != 1) {
            logger.warning(String.format(" Improper number of arguments (must be 1) for -m mode %d", refineMode));
            return cli.usage();
        }
        /*if (matchDataFiles) {
            logger.warning(String.format(" dataType argument invalid for refine modes 0 or 1."));
            return cli.usage();
        }*/
        break;
    case 2:
    case 3:
        if (arguments.size() > 2) {
            logger.warning(String.format(" Improper number of arguments (must be 1 or 2) for -m mode %d", refineMode));
            return cli.usage();
        }
        break;
    default:
        logger.warning(" Improper use of -m flag: must be 0, 1, 2, or 3");
        return cli.usage();
    }
}

Path pwdPath = Paths.get(new File("").getCanonicalPath());

if (options.s) {
    suffix = options.s;
}

if (options.t) {
    acceptThreshold = Double.parseDouble(options.t);
    if (acceptThreshold < 0.0) {
        logger.warning(" Accept threshold should be > 0. Defaulting to numToAccept.");
    }
}

if (options.pm) {
    printModels = Boolean.parseBoolean(options.pm);
}

if (options.n) {
    numToAccept = Integer.parseInt(options.n);
}

if (options.ir) {
    includeRejected = Boolean.parseBoolean(options.ir);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.i) {
    maxiter = Integer.parseInt(options.i);
}

if (options.r) {
    rescoreFileName = options.r;
}
resultsFile = new File(rescoreFileName);
if (resultsFile.exists()) {
    for (int i = 2; i <= 1000; i++) {
        resultsFile = new File(rescoreFileName + "_" + i);
        if (!resultsFile.exists()) {
            break;
        }
    }
    if (resultsFile.exists()) {
        logger.warning(String.format(" Versioning failed: appending rescore file to end of file %s", resultsFile.getName()));
        versioned = false;
    }
}

if (/*!matchDataFiles && */arguments.size() > 1) {
    if (refineMode == 2) {
        DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false);
        diffractionfiles.add(diffractionfile);
    } else if (refineMode == 3) {
        RealSpaceFile realSpaceFile = new RealSpaceFile(arguments.get(1), 1.0);
        mapFiles.add(realSpaceFile);
    }
}
if (/*!matchDataFiles && */options.xd) {
    if (options.rd) { // At present, double-minimization not allowed.
        logger.warning(String.format(" Only one refinement mode (%d) will be used.", refineMode));
        switch (refineMode) {
        case 0:
            logger.info(" No minimization will be applied.");
            break;
        case 1:
            logger.info(" Force field minimization will be applied.");
            break;
        case 2:
            logger.info(" X-ray refinement will be applied.");
            break;
        case 3:
            logger.info(" Real space refinement will be applied.");
            break;
        default:
            throw new IllegalArgumentException(" Program should already have returned usage due to invalid -m flag.");
        }
    }
    switch (refineMode) {
    case 0:
    case 1:
        logger.warning(" Data file is only used for X-ray and real-space refinements.");
        break;
    case 2:
        for (int i = 0; i < options.xds.size(); i += 3) {
            double wA = Double.parseDouble(options.xds[i+1]);
            boolean neutron = Boolean.parseBoolean(options.xds[i+2]);
            DiffractionFile diffractionFile = new DiffractionFile(options.xds[i], wA, neutron);
            diffractionfiles.add(diffractionFile);
        }
        break;
    case 3:
        logger.warning("The xrayData flag is only for use with x-ray refinement.");
        break;
    }
}

if (/*!matchDataFiles && */options.rd) {
    switch (refineMode) {
    case 0:
    case 1:
        logger.warning(" Data file is only used for X-ray and real-space refinements.");
        break;
    case 2: 
        logger.warning("The realspaceData flag is only for use with realspace refinement.");
        break;
    case 3:
        for (int i = 0; i < options.rds.size(); i += 2) {
            double wA = Double.parseDouble(options.rds[i+1]);
            RealSpaceFile realSpaceFile = new RealSpaceFile(options.rds[i], wA);
            mapFiles.add(realSpaceFile);
        }
        break;
    }
}

try {
    File sourceFile = new File(arguments.get(0));
    CoordinateFileFilter filter = new CoordinateFileFilter();
    if (!sourceFile.exists()) {
        throw new FileNotFoundException("Source file or directory does not exist");
    } else if (sourceFile.isFile()) {
        BufferedReader sfReader = new BufferedReader(new FileReader(sourceFile));
        try {
            String line = sfReader.readLine();
            while (line != null) {
                if (line == null) {
                    break;
                }
                File model = new File(line);
                if (!model.exists()) {
                    logger.warning(String.format(" File %s does not exist", line));
                } else if (filter.acceptDeep(model)) {
                    modelFiles.add(model);
                } else {
                    logger.info(String.format(" File %s not recognized as valid PDB, XYZ, INT, or ARC file.", line));
                }
                line = sfReader.readLine();
            }
        } finally {
            sfReader.close();
        }
        if (modelFiles.isEmpty()) {
            throw new IOException(" Source file contained no valid file names.");
        }
    } else if (sourceFile.isDirectory()) {
        File[] allFiles = sourceFile.listFiles();
        for (File file : allFiles) {
            if (filter.acceptDeep(file)) {
                modelFiles.add(file);
            }
        }
        if (modelFiles.isEmpty()) {
            throw new IOException(" Source directory contained no files.");
        }
    }
} catch (IOException ex) {
    logger.severe(String.format(" Exception in generating file list: %s \n %s", ex.toString(), ex.printStackTrace()));
}

rescoredFiles = new File[modelFiles.size()];

if (options.dr) {
    directory = new File(options.dr);
    if (directory.exists()) {
        if (directory.isFile()) {
            logger.severe(String.format(" FFX shutting down: Pathname %s provided for directory points to a file.", directory.toString()));
        }
    } else {
        boolean success = directory.mkdirs();
        if (!success) {
            logger.severe(String.format(" FFX shutting down: could not create directory %s", directory.toString()));
        }
    }
    rescoreFileName = FilenameUtils.concat(directory.getPath(), rescoreFileName);
    useDirectory = true;
    directoryPath = pwdPath.relativize(Paths.get(directory.getCanonicalPath()));
}

/*if (matchDataFiles) {
    
}*/

int totalFiles = modelFiles.size();
int counter = 1;
DoubleIndexPair[] energies = new DoubleIndexPair[totalFiles];
// Energies contains the index of the file in modelFiles and printedFiles, plus the calculated energies.

for (int i = 0; i < totalFiles; i++) {
    File file = modelFiles.get(i);
    Path filepath = Paths.get(file.getCanonicalPath());
    String filename = pwdPath.relativize(filepath).toString();
    logger.info(String.format("\n\n Rescoring file %d of %d: %s", counter++, totalFiles, filename));
    try {
        systems = open(filename);
        switch (refineMode) {
        case 0:
            rescoredFiles[i] = file;
            break;
        case 1:
            logger.info("\n Running minimize on " + filename);
            logger.info(" RMS gradient convergence criteria: " + eps);
            energy();
            
            // Do the minimization
            minimize(eps);
            
            String ext = FilenameUtils.getExtension(filename);
            ext = ".".concat(ext);
            
            if (useDirectory) {
                filename = FilenameUtils.getBaseName(filename);
                filename = FilenameUtils.concat(directoryPath.toString(), filename);
                filename = filename.concat(suffix).concat(ext);
            } else {
                filename = FilenameUtils.removeExtension(filename);
                filename = filename.concat(suffix).concat(ext);
            }
            rescoredFiles[i] = new File(filename);
            if (ext.toUpperCase().contains("XYZ")) {
                saveAsXYZ(new File(rescoredFiles[i]));
            } else {
                saveAsPDB(systems, new File(rescoredFiles[i]));
            }
            break;
        case 2:
            logger.info("\n Running x-ray minimize on " + filename);
            
            DiffractionFile diffractionfile = null;
            if (diffractionfiles.size() == 0) {
                diffractionfile = new DiffractionFile(systems, 1.0, false);
                diffractionfiles.add(diffractionfile);
            }
            DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));
            
            diffractiondata.scaleBulkFit();
            diffractiondata.printStats();
            energy();
        
            RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode);
            if (eps < 0.0) {
                eps = refinementMinimize.getEps();
            }
            logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
            refinementMinimize.minimize(eps, maxiter);
            diffractiondata.scaleBulkFit();
            diffractiondata.printStats();
            
            String ext = FilenameUtils.getExtension(filename);
            ext = ".".concat(ext);
        
            if (useDirectory) {
                filename = FilenameUtils.getBaseName(filename);
                filename = FilenameUtils.concat(directoryPath.toString(), filename);
                filename = filename.concat(suffix);
                diffractiondata.writeData(filename + ".mtz");
                filename = filename.concat(ext);
                diffractiondata.writeModel(filename);
            } else {
                filename = FilenameUtils.removeExtension(filename);
                filename = filename.concat(suffix);
                diffractiondata.writeData(filename + ".mtz");
                filename = filename.concat(ext);
                diffractiondata.writeModel(filename);
            }
            rescoredFiles[i] = new File(filename);
            if (diffractionfile != null) {
                try {
                    diffractionfiles.remove(diffractionfile);
                } catch (UnsupportedOperationException ex) {
                    // This should never occur, because diffractionfiles should be of a List type supporting remove(object).
                    diffractionfiles = new ArrayList<>();
                }
            }
            break;
        case 3:
            logger.info("\n Running real-space minimize on " + filename);
            
            RealSpaceFile realspacefile = null;
            if (mapFiles.size() == 0) {
                realspacefile = new RealSpaceFile(systems);
                mapFiles.add(realspacefile);
            }
            RealSpaceData realspacedata = new RealSpaceData(systems, systems[0].getProperties(), mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
            energy();
            
            RefinementMinimize refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES);
            
            if (eps < 0.0) {
                eps = 1.0;
            }
            logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
            refinementMinimize.minimize(eps, maxiter);
            
            String ext = FilenameUtils.getExtension(filename);
            ext = ".".concat(ext);
            
            if (useDirectory) {
                filename = FilenameUtils.getBaseName(filename);
                filename = FilenameUtils.concat(directoryPath.toString(), filename);
                filename = filename.concat(suffix).concat(ext);
            } else {
                filename = FilenameUtils.removeExtension(filename);
                filename = filename.concat(suffix).concat(ext);
            }
            rescoredFiles[i] = new File(filename);
            if (ext.toUpperCase().contains("XYZ")) {
                saveAsXYZ(rescoredFiles[i]);
            } else {
                saveAsPDB(systems, rescoredFiles[i]);
            }
            
            if (realspacefile != null) {
                try {
                    mapFiles.remove(realspacefile);
                } catch (UnsupportedOperationException ex) {
                    // This should never occur, because diffractionfiles should be of a List type supporting remove(object).
                    mapFiles = new ArrayList<>();
                }
            }
            break;
        }
        double e = returnEnergy();
        energies[i] = new DoubleIndexPair(i, e);
        close();
    } catch (Exception ex) {
        logger.warning(String.format(" Exception rescoring on file %s", filename));
        logger.info(ex.toString());
    }
}
/*Arrays.sort(energies);
ArrayList<Integer> acceptedFileIndices = new ArrayList<>();

if (acceptThreshold > 0.0) {
    int numAccepted = totalFiles; // Default condition: all accepted.
    for (int i = 0; i < totalFiles; i++) {
        DoubleIndexPair currentPair = energies.get(i);
        double e = currentPair.getDoubleValue();
        if (e - minEnergy > acceptThreshold) {
            numAccepted = i; // I am pretty sure this is the correct indexing for the first non-accepted pair.
            break;
        } else {
            acceptedFileIndices.add(currentPair.getIndex());
        }
    }
} else if (numToAccept > 0) {
    numToAccept = Math.min(numToAccept, totalFiles);
    for (int i = 0; i < numToAccept; i++) {
        acceptedFileIndices.add(energies.get(i).getIndex());
    }
}*/

Arrays.sort(energies);

double minEnergy = energies[0].getDoubleValue();
int numAccepted = 0;
if (acceptThreshold > 0.0) {
    for (int i = 0; i < totalFiles; i++) {
        if (energies[i].getDoubleValue() > (minEnergy + acceptThreshold)) {
            break;
        } else {
            ++numAccepted;
        }
    }
} else {
    numAccepted = numToAccept < totalFiles ? numToAccept : totalFiles; //Minimum of numToAccept or totalFiles
}

/*for (int i = 0; i < numAccepted; i++) {
    File filei = rescoredFiles[energies[i].getIndex()];
    Path pathi = Paths.get(filei.getCanonicalPath());
    String relPath = pwdPath.relativize(pathi).toString();
    logger.info(String.format(" Accepted: %s at %f kcal/mol", relPath, energies[i].getDoubleValue()));
}
for (int i = numAccepted; i < totalFiles; i++) { // Rejected files
    File filei = rescoredFiles[energies[i].getIndex()];
    Path pathi = Paths.get(filei.getCanonicalPath());
    String relPath = pwdPath.relativize(pathi).toString();
    logger.info(String.format(" Rejected: %s at %f kcal/mol", relPath, energies[i].getDoubleValue()));
}*/

BufferedWriter bw;
try {
    bw = new BufferedWriter(new FileWriter(resultsFile, true));
    Path rscFilePath = Paths.get(resultsFile.getCanonicalPath());
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
    logger.info(" " + message);
    for (int i = 0; i < numAccepted; i++) {
        int fileIndex = energies[i].getIndex();
        File pointedFile = rescoredFiles[fileIndex];
        Path pointedPath = Paths.get(pointedFile.getCanonicalPath());
        String relPath = pwdPath.relativize(pointedPath).toString();
        double thisEnergy = energies[i].getDoubleValue();
        
        logger.info(String.format(" Accepted file %d energy %9.3f < %9.3f kcal/mol", (i+1), thisEnergy, minEnergy + acceptThreshold));
        logger.info(String.format(" %s", relPath));
        try {
            bw.write(String.format("Accepted file: %s rank %d energy %f\n", relPath, (i+1), thisEnergy));
            if (printModels) {
                bw.newLine();
                BufferedReader br;
                try {
                    br = new BufferedReader(new FileReader(pointedFile));
                    String line = br.readLine();
                    while (line != null) {
                        bw.write(line);
                        bw.newLine();
                        line = br.readLine();
                    }
                } finally {
                    if (br != null) br.close();
                }
                bw.newLine();
            }
        } catch (IOException ex) {
            logger.info(String.format(" Exception printing rescore file for file %s", relPath));
        }
    }
    
    logger.info(String.format("\n Number of files not accepted: %d", totalFiles - numAccepted));
    bw.newLine();
    bw.write(String.format("Number of files not accepted: %d", totalFiles - numAccepted));
    bw.newLine();
    if (includeRejected) {
        for (int i = numAccepted; i < totalFiles; i++) {
            int fileIndex = energies[i].getIndex();
            File pointedFile = rescoredFiles[fileIndex];
            Path pointedPath = Paths.get(pointedFile.getCanonicalPath());
            String relPath = pwdPath.relativize(pointedPath).toString();
            double thisEnergy = energies[i].getDoubleValue();
            
            logger.info(String.format(" Non-accepted file %d energy %9.3f > %9.3f kcal/mol", (i+1), thisEnergy, minEnergy + acceptThreshold));
        logger.info(String.format(" %s", relPath));
            try {
                bw.write(String.format("Non-accepted file: %s rank %d energy %f", relPath, (i+1), thisEnergy));
                bw.newLine();
            } catch (IOException ex) {
                logger.info(String.format(" Exception printing rescore file for file %s", relPath));
            }
        }
    }
} finally {
    if (bw != null) bw.close();
}

