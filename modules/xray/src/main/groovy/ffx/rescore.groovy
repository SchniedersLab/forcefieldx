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

package ffx.xray

import java.nio.file.Path
import java.nio.file.Paths

import org.apache.commons.io.FilenameUtils

import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.potential.parsers.CoordinateFileFilter
import ffx.realspace.RealSpaceFile
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.Rescore.RescoreStrategy
import ffx.xray.parsers.DiffractionFile

// RMS gradient per atom convergence criteria
double eps = -1.0;

// maximum number of refinement cycles
int maxiter = 1000;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;

String suffix = "_rsc";
boolean addResultsFile = true;
int refineMode = 0;
List<DiffractionFile> diffractionFiles = new ArrayList<>();
List<RealSpaceFile> mapFiles = new ArrayList<>();
List<File> modelFiles = new ArrayList<>();
double acceptThreshold = 0.0;
int numToAccept = 10;
File resultsFile;
String resultsFileName = "ffx_rescore.txt";
boolean versioned = true;
File resultDirectory = null;
boolean printModels = false;
boolean includeRejected = true;
Path pwdPath;
RescoreStrategy rscType = RescoreStrategy.NO_RESCORE;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rescore [options] <modelfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.m(longOpt: 'minimize', args: 1, argName:'1', 'Rescore on energy evaluation (1), minimization (2), minimization with X-ray target (3), or minimization with real-space target (4).');
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

try {
    pwdPath = Paths.get(new File("").getCanonicalPath());
} catch (IOException ex) {
    logger.severe(" FFX could not establish path to current working directory. FFX will not continue.");
}

if (options.m) {
    refineMode = Integer.parseInt(options.m);
    switch (refineMode) {
    case 1:
        rscType = RescoreStrategy.ENERGY_EVAL;
        if (arguments.size() != 1) {
            logger.warning(String.format(" Improper number of arguments (must be 1) for -m mode %d", refineMode));
            return cli.usage();
        }
        break;
    case 2:
        rscType = RescoreStrategy.MINIMIZE;
        if (arguments.size() != 1) {
            logger.warning(String.format(" Improper number of arguments (must be 1) for -m mode %d", refineMode));
            return cli.usage();
        }
        break;
    case 3:
        rscType = RescoreStrategy.XRAY_MIN;
        if (arguments.size() > 2) {
            logger.warning(String.format(" Improper number of arguments (must be 1 or 2) for -m mode %d", refineMode));
            return cli.usage();
        }
        break;
    case 4:
        rscType = RescoreStrategy.RS_MIN;
        if (arguments.size() > 2) {
            logger.warning(String.format(" Improper number of arguments (must be 1 or 2) for -m mode %d", refineMode));
            return cli.usage();
        }
        break;
    default:
        logger.warning(" Improper use of -m flag: must be 1, 2, 3, or 4.");
        return cli.usage();
    }
}

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
    resultsFileName = options.r;
} // Results file generated later.

if (/*!matchDataFiles && */arguments.size() > 1) {
    if (refineMode == 2) {
        DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false);
        diffractionFiles.add(diffractionfile);
    } else if (refineMode == 3) {
        RealSpaceFile realSpaceFile = new RealSpaceFile(arguments.get(1), 1.0);
        mapFiles.add(realSpaceFile);
    }
}
if (/*!matchDataFiles && */options.xd) {
    if (options.rd) { // At present, double-minimization not allowed.
        logger.warning(String.format(" Only one refinement mode (%d) will be used.", refineMode));
        switch (refineMode) {
        case 1:
            logger.info(" No minimization will be applied.");
            break;
        case 2:
            logger.info(" Force field minimization will be applied.");
            break;
        case 3:
            logger.info(" X-ray refinement will be applied.");
            break;
        case 4:
            logger.info(" Real space refinement will be applied.");
            break;
        default:
            logger.severe(String.format(" Invalid -m flag %d", refineMode));
        }
    }
    switch (refineMode) {
    case 0:
    case 1:
    case 2:
        logger.warning(" Data file is only used for X-ray and real-space refinements.");
        break;
    case 3:
        for (int i = 0; i < options.xds.size(); i += 3) {
            double wA = Double.parseDouble(options.xds[i+1]);
            boolean neutron = Boolean.parseBoolean(options.xds[i+2]);
            DiffractionFile diffractionFile = new DiffractionFile(options.xds[i], wA, neutron);
            diffractionFiles.add(diffractionFile);
        }
        break;
    case 4:
        logger.warning("The xrayData flag is only for use with x-ray refinement.");
        break;
    }
}

if (/*!matchDataFiles && */options.rd) {
    switch (refineMode) {
    case 0:
    case 1:
    case 2:
        logger.warning(" Data file is only used for X-ray and real-space refinements.");
        break;
    case 3:
        logger.warning("The realspaceData flag is only for use with realspace refinement.");
        break;
    case 4:
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

if (options.dr) {
    resultDirectory = new File(options.dr);
    if (resultDirectory.exists()) {
        if (resultDirectory.isFile()) {
            logger.severe(String.format(" FFX shutting down: Pathname %s provided for directory points to a file.", resultDirectory.toString()));
        }
    } else {
        boolean success = resultDirectory.mkdirs();
        if (!success) {
            logger.severe(String.format(" FFX shutting down: could not create directory %s", resultDirectory.toString()));
        }
    }
    resultsFileName = FilenameUtils.concat(resultDirectory.getPath(), resultsFileName);
    directoryPath = pwdPath.relativize(Paths.get(resultDirectory.getCanonicalPath()));
}

resultsFile = new File(resultsFileName);
if (resultsFile.exists()) {
    for (int i = 2; i <= 1000; i++) {
        resultsFile = new File(resultsFileName + "_" + i);
        if (!resultsFile.exists()) {
            break;
        }
    }
    if (resultsFile.exists()) {
        logger.warning(String.format(" Versioning failed: appending rescore file to end of file %s", resultsFile.getName()));
        versioned = false;
    }
}

AlgorithmFunctions utils;
try {
    utils = getAlgorithmUtils();
} catch (MissingMethodException ex) {
    utils = new AlgorithmUtils();
}

File[] modelFileArray = new File[modelFiles.size()];
modelFiles.toArray(modelFileArray);
RescoreAndCluster rescorer = new RescoreAndCluster(utils);
rescorer.setAcceptThreshold(acceptThreshold);
rescorer.setDiffractionFiles(diffractionFiles);
rescorer.setFileSuffix(suffix);
rescorer.setIncludeRejected(includeRejected);
rescorer.setMapFiles(mapFiles);
rescorer.setMaxIter(maxiter);
rescorer.setNumToAccept(numToAccept);
rescorer.setPrintModels(printModels);
rescorer.setRefineEps(eps);
rescorer.setRescoreStrategy(rscType);
rescorer.setResultsFile(resultsFile);
try {
    rescorer.setResultDirectory(resultDirectory);
} catch (IOException ex) {
    logger.severe(ex.toString());
}

rescorer.runRsc(modelFileArray);
