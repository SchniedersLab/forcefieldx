/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// ENERGY
import ffx.potential.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3;
import ffx.potential.bonded.Residue.ResidueType;

// X-Ray imports
import ffx.xray.RealSpaceData;
import ffx.xray.RefinementEnergy;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.RotamerOptimization
import ffx.algorithms.RotamerOptimization.Direction;
import edu.rit.pj.Comm
import java.util.Scanner;

// Things below this line normally do not need to be changed.
// ===============================================================================================
int library = 2;
int startResID = -1;
int allStartResID = 1;
int finalResID = -1;
int windowSize = 7;
int increment = 3;
double distance = 2.0;
Direction direction = Direction.FORWARD;
boolean revert = false;
int algorithm = 1;
boolean min = false;
double eps = 0.01;
boolean threeBodyTerm = true;
boolean useGoldstein = true;
String chain = "A";
int pruningType = 2;
int ensemble = 1;
double buffer = 0.0;
double ensembleTarget = 0.0;
int counter = 1;
double nucleicCorrectionThreshold = 0;
int minimumNumberAcceptedNARotamers = 10;
double pruningFactor = 1.0;
double singletonNAPruningFactor = 1.5;
boolean verboseEnergies = true;
boolean useOrigCoordsRotamer = true;
//double threeBodyCutoffDist = 9.0;
boolean decomposeOriginal = false;
boolean parallelEnergies = true;
boolean useEnergyRestart = false;
File energyRestartFile;
int[] numXYZBoxes = new int[3];
double boxBorderSize = 3.0;
double approxBoxLength = 0;
int boxInclusionCriterion = 1;
int boxStart = -1;
int boxEnd = -1;
/*boolean useBoxDimensions = false;
double[] boxDimensions = new double[6];
double superboxBuffer = 2.0; // Actual default behavior is 8.0, but if useBoxDimensions, this should be the old value of 2A.*/
int forceResiduesStart = -1;
int forceResiduesEnd = -1;
double superpositionThreshold = 0.1;
double singletonClashThreshold = 20.0;
double pairClashThreshold = 50.0;

// prototype
boolean video_writeVideo = false;
boolean video_ignoreInactiveAtoms = false;
boolean video_skipEnergies = false;

// maximum number of refinement cycles
int maxiter = 1000;

// suffix to append to output data
String suffix = "_rsrefine";

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc realSpace.rotamer [options] <filename> [data]');
cli.h(longOpt:'help', 'Print this help message.');
cli.l(longOpt:'library', args:1, argName:'2', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.a(longOpt:'algorithm', args:1, argName:'1', 'Choices are independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5).');
cli.w(longOpt:'window', args:1, argName:'7', 'Size of the sliding window with respect to adjacent residues');
cli.r(longOpt:'cutoff', args:1, argName:'2.0', 'The sliding window cutoff radius (Angstroms).')
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting residue to perform the rotamer search on (-1 exits). For box optimization, first box to optimize.');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Final residue to perform the rotamer search on (-1 exits). For box optimization, last box to optimize.');
cli.c(longOpt:'chain', args:1, argName:'A', 'Chain the residue selection belongs to.');
cli.m(longOpt:'minimize', args:1, argName:'0.01', 'Minimize the final structure to the given RMS gradient (Kcal/mole/A).');
cli.t(longOpt:'threeBody', args:1, argName:'true', 'Include 3-Body interactions in elimination criteria.');
cli.p(longOpt:'pruning', args:1, argName:'2', 'Prune no clashes (0), only single clashes (1), or all clashes (2).');
cli.x(longOpt:'all', args:1, argName:'1', 'Optimize all residues in the system beginning from the passed residue number (overrides other options); for box optimization, optimizes all boxes from the passed index.');
cli.v(longOpt:'verboseEnergies', args:1, argName: 'true', 'Calculates and prints beginning and default-conformation energies not necessary for the DEE algorithm');
cli.dt(longOpt:'data', args:2, valueSeparator:',', argName:'data.map,1.0', 'specify input data filename (or simply provide the datafilename argument after the PDB file) and weight to apply to the data (wA)');
cli.I(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed refinement iterations');
cli.X(longOpt:'suffix', args:1, argName:'_rsrefine', 'output suffix');
cli.o(longOpt:'includeOriginal', args:1, argName:'true', 'Include starting coordinates as their own rotamer.');
cli.dO(longOpt:'decomposeOriginal', args: 1, argName:'false', 'Only print energy decomposition of original-coordinates rotamers; overrides other flags!');
cli.eR(longOpt:'energyRestart', args: 1, argName:'filename', 'Load energy restart file from a previous run. Ensure that all parameters are the same!');
cli.nB(longOpt:'numXYZBoxes', args:1, argName: '3,3,3', 'Specify number of boxes along X, Y, and Z');
cli.bB(longOpt:'boxBorderSize', args: 1, argName: '3.0', 'Extent of overlap between optimization boxes (Angstroms).');
cli.bL(longOpt:'approxBoxLength', args: 1, argName: '0.0', 'Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes). Box sizes rounded up to make whole number of boxes along each axis. 0 disables this function.');
//cli.bO(longOpt: 'boxOrder', args: 1, argName: 'IncreasingZYX', 'Order of boxes optimized. Options: Random, IncreasingZYX, DecreasingZYX, IncreasingXYZ, DecreasingXYZ, File (specify file name). Presently, always uses default behavior.');
cli.bC(longOpt:'boxInclusionCriterion', args: 1, argName: '1', 'Criterion to use for adding residues to boxes. (1) uses C alpha only (N1/9 for nucleic acids). (2) uses any atom. (3) uses any rotamer');
cli.fR(longOpt:'forceResidues', args: 1, argName: '-1,-1', 'Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.');
cli.lR(longOpt:'listResidues', args: 1, argName: '-1', 'Choose a list of individual residues to optimize (eg. A11,A24,B40).');
cli.vw(longOpt:'videoWriter', args: 0, 'Prototype video snapshot output; skips energy calculation.');
cli.sO(longOpt:'sequenceOptimization', args:1, argName: '-1', 'Choose a list of individual residues to sequence optimize.');
cli.tO(longOpt:'titrationOptimization', args:1, argName: '-1', 'Choose a list of individual residues to titrate (protonation state optimization).');
cli.nt(longOpt:'nucleicCorrectionThreshold', args:1, argName: '0', 'Nucleic acid Rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).');
cli.mn(longOpt:'minimumAcceptedNARotamers', args:1, argName: '10', 'Minimum number of NA rotamers to be accepted if a threshold distance is enabled.');
/**
 * Now handled by system keys.
// bD is used as a single string argument to allow for easier argument size checking.
cli.bD(longOpt:'boxDimensions', args:1, argName:'buffer,xmin,xmax,ymin,ymax,zmin,zmax', 'If set, box optimization only uses supplied coordinates plus buffer.');
cli.i(longOpt:'increment', args:1, argName:'3', 'Distance sliding window shifts with respect to adjacent residues');
cli.d(longOpt:'direction', args:1, argName:'Forward', 'Direction of the sliding window or box optimization (boxes indexed by increasing Z,Y,X): [Forward / Backward]');
cli.z(longOpt:'undo', args:1, argName:'false', 'Window optimizations that do not lower the energy are discarded.');
cli.g(longOpt:'Goldstein', args:1, argName:'true', 'True to use Goldstein Criteria, False to use DEE');
cli.pE(longOpt:'parallelEnergies', args: 1, argName:'true', 'Compute rotamer energies in parallel.');
cli.sT(longOpt:'superpositionThreshold', args: 1, argName: '0.1', 'Sets the maximum atom-atom distance (Angstroms) which will cause a pair or triple energy to be defaulted to 1.0E100 kcal/mol.');
cli.e(longOpt:'ensemble', args:1, argName:'1', 'Produce an ensemble of this many of the most favorable structures.');
cli.eT(longOpt: 'ensembleTarget', args:1, argName:'0.0', 'Produces an ensemble of structures with energies between GMEC and GMEC + value (kcal/mol).');
cli.b(longOpt:'buffer', args:1, argName:'0.0/5.0', 'Sets a starting energy buffer value for use with ensemble search.');
cli.td(longOpt:'threeBodyCutoffDist', args:1, argName: '9.0', 'Angstrom distance beyond which three-body interactions will be truncated (-1 for no cutoff).');
cli.pf(longOpt:'pruningFactor', args:1, argName: '1.0', 'Multiplier of pruning constraints for nucleic acids.');
cli.sf(longOpt:'nucleicSinglesPruningFactor', args:1, argName: '1.5', 'Constant multiplier of singleton pruning constraint for nucleic acids');
cli.sC(longOpt:'singletonClashThreshold', args: 1, argName: '20.0', 'Sets the threshold for singleton pruning.');
cli.pC(longOpt:'pairClashThreshold', args: 1, argName: '50.0', 'Sets the threshold for pair pruning');
*/

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

List<String> resList = new ArrayList<>();
if (options.lR) {
    def tok = (options.lR).tokenize('.');
    for (String t : tok) {
        logger.info("Adding " + t);
        resList.add(t);
    }
}

List<String> sequenceOptimizationList = new ArrayList<>();
if (options.sO) {
    def tok = (options.sO).tokenize('.');
    for (String t : tok) {
        logger.info(" Sequence optimizing " + t);
        sequenceOptimizationList.add(t);
    }
    if (System.getProperty("RELATIVE_SOLVATION") == null) {
        System.setProperty("RELATIVE_SOLVATION", "AUTO");
    }
}

List<String> titrationOptimizationList = new ArrayList<>();
if (options.tO) {
    def tok = (options.tO).tokenize('.');
    for (String t : tok) {
        logger.info(" Protonation state optimizing " + t);
        titrationOptimizationList.add(t);
    }
}

if (options.vw) {
    video_writeVideo = true;
    video_ignoreInactiveAtoms = true;
    video_skipEnergies = true;
}

// Ensemble.
/*if (options.pE) {
    parallelEnergies = Boolean.parseBoolean(options.pE);
}*/

if (options.e) {
    ensemble = Integer.parseInt(options.e);
    buffer = 5.0;
}

// Buffer.
if (options.b) {
    buffer = Double.parseDouble(options.b);
}

// Rotamer Library.
if (options.l) {
    library = Integer.parseInt(options.l);
}

// Algorithm.
if (options.a) {
    algorithm = Integer.parseInt(options.a);
}

// Implements different start/stop behavior for box, non-box optimizations.
if (algorithm == 5) {
    if (options.x) {
        boxStart = Integer.parseInt(options.x);
        --boxStart; // Internal machinery indexed 0 to (n-1)
        if (boxStart < 0) {
            logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
            return;
        }
    } else {
        if (options.s) {
            boxStart = Integer.parseInt(options.s);
            --boxStart; // Internal machinery indexed 0 to (n-1)
        }
        if (options.f) {
            boxEnd = Integer.parseInt(options.f);
            --boxEnd; // Internal machinery indexed 0 to (n-1)
        }
        if ((options.f && boxEnd < boxStart) || boxStart < 0 || (options.f && boxEnd < 0)) {
            logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
            return;
        }
    }
    // boxEnd is set to -1 by default as a flag for "no end box set".
} else {
    // Starting residue.
    if (options.s) {
        startResID = Integer.parseInt(options.s);
    }
    // Load the number iterations.
    if (options.f) {
        finalResID = Integer.parseInt(options.f);
    }
    if (options.x) {
        allStartResID = Integer.parseInt(options.x);
    } else if (!options.lR) {
        if (finalResID < startResID || startResID < 0 || finalResID < 0) {
            return;
        }
    }
}

// Sliding window size.
if (options.w) {
    windowSize = Integer.parseInt(options.w);
}

// Distance sliding window moves.
if (options.i) {
    increment = Integer.parseInt(options.i);
}

// Amount of space distance matrix checks for.
if (options.r) {
    distance = Double.parseDouble(options.r);
}

// Direction of sliding window.
if (options.d) {
    try {
        direction = Direction.valueOf(options.d.toUpperCase());
    } catch (Exception e) {
        direction = null;
    }
}

// Undo an unfavorable change.
if (options.z) {
    revert = Boolean.parseBoolean(options.z);
}

if (options.m) {
    min = true;
    eps = Double.parseDouble(options.m);
}

if (options.t) {
    threeBodyTerm = Boolean.parseBoolean(options.t);
}

/*if (options.td) {
    threeBodyCutoffDist = Double.parseDouble(options.td);
}*/

if (options.g) {
    useGoldstein = Boolean.parseBoolean(options.g);
}

if (options.p) {
    pruningType = Integer.parseInt(options.p);
}

// Read in command line.
String filename = arguments.get(0);

systems = open(filename);

// set up real space map data (can be multiple files)
List mapFiles = new ArrayList();
int nDiffractionData = 0;
if (arguments.size() > 1) {
    String dataFileName = arguments.get(1);
    if (FilenameUtils.isExtension(dataFileName, "map")) {
        RealSpaceFile realspacefile = new RealSpaceFile(dataFileName, 1.0);
        mapFiles.add(realspacefile);
    } else {
        DiffractionFile diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
        DiffractionData diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
        diffractionData.scaleBulkFit();
        diffractionData.printStats();
        String mapFileName = String.format("%s_ffx_%d", FilenameUtils.removeExtension(dataFileName), ++nDiffractionData);
        diffractionData.writeMaps(mapFileName);
        mapFiles.add(new RealSpaceFile(mapFileName + "_2fofc.map", 1.0));
    }
}

if (options.dt) {
    for (int i=0; i<options.dts.size(); i+=2) {
	double wA = Double.parseDouble(options.dts[i+1]);
        String dataFileName = options.dts[i];
        if (FilenameUtils.isExtension(dataFileName, "map")) {
            RealSpaceFile realspacefile = new RealSpaceFile(dataFileName, wA);
            mapFiles.add(realspacefile);
        } else {
            DiffractionFile diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
            DiffractionData diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
            diffractionData.scaleBulkFit();
            diffractionData.printStats();
            String mapFileName = String.format("%s_ffx_%d", FilenameUtils.removeExtension(dataFileName), ++nDiffractionData);
            diffractionData.writeMaps(mapFileName);
            mapFiles.add(new RealSpaceFile(mapFileName + "_2fofc.map", wA));
        }
    }
}

if (options.I) {
    maxiter = Integer.parseInt(options.I);
}

if (options.X) {
    suffix = options.X;
}

if (options.nt) {
    nucleicCorrectionThreshold = Double.parseDouble(options.nt);
}

if (options.mn) {
    minimumNumberAcceptedNARotamers = Integer.parseInt(options.mn);
}

if (options.pf) {
    pruningFactor = Double.parseDouble(options.pf);
}

if (options.sf) {
    singletonNAPruningFactor = Double.parseDouble(options.sf);
}

if (options.sT) {
    superpositionThreshold = Double.parseDouble(options.sT);
}

if (options.sC) {
    singletonClashThreshold = Double.parseDouble(options.sC);
}

if (options.pP) {
    pairClashThreshold = Double.parseDouble(options.pP);
}

if (options.v) {
    verboseEnergies = Boolean.parseBoolean(options.v);
}

if (options.o) {
    useOrigCoordsRotamer = Boolean.parseBoolean(options.o);
}

if (options.dO) {
    decomposeOriginal = Boolean.parseBoolean(options.dO);
    if (decomposeOriginal) {
        algorithm = 0;
        useOrigCoordsRotamer = true;
    }
}

if (options.eR) {
    /*if (!parallelEnergies || algorithm == 4) {
        logger.severe(" FFX shutting down: energy restart only implemented for parallelized global optimizations.");
    }*/
    useEnergyRestart = true;
    energyRestartFile = new File(options.eR);
}

if (algorithm == 5 && options.nB) {
    String input = options.nB;
    //logger.info("Input: " + input);
    Scanner boxNumInput = new Scanner(input);
    boxNumInput.useDelimiter(",");
    int inputLoopCounter = 0;
    numXYZBoxes[0] = 3; // Default
    while (inputLoopCounter < numXYZBoxes.length) {
        if (boxNumInput.hasNextInt()) {
            numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt();
            inputLoopCounter++;
        } else if (boxNumInput.hasNextDouble()) {
            numXYZBoxes[inputLoopCounter] = boxNumInput.nextDouble();
            inputLoopCounter++;
            logger.info("Double input to nB truncated to integer.");
        } else if (boxNumInput.hasNext()) {
            logger.info("Non-numeric input to nB discarded");
            boxNumInput.next();
        } else {
            logger.info("Insufficient input to nB. Non-input values assumed either equal to X or default to 3");
            break;
        }
    }
    boxNumInput.close();
    for (int i = inputLoopCounter; i < numXYZBoxes.length; i++) {
        numXYZBoxes[i] = numXYZBoxes[0];
    }
    for (int i = 0; i < numXYZBoxes.length; i++) {
        if (numXYZBoxes[i] == 0) {
            numXYZBoxes[i] = 3;
            logger.info("Input of zero to nB reset to default of three.");
        } else if (numXYZBoxes[i] < 0) {
            numXYZBoxes[i] = -1*numXYZBoxes[i];
            logger.info("Input of negative number to nB reset to positive number");
        }
    }
} else {
    for (int i = 0; i < numXYZBoxes.length; i++) {
        // Default case
        numXYZBoxes[i] = 3;
    }
}

/*if (options.bD) {
    try {
        String[] bdTokens = options.bD.split(",+");
        if (bdTokens.length != 7) {
            logger.warning(" Improper number of arguments to boxDimensions; default settings used.");
        } else {
            for (int i = 1; i < 7; i+=2) {
                boxDimensions[i-1] = Double.parseDouble(bdTokens[i]);
                boxDimensions[i] = Double.parseDouble(bdTokens[i+1]);
                if (boxDimensions[i] < boxDimensions[i-1]) {
                    logger.info(String.format(" Improper dimension min %8.5f > max %8.5f; max/min reversed.", boxDimensions[i-1], boxDimensions[i]));
                    double temp = boxDimensions[i];
                    boxDimensions[i] = boxDimensions[i-1];
                    boxDimensions[i-1] = temp;
                }
            }
            superboxBuffer = Double.parseDouble(bdTokens[0]);
            useBoxDimensions = true;
        }
    } catch (Exception ex) {
        logger.warning(String.format(" Error in parsing box dimensions: input discarded and defaults used: %s.", ex));
        useBoxDimensions = false;
    }
}*/

if (options.bC) {
    boxInclusionCriterion = Integer.parseInt(options.bC);
}

if (options.bB) {
    boxBorderSize = Double.parseDouble(options.bB);
}
if (options.bL) {
    approxBoxLength = Double.parseDouble(options.bL);
    if (approxBoxLength < 0) {
        logger.info(" Negative box length value inputted; changing to -1 * input.")
        approxBoxLength *= -1;
    }
}

// Only makes sense for sliding window; otherwise don't touch anything.
if (options.fR && algorithm == 4) {
    String input = options.fR;
    Scanner frScan = new Scanner(input);
    frScan.useDelimiter(",");
    try {
        if (!frScan.hasNextInt()) {
            frScan.next(); // Discards extra input to indicate a negative value of frStart.
        }
        forceResiduesStart = frScan.nextInt();
        forceResiduesEnd = frScan.nextInt();
    } catch (Exception ex) {
        logger.severe(String.format(" FFX shutting down: input to -fR could not be parsed as a pair of integers: %s", input));
    }
    if (forceResiduesStart > forceResiduesEnd) {
        logger.info(" Start of range higher than ending: start flipped with end.");
        int temp = forceResiduesStart;
        forceResiduesStart = forceResiduesEnd;
        forceResiduesEnd = temp;
    }
    if (forceResiduesEnd < 1) {
        logger.severe(String.format(" FFX shutting down: end range for -fR must be at least 1; input range %d to %d", forceResiduesStart, forceResiduesEnd));
    }
}


if (algorithm != 5) {
    if (options.lR) {
        String info = "\n Evaluating rotamers for residues ";
        for (String i : resList) {
            info += String.format("%s, ", i);
        }
        logger.info(info);
    } else if (!options.x) {
        logger.info("\n Evaluating rotamers for residues " + startResID + " to " + finalResID);
    } else {
        logger.info("\n Evaluating rotamers for all residues beginning at " + allStartResID);
    }
} else {
    if (!options.x) {
        logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1));
    } else {
        logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1));
    }
}

if (nucleicCorrectionThreshold < 0) {
    logger.warning("\n Correction threshold must be >= 0.  Setting to default of 0 (threshold inactive).\n");
    nucleicCorrectionThreshold = 0;
}

if (minimumNumberAcceptedNARotamers < 1) {
    logger.warning("\n Minimum number of accepted NA rotamers must be a positive integer.\n Setting to default value 10.\n");
    minimumNumberAcceptedNARotamers = 10;
}
/**
 * Now handled by system keys.
 * 
if (pruningFactor < 0) {
    logger.warning("\n Pruning factor must be >= 0.  Setting to default of 1.0.\n");
    pruningFactor = 1;
}

if (singletonNAPruningFactor < 0) {
    logger.warning("\n Pruning factor must be >= 0.  Setting to default of 1.5.\n");
    singletonNAPruningFactor = 1.5;
}
*/

if (mapFiles.size() == 0) {
    RealSpaceFile realspacefile = new RealSpaceFile(systems);
    mapFiles.add(realspacefile);
}

Atom[] atoms = systems[0].getAtomArray();
RealSpaceData realSpaceData = new RealSpaceData(systems, systems[0].getProperties(), mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
RefinementEnergy refinementEnergy = RefinementEnergy.refinementEnergyFactory(realSpaceData, RefinementMode.COORDINATES, null);

double [] x = new double[atoms.length*3];
refinementEnergy.getCoordinates(x);
refinementEnergy.energy(x);

RotamerOptimization rotamerOptimization = new RotamerOptimization(active, refinementEnergy, sh);
rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
rotamerOptimization.setPruning(pruningType);
rotamerOptimization.setWindowSize(windowSize);
rotamerOptimization.setDistanceCutoff(distance);
rotamerOptimization.setNucleicCorrectionThreshold(nucleicCorrectionThreshold);
rotamerOptimization.setMinimumNumberAcceptedNARotamers(minimumNumberAcceptedNARotamers);
/**
 * Now handled by system keys.
rotamerOptimization.setIncrement(increment);
rotamerOptimization.setDirection(direction);
rotamerOptimization.setRevert(revert);
rotamerOptimization.setGoldstein(useGoldstein);
rotamerOptimization.setEnsemble(ensemble, buffer);
rotamerOptimization.setEnsembleTarget(ensembleTarget);
rotamerOptimization.setThreeBodyCutoffDist(threeBodyCutoffDist);
rotamerOptimization.setPruningFactor(pruningFactor);
rotamerOptimization.setSingletonNAPruningFactor(singletonNAPruningFactor);
rotamerOptimization.setSingletonClashThreshold(singletonClashThreshold);
rotamerOptimization.setPairClashThreshold(pairClashThreshold);
rotamerOptimization.setParallelEnergies(parallelEnergies);
rotamerOptimization.setSuperpositionThreshold(superpositionThreshold);
*/
rotamerOptimization.setVerboseEnergies(verboseEnergies);
rotamerOptimization.setBoxBorderSize(boxBorderSize);
rotamerOptimization.setApproxBoxLength(approxBoxLength);
rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
if (useEnergyRestart) {
    rotamerOptimization.setEnergyRestartFile(energyRestartFile);
}
/*if (useBoxDimensions) {
    rotamerOptimization.setBoxDimensions(boxDimensions, superboxBuffer);
}*/
if (options.vw) {
    rotamerOptimization.setVideoWriter(video_writeVideo, video_ignoreInactiveAtoms, video_skipEnergies);
}

/*if (useBoxFile) {
rotamerOptimization.setBoxOrder(boxFile);
} else {
boxOrder = boxOrder.toLowerCase();
rotamerOptimization.setBoxOrder(boxOrder);
}*/

if (library == 1) {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
} else {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
}

if (useOrigCoordsRotamer) {
    RotamerLibrary.setUseOrigCoordsRotamer(true);
}

if (algorithm != 5) {
    if (options.x) {
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = active.getPolymers();
        int nPolymers = polymers.length;
        for (int p=0; p<nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            int nResidues = residues.size();
            for (int i=0; i<nResidues; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers();
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    } else if (nrot > 1) {
                        if (counter >= allStartResID) {
                            residueList.add(residue);
                        }
                    }
                } else if (options.fR) {
                    if (counter >= allStartResID && counter >= forceResiduesStart
                        && counter <= forceResiduesEnd) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }
        rotamerOptimization.setResidues(residueList);
    } else if (options.lR) {
        ArrayList<Residue> residueList = new ArrayList<>();
        Polymer[] polymers = active.getPolymers();
        int n = 0;
        for (String s : resList) {
            Character chainID = s.charAt(0);
            int i = Integer.parseInt(s.substring(1));
            for (Polymer p : polymers) {
                if (p.getChainIDChar() == chainID) {
                    List<Residue> rs = p.getResidues();
                    for (Residue r : rs) {
                        if (r.getResidueIndex() == i) {
                            residueList.add(r);
                            Rotamer[] rotamers = r.getRotamers();
                            if (rotamers != null && rotamers.size() > 1) {
                                n++;
                            }
                        }
                    }
                }
            }
        }
        rotamerOptimization.setResiduesIgnoreNull(residueList);
        if (n < 1) {
            return;
        }
    } else if (options.c) {
        rotamerOptimization.setResidues(options.c, startResID, finalResID);
    } else {
        rotamerOptimization.setResidues(startResID, finalResID);
    }
} else {
    boolean ignoreNA = false;
    String ignoreNAProp = System.getProperty("ignoreNA");
    if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
        ignoreNA = true;
    }
    ArrayList<Residue> residueList = new ArrayList<Residue>();
    Polymer[] polymers = active.getPolymers();
    int nPolymers = polymers.length;
    for (int p=0; p<nPolymers; p++) {
        Polymer polymer = polymers[p];
        ArrayList<Residue> residues = polymer.getResidues();
        int nResidues = residues.size();
        for (int i=0; i<nResidues; i++) {
            Residue residue = residues.get(i);
            if (ignoreNA && residue.getResidueType() == ResidueType.NA) {
                continue;
            }
            Rotamer[] rotamers = residue.getRotamers();
            if (rotamers != null) {
                int nrot = rotamers.length;
                if (nrot == 1) {
                    RotamerLibrary.applyRotamer(residue, rotamers[0]);
                } else if (nrot > 1) {
                    residueList.add(residue);
                }
            }
            counter++;
        }
    }
    rotamerOptimization.setResidues(residueList);
    rotamerOptimization.setBoxStart(boxStart);
    if (!options.x) {
        rotamerOptimization.setBoxEnd(boxEnd);
    }
}

ArrayList<Residue> residueList = rotamerOptimization.getResidues();

if (options.sO) {
    for (String s : sequenceOptimizationList) {
        Character chainID = s.charAt(0);
        int num = Integer.parseInt(s.substring(1));
        for (int i = 0; i < residueList.size(); i++) {
            Residue res = residueList.get(i);
            if (res.getChainID() == chainID && res.getResidueIndex() == num) {
                MultiResidue multiRes = new MultiResidue(res, active.getForceField(), active.getPotentialEnergy());
                for (Polymer polymer : active.getPolymers()) {
                    if (polymer.getChainIDChar() == chainID) {
                        logger.info(String.format(" Adding multiresidue %s to chain %c.", multiRes, chainID));
                        polymer.addMultiResidue(multiRes);
                    }
                }
                for (CommonAminoAcid3 aa : CommonAminoAcid3.values()) {
                    if (aa.toString().equals("PRO") || aa.toString().equals("GLY")) {
                        continue;
                    }
                    if (!aa.toString().equalsIgnoreCase(res.getName())) {
                        logger.info(String.format(" Adding %s to residue %s.", aa.toString(), multiRes.toString()));
                        multiRes.addResidue(new Residue(aa.toString(), res.getResidueIndex(), ResidueType.AA));
                    }
                }
                multiRes.finalize();
                multiRes.requestSetActiveResidue(ResidueEnumerations.AminoAcid3.valueOf(res.getName()));
                active.getPotentialEnergy().reInit();
                residueList.remove(i);
                residueList.add(i, multiRes);
            }
        }
    }
}

if (options.tO) {
    ArrayList<Residue> titrating = new ArrayList<>();
    for (String s : titrationOptimizationList) {
        Character chainID = s.charAt(0);
        int num = Integer.parseInt(s.substring(1));
        for (int i = 0; i < residueList.size(); i++) {
            Residue res = residueList.get(i);
            if (res.getChainID() == chainID && res.getResidueIndex() == num) {
                titrating.add(res);
            }
        }
    }
    rotamerOptimization.titrationSetResidues(titrating);
}

boolean master = true;
if (Comm.world().size() > 1) {
    int rank = Comm.world().rank();
    if (rank != 0) {
        master = false;
    }
}

if (nDiffractionData != 0) {
    logger.info("\n Diffraction files treated as X-ray data; neutron data currently requires generation of maps via xray.scaleBulk");
}

energy();
RotamerLibrary.measureRotamers(residueList, false);

if (decomposeOriginal) {
    RotamerLibrary.setUseOrigCoordsRotamer(true);
    boolean doQuadsInParallel = true;
    if (!doQuadsInParallel) {
        String quadsProp = System.getProperty("evalQuad");
        if (quadsProp != null && quadsProp.equalsIgnoreCase("true")) {
            Residue[] residueArray = residueList.toArray(new Residue[residueList.size()]);

            String quadsCutoffProp = System.getProperty("quadCutoff");
            double quadsCutoff = 5.0;
            if (quadsCutoffProp != null) {
                try {
                    quadsCutoff = Double.parseDouble(quadsCutoffProp);
                } catch (Exception ex) {
                    logger.warning(String.format(" Exception in parsing quads cutoff: %s", ex));
                }
                quadsCutoff = quadsCutoff <= 0 ? 5.0 : quadsCutoff;
            }

            String quadsEvalProperty = System.getProperty("numQuads");
            int numQuads = 1000000;
            if (quadsEvalProperty != null) {
                try {
                    numQuads = Integer.parseInt(System.getProperty("numQuads"));
                } catch (Exception ex) {
                    logger.warning(String.format(" Exception in parsing number of quads to evaluate: %s", ex));
                }
                numQuads = numQuads <= 0 ? 1000000 : numQuads;
            }

            rotamerOptimization.decomposeOriginalQuads(quadsCutoff, numQuads);
        }
    } else if (options.x) {
        rotamerOptimization.decomposeOriginalParallel();
    } else {
        Residue[] residueArray = residueList.toArray(new Residue[residueList.size()]);
        rotamerOptimization.decomposeOriginal(residueArray);
    }
    if (master) {
        logger.info(String.format("\n"));
        energy();
    }
    return;
}

if (algorithm == 1) {
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT);
} else if (algorithm == 2) {
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.GLOBAL_DEE);
} else if (algorithm == 3) {
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.GLOBAL);
} else if (algorithm == 4) {
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.SLIDING_WINDOW);
} else if (algorithm == 5) {
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX_OPTIMIZATION);
}

if (master) {
    if (min) {
        minimize(eps);
    }
    logger.info(" Final Minimum Energy");
    energy();
    String ext = FilenameUtils.getExtension(filename);
    filename = FilenameUtils.removeExtension(filename);
    if (ext.toUpperCase().contains("XYZ")) {
        saveAsXYZ(new File(filename + ".xyz"));
    } else {
        saveAsPDB(new File(filename + ".pdb"));
    }
}
