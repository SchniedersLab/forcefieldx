
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

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// ENERGY
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;

import static java.lang.Math.exp;
import static java.lang.Math.random;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.RotamerOptimization;
import ffx.algorithms.RotamerOptimization.Direction;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import edu.rit.pj.Comm;
import ffx.algorithms.MolecularDynamics;
import java.util.Scanner;

// Things below this line normally do not need to be changed.
// ===============================================================================================
int library = 1;
int startResID = -1;
int allStartResID = 1;
int finalResID = -1;
int windowSize = 7;
int increment = 3;
double distance = 2.0;
Direction direction = Direction.FORWARD;
boolean revert = false;
int algorithm = 1;
int numberOfIterations = 1;
boolean min = false;
double eps = 0.01;
boolean threeBodyTerm = true;
boolean useGoldstein = true;
String chain = "A";
int pruningType = 2;
int ensemble = 1;
double buffer = 0.0;
int counter = 1;
double nucleicCorrectionThreshold = 0;
int minimumNumberAcceptedNARotamers = 10;
double pruningFactor = 1.0;
double singletonNAPruningFactor = 1.5;
boolean verboseEnergies = true;
boolean useOrigCoordsRotamer = true;
double threeBodyCutoffDist = 9.0;
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
int forceResiduesStart = -1;
int forceResiduesEnd = -1;
double superpositionThreshold = 0.1;
double singletonClashThreshold = 20.0;
double pairClashThreshold = 50.0;
Thermostats thermostat = null;
Integrators integrator = null;
String fileType = "PDB";
double restartFrequency = 10000;
int nSteps = 1000;
double timeStep = 1;
double printInterval = 0.01;
double saveInterval = 10;
boolean initVelocities = true;
double temperature = 298.15;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rotamer.md [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.l(longOpt:'library', args:1, argName:'1', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.a(longOpt:'algorithm', args:1, argName:'1', 'Algorithm for rotamer optimization; choices are independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5).');
cli.iT(longOpt:'iterations', args:1, argName:'1', 'Number of total moves.');
cli.tH(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.iN(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic]');
cli.w(longOpt:'window', args:1, argName:'7', 'Size of the sliding window with respect to adjacent residues');
cli.i(longOpt:'increment', args:1, argName:'3', 'Distance sliding window shifts with respect to adjacent residues');
cli.r(longOpt:'cutoff', args:1, argName:'2.0', 'The sliding window cutoff radius (Angstroms).')
cli.d(longOpt:'direction', args:1, argName:'Forward', 'Direction of the sliding window or box optimization (boxes indexed by increasing Z,Y,X): [Forward / Backward]');
cli.z(longOpt:'undo', args:1, argName:'false', 'Window optimizations that do not lower the energy are discarded.');
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting residue to perform the rotamer search on (-1 exits). For box optimization, first box to optimize.');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Final residue to perform the rotamer search on (-1 exits). For box optimization, last box to optimize.');
cli.c(longOpt:'chain', args:1, argName:'A', 'Chain the residue selection belongs to.');
cli.m(longOpt:'minimize', args:1, argName:'0.01', 'Minimize the final structure to the given RMS gradient (Kcal/mole/A).');
cli.t(longOpt:'threeBody', args:1, argName:'true', 'Include 3-Body interactions in elimination criteria.');
cli.g(longOpt:'Goldstein', args:1, argName:'true', 'True to use Goldstein Criteria, False to use DEE');
cli.p(longOpt:'pruning', args:1, argName:'2', 'Prune no clashes (0), only single clashes (1), or all clashes (2).');
cli.e(longOpt:'ensemble', args:1, argName:'1', 'Produce an ensemble of this many of the most favorable structures.');
cli.b(longOpt:'buffer', args:1, argName:'0.0/5.0', 'Sets a starting energy buffer value for use with ensemble search.');
cli.x(longOpt:'all', args:1, argName:'1', 'Optimize all residues in the system beginning from the passed residue number (overrides other options); for box optimization, optimizes all boxes from the passed index.');
cli.nt(longOpt:'nucleicCorrectionThreshold', args:1, argName: '0', 'Nucleic acid Rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).');
cli.mn(longOpt: 'minimumNumberAcceptedNARotamers', args:1, argName: '10', 'Minimum number of NA rotamers to be accepted if a threshold distance is enabled.');
cli.pf(longOpt: 'pruningFactor', args:1, argName: '1.0', 'Multiplier of pruning constraints for nucleic acids.');
cli.sf(longOpt: 'nucleicSinglesPruningFactor', args:1, argName: '1.5', 'Constant multiplier of singleton pruning constraint for nucleic acids');
cli.v(longOpt: 'verboseEnergies', args:1, argName: 'true', 'Calculates and prints beginning and default-conformation energies not necessary for the DEE algorithm');
cli.o(longOpt:'includeOriginal', args:1, argName:'true', 'Include starting coordinates as their own rotamer.');
cli.td(longOpt: 'threeBodyCutoffDist', args:1, argName: '9.0', 'Angstrom distance beyond which three-body interactions will be truncated (-1 for no cutoff).');
cli.dO(longOpt: 'decomposeOriginal', args: 1, argName:'false', 'Only print energy decomposition of original-coordinates rotamers; overrides other flags!');
cli.pE(longOpt: 'parallelEnergies', args: 1, argName:'true', 'Compute rotamer energies in parallel.');
cli.eR(longOpt: 'energyRestart', args: 1, argName:'filename', 'Load energy restart file from a previous run. Ensure that all parameters are the same!');
cli.nB(longOpt: 'numXYZBoxes', args:1, argName: '3,3,3', 'Specify number of boxes along X, Y, and Z');
cli.bB(longOpt: 'boxBorderSize', args: 1, argName: '3.0', 'Extent of overlap between optimization boxes (Angstroms).');
cli.bL(longOpt: 'approxBoxLength', args: 1, argName: '0.0', 'Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes). Box sizes rounded up to make whole number of boxes along each axis. 0 disables this function.');
//cli.bO(longOpt: 'boxOrder', args: 1, argName: 'IncreasingZYX', 'Order of boxes optimized. Options: Random, IncreasingZYX, DecreasingZYX, IncreasingXYZ, DecreasingXYZ, File (specify file name). Presently, always uses default behavior.');
cli.bC(longOpt: 'boxInclusionCriterion', args: 1, argName: '1', 'Criterion to use for adding residues to boxes. (1) uses C alpha only (N1/9 for nucleic acids). (2) uses any atom. (3) uses any rotamer');
cli.fR(longOpt: 'forceResidues', args: 1, argName: '-1,-1', 'Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.');
cli.sT(longOpt: 'superpositionThreshold', args: 1, argName: '0.1', 'Sets the maximum atom-atom distance (Angstroms) which will cause a pair or triple energy to be defaulted to 1.0E100 kcal/mol.');
cli.sC(longOpt: 'singletonClashThreshold', args: 1, argName: '20.0', 'Sets the threshold for singleton pruning.');
cli.pC(longOpt: 'pairClashThreshold', args: 1, argName: '50.0', 'Sets the threshold for pair pruning');
cli.mB(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.mD(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
cli.mI(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic]');
cli.mL(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.mN(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.mT(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.mW(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.mS(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.mF(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Ensemble.
if (options.pE) {
    parallelEnergies = Boolean.parseBoolean(options.pE);
}

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

// Number of iterations.
if (options.iT) {
    numberOfIterations = Integer.parseInt(options.iT);
}
// Load the number of molecular dynamics steps.
if (options.mN) {
    nSteps = Integer.parseInt(options.mN);
}

// Write dyn interval in picoseconds
if (options.mS) {
    restartFrequency = Double.parseDouble(options.mS);
}

if (options.mF) {
    fileType = options.mF.toUpperCase();
}
// Load the time steps in femtoseconds.
if (options.mD) {
    timeStep = Double.parseDouble(options.mD);
}

// Report interval in picoseconds.
if (options.mL) {
    printInterval = Double.parseDouble(options.mL);
}

// Write snapshot interval in picoseconds.
if (options.mW) {
    saveInterval = Double.parseDouble(options.mW);
}

// Temperature in degrees Kelvin.
if (options.mT) {
    temperature = Double.parseDouble(options.mT);
}

// Thermostat.
if (options.mB) {
    try {
        thermostat = Thermostats.valueOf(options.mB.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

// Integrator.
if (options.mI) {
    try {
        integrator = Integrators.valueOf(options.mI.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
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
    } else {
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

if (options.td) {
    threeBodyCutoffDist = Double.parseDouble(options.td);
}

if (options.g) {
    useGoldstein = Boolean.parseBoolean(options.g);
}

if (options.p) {
    pruningType = Integer.parseInt(options.p);
}

// Read in command line.
String filename = arguments.get(0);

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

if (options.sP) {
    singletonClashThreshold = Double.parseDouble(options.sP);
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
    if (!parallelEnergies || algorithm == 4 || algorithm == 5) {
        logger.severe(" FFX shutting down: energy restart only implemented for parallelized global optimizations.");
    }
    useEnergyRestart = true;
    energyRestartFile = new File(options.eR);
}

if (options.nB) {
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
}

if (options.bC) {
    boxInclusionCriterion = Integer.parseInt(options.bC);
}

if (options.bS) {
    boxBorderSize = Double.parseDouble(options.bS);
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
    if (!options.x) {
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

if (pruningFactor < 0) {
    logger.warning("\n Pruning factor must be >= 0.  Setting to default of 1.0.\n");
    pruningFactor = 1;
}

if (singletonNAPruningFactor < 0) {
    logger.warning("\n Pruning factor must be >= 0.  Setting to default of 1.5.\n");
    singletonNAPruningFactor = 1.5;
}

open(filename);

Potential potential;
RotamerOptimization rotamerOptimization = new RotamerOptimization(active, active.getPotentialEnergy(), sh);
rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
rotamerOptimization.setThreeBodyCutoffDist(threeBodyCutoffDist);
rotamerOptimization.setGoldstein(useGoldstein);
rotamerOptimization.setPruning(pruningType);
rotamerOptimization.setEnsemble(ensemble, buffer);
rotamerOptimization.setWindowSize(windowSize);
rotamerOptimization.setDirection(direction);
rotamerOptimization.setDistanceCutoff(distance);
rotamerOptimization.setIncrement(increment);
rotamerOptimization.setRevert(revert);
rotamerOptimization.setNucleicCorrectionThreshold(nucleicCorrectionThreshold);
rotamerOptimization.setMinimumNumberAcceptedNARotamers(minimumNumberAcceptedNARotamers);
rotamerOptimization.setPruningFactor(pruningFactor);
rotamerOptimization.setSingletonNAPruningFactor(singletonNAPruningFactor);
rotamerOptimization.setVerboseEnergies(verboseEnergies);
rotamerOptimization.setParallelEnergies(parallelEnergies);
rotamerOptimization.setBoxBorderSize(boxBorderSize);
rotamerOptimization.setApproxBoxLength(approxBoxLength);
rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
rotamerOptimization.setSuperpositionThreshold(superpositionThreshold);
rotamerOptimization.setSingletonClashThreshold(singletonClashThreshold);
rotamerOptimization.setPairClashThreshold(pairClashThreshold);
if (useEnergyRestart) {
    rotamerOptimization.setEnergyRestartFile(energyRestartFile);
}
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
        Polymer[] polymers = active.getChains();
        int nPolymers = polymers.length;
        for (int p=0; p<nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            int nResidues = residues.size();
            for (int i=0; i<nResidues; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
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
    } else if (options.c) {
        rotamerOptimization.setResidues(options.c, startResID, finalResID);
    } else {
        rotamerOptimization.setResidues(startResID, finalResID);
    }
} else {
    ArrayList<Residue> residueList = new ArrayList<Residue>();
    Polymer[] polymers = active.getChains();
    int nPolymers = polymers.length;
    for (int p=0; p<nPolymers; p++) {
        Polymer polymer = polymers[p];
        ArrayList<Residue> residues = polymer.getResidues();
        int nResidues = residues.size();
        for (int i=0; i<nResidues; i++) {
            Residue residue = residues.get(i);
            Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
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

boolean master = true;
if (Comm.world().size() > 1) {
    int rank = Comm.world().rank();
    if (rank != 0) {
        master = false;
    }
}

energy();
RotamerLibrary.measureRotamers(residueList, false);

if (decomposeOriginal && master) {
    RotamerLibrary.setUseOrigCoordsRotamer(true);
    if (options.x) {
        rotamerOptimization.decomposeOriginal();
    } else {
        Residue[] residueArray = residueList.toArray(new Residue[residueList.size()]);
        rotamerOptimization.decomposeOriginal(residueArray);
    }
    logger.info(String.format("\n"));
    energy();
    return;
}

double[][][] original;
double[][][] move;
boolean accepted = false;
// Store the original coordinates.
original = rotamerOptimization.storeCoordinates(residueList);
// Do the first rotamer optimization.
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
e1 = returnEnergy();
rotamerOptimization.revertCoordinates(residueList, original);
for (int i = 0; i < numberOfIterations; i++) {
    logger.info(String.format("\n Iteration %d of %d.", i+1, numberOfIterations));
    if (accepted) {
        original = rotamerOptimization.storeCoordinates(residueList);
    }
    accepted = true;
    // Do the molecular dynamics run.
    File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn");
    if (!dyn.exists()) {
        dyn = null;
    }
    MolecularDynamics molecularDynamics = new MolecularDynamics(active, active.getPotentialEnergy(), active.getProperties(), sh, thermostat, integrator);
    molecularDynamics.setFileType(fileType);
    molecularDynamics.setRestartFrequency(restartFrequency);
    molecularDynamics.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
    move = rotamerOptimization.storeCoordinates(residueList);
    // Do the second rotamer optimization.
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
    e2 = returnEnergy();
    if (e1 >= e2) {
        logger.info(String.format(" %16.8f <= %16.8f\n\n Move accepted.\n",e2,e1));
        // Check if it is the last iteration.
        if (i != numberOfIterations-1) {
            rotamerOptimization.revertCoordinates(residueList, move);
            e1 = returnEnergy();
        }
    } else if (e1 < e2) {
        // Metropolis criterion.
        double k = 1.98720415e-3;
        double deltaE = e2 - e1;
        double metropolis = random();
        double boltzmann = exp(-deltaE/(k*temperature));
        if (metropolis < boltzmann) {
            logger.info(String.format("\n%16.8f < %16.8f Move accepted.\n",metropolis,boltzmann));
            if (i != numberOfIterations-1) {
                rotamerOptimization.revertCoordinates(residueList, move);
                e1 = returnEnergy();
            }
        } else {
            logger.info(String.format("\n%16.8f >= %16.8f Move rejected.\n",metropolis,boltzmann));
            if (i != numberOfIterations-1) {
                rotamerOptimization.revertCoordinates(residueList, original);
                e1 = returnEnergy();
            }
            accepted = false;
        }
    }
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