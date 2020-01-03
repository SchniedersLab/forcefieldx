//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.xray.groovy.test

import org.apache.commons.io.FilenameUtils

import groovy.cli.picocli.CliBuilder

import edu.rit.mp.DoubleBuf
import edu.rit.pj.Comm

import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.integrators.IntegratorEnum
import ffx.algorithms.dynamics.thermostats.ThermostatEnum
import ffx.algorithms.mc.MCLoop
import ffx.algorithms.optimize.anneal.SimulatedAnnealing
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.OptimizationParameters
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Angle
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.xray.CrystalReciprocalSpace.SolventModel
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.parsers.DiffractionFile

// Default convergence criteria.
double eps = 0.1;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100;

// Frequency to write out restart information in picoseconds.
double restartInterval = 1.0;

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 50000;

// ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
ThermostatEnum thermostat = ThermostatEnum.BERENDSEN;

// IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC ]
IntegratorEnum integrator = IntegratorEnum.RESPA;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// File type of coordinate snapshots to write out.
String fileType = "PDB";

// Value of Lambda.
double lambda = 0.0;

// Monte-Carlo step frequencies for loop moves.
int mcStepFrequency = 1000;

// Rotamer Optimization
boolean runRotamer = false;

// Simulated Annealing
boolean runSimulatedAnnealing = false;

// OST
boolean runOST = true;

// Monte Carlo with KIC
boolean runMCLoop = false;

// Local minimization mode
boolean localMin = false;

// Check if OST found any loop
boolean loopBuildError = false;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc loopBuilder [options] <pdbFile> <diffractionDataFile>');
cli.h(longOpt: 'help', 'Print this help message.');
cli.e(longOpt: 'eps', args: 1, argName: '1.0', 'RMS gradient convergence criteria');
cli.n(longOpt: 'steps', args: 1, argName: '10000', 'Number of molecular dynamics steps.');
cli.d(longOpt: 'dt', args: 1, argName: '1.0', 'Time discretization step (fsec).');
cli.m(longOpt: 'minimize', 'Local minimization of loop residues (need -s and -f flags).');
cli.r(longOpt: 'report', args: 1, argName: '0.01', 'Interval to report thermodyanamics (psec).');
cli.w(longOpt: 'write', args: 1, argName: '100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt: 'temperature', args: 1, argName: '298.15', 'Temperature in degrees Kelvin.');
cli.g(longOpt: 'bias', args: 1, argName: '0.01', 'Gaussian bias magnitude (kcal/mol).');
cli.ost(longOpt: 'OST', 'Run OST.');
cli.sa(longOpt: 'simulatedAnnealing', 'Run simulated annealing.');
cli.rot(longOpt: 'rotamer', 'Run rotamer optimization.');
cli.mc(longOpt: 'mcLoop', 'Run Monte Carlo KIC');
cli.a(longOpt: 'all', 'Run optimal pipeline of algorithms.');
cli.s(longOpt: 'start', args: 1, argName: '1', 'Starting residue of existing loop.');
cli.f(longOpt: 'final', args: 1, argName: '-1', 'Final residue of an existing loop.');
cli.mcn(longOpt: 'mcStepFreq', args: 1, argName: '10', 'Number of MD steps between Monte-Carlo protonation changes.')

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Starting and ending loop atoms.
if (options.s && options.f) {
    loopStart = Integer.parseInt(options.s);
    loopStop = Integer.parseInt(options.f);
} else if (options.s || options.f) {
    logger.info("Starting atom and final atom numbers are need to use this option.");
}

// Set Monte Carlo step frequency
if (options.mcn) {
    mcStepFrequency = Integer.parseInt(options.mcn);
}

// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

// Report interval in picoseconds.
if (options.r) {
    printInterval = Double.parseDouble(options.r);
}

// Write interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);
}

// Run OST
if (options.ost) {
    runOST = true;
}

// Gaussian bias magnitude (kcal/mol).
if (options.g) {
    biasMag = Double.parseDouble(options.g);
}

// Local minimization mode
if (options.m) {
    localMin = true;
}

// Run Simulated Annealing
if (options.sa) {
    runSimulatedAnnealing = true;
}

// Run Rotamer Optimization
if (options.rot) {
    runRotamer = true;
}

// Default
if (!(options.ost && options.sa)) {
    runOST = true;
}

// Run MC Loop Optimization
if (options.mc) {
    runMCLoop = true;
    MCLoop mcLoop;
}

// Robust Default
if (options.a) {
    runOST = true;
    runRotamer = true;
}
// Build loop with PDBFilter if an existing loop is not provided
if (!(options.s && options.f)) {
    System.setProperty("buildLoops", "true");
}

// Ideally the Initial force field will not include non-bonded iteractions,
// however, the real space term depends on vdW radii.
System.setProperty("vdwterm", "true");
System.setProperty("mpoleterm", "false");

List<String> arguments = options.arguments();
String filename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    filename = arguments.get(0);
} else {
    return cli.usage();
}

File structureFile = new File(FilenameUtils.normalize(filename));
structureFile = new File(structureFile.getAbsolutePath());
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File histogramRestart = new File(baseFilename + ".his");
File lambdaRestart = new File(baseFilename + ".lam");
File dyn = new File(baseFilename + ".dyn");

Comm world = Comm.world();
int size = world.size();
int rank = 0;
double[] energyArray = new double[world.size()];
for (int i = 0; i < world.size(); i++) {
    energyArray[i] = Double.MAX_VALUE;
}

// For a multi-process job, try to get the restart files from rank sub-directories.
if (size > 1) {
    rank = world.rank();
    File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank));
    if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
    }
    lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
    dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
    structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
}

if (!dyn.exists()) {
    dyn = null;
}

systems = open(filename);

// Set built atoms active/use flags to true (false for other atoms).
Atom[] atoms = active.getAtomArray();

/**
 * If existing loop is used, set loop atoms to match atoms built with PDBFilter.
 */
if (options.s && options.f) {
    for (int i = 0; i < atoms.length; i++) {
        Atom ai = atoms[i];
        if (!options.c || chain == ai.getChainID()) {
            if (ai.getResidueNumber() >= loopStart && ai.getResidueNumber() <= loopStop) {
                ai.setBuilt(true);
            }
        }
    }
} else {
    //create array of built residues
    ArrayList<Residue> loopResidues = new ArrayList<>();
    for (int i = 0; i < active.getChains().size(); i++) {
        ArrayList<Residue> allResidues = active.getChains()[i].getResidues();

        for (int j = 0; j < allResidues.size(); j++) {
            Residue temp = allResidues[j];
            if (temp.getBackboneAtoms().get(0).getBuilt()) {
                loopResidues.add(allResidues[j]);
            }
        }
    }
    loopStart = loopResidues.get(0).getResidueNumber();
    loopStop = loopResidues.get(0).getResidueNumber();
    for (int i = 0; i < loopResidues.size(); i++) {
        if (loopResidues.get(i).getChainID() == loopResidues.get(0).getChainID()) {
            if (loopStop + 1 == loopResidues.get(i).getResidueNumber()) {
                loopStop = loopResidues.get(i).getResidueNumber();
            }
        }
    }
}

// Get a reference to the first system's ForceFieldEnergy.
ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
forceFieldEnergy.setPrintOnFailure(false, false);

// Configure built atoms.
logger.info("\n Built Atoms:");
for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setActive(true);
        ai.setUse(true);
        logger.info(String.format(" %d %s", ai.getResidueNumber(), ai.getAtomType().toString()));
    } else {
        ai.setActive(false);
        ai.setUse(true);
    }
}

// If this is a multi-process job, set the structure file to come from the subdirectory.
if (size > 1) {
    active.setFile(structureFile);
}

DiffractionFile diffractionFile = null;
DiffractionData diffractionData = null;
// Set up real space map data (can be multiple files)
List mapFiles = new ArrayList();
int nDiffractionData = 0;
if (arguments.size() > 1) {
    String dataFileName = arguments.get(1);
    diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
    diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();
} else {
    logger.info(" Pleases specify diffraction data.");
    return;
}

logger.info("\n Running minimize on built atoms of " + active.getName());
logger.info(" RMS gradient convergence criteria: " + eps);

RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData, RefinementMode.COORDINATES);
if (localMin) {
    runOST = false;
} else {
    // Initial minimization before OST.
    refinementMinimize.minimize(eps);
    energy();
}

if (runOST) {

    // Run OST.
    System.setProperty("vdwterm", "true");
    System.setProperty("vdw-cutoff", "7.0");
    System.setProperty("mpoleterm", "true");
    System.setProperty("polarization", "none");
    System.setProperty("intramolecular-softcore", "true");
    System.setProperty("intermolecular-softcore", "true");
    System.setProperty("lambdaterm", "true");
    System.setProperty("torsion-lambdaterm", "true");
    System.setProperty("ligand-vapor-elec", "false");
    System.setProperty("lambda-bias-cutoff", "3");
    if (options.g) {
        System.setProperty("bias-gaussian-mag", String.format("%f", biasMag));
    } else {
        System.setProperty("bias-gaussian-mag", "0.002");
    }
    System.setProperty("lambda-bin-width", "0.01");

    // Set the thermostat time constant (in psec) to the the time step (i.e. to give velocity rescaling).
    System.setProperty("tau-temperature", String.format("%f", timeStep * 1.0e-3));
    System.setProperty("integrate", "respa");

    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        if (ai.getBuilt()) {
            ai.setApplyLambda(true);
        } else {
            ai.setApplyLambda(false);
        }
    }

    // Set atoms that are 1-2 and 1-3 to Built atoms to be active.
    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        if (ai.getBuilt()) {
            ArrayList<Angle> angles = ai.getAngles();
            for (int j = 0; j < angles.size; j++) {
                Atom[] angleAtoms = angles[j].getAtomArray();
                for (int k = 0; k < angleAtoms.length; k++) {
                    angleAtoms[k].setActive(true);
                }
            }
        }
    }

    forceFieldEnergy = ForceFieldEnergy.energyFactory(active);
    forceFieldEnergy.setPrintOnFailure(false, false);
    forceFieldEnergy.setLambda(lambda);
    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    forceFieldEnergy.getCrystal().setSpecialPositionCutoff(0.0);

    diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();

    RefinementEnergy refinementEnergy = new RefinementEnergy(diffractionData, RefinementMode.COORDINATES);
    refinementEnergy.setLambda(lambda);

    energy();

    boolean asynchronous = true;
    Potential ost = new OrthogonalSpaceTempering(refinementEnergy, refinementEnergy,
                lambdaRestart, histogramRestart, active.getProperties(),
                (temperature), timeStep, printInterval, saveInterval, asynchronous, sh);
    ost.setLambda(lambda);
    ost.setThetaMass(5.0e-19);
    ost.getOptimizationParameters().setOptimization(true, active);
    // Create the MolecularDynamics instance.
    MolecularDynamics molDyn = new MolecularDynamics(active, ost, active.getProperties(),
            null, thermostat, integrator);

    if (runMCLoop) {
        mcLoop = new MCLoop(active, mcStepFrequency, molDyn.getThermostat(), loopStart, loopStop);
        molDyn.addMCListener(mcLoop);
        mcLoop.addMolDyn(molDyn);
        mcLoop.addLambdaInterface(ost.getLambdaInterface());
        mcLoop.setIterations(20);
    }

    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
            fileType, restartInterval, dyn);

    logger.info("Obtaining low energy coordinates");
    OptimizationParameters opt = ost.getOptimizationParameters();
    double[] lowEnergyCoordinates = opt.getOptimumCoordinates();
    double currentOSTOptimum = opt.getOptimumEnergy();
    if (lowEnergyCoordinates != null) {
        forceFieldEnergy.setCoordinates(lowEnergyCoordinates);
    } else {
        logger.info("OST stage did not succeed in finding a loop.");
        loopBuildError = true;
    }
}

if (runSimulatedAnnealing) {
    logger.severe(" Script must be re-factored for new simulated annealing API (and likely also to new Groovy script style!)")
    // Minimize with vdW.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "false");

    energy = ForceFieldEnergy.energyFactory(active);
    energy.setPrintOnFailure(false, false);

    diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();

    refinementMinimize = new RefinementMinimize(diffractionData, RefinementMode.COORDINATES);
    refinementMinimize.minimize(eps);

    // SA with vdW.
    logger.info("\n Running simulated annealing on " + active.getName());
    double[] heatUpTemperatures = [150, 250, 400, 700, 1000];
    // Number of molecular dynamics steps at each temperature.
    int steps = 267; //267 at 3
    // Time step in femtoseconds.
    timeStep = 3.0;
    // ThermostatEnum[ ADIABATIC, BERENDSEN, BUSSI ]
    thermostat = ThermostatEnum.ADIABATIC;
    // IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC]
    integrator = IntegratorEnum.STOCHASTIC;

    refinementEnergy = new RefinementEnergy(realSpaceData, RefinementMode.COORDINATES, null);
    SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, refinementEnergy,
            active.getProperties(), null, thermostat, integrator);
    simulatedAnnealing.annealToTargetValues(heatUpTemperatures, steps, timeStep);

    double[] annealingTargetTemperatures = [1000, 800, 600, 500, 400, 300];
    steps = 800; //800 at 3
    simulatedAnnealing.annealToTargetValues(annealingTargetTemperatures, steps, timeStep);
}

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    ai.setUse(true);
    ai.setApplyLambda(false);
}

if (!loopBuildError) {
    // Optimize with the full AMOEBA potential energy.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "true");
    System.setProperty("polarization", "direct");
    System.setProperty("intramolecularSoftcore", "false");
    System.setProperty("intermolecularSoftcore", "false");
    System.setProperty("lambdaterm", "false");
    System.setProperty("torsion-lambdaterm", "false");

    forceFieldEnergy = ForceFieldEnergy.energyFactory(active);
    forceFieldEnergy.setPrintOnFailure(false, false);
    diffractionData = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionFile);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();

    refinementMinimize = new RefinementMinimize(diffractionData, RefinementMode.COORDINATES);
    refinementMinimize.minimize(eps);

    energy();
    if (size > 1) {
        structureFile = new File("postOST." + String.format("%d", world.rank()) + "." + structureFile.getName());
    } else {
        structureFile = new File("postOST." + structureFile.getName());
    }
    saveAsPDB(structureFile);
}

if (runOST && size > 1) {

    DoubleBuf receiveBuffer = DoubleBuf.buffer(energyArray);
    if (!(world.rank() == 0)) {
        world.receive(world.rank() - 1, receiveBuffer);
        energyArray[world.rank() - 1] = receiveBuffer.get(world.rank() - 1);
    }
    if (!loopBuildError) {
        energyArray[world.rank()] = active.getPotentialEnergy().getTotalEnergy();
    } else {
        energyArray[world.rank()] = Double.MAX_VALUE;
    }
    if (world.rank() < world.size() - 1) {
        DoubleBuf sendBuffer = DoubleBuf.buffer(energyArray);
        world.send(world.rank() + 1, sendBuffer);
    }
    world.barrier();
    if (world.rank() == world.size() - 1) {
        for (int i = 0; i < world.size(); i++) {
            String resultFileName = "Loop.txt";
            File rankAndEnergyFile = new File(resultFileName);
            BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(rankAndEnergyFile, true));

            String rankString = Integer.toString(i);
            String energyString = Double.toString(energyArray[i]);

            bufferedWriter.write(rankString + ":" + energyString);
            bufferedWriter.newLine();
            bufferedWriter.flush();
        }
    }
    saveAsPDB(structureFile);
}
