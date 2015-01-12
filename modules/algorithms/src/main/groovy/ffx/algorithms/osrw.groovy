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

// ORTHOGONAL SPACE RANDOM WALK

// Apache Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Paralle Java Imports
import edu.rit.pj.Comm;

// Force Field X Imports
//import ffx.algorithms.Barostat;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// Asychronous communication between walkers.
boolean asynchronous = false;

// Well-Tempered Metadynamics.
boolean wellTempered = false;

// First atom of the ligand.
int ligandStart = 1;

// Last atom of the ligand.
int ligandStop = -1;

// First atom of the ligand for the 2nd topology.
int ligandStart2 = 1;

// Last atom of the ligand for the 2nd topology.
int ligandStop2 = -1;

// No electrostatics for ligand 1
boolean noElec = false;

// No electrostatics for ligand 2
boolean noElec2 = false;

// Initial lambda value (0 is ligand in vacuum; 1 is ligand in PBC).
double lambda = 0.0;

// Friction coefficient on the Lambda particle.
double friction = 1.0e-18;

// Mass of the Lambda particle.
double mass = 1.0e-18;

// Interval between OSRW recursion kernel counts.
int count = 10;

// Magnitude of each Gaussian bias (kcal/mole).
double biasMag = 0.002;

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 100000000;

// Number of equilibration steps: default is 1 picosecond.
int eSteps = 1000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 1.0;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100.0;

// File type of coordinate snapshots to write out.
String fileType = "XYZ";

// Frequency to write out restart information in picoseconds.
double restartInterval = 0.1;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC ]
Integrators integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Reset OSRW statistics
boolean resetStatistics = false;

// Constant pressure
boolean NPT = false;
int meanInterval = 10;
double minDensity = 0.5;
double maxDensity = 1.5;
double maxSideMove = 0.25;
double maxAngleMove = 1.0;

// Write traversal snapshots
boolean writeTraversals = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc osrw [options] <filename> [filename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'async', args:0, 'Walker communication is asynchronous.');
cli.wt(longOpt:'wellTempered', args:0, 'Use Well-Tempered Metadynamics.');
cli.n(longOpt:'steps', args:1, argName:'10000000', 'Number of molecular dynamics steps.');
cli.q(longOpt:'equilibrate', args:1, argName:'1000', 'Equilibration steps prior to OSRW counts.');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization step (fsec).');
cli.r(longOpt:'report', args:1, argName:'1.0', 'Interval to report thermodyanamics (psec).');
cli.w(longOpt:'write', args:1, argName:'100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.i(longOpt:'integrator', args:1, argName:'Beeman', 'Integrator: [Beeman / Respa / Stochastic]');
cli.s(longOpt:'start', args:1, argName:'1', 'Starting ligand atom.');
cli.s2(longOpt:'start2', args:1, argName:'1', 'Starting ligand atom for the 2nd topology.');
cli.f(longOpt:'final', args:1, argName:'-1', 'Final ligand atom.');
cli.f2(longOpt:'final2', args:1, argName:'-1', 'Final ligand atom for the 2nd topology.');
cli.l(longOpt:'lambda', args:1, argName:'0.0', 'Initial lambda value (> 1.0 distributes lambda across walkers)');
cli.c(longOpt:'count', args:1, argName:'10', 'Time steps between OSRW counts.');
cli.g(longOpt:'bias', args:1, argName:'0.002', 'Gaussian bias magnitude (kcal/mol).');
cli.m(longOpt:'mass', args:1, argName:'1e-18', 'Lambda particle mass.');
//cli.p(longOpt:'npt', args:0, 'Constant pressure MD (1 atm).');
cli.e(longOpt:'elec', args:0, 'No electrostatics on ligand 1.');
cli.e2(longOpt:'elec2', args:0, 'No electrostatics on ligand 2.');
cli.x(longOpt:'friction', args:1, argName:'1e-18', 'Lambda particle friction.');
cli.W(longOpt:'notraversals', args:0, 'Don\'t write out lambda-traversal snapshots.');
//cli.ld(longOpt:'minDensity', args:1, argName:'0.5', 'Minimum density allowed by the barostat.');
//cli.hd(longOpt:'maxDensity', args:1, argName:'1.5', 'Maximum density allowed by the barostat.');
//cli.sm(longOpt:'maxSideMove', args:1, argName:'0.25', 'Maximum side move allowed by the barostat.');
//cli.am(longOpt:'maxAngleMove', args:1, argName:'1.0', 'Maximum angle move allowed by the barostat.');
//cli.mi(longOpt:'meanInterval', args:1, argName:'10', 'Mean number of MD steps between applications of the barostat.');
cli.rt(longOpt:'reset', args:0, 'Reset OSRW histogram once, when lambda reaches 0.99.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() > 2) {
    return cli.usage();
}

// Read in command line file.
String filename = arguments.get(0);

// Asynchronous?
if (options.a) {
    asynchronous = true;
}

// Well-Tempered?
if (options.wt) {
    wellTempered = true;
}

// Constant Pressue?
if (options.p) {
    NPT = true;
}

// Reset OSRW statistics when L=0.99 is reached the first time.
if (options.rt) {
    resetStatistics = true;
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load the number of molecular dynamics steps.
if (options.q) {
    eSteps = Integer.parseInt(options.q);
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

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

// Starting ligand atom.
if (options.s) {
    ligandStart = Integer.parseInt(options.s);
}

// Final ligand atom.
if (options.f) {
    ligandStop = Integer.parseInt(options.f);
}

// No electrostatics on ligand 1.
if (options.e) {
    noElec = true;
}

// Starting ligand atom for the 2nd topology.
if (options.s2) {
    ligandStart2 = Integer.parseInt(options.s2);
}

// Final ligand atom for the 2nd topology.
if (options.f2) {
    ligandStop2 = Integer.parseInt(options.f2);
}

// No electrostatics on ligand 2.
if (options.e2) {
    noElec2 = true;
}

// Starting lambda value.
if (options.l) {
    lambda = Double.parseDouble(options.l);
}

// Number of time steps between counts.
if (options.c) {
    count = Integer.parseInt(options.c);
}

// Gaussian bias magnitude (kcal/mol).
if (options.g) {
    biasMag = Double.parseDouble(options.g);
}

// Lambda particle mass.
if (options.m) {
    mass = Double.parseDouble(options.m);
}

// Lambda particle friction.
if (options.x) {
    friction = Double.parseDouble(options.x);
}

// Minimum density
if (options.ld) {
    minDensity = Double.parseDouble(options.ld);
}
// Maximum density
if (options.hd) {
    maxDensity = Double.parseDouble(options.hd);
}
// Max side move
if (options.sm) {
    maxSideMove = Double.parseDouble(options.sm);
}
// Max angle move
if (options.am) {
    maxAngleMove = Double.parseDouble(options.am);
}
// Max angle move
if (options.mi) {
    meanInterval = Integer.parseInt(options.mi);
}
// Traversal snapshots
if (options.W) {
    writeTraversals = false;
}

println("\n Running Orthogonal Space Random Walk on " + filename);

File structureFile = new File(FilenameUtils.normalize(filename));
structureFile = new File(structureFile.getAbsolutePath());
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File histogramRestart = new File(baseFilename + ".his");
File lambdaOneFile = null;
File lambdaZeroFile = null;
File lambdaRestart = null;
File dyn = null;

Comm world = Comm.world();
int size = world.size();
int rank = 0;

// For a multi-process job, try to get the restart files from rank sub-directories.
if (size > 1) {
    rank = world.rank();
    File rankDirectory = new File(structureFile.getParent() + File.separator
        + Integer.toString(rank));
    if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
    }
    lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
    dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
    lambdaOneFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam1");
    lambdaZeroFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam0");
    structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
} else {
    // For a single process job, try to get the restart files from the current directory.
    lambdaRestart = new File(baseFilename + ".lam");
    dyn = new File(baseFilename + ".dyn");
    lambdaOneFile = new File(baseFilename + ".lam1");
    lambdaZeroFile = new File(baseFilename + ".lam0");
}

if (!dyn.exists()) {
    dyn = null;
}

// Turn on computation of lambda derivatives.
System.setProperty("lambdaterm", "true");

// Relative free energies via the DualTopologyEnergy class require different
// default OSRW parameters than absolute free energies.
if (arguments.size() == 2) {
    // Condensed phase polarization is evaluated over the entire range.
    System.setProperty("polarization-lambda-start","0.0");
    // Polarization energy is not scaled individually by lambda, but
    // along with the overall potential energy of a topology.
    System.setProperty("polarization-lambda-exponent","0.0");
    // Ligand vapor electrostatics are not calculated. This cancels when the
    // difference between protein and water environments is considered.
    System.setProperty("ligand-vapor-elec","false");
    // Condensed phase polarization, without the ligand present, is unecessary.
    System.setProperty("no-ligand-condensed-scf","false");
}

// Open the first system
open(filename);
// If this is a multi-process job, set the structure file to come from the subdirectory.
if (size > 1) {
    active.setFile(structureFile);
}

// Get a reference to the first system's ForceFieldEnergy and atom array.
ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
    ai.print();
}

// Turn off checks for overlapping atoms, which is expected for lambda=0.
energy.getCrystal().setSpecialPositionCutoff(0.0);
// Turn off ligand electrostatics if requested.
energy.setNoSoftCoreElectrostatics(noElec);
// OSRW will be configured for either single or dual topology.
OSRW osrw = null;
// Save a reference to the first topology.
topology1 = active;

if (arguments.size() == 1) {
    // Check for constant pressure
    if (NPT) {
        //        // Create a barostat.
        //        Barostat barostat = new Barostat(active);
        //        barostat.setMaxdUdL(1000.0);
        //        barostat.setMaxDensity(maxDensity);
        //        barostat.setMinDensity(minDensity);
        //        barostat.setMaxSideMove(maxSideMove);
        //        barostat.setMaxAngleMove(maxAngleMove);
        //        barostat.setMeanBarostatInterval(meanInterval);
        //
        //        // Create the OSRW instance.
        //        osrw = new OSRW(energy, barostat, lambdaRestart, histogramRestart, active.getProperties(),
        //            temperature, timeStep, printInterval, saveInterval, asynchronous, sh, wellTempered);
        //        osrw.setResetStatistics(resetStatistics);
        if (writeTraversals) {
            osrw.setTraversalOutput(lambdaOneFile, topology1, lambdaZeroFile, topology1);
        }
    } else {
        // Wrap the single topology ForceFieldEnergy inside an OSRW instance.
        osrw = new OSRW(energy, energy, lambdaRestart, histogramRestart, active.getProperties(),
            temperature, timeStep, printInterval, saveInterval, asynchronous, sh, wellTempered);
        osrw.setResetStatistics(resetStatistics);
        if (writeTraversals) {
            osrw.setTraversalOutput(lambdaOneFile, topology1, lambdaZeroFile, topology1);
        }
    }
} else {
    // Open the 2nd topology.
    filename = arguments.get(1);
    open(filename);
    // If this is a multi-process job, set the structure file to come from the subdirectory.
    if (size > 1) {
        active.setFile(structureFile);
    }
    energy = active.getPotentialEnergy();
    atoms = active.getAtomArray();
    // Apply the ligand atom selection for the 2nd topology.
    for (int i = ligandStart2; i <= ligandStop2; i++) {
        Atom ai = atoms[i - 1];
        ai.setApplyLambda(true);
        ai.print();
    }
    // Save a reference to the second topology.
    topology2 = active;
    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff(0.0);
    // Turn off ligand electrostatics if requested.
    energy.setNoSoftCoreElectrostatics(noElec2);
    // Create the DualTopology potential energy.
    DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topology1, active);
    // Wrap the DualTopology potential energy inside an OSRW instance.
    osrw = new OSRW(dualTopologyEnergy, dualTopologyEnergy, lambdaRestart, histogramRestart, active.getProperties(),
        temperature, timeStep, printInterval, saveInterval, asynchronous, sh, wellTempered);
    osrw.setResetStatistics(resetStatistics);
    if (writeTraversals) {
        osrw.setTraversalOutput(lambdaOneFile, topology1, lambdaZeroFile, topology2);
    }
}

// Apply the command line lambda value if a lambda restart file does not exist.
if (!lambdaRestart.exists()) {
    if (lambda >= 0.0 && lambda <= 1.0) {
        osrw.setLambda(lambda);
        logger.info(String.format(" Setting lambda to %5.3f.", lambda));
    } else {
        if (size > 1) {
            dL = 1.0 / (size - 1.0);
            lambda = rank * dL;
            if (lambda > 1.0) {
                lambda = 1.0;
            }
            if (lambda < 0.0) {
                lambda = 0.0;
            }
            logger.info(String.format(" Setting lambda to %5.3f.", lambda));
            osrw.setLambda(lambda);
        } else {
            lambda = 0.5;
            logger.info(String.format(" Setting lambda to %5.3f", lambda));
            osrw.setLambda(lambda);
        }
    }
}

// Apply the command line OSRW values if a histogram restart file does not exist.
if (!histogramRestart.exists()) {
    osrw.setThetaFrication(friction);
    osrw.setThetaMass(mass);
    osrw.setCountInterval(count);
    osrw.setBiasMagnitude(biasMag);
}

// Create the MolecularDynamics instance.
MolecularDynamics molDyn = new MolecularDynamics(topology1, osrw, topology1.getProperties(), null, thermostat, integrator);

// Start sampling.
if (eSteps > 0) {
    logger.info(" Beginning equilibration");
    osrw.setPropagateLambda(false);
    molDyn.dynamic(eSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
    logger.info(" Beginning OSRW sampling");
    osrw.setPropagateLambda(true);
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, false,
        fileType, restartInterval, dyn);
} else {
    logger.info(" Beginning OSRW sampling without equilibration");
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
        fileType, restartInterval, dyn);
}

