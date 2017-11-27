/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ffx.algorithms

//Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

//FFX Imports
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.OpenMMMolecularDynamics;
import ffx.algorithms.Barostat;
import ffx.algorithms.MolecularDynamics;
import ffx.potential.OpenMMForceFieldEnergy;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.MolecularAssembly
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.ForceFieldEnergy;
import ffx.numerics.Potential;
import ffx.crystal.CrystalPotential

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double timeStep = 1.0;

int intervalSteps = 100;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 100.0;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 2.0;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ANDERSON ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET] for regular MolecularDynamics
Integrators integrator = null;

// Integrator string, sets integrator for OpenMMMolecularDynamics
//String stringIntegrator = "VERLET";

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 2.0;

// Coefficient of Friction for Langevin itegrator
double frictionCoeff = 91.0;

// Collision Frequency for Anderson Thermostat used in conjucntion with Verlet Integrator
double collisionFreq = 0.01;

List<CoordRestraint> restraints = null;

// Pressure in atm; invalid (such as negative) values mean do not run NPT.
double pressure = -1.0;

// Set min and max density and set the maximum MC move size
int meanInterval = 10;
double minDensity = 0.5;
double maxDensity = 1.5;
double maxSideMove = 0.25;
double maxAngleMove = 1.0;

int numThreads = 1;

// File type of snapshots.
String fileType = "PDB";

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc md [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Anderson]');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
cli.z(longOpt: 'intervalSteps', args: 1, argName:'100', 'Number of steps for each MD trajectory (fsec)');
//cli.si(longOpt:'integrate', args:1, argName:'Verlet', 'Integrator: [Langevin / Verlet]');
cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.a(longOpt:'pressure', args:1, argName:'Disabled', 'If numeric argument, then sets NPT pressure, else do not run NPT.')
cli.ld(longOpt:'minDensity', args:1, argName:'0.5', 'Minimum density allowed by the barostat.');
cli.hd(longOpt:'maxDensity', args:1, argName:'1.5', 'Maximum density allowed by the barostat.');
cli.sm(longOpt:'maxSideMove', args:1, argName:'0.25', 'Maximum side move allowed by the barostat.');
cli.am(longOpt:'maxAngleMove', args:1, argName:'1.0', 'Maximum angle move allowed by the barostat.');
cli.mi(longOpt:'meanInterval', args:1, argName:'10', 'Mean number of MD steps between applications of the barostat.');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
cli.cf(longOpt:'Coefficient of friction', args:1, argName:'1.0', 'Coefficient of friction to be used with Langevin and Brownian integrators')
cli.q(longOpt:'Collision frequency in 1/ps', args:1, argName:'10.0', 'Collision frequency to be set when Anderson Thermostat is created: Verlet integrator requires this')
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

if (options.z) {
    intervalSteps = Integer.parseInt(options.z);
}

// Integrator for OpenMMMolecularDynamics
if (options.si) {
    stringIntegrator = options.i.toUpperCase();
}

// Regular MolecularDynamics integrator option, not used for OpenMMMolecularDynamics
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

// Report interval in picoseconds.
if (options.l) {
    printInterval = Double.parseDouble(options.l);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Pressure in atm.
if (options.a) {
    pressure = Double.parseDouble(options.a);
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

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

// Write snapshot interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Write dyn interval in picoseconds
if (options.s) {
    restartFrequency = Double.parseDouble(options.s);
}

//
if (options.f) {
    fileType = options.f.toUpperCase();
}

// Coefficient of friction to be used with Langevin and Brownian integrators
if (options.cf) {
    frictionCoeff = Double.parseDouble(options.cf);
}

//Collision frequency (in 1/ps) set when using Anderson Thermostat: Verlet integrator requires this
if (options.q) {
    collisionFreq = Double.parseDouble(options.q);
}


List<String> arguments = options.arguments();
String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    modelfilename = arguments.get(0);
    open(modelfilename);
} else if (active == null) {
    return cli.usage();
} else {
    modelfilename = active.getFile();
}

ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
logger.info(" Starting energy (before .dyn restart loaded):");
boolean updatesDisabled = active.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
if (updatesDisabled) {
    logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
}
double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
forceFieldEnergy.getCoordinates(x);
forceFieldEnergy.energy(x, true);

if (pressure > 0) {
    if (potential instanceof OpenMMForceFieldEnergy) {
        logger.warning(" NPT with OpenMM acceleration is still experimental and may not function correctly.");
    }
    logger.info(String.format(" Running NPT dynamics at pressure %7.4g", pressure));
    MolecularAssembly molecAssem = active; // Why this is needed is utterly inexplicable to me.
    CrystalPotential cpot = (CrystalPotential) potential;
    Barostat barostat = new Barostat(molecAssem, cpot);
    barostat.setPressure(pressure);
    barostat.setMaxDensity(maxDensity);
    barostat.setMinDensity(minDensity);
    barostat.setMaxSideMove(maxSideMove);
    barostat.setMaxAngleMove(maxAngleMove);
    barostat.setMeanBarostatInterval(meanInterval);
    potential = barostat;
}

logger.info("\n Running molecular dynmaics on " + modelfilename);

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

logger.info(" about to create dynamics object");

//if (forceFieldEnergy instanceof OpenMMForceFieldEnergy) {
MolecularDynamics moldyn = MolecularDynamics.dynamicsFactory(active, forceFieldEnergy, active.getProperties(), sh, thermostat, integrator)
if (moldyn instanceof OpenMMMolecularDynamics){
    moldyn.setRestartFrequency(restartFrequency);
    //moldyn.setIntegratorString(stringIntegrator);
    //moldyn.setFrictionCoefficient(frictionCoeff);
    //moldyn.setCollisionFrequency(collisionFreq);
    moldyn.setIntervalSteps(intervalSteps);
    moldyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
}
else{
    logger.info(" in set up for java md");
    MolecularDynamics molDyn_java = new MolecularDynamics(active, forceFieldEnergy, active.getProperties(), sh, thermostat, integrator);
    molDyn_java.setFileType(fileType);
    molDyn_java.setRestartFrequency(restartFrequency);
    molDyn_java.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
}
//}
