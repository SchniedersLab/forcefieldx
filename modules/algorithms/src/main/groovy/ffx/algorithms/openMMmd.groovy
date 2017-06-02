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
import ffx.potential.OpenMMForceFieldEnergy;
import ffx.potential.nonbonded.CoordRestraint;

import ffx.numerics.Potential;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.ForceFieldEnergy;
import ffx.algorithms.MolecularDynamics;

// Number of molecular dynamics steps
int nSteps = 2000000;

// Time step in femtoseconds.
double timeStep = 1.0;

int intervalSteps = 100;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.025;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 2000.0;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ANDERSON ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET] for regular MolecularDynamics
Integrators integrator = null;

// Integrator string, sets integrator for OpenMMMolecularDynamics
String stringIntegrator = "VERLET";

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 2000.0;

double frictionCoeff = 91.0;

double collisionFreq = 0.1;

List<CoordRestraint> restraints = null;

int numThreads = 1;

// File type of snapshots.
String fileType = "PDB";

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc md [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Anderson]');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
cli.si(longOpt:'integrate', args:1, argName:'Verlet', 'Integrator: [Brownian / Langevin / Verlet / Custom / Compound]');
cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
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

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Write dyn interval in picoseconds
if (options.s) {
    restartFrequency = Double.parseDouble(options.s);
}

//
if (options.f) {
    fileType = options.f.toUpperCase();
}
// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

// Report interval in picoseconds.
if (options.l) {
    printInterval = Double.parseDouble(options.l);
}

// Write snapshot interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
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

Potential potential = active.getPotentialEnergy();
logger.info(" Starting energy (before .dyn restart loaded):");
boolean updatesDisabled = active.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
if (updatesDisabled) {
    logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
}
double[] x = new double[potential.getNumberOfVariables()];
potential.getCoordinates(x);

potential.energy(x, true);

logger.info("\n Running molecular dynmaics on " + modelfilename);

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}


ForceFieldEnergy forceFieldEnergy = ForceFieldEnergy.energyFactory(active, restraints, numThreads);


if (forceFieldEnergy instanceof OpenMMForceFieldEnergy) {
    
    MolecularDynamics moldyn = MolecularDynamics.dynamicsFactory(active, forceFieldEnergy, active.getProperties(), sh, thermostat, integrator)

    if (moldyn instanceof OpenMMMolecularDynamics){
        moldyn.setRestartFrequency(restartFrequency);
        moldyn.setIntegratorString(stringIntegrator);
        moldyn.setFrictionCoefficient(frictionCoeff);
        moldyn.setCollisionFrequency(collisionFreq);
        moldyn.start(nSteps, intervalSteps, temperature, saveInterval, timeStep, dyn, initVelocities);
    }
} 