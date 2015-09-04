
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

// LOOP BUILDER

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports

import ffx.algorithms.Minimize;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.SimulatedAnnealing;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// Default convergence criteria.
double eps = 1.0;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Time step in femtoseconds.
double timeStep = 2.5;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100.0;

// Frequency to write out restart information in picoseconds.
double restartInterval = 1.0;

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 10000;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;

// Integrators [ BEEMAN, RESPA, STOCHASTIC ]
Integrators integrator = Integrators.RESPA;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// File type of coordinate snapshots to write out.
String fileType = "PDB";

// Value of Lambda.
double lambda = 0.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc loopBuilder [options] <filename1>');
cli.h(longOpt:'help', 'Print this help message.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.n(longOpt:'steps', args:1, argName:'10000', 'Number of molecular dynamics steps.');
cli.d(longOpt:'dt', args:1, argName:'2.5', 'Time discretization step (fsec).');
cli.r(longOpt:'report', args:1, argName:'0.01', 'Interval to report thermodyanamics (psec).');
cli.w(longOpt:'write', args:1, argName:'100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
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

System.setProperty("buildLoops", "true");
System.setProperty("vdwterm", "false");

List<String> arguments = options.arguments();
String filename = null;
MolecularAssembly[] systems = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    filename = arguments.get(0);
    systems = open(filename);
} else {
    return cli.usage();
}

// Get a reference to the first system's ForceFieldEnergy.
ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
// Set built atoms active/use flags to true (false for other atoms).
Atom[] atoms = active.getAtomArray();
for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setActive(true);
        ai.setUse(true);
    } else {
        ai.setActive(false);
        ai.setUse(true);
    }
}

logger.info("\n Running minimize on built atoms of " + active.getName());
logger.info(" RMS gradient convergence criteria: " + eps);

// Minimization without vdW.
e = minimize(eps);

energy();

// Run OSRW with a vdW potential.
System.setProperty("vdwterm", "true");
System.setProperty("mpoleterm", "true");
System.setProperty("polarization", "none");
System.setProperty("intramolecular-softcore", "true");
System.setProperty("intermolecular-softcore", "true");
System.setProperty("lambdaterm", "true");
System.setProperty("ligand-vapor-elec","false");
System.setProperty("ewald-alpha","0.0");
System.setProperty("lambda-bias-cutoff", "3");
System.setProperty("bias-gaussian-mag", "0.01");
System.setProperty("lambda-bin-width", "0.01");

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setApplyLambda(true);
    } else {
        ai.setApplyLambda(false);
    }
}

forceFieldEnergy= new ForceFieldEnergy(active);
forceFieldEnergy.setLambda(lambda);

energy();

// Turn off checks for overlapping atoms, which is expected for lambda=0.
forceFieldEnergy.getCrystal().setSpecialPositionCutoff(0.0);
// OSRW will be configured for a single topology.
File lambdaRestart = null;
File histogramRestart = null;
boolean asynchronous = false;
boolean wellTempered = false;
OSRW osrw =  new OSRW(forceFieldEnergy, forceFieldEnergy, lambdaRestart, histogramRestart, active.getProperties(),
    temperature, timeStep, printInterval, saveInterval, asynchronous, sh, wellTempered);
osrw.setLambda(lambda);
// Create the MolecularDynamics instance.
MolecularDynamics molDyn = new MolecularDynamics(active, osrw, active.getProperties(),
    null, thermostat, integrator);
File dyn = null;
molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
    fileType, restartInterval, dyn);

energy();

if (false) {
    e = minimize(eps);

    // SA with vdW.
    logger.info("\n Running simulated annealing on " + active.getName());

    // High temperature starting point.
    double high = 1000.0;
    // Low temperature end point.
    double low = 10.0;
    // Number of annealing steps.
    int windows = 10;
    // Number of molecular dynamics steps at each temperature.
    int steps = 100;
    // Time step in femtoseconds.
    timeStep = 1.0;
    // Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
    thermostat = Thermostats.BERENDSEN;
    // Integrators [ BEEMAN, RESPA, STOCHASTIC]
    integrator = null;

    SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, forceFieldEnergy,
        active.getProperties(), null, thermostat, integrator);
    simulatedAnnealing.anneal(high, low, windows, steps, timeStep);

    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        ai.setUse(true);
    }
}

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    ai.setUse(true);
    ai.setApplyLambda(false);
}

// Optimize with the full AMOEBA potential energy.
System.setProperty("vdwterm", "true");
System.setProperty("mpoleterm", "true");
System.setProperty("polarization", "direct");
System.setProperty("intramolecularSoftcore", "false");
System.setProperty("intermolecularSoftcore", "false");
System.setProperty("lambdaterm", "false");

forceFieldEnergy = new ForceFieldEnergy(active);
e = minimize(eps);

if (false) {
    // MD, SA, OSRW and/or Side-Chain DEE
    simulatedAnnealing = new SimulatedAnnealing(active, forceFieldEnergy, active.getProperties(),
        null, thermostat, integrator);
    simulatedAnnealing.anneal(high, low, windows, steps, timeStep);
    e = minimize(0.1);
}


String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);
saveAsPDB(systems, new File(filename + ".pdb"));


