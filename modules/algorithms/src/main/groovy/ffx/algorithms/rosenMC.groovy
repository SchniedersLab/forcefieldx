/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

// MOLECULAR & STOCHASTIC DYNAMICS

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.MonteCarloListener;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.mc.RosenbluthOBMC;
import ffx.algorithms.mc.RosenbluthCBMC;
import ffx.potential.bonded.Residue;

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET]
Integrators integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 1000;

// File type of snapshots.
String fileType = "PDB";

// Parameters to RRMC.
List<Residue> targets = new ArrayList<>();
int resNum = 1;
int mcFrequency = 50;
int trialSetSize = 10;
boolean writeSnapshots = false;
boolean dynamics = false;
boolean bias = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc md [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
cli.rmcR(longOpt:'residue', args:1, argName:'1', 'RRMC target residue number.');
cli.rmcF(longOpt:'mcFreq', args:1, argName:'50', 'RRMC frequency.');
cli.rmcK(longOpt:'trialSetSize', args:1, argName:'10', 'Size of RRMC trial sets.');
cli.rmcW(longOpt:'writeSnapshots', args:1, argName:'false', 'Output PDBs of trial sets and orig/proposed conformations.');
cli.rmcD(longOpt:'dynamics', args:1, argName:'false', 'Skip molecular dynamics; do only Monte Carlo moves.')
cli.rmcB(longOpt:'bias', args:1, argName:'true', 'For validation. Skips biasing.');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

if (options.rmcB) {
    bias = Boolean.parseBoolean(options.rmcB);
}

if (options.rmcD) {
    dynamics = Boolean.parseBoolean(options.rmcD);
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

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

if (options.rmcR) {
    resNum = Integer.parseInt(options.rmcR);
}

if (options.rmcF) {
    mcFrequency = Integer.parseInt(options.rmcF);
}

if (options.rmcK) {
    trialSetSize = Integer.parseInt(options.rmcK);
}

if (options.rmcW) {
    writeSnapshots = Boolean.parseBoolean(options.rmcW);
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

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

if (!dynamics) {
    int numAccepted = 0;
    logger.info("\n Running CBMC on " + modelfilename);
    System.setProperty("forcefield","AMOEBA_PROTEIN_2013");
    mcFrequency = 1;
    targets.add(active.getChains()[0].getResidues().get(resNum));
    MonteCarloListener rrmc = new RosenbluthCBMC(active, active.getPotentialEnergy(), null,
        targets, mcFrequency, trialSetSize, writeSnapshots);
    while (numAccepted < nSteps) {
        if (bias) {
            boolean accepted = rrmc.cbmcStep();
            if (accepted) {
                numAccepted++;
            }
        } else {
            boolean accepted = rrmc.controlStep();
            if (accepted) {
                numAccepted++;
            }
        }
    }
    return;
}

logger.info("\n Running dynamics with CBMC on " + modelfilename);
MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), active.getProperties(), sh, thermostat, integrator);
molDyn.setFileType(fileType);
molDyn.setRestartFrequency(restartFrequency);

targets.add(active.getChains()[0].getResidues().get(resNum));
MonteCarloListener rrmc = new RosenbluthCBMC(active, active.getPotentialEnergy(), molDyn.getThermostat(),
    targets, mcFrequency, trialSetSize, writeSnapshots);
molDyn.addMCListener(rrmc);

molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
