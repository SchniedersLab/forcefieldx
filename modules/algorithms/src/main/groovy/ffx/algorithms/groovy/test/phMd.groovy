//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.groovy.test

import org.apache.commons.io.FilenameUtils

import groovy.cli.picocli.CliBuilder

import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.integrators.Integrator
import ffx.algorithms.dynamics.integrators.IntegratorEnum
import ffx.algorithms.dynamics.thermostats.Thermostat
import ffx.algorithms.dynamics.thermostats.ThermostatEnum
import ffx.algorithms.ph.PhMD
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3
import ffx.potential.extended.TitrationUtils

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

// ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
ThermostatEnum thermostat = null;

// IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET]
IntegratorEnum integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 1000;

// File type of snapshots.
String fileType = "PDB";

// Monte-Carlo step frequencies for titration and rotamer moves.
int mcStepFrequency = 10;
int rotamerStepFrequency = 0;

// Simulation pH
double pH = 7.4;

// Single-residue titration option.
Character chainID = ' ';
int resID = -1;
List<String> resList = new ArrayList<>();
double window = 2.0;
boolean dynamics = true;
boolean titrateTermini = false;
boolean discreteOnly = true;
PhMD.Distribution distribution = PhMD.Distribution.DISCRETE;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
// Holy mother of god.
def cli = new CliBuilder(usage: ' ffxc test.phMD [options] <filename>');
cli.h(longOpt: 'help', 'Print this message.');
cli.d(longOpt: 'distribution', args: 1, argName: 'Discrete', 'Distribution: [DISCRETE / CONTINUOUS]');
cli.b(longOpt: 'thermostat', args: 1, argName: 'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.dt(longOpt: 'dt', args: 1, argName: '1.0', 'Time discretization (fsec).');
cli.i(longOpt: 'integrate', args: 1, argName: 'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt: 'log', args: 1, argName: '0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt: 'steps', args: 1, argName: '1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt: 'polarization', args: 1, argName: 'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt: 'temperature', args: 1, argName: '298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt: 'save', args: 1, argName: '0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt: 'restart', args: 1, argName: '0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt: 'file', args: 1, argName: 'PDB', 'Choose file type to write to [PDB/XYZ]');
cli.ra(longOpt: 'resAll', 'Titrate all residues.');
cli.rl(longOpt: 'resList', args: 1, 'Titrate a list of residues (eg A4.A8.B2.B34)');
cli.rn(longOpt: 'resName', args: 1, 'Titrate a list of residue names (eg "LYS,TYR,HIS")');
cli.rw(longOpt: 'resWindow', args: 1, 'Titrate all residues with intrinsic pKa within [arg] units of simulation pH.');
cli.pH(longOpt: 'pH', args: 1, argName: '7.4', 'Constant simulation pH.');
cli.mc(longOpt: 'mc-frequency', args: 1, argName: '10', '[DISCRETE only] Number of MD steps between Monte-Carlo attempts.')
cli.nomd(longOpt: 'no-dynamics', 'Testing; no movement.');
cli.noref(longOpt: 'no-references', 'Testing; zero out reference energies.')

def options = cli.parse(args);
if (options.h) {
    return cli.usage();
}

/* Check for missing or inconsistent flags. */
if ((options.rw && (options.ra || options.rl)) || (options.ra && options.rl)) {
    return cli.usage();
    logger.info(" Must specify one of the following: -ra, -rl, or -rw.");
}
if (!options.ra && !options.rl && !options.rw && !options.rn) {
    return cli.usage();
    logger.info(" Must specify one of the following: -ra, -rl, -rn, or -rw.");
}
if (!options.pH) {
    return cli.usage();
    logger.info(" Must specify a solution pH.");
}

/* Load command-line parameters. */
if (options.d) {
    if (options.d.toUpperCase().startsWith("C")) {
        distribution = PhMD.Distribution.CONTINUOUS;
    }
}
if (options.noref) {
    System.setProperty("cphmd-zeroReferences", "true");
}
if (options.tt) {
    titrateTermini = true;
    System.setProperty("cphmd-termini", "true");
}
if (options.nomd) {
    dynamics = false;
}
if (options.rw) {
    window = Double.parseDouble(options.rw);
}
if (options.mc) {
    mcStepFrequency = Integer.parseInt(options.mc);
}
if (options.mcr) {
    rotamerStepFrequency = Integer.parseInt(options.mcr);
}
if (options.pH) {
    pH = Double.parseDouble(options.pH);
}
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}
if (options.s) {
    restartFrequency = Double.parseDouble(options.s);
}
if (options.f) {
    fileType = options.f.toUpperCase();
}
if (options.dt) {
    timeStep = Double.parseDouble(options.dt);
}
if (options.l) {
    printInterval = Double.parseDouble(options.l);
}
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}
if (options.t) {
    temperature = Double.parseDouble(options.t);
}
if (options.p) {
    System.setProperty("polarization", options.p);
}

if (options.b) {
    thermostat = Thermostat.parseThermostat(options.b)
}

if (options.i) {
    integrator = Integrator.parseIntegrator(options.i)
}

List<String> arguments = options.arguments();
if (arguments == null || arguments.size() != 1) {
    return usage();
}
String filename = arguments.get(0);
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

/* Setup system properties before loading input file. */
TitrationUtils.initDiscountPreloadProperties();
MolecularAssembly mola = TitrationUtils.openFullyProtonated(new File(filename));

List<Residue> titrating = null;
if (options.ra) {
    titrating = TitrationUtils.chooseTitratables(mola);
} else if (options.rl) {
    titrating = TitrationUtils.chooseTitratables(options.rl, mola);
} else if (options.rw) {
    titrating = TitrationUtils.chooseTitratables(pH, windows, mola);
} else if (options.rn) {
    titrating = TitrationUtils.chooseTitratables(AminoAcid3.valueOf(options.rn), mola);
}

MolecularDynamics molDyn = new MolecularDynamics(mola, mola.getPotentialEnergy(), mola.getProperties(), sh, thermostat, integrator);
molDyn.setFileType(fileType);
molDyn.setRestartFrequency(restartFrequency);

/* Pass MolecularDynamics to the PhMD constructor, which will attach itself as a MonteCarloListener if necessary. */
PhMD phmd = new PhMD(distribution, mola, molDyn, titrating, pH, mcStepFrequency);

if (!dynamics) {    // For testing.
    for (int i = 0; i < nSteps; i++) {
        phmd.mcUpdate();
    }
    return;
}

// Ready!
molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
