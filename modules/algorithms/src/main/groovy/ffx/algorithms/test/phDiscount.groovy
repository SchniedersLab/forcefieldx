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

// MOLECULAR & STOCHASTIC DYNAMICS

import org.apache.commons.io.FilenameUtils

import groovy.cli.picocli.CliBuilder

import ffx.algorithms.MolecularDynamics
import ffx.algorithms.PhDiscount
import ffx.algorithms.integrators.Integrator
import ffx.algorithms.integrators.IntegratorEnum
import ffx.algorithms.thermostats.Thermostat
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.potential.MolecularAssembly
import ffx.potential.extended.ExtendedSystem
import ffx.potential.extended.TitrationUtils

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double dt = 1.0;

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
int moveFrequency = 10;
int stepsPerMove = 100;
int rotamerMoveRatio = 0;

// Simulation pH
Double pH;

// Titrating residue list.
List<String> rlTokens = new ArrayList<>();
double cutoffs = 10.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc test.phDiscount [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
cli.rl(longOpt:'resList', args:1, 'Titrate a list of residues (eg A4.A8.B2.B34)');
cli.pH(longOpt:'pH', args:1, argName:'7.4', 'Constant simulation pH.');
cli.mc(longOpt:'titrationFrequency', args:1, argName:'100', 'Number of steps between Monte-Carlo proton attempts.')
cli.mcd(longOpt:'titrationDuration', args:1, argName:'100', 'Number of steps for which to run continuous proton dynamics during MC move.');
cli.mcr(longOpt:'rotamerMoveRatio', args:1, argName:'0', 'Number of steps between Monte-Carlo rotamer attempts.')
cli.cut(longOpt:'cutoff', args:1, argName:'1000', 'Value of vdw-cutoff and pme-cutoff.');

//cli.ra(longOpt:'resAll', 'Titrate all residues.');
//cli.rn(longOpt:'resName', args:1, 'Titrate a list of residue names (eg "LYS,TYR,HIS")');
//cli.rw(longOpt:'resWindow', args:1, 'Titrate all residues with intrinsic pKa within [arg] units of simulation pH.');
//cli.tt(longOpt:'titrateTermini', args:1, argName:'false', 'Titrate amino acid chain ends.');

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Required Options
if (!options.rl) {
    logger.warning("Must list titrating residues (-rl).");
    return;
} else {
    def tok = (options.rl).tokenize('.');
    for (String t : tok) {
        rlTokens.add(t);
    }
}

if (options.cut) {
    cutoffs = Double.parseDouble(options.cut);
}

// Suggested Options
if (!options.pH) {
    /* DISCOUNT interprets null as 7.4, then issues a warning (in a more readable location). */
    pH = null;
} else {
    pH = Double.parseDouble(options.pH);
}

if (options.rw) {
    window = Double.parseDouble(options.rw);
}

if (options.mc) {
    moveFrequency = Integer.parseInt(options.mc);
}
if (options.mcd) {
    stepsPerMove = Integer.parseInt(options.mcd);
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
if (options.d) {
    dt = Double.parseDouble(options.d);
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

if (options.mcmd) {
    mcRunTime = Integer.parseInt(options.mcmd);
}
if (options.i) {
    integrator = Integrator.parseIntegrator(options.i)
}

List<String> arguments = options.arguments();
if (arguments == null || arguments.size() != 1) {
    logger.warning("No input file supplied.");
    return;
}
String filename = arguments.get(0);
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

TitrationUtils.initDiscountPreloadProperties(cutoffs);
MolecularAssembly mola = TitrationUtils.openFullyProtonated(filename);

ExtendedSystem esvSystem = new ExtendedSystem(mola, pH);
esvSystem.populate(options.rl);

// Create the MolecularDynamics.
MolecularDynamics molDyn = new MolecularDynamics(mola, mola.getPotentialEnergy(), mola.getProperties(), sh, thermostat, integrator);
molDyn.setFileType(fileType);
molDyn.setRestartFrequency(restartFrequency);

// Create the DISCOuNT object, linking it to the MD.
PhDiscount phmd = new PhDiscount(mola, esvSystem, molDyn,
    dt, printInterval, saveInterval, initVelocities,
    fileType, restartFrequency, dyn);

// Launch dynamics through the DISCOuNT controller.
phmd.dynamic(nSteps, pH, temperature, moveFrequency, stepsPerMove);
