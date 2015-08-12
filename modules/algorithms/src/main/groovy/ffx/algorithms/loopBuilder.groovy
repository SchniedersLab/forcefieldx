
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

// MINIMIZE

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports

import ffx.algorithms.Minimize;
import ffx.algorithms.SimulatedAnnealing;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// Default convergence criteria.
double eps = 1.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================
for (String arg : args) {
    logger.info(arg);
}

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc loopBuilder [options] <filename1>');
cli.h(longOpt:'help', 'Print this help message.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);
}

System.setProperty("buildLoops", "true");
System.setProperty("polarization", "direct");
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
ForceFieldEnergy energy = active.getPotentialEnergy();
// Set built atoms active/use flags to true (false for other atoms).
Atom[] atoms = active.getAtomArray();
for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setActive(true);
        ai.setUse(true);
    } else {
        ai.setActive(false);
        ai.setUse(false);
    }
}

logger.info("\n Running minimize on built atoms of " + active.getName());
logger.info(" RMS gradient convergence criteria: " + eps);

// Do the minimization without vdW.
// energy.setVDW(false);

e = minimize(eps);
// Run SA with vdW.

logger.info("\n Running simulated annealing on " + active.getName());
// energy.setVDW(true);

// High temperature starting point.
double high = 1000.0;
// Low temperature end point.
double low = 10.0;
// Number of annealing steps.
int windows = 10;
// Number of molecular dynamics steps at each temperature.
int steps = 1000;
// Time step in femtoseconds.
double timeStep = 1.0;
// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;
// Integrators [ BEEMAN, RESPA, STOCHASTIC]
Integrators integrator = null;

SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, energy, active.getProperties(), null, thermostat, integrator);
simulatedAnnealing.anneal(high, low, windows, steps, timeStep);

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    ai.setUse(true);
}

// Do the minimization
e = minimize(0.1);

// MD, SA, OSRW and/or Side-Chain DEE

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);
saveAsPDB(systems, new File(filename + ".pdb"));


