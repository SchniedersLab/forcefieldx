/**
 * Title: Force Field X Description: Force Field X - Software for Molecular
 * Biophysics Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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

// SIMULATED ANNEALING

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.SimulatedAnnealing;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;

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

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc anneal [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.p(longOpt:'polarization', args:1, 'polarization model: [none / direct / mutual]');
cli.n(longOpt:'steps', args:1, argName:'1000', 'Number of molecular dynamics steps per annealing window.');
cli.f(longOpt:'dt', args:1, argName:'1.0', 'Time step in femtoseconds.');
cli.w(longOpt:'windows', args:1, argName:'10', 'Number of annealing windows.');
cli.l(longOpt:'low', args:1, argName:'10.0', 'Low temperature limit in degrees Kelvin.');
cli.t(longOpt:'high', args:1, argName:'1000.0', 'High temperature limit in degrees Kelvin.');
cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic]');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic/Berendsen/Bussi]');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    steps = Integer.parseInt(options.n);
}

// Load the number of annealing windows.
if (options.w) {
    windows = Integer.parseInt(options.w);
}

// Low temperature in degrees Kelvin.
if (options.l) {
    low =  Double.parseDouble(options.l);
}

// High temperature in degrees Kelvin.
if (options.t) {
    high =  Double.parseDouble(options.t);
}

// Load the time steps in femtoseconds.
if (options.f) {
    timeStep = Double.parseDouble(options.f);
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

logger.info("\n Running simulated annealing on " + filename);
systems = open(filename);
SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, active.getPotentialEnergy(),
    active.getProperties(), null, thermostat, integrator);
simulatedAnnealing.anneal(high, low, windows, steps, timeStep);

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}
