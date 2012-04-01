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
simulatedAnnealing.anneal(high, low, windows, steps);

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}
