// SIMULATED ANNEALING

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.SimulatedAnnealing;

// High temperature starting point.
double high = 1000.0;

// Low temperature end point.
double low = 10.0;

// Number of annealing steps.
int windows = 10;

// Number of molecular dynamics steps at each temperature.
int steps = 1000;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc anneal [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.n(longOpt:'steps', args:1, argName:'1000', 'Number of molecular dynamics steps per annealing window.');
cli.w(longOpt:'windows', args:1, argName:'10', 'Number of annealing windows.');
cli.l(longOpt:'low', args:1, argName:'10.0', 'Low temperature limit in degrees Kelvin.');
cli.t(longOpt:'high', args:1, argName:'1000.0', 'High temperature limit in degrees Kelvin.');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

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

logger.info("\n Running simulated annealing on " + filename);
open(filename);
SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, active.getPotentialEnergy(),
    active.getProperties(), null);
simulatedAnnealing.anneal(high, low, windows, steps);

