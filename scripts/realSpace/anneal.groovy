// REALSPACE SIMULATED ANNEALING

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.algorithms.SimulatedAnnealing;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementEnergy;
import ffx.xray.RefinementMinimize.RefinementMode;

// suffix to append to output data
String suffix = "_anneal";

// include SCF/polarization?
boolean noscf = false;

// starting temp
double highTemperature = 1000.0;

// ending temp
double lowTemperature = 100.0;

// number of steps to take between high and low temps
int annealingSteps = 10;

// number of MD steps at each annealing step
int mdSteps = 200;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;


// Things below this line normally do not need to be changed.
// ===============================================================================================

def today = new Date();
logger.info(" " + today);
logger.info(" command line variables:");
logger.info(" " + args + "\n");

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc realspace.anneal [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:2, valueSeparator:',', argName:'data.map,1.0', 'specify input data filename (or simply provide the datafilename argument after the PDB file) and weight to apply to the data (wA)');
cli.s(longOpt:'suffix', args:1, argName:'_anneal', 'output suffix');
cli.S(longOpt:'scf', 'set to turn off SCF/polarization');
cli.H(longOpt:'hightemp', args:1, argName:'1000.0', 'starting temperature');
cli.L(longOpt:'lowtemp', args:1, argName:'100.0', 'ending temperature');
cli.N(longOpt:'annealsteps', args:1, argName:'10', 'Number of steps between high and low temperature');
cli.n(longOpt:'mdsteps', args:1, argName:'200', 'Number of molecular dynamics steps at each temperature.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Name of the file (PDB or XYZ).
String modelfilename = arguments.get(0);

// set up real space map data (can be multiple files)
List mapfiles = new ArrayList();
if (arguments.size() > 1) {
    RealSpaceFile realspacefile = new RealSpaceFile(arguments.get(1), 1.0);
    mapfiles.add(realspacefile);
}
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=2) {
	double wA = Double.parseDouble(options.ds[i+1]);
	RealSpaceFile realspacefile = new RealSpaceFile(options.ds[i], wA);
	mapfiles.add(realspacefile);
    }
}

if (options.s) {
    suffix = options.s;
}

if (options.S) {
    noscf = true;
    logger.info(" setting polarization to direct (turning off SCF)!");
    System.setProperty("polarization","direct");
    System.setProperty("tau-temperature","0.001");
}

if (options.H) {
    highTemperature = Double.parseDouble(options.H);
}

if (options.L) {
    lowTemperature = Double.parseDouble(options.L);
}

if (options.N) {
    annealingSteps = Integer.parseInteger(options.N);
}

if (options.n) {
    mdSteps = Integer.parseInteger(options.n);
}

logger.info("\n Running simulated annealing on " + modelfilename);
open(modelfilename);

if (mapfiles.size() == 0) {
    RealSpaceFile realspacefile = new RealSpaceFile(active, 1.0);
    mapfiles.add(realspacefile);
}

RealSpaceData realspacedata = new RealSpaceData(active, active.getProperties(), mapfiles.toArray(new RealSpaceFile[mapfiles.size()]));

energy();

RefinementEnergy refinementEnergy = new RefinementEnergy(realspacedata, RefinementMode.COORDINATES);
SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, refinementEnergy, active.getProperties(), refinementEnergy);
simulatedAnnealing.anneal(highTemperature, lowTemperature, annealingSteps, mdSteps);
energy();

saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb"));
