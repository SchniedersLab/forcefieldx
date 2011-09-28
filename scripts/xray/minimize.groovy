// XRAY MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// RMS gradient per atom convergence criteria
double eps = -1.0;

// RMS gradient convergence criteria for three stage refinement
double coordeps = -1.0;
double beps = -1.0;
double occeps = -1.0;

// maximum number of refinement cycles
int maxiter = 1000;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;

// do 3 stage refinement (coordinates, then B, then occupancies)?
boolean threestage = false;

// suffix to append to output data
String suffix = "_refine";

// include SCF/polarization?
boolean noscf = false;


// Things below this line normally do not need to be changed.
// ===============================================================================================

def today = new Date();
logger.info(" " + today);
logger.info(" command line variables:");
logger.info(" " + args + "\n");

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.minimize [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.D(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.e(longOpt:'eps', args:1, argName:'-1.0', 'RMS gradient convergence criteria (negative: automatically determine based on refinement type)');
cli.f(longOpt:'threeeps', args:3, valueSeparator:',', argName:'-1.0,-1.0,-1.0', 'RMS gradient convergence criteria for three stage refinement (negative: automatically determine for each stage)');
cli.m(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed refinement iterations');
cli.r(longOpt:'mode', args:1, argName:'coordinates', 'type of refinement: [coordinates / bfactors / coordinates_and_bfactors / occupancies / bfactors_and_occupancies / coordinates_and_occupancies / coordinates_and_bfactors_and_occupancies]');
cli.t(longOpt:'threestage', 'set to perform refinement in 3 stages: coordinates, bfactors, then occupancies - overrides mode setting if true');
cli.s(longOpt:'suffix', args:1, argName:'_refine', 'output suffix');
cli.S(longOpt:'scf', 'set to turn off SCF/polarization');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Name of the file (PDB or XYZ).
String modelfilename = arguments.get(0);

// set up diffraction data (can be multiple files)
List diffractionfiles = new ArrayList();
if (arguments.size() > 1) {
    DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false);
    diffractionfiles.add(diffractionfile);
}
if (options.D) {
    for (int i=0; i<options.Ds.size(); i+=3) {
	double wA = Double.parseDouble(options.Ds[i+1]);
	boolean neutron = Boolean.parseBoolean(options.Ds[i+2]);
	DiffractionFile diffractionfile = new DiffractionFile(options.Ds[i], wA, neutron);
	diffractionfiles.add(diffractionfile);
    }
}

if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.f) {
    coordeps = Double.parseDouble(options.fs[0]);
    beps = Double.parseDouble(options.fs[1]);
    occeps = Double.parseDouble(options.fs[2]);
}

if (options.m) {
    maxiter = Integer.parseInt(options.m);
}

if (options.r) {
    try {
	refinementmode = RefinementMode.valueOf(options.r.toUpperCase());
    } catch (Exception e) {
	refinementmode = RefinementMode.COORDINATES;
    }
}

if (options.t) {
    threestage = true;
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

logger.info("\n Running x-ray minimize on " + modelfilename);
systems = open(modelfilename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

if (threestage) {
    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES);
    if (coordeps < 0.0) {
	coordeps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + coordeps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(coordeps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();
    energy();

    refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.BFACTORS);
    if (beps < 0.0) {
	beps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + beps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(beps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();

    if (diffractiondata.getAltResidues().size() > 0
	|| diffractiondata.getAltMolecules().size() > 0){
	refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.OCCUPANCIES);
	if (occeps < 0.0){
	    occeps = refinementMinimize.getEps();
	}
	logger.info("\n RMS gradient convergence criteria: " + occeps + " max number of iterations: " + maxiter);
	refinementMinimize.minimize(occeps, maxiter);
	diffractiondata.scaleBulkFit();
	diffractiondata.printStats();
    } else {
	logger.info("Occupancy refinement not necessary, skipping");
    }
} else {
    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode);
    if (eps < 0.0) {
	eps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(eps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();
}

energy();

diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb");
diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + suffix + ".mtz");
