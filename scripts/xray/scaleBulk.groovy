// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

// include SCF/polarization?
boolean noscf = false;

boolean writemaps = false;

boolean writemtz = false;

boolean timings = false;


// Things below this line normally do not need to be changed.
// ===============================================================================================

def today = new Date();
logger.info(" " + today);
logger.info(" command line variables:");
logger.info(" " + args + "\n");

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.scaleBulk [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.S(longOpt:'scf', 'set to turn off SCF/polarization');
cli.m(longOpt:'maps', 'set to output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps');
cli.t(longOpt:'timings', 'set to perform FFT test timings');
cli.w(longOpt:'mtz', 'write out MTZ containing structure factor coefficients');
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
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=3) {
	double wA = Double.parseDouble(options.ds[i+1]);
	boolean neutron = Boolean.parseBoolean(options.ds[i+2]);
	DiffractionFile diffractionfile = new DiffractionFile(options.ds[i], wA, neutron);
	diffractionfiles.add(diffractionfile);
    }
}

if (options.S) {
    noscf = true;
    logger.info(" setting polarization to direct (turning off SCF)!");
    System.setProperty("polarization","direct");
    System.setProperty("tau-temperature","0.001");
}

if (options.m) {
    writemaps = true;
}

if (options.t) {
    timings = true;
}

if (options.w) {
    writemtz = true;
}

systems = open(modelfilename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

if (writemtz) {
    diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_ffx.mtz");
}

if (writemaps) {
    diffractiondata.writeMaps(FilenameUtils.removeExtension(modelfilename) + "_ffx");
}

if (timings) {
    diffractiondata.timings();
}
