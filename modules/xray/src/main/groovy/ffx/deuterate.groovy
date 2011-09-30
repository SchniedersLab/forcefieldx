// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.MSNode;

// Things below this line normally do not need to be changed.
// ===============================================================================================

def today = new Date();
logger.info(" " + today);
logger.info(" command line variables:");
logger.info(" " + args + "\n");

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc deuterate [options] <pdbfilename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Name of the PDB with crystal header information
String modelfilename = arguments.get(0);

systems = open(modelfilename);

for (int i=0; i<systems.length; i++) {
    Atom[] atoms = systems[i].getAtomArray();
    for (Atom a : atoms) {
	if (a.getAtomicNumber() == 1) {
	    Atom b = a.getBonds().get(0).get1_2(a);

	    // criteria for converting H to D
	    if (b.getAtomicNumber() == 7
		|| b.getAtomicNumber() == 8) {
		String name = a.getName().replaceFirst("H", "D");
		a.setName(name);
	    }
	}
    }

    ArrayList<MSNode> waters = systems[i].getWaters();
    for (MSNode node : waters) {
    	Molecule water = (Molecule) node;
	water.setName("DOD");
    }
}

saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + "_deuterate.pdb"));
