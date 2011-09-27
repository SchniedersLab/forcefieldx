// CONVERT FROM CARTESIAN TO FRACTIONAL COORDINATES

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc cart2frac [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Converting from Cartesian to fractional coordinates for " + filename);

systems = open(filename);

// Loop over each system.
for (int i=0; i<systems.length; i++) {
    system = systems[i];
    Crystal crystal = system.getCrystal().getUnitCell();

    List<Atom> atoms = system.getAtomList();
    int nAtoms = atoms.size();
    double[] frac = new double[3];
    double[] cart = new double[3];

    for (Atom atom in atoms) {
        atom.getXYZ(cart);
        crystal.toFractionalCoordinates(cart, frac);
        atom.moveTo(frac);
    }
}

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}