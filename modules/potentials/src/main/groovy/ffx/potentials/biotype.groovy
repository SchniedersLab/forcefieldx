// BIOTYPE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc biotype [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

systems = open(xyzname);
energy();

String mol = FilenameUtils.getBaseName(xyzname);
List atoms = active.getAtomList();
int index = 1;
for (Atom atom : atoms) {
    logger.info(String.format(" biotype %3d %4s \"%s\" %3d", index++, atom.getName(), mol, atom.getAtomType().type));
    List bonds = atom.getBonds();
    if (bonds != null) {
        for (Bond bond : bonds) {
            print String.format(" %4s", bond.get1_2(atom).getName());
        }
    }
    logger.info("");
}

return;

