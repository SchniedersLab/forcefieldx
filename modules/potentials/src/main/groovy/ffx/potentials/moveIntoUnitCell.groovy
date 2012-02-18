// MOVE MOLECULES INTO THE UNIT CELL

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Polymer;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc moveIntoUnitCell [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Moving molecular centers of mass into the unit cell for " + filename);

systems = open(filename);

// Loop over each system.
for (int i=0; i<systems.length; i++) {
    system = systems[i];
    Crystal crystal = system.getCrystal().getUnitCell();

    int nAtoms = 0;
    double[] com = new double[3];
    double[] translate = new double[3];

    // Move the polymers together
    List<Polymer> polymers = system.getChains();
    if (polymers != null && polymers.size() > 0) {

        // Find the center of mass
        for (polymer in polymers) {
            List<Atom> atoms = polymer.getAtomList();
            nAtoms = nAtoms + atoms.size();
            for (atom in atoms) {
                com[0] = com[0] + atom.getX();
                com[1] = com[1] + atom.getY();
                com[2] = com[2] + atom.getZ();
            }
        }
        com[0] = com[0] / nAtoms;
        com[1] = com[1] / nAtoms;
        com[2] = com[2] / nAtoms;

        // Calculate the translation vector for the center of mass
        crystal.toPrimaryCell(com, translate);
        translate[0] = translate[0] - com[0];
        translate[1] = translate[1] - com[1];
        translate[2] = translate[2] - com[2];

        // Move each atom
        for (polymer in polymers) {
            List<Atom> atoms = polymer.getAtomList();
            for (atom in atoms) {
                atom.move(translate);
            }
        }
    }

    // Loop over each molecule
    List<Molecule> molecules = system.getMolecules();
    for (molecule in molecules) {
        List<Atom> atoms = molecule.getAtomList();
        // Find the center of mass
        com[0] = 0.0;
        com[1] = 0.0;
        com[2] = 0.0;
        for (atom in atoms) {
            com[0] = com[0] + atom.getX();
            com[1] = com[1] + atom.getY();
            com[2] = com[2] + atom.getZ();
        }

        nAtoms = atoms.size();
        com[0] = com[0] / nAtoms;
        com[1] = com[1] / nAtoms;
        com[2] = com[2] / nAtoms;

        // Calculate the translation vector for the center of mass
        crystal.toPrimaryCell(com, translate);
        translate[0] = translate[0] - com[0];
        translate[1] = translate[1] - com[1];
        translate[2] = translate[2] - com[2];

        // Move each atom
        for (atom in atoms) {
            atom.move(translate);
        }
    }

    // Loop over each water
    List<MSNode> waters = system.getWaters();
    for (water in waters) {
        List<Atom> atoms = water.getAtomList();
        // Find the center of mass
        com[0] = 0.0;
        com[1] = 0.0;
        com[2] = 0.0;
        for (atom in atoms) {
            com[0] = com[0] + atom.getX();
            com[1] = com[1] + atom.getY();
            com[2] = com[2] + atom.getZ();
        }

        nAtoms = atoms.size();
        com[0] = com[0] / nAtoms;
        com[1] = com[1] / nAtoms;
        com[2] = com[2] / nAtoms;

        // Calculate the translation vector for the center of mass
        crystal.toPrimaryCell(com, translate);
        translate[0] = translate[0] - com[0];
        translate[1] = translate[1] - com[1];
        translate[2] = translate[2] - com[2];

        // Move each atom
        for (atom in atoms) {
            atom.move(translate);
        }
    }

    // Loop over each ion
    List<MSNode> ions = system.getIons();
    for (ion in ions) {
        List<Atom> atoms = ion.getAtomList();
        // Find the center of mass
        com[0] = 0.0;
        com[1] = 0.0;
        com[2] = 0.0;
        for (atom in atoms) {
            com[0] = com[0] + atom.getX();
            com[1] = com[1] + atom.getY();
            com[2] = com[2] + atom.getZ();
        }

        nAtoms = atoms.size();
        com[0] = com[0] / nAtoms;
        com[1] = com[1] / nAtoms;
        com[2] = com[2] / nAtoms;

        // Calculate the translation vector for the center of mass
        crystal.toPrimaryCell(com, translate);
        translate[0] = translate[0] - com[0];
        translate[1] = translate[1] - com[1];
        translate[2] = translate[2] - com[2];

        // Move each atom
        for (atom in atoms) {
            atom.move(translate);
        }
    }
}

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}