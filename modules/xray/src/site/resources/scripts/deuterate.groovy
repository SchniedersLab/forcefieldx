// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.MSNode;

// input PDB file
String filename = args[0];
if (filename == null){
  filename = "examples/1N7S.pdb";
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

systems = open(filename);

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

saveAsPDB(systems, new File(FilenameUtils.removeExtension(filename) + "_deuterate.pdb"));
