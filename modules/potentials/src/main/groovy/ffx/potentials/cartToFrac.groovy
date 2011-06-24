// CONVERT FROM CARTESIAN TO FRACTIONAL COORDINATES

import org.apache.commons.io.FilenameUtils;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
    println("Usage: ffxc cartToFrac filname");
    return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Converting from Cartesian to fractional coordinates for " + filename);
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