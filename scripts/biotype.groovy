// BIOTYPE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// FFX Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;

// Input PDB file
String xyzname = args[0];
if (xyzname == null) {
   println "Usage: ffxc biotype.ffx poltype.xyz"
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

systems = open(xyzname);
energy();

String mol = FilenameUtils.getBaseName(xyzname);
List atoms = active.getAtomList();
int index = 1;
for (Atom atom : atoms) {
   print String.format(" biotype %3d %4s \"%s\" %3d", index++, atom.getName(), mol, atom.getAtomType().type); 
   List bonds = atom.getBonds();
   if (bonds != null) {
      for (Bond bond : bonds) {
         print String.format(" %4s", bond.get1_2(atom).getName());
      }  
   }
   println "";
}

return;

