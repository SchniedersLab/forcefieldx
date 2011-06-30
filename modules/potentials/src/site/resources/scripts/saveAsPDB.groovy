// SAVE AS PDB

import org.apache.commons.io.FilenameUtils;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("Usage: ffxc saveAsPDB filename");
   return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Writing out PDB for " + filename);

systems = open(filename);
filename = FilenameUtils.removeExtension(filename) + ".pdb";
saveAsPDB(systems, new File(filename));

