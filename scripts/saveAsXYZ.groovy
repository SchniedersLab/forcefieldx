// SAVE AS XYZ

import org.apache.commons.io.FilenameUtils;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("Usage: ffxc saveAsXYZ filename");
   return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Writing out XYZ for " + filename);

open(filename);
filename = FilenameUtils.removeExtension(filename) + ".xyz";
save(new File(filename));
