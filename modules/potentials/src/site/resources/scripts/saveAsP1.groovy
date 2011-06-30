// SAVE AS P1

import org.apache.commons.io.FilenameUtils;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("Usage: ffxc saveAsP1 filename");
   return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Expanding to P1 for " + filename);

open(filename);

saveAsP1(new File(filename));
