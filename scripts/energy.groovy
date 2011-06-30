// ENERGY

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("\n Usage: ffxc energy filename");
   return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Running energy on " + filename);

open(filename); 
energy();
