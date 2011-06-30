// MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("\n Usage: ffxc minimize filename");
   return;
}

// Set the RMS gradient per atom convergence criteria
String epsString = args[1];
double eps = 1.0;
if (epsString != null) {
   eps = Double.parseDouble(epsString);
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Running minimize on " + filename);
open(filename);

println("\n RMS gradient convergence criteria: " + eps);
e = minimize(eps);

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}

