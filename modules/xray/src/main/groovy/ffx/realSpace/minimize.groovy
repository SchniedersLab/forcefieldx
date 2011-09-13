// Real Space MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// Name of the file (PDB or XYZ).
String modelfilename = args[0];

// input CCP4 map data (optional - if not given, data must be present as pdbfilename.[map]
String datafilename = args[1];

// data weight
double wA = 1.0;

// Set the RMS gradient per atom convergence criteria (optional)
String epsString = args[2];
//default if epsString is not given on the command line
double eps = 1.0;

// set the maximum number of refinement cycles
int maxiter = 50000;


// Things below this line normally do not need to be changed.
// ===============================================================================================

if (epsString != null) {
   eps = Double.parseDouble(epsString);
}

logger.info("\n Running real space minimization on " + modelfilename);
systems = open(modelfilename);

RealSpaceFile mapfile = null;
if (datafilename != null) {
  mapfile = new RealSpaceFile(datafilename, wA);
} else {
  mapfile = new RealSpaceFile(systems, wA);
}

RealSpaceData realspacedata = new RealSpaceData(systems, systems[0].getProperties(), mapfile);

energy();

RefinementMinimize refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES);

logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(eps, maxiter);

energy();

saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + "_rsrefine.pdb"));
