// XRAY MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// Name of the file (PDB or XYZ).
String modelfilename = args[0];

// input MTZ/CIF/CNS data (optional - if not given, data must be present as pdbfilename.[mtz/cif/ent/cns]
String datafilename = args[1];

// data weight
double wA = 1.0;

// is this neutron data?
boolean neutron = false;

// Set the RMS gradient per atom convergence criteria (optional)
String epsString = args[2];

// set the maximum number of refinement cycles
int maxiter = 50000;

// type of refinement
RefinementMode refinementmode = RefinementMode.BFACTORS;


// Things below this line normally do not need to be changed.
// ===============================================================================================

double eps = -1.0;
if (epsString != null) {
   eps = Double.parseDouble(epsString);
}

println("\n Running x-ray minimize on " + modelfilename);
systems = open(modelfilename);

DiffractionFile diffractionfile = null;
if (datafilename != null) {
  diffractionfile = new DiffractionFile(datafilename, wA, neutron);
} else {
  diffractionfile = new DiffractionFile(systems, wA, neutron);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfile);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode);
if (eps < 0.0) {
    eps = refinementMinimize.getEps();
}
println("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(eps, maxiter);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + "_b.pdb");
diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_b.mtz");
