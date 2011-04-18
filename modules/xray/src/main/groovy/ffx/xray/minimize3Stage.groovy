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
String coordepsString = args[2];

// same, but for B factors
String bepsString = args[3];

// and occupancies (if necessary)
String occepsString = args[4];

// set the maximum number of refinement cycles
int maxiter = 50000;


// Things below this line normally do not need to be changed.
// ===============================================================================================

if (modelfilename == null) {
    modelfilename = "examples/1n7s.P212121.xyz";
}

double coordeps = -1.0;
if (coordepsString != null) {
    coordeps = Double.parseDouble(coordepsString);
}

double beps = -1.0;
if (bepsString != null) {
    beps = Double.parseDouble(bepsString);
}

double occeps = -1.0;
if (occepsString != null) {
    occeps = Double.parseDouble(occepsString);
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

diffractiondata.scalebulkfit();
diffractiondata.printstats();
energy();

RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES);
if (coordeps < 0.0) {
    coordeps = refinementMinimize.getEps();
}
println("\n RMS gradient convergence criteria: " + coordeps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(coordeps, maxiter);
diffractiondata.scalebulkfit();
diffractiondata.printstats();
energy();

refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.BFACTORS);
if (beps < 0.0) {
    beps = refinementMinimize.getEps();
}
println("\n RMS gradient convergence criteria: " + beps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(beps, maxiter);
diffractiondata.scalebulkfit();
diffractiondata.printstats();

if (diffractiondata.getAltResidues().size() > 0
    || diffractiondata.getAltMolecules().size() > 0){
    refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.OCCUPANCIES);
    if (occeps < 0.0){
        occeps = refinementMinimize.getEps();
    }
    println("\n RMS gradient convergence criteria: " + occeps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(occeps, maxiter);
    diffractiondata.scalebulkfit();
    diffractiondata.printstats();
} else {
    println("Occupancy refinement not necessary, skipping");
}

energy();

diffractiondata.writedata(FilenameUtils.removeExtension(modelfilename) + "_refine.mtz");
saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + "_refine.pdb"));
