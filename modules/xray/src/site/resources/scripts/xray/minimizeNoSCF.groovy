import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.SimulatedAnnealing;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.RefinementEnergy;
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
//default if epsString is not given on the command line
double eps = 5.0;

// set the maximum number of refinement cycles
int maxiter = 50000;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;


// Things below this line normally do not need to be changed.
// ===============================================================================================

if (epsString != null) {
   eps = Double.parseDouble(epsString);
}

logger.info("\n Running x-ray minimize without an SCF on " + modelfilename);
System.setProperty("polarization","direct");
System.setProperty("tau-temperature","0.001");
systems = open(modelfilename);
energy();

DiffractionFile diffractionfile = null;
if (datafilename != null) {
  diffractionfile = new DiffractionFile(datafilename, wA, neutron);
} else {
  diffractionfile = new DiffractionFile(systems, wA, neutron);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), diffractionfile);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

// Do an initial loose optimization without an SCF.
RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode);

logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(eps, maxiter);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + "_refine_noscf.pdb");
diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_refine_noscf.mtz");
