import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.SimulatedAnnealing;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
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
String coordepsString = args[2];

// same, but for B factors
String bepsString = args[3];

// set the maximum number of refinement cycles
int maxiter = 50000;

// simulated annealing settings:
// type of simulated annealing refinement
RefinementMode refinementmode = RefinementMode.COORDINATES_AND_BFACTORS;
// starting temp
double highTemperature = 1000.0;
// ending temp
double lowTemperature = 100.0;
// number of steps to take between high and low temps
int annealingSteps = 10;
// number of MD steps at each annealing step
int mdSteps = 200;

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

System.setProperty("polarization","direct");
System.setProperty("tau-temperature","0.001");
open(modelfilename);
energy();

DiffractionFile diffractionfile = null;
if (datafilename != null) {
  diffractionfile = new DiffractionFile(datafilename, 1.0, false);
} else {
  diffractionfile = new DiffractionFile(active);
}

DiffractionData diffractiondata = new DiffractionData(active, active.getProperties(), SolventModel.POLYNOMIAL, diffractionfile);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();

// Do an initial loose optimization without an SCF.
RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES);
refinementMinimize.minimize(5.0, maxiter);

// Create a new ForceFieldEnergy instance with an SCF.
ForceField forceField = active.getForceField();
forceField.addForceFieldString(ForceFieldString.POLARIZATION, "mutual");
forceField.addForceFieldDouble(ForceFieldDouble.POLAR_EPS, 0.01);
ForceFieldEnergy forceFieldEnergy = new ForceFieldEnergy(active);

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

RefinementEnergy refinementEnergy = new RefinementEnergy(diffractiondata, refinementmode);
SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, refinementEnergy, active.getProperties(), refinementEnergy);
simulatedAnnealing.anneal(highTemperature, lowTemperature, annealingSteps, mdSteps);
diffractiondata.scaleBulkFit();
diffractiondata.printStats();

refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES);
if (coordeps < 0.0) {
    coordeps = refinementMinimize.getEps();
}
logger.info("\n RMS gradient convergence criteria: " + coordeps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(coordeps, maxiter);
diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.BFACTORS);
if (beps < 0.0) {
    beps = refinementMinimize.getEps();
}
logger.info("\n RMS gradient convergence criteria: " + beps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(beps, maxiter);
diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + "_refine.pdb");
diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_refine.mtz");
