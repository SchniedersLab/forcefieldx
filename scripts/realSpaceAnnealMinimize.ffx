import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.SimulatedAnnealing;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementEnergy;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;


// Name of the file (PDB or XYZ).
String modelfilename = args[0];

// input CCP4 map data (optional - if not given, data must be present as pdbfilename.[map]
String datafilename = args[1];

// data weight
double wA = 1.0;

// Set the RMS gradient per atom convergence criteria (optional)
String coordepsString = args[2];
// default if epsString not given on the command line
double coordeps = 1.0;

// set the maximum number of refinement cycles
int maxiter = 50000;

// simulated annealing settings:
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

if (coordepsString != null) {
   coordeps = Double.parseDouble(coordepsString);
}

System.setProperty("polarization","direct");
System.setProperty("tau-temperature","0.001");
open(modelfilename);
energy();

RealSpaceFile mapfile = null;
if (datafilename != null) {
  mapfile = new RealSpaceFile(datafilename, wA);
} else {
  mapfile = new RealSpaceFile(active, wA);
}

RealSpaceData realspacedata = new RealSpaceData(active, active.getProperties(), mapfile);

// Do an initial loose optimization without an SCF.
RefinementMinimize refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES);
refinementMinimize.minimize(5.0, maxiter);

// Create a new ForceFieldEnergy instance with an SCF.
ForceField forceField = active.getForceField();
forceField.addForceFieldString(ForceFieldString.POLARIZATION, "mutual");
forceField.addForceFieldDouble(ForceFieldDouble.POLAR_EPS, 0.01);
ForceFieldEnergy forceFieldEnergy = new ForceFieldEnergy(active);

energy();

RefinementEnergy refinementEnergy = new RefinementEnergy(realspacedata, RefinementMode.COORDINATES);
SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, refinementEnergy, active.getProperties(), refinementEnergy);
simulatedAnnealing.anneal(highTemperature, lowTemperature, annealingSteps, mdSteps);

refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES);
refinementMinimize.minimize(coordeps, maxiter);

energy();

saveAsPDB(new File(FilenameUtils.removeExtension(modelfilename) + "_rsrefine.pdb"));
