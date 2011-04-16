// Lambda DYNAMICS

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("\n Usage: ffxc lambda filename lambda fiction ligandStart ligandStop");
   return;
}

double initialLambda = 1.0;
if (args.size() > 1) {
    initialLambda = Double.parseDouble(args[1]);
}

double initialFriction = 60.0;
if (args.size() > 2) {
    initialFriction = Double.parseDouble(args[2]);
}

int ligandStart = 1;
if (args.size() > 3) {
    ligandStart = Integer.parseInt(args[3]);
}

int ligandStop = ligandStart;
if (args.size() > 4) {
    ligandStop = Integer.parseInt(args[4]);
}

// Restart File
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn"); 
if (!dyn.exists()) {
   dyn = null;
} 

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 300.0;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Running molecular dynmaics on " + filename);

System.setProperty("mpoleterm","false");
System.setProperty("polarizeterm","false");

open(filename);

ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
int n = atoms.length;

// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
}

MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), 
    active.getProperties(), null, thermostat);

molDyn.doLambdaDynamics(true);
molDyn.setLambda(initialLambda);
molDyn.setFriction(initialFriction);

molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
