// MOLECULAR DYNAMICS

// Apache Imports
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Thermostat.Thermostats;

// Name of the file (PDB or XYZ).
String filename = args[0];

// If the filename is null or "-help" is specified print instructions.
if (filename == null || filename.equalsIgnoreCase("-help")) {
    println("\n Usage: ffxc md filename"
    + "\n\n Argument        Default"
    + "\n steps           1000000"
    + "\n time-step          1.00 (femtoseconds)" 
    + "\n print-interval     0.01 (picoseconds)"
    + "\n save-interval      0.10 (picoseconds)"
    + "\n temperature      298.14 (Kelvin)"
    + "\n thermostat    Berendsen [Adiabatic/Berendsen/Bussi]");
    return;
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
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Running molecular dynmaics on " + filename);

open(filename);

CompositeConfiguration properties = active.getProperties();

// Overwrite default script arguments
nSteps = properties.getInt("steps", nSteps);
timeStep = properties.getInt("time-step", timeSteps);
printInterval = properties.getInt("print-interval", printInterval);
saveInterval = properties.getInt("save-interval", saveInterval);
temperature = properties.getDouble("temperature", temperature);
String newThermostat = properties.getString("thermostat", thermostat.toString());
thermostat = Thermostats.valueOf(newThermostat.toUpperCase());

// Restart File
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn"); 
if (!dyn.exists()) {
    dyn = null;
}

MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), properties, null, thermostat);

molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
