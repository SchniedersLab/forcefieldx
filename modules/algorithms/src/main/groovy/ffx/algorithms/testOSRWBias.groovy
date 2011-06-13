// TEST LAMBDA GRADIENTS

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
    println("\n Usage: ffxc testLambdaBias filename ligandStart ligandStop lambda fiction mass ");
    return; 
}

int ligandStart = 1;
if (args.size() > 1) {
    ligandStart = Integer.parseInt(args[1]);
}

int ligandStop = ligandStart;
if (args.size() > 2) {
    ligandStop = Integer.parseInt(args[2]);
}

double initialLambda = 0.5;
if (args.size() > 3) {
    initialLambda = Double.parseDouble(args[3]);
}

double friction = 1.0e-16;
if (args.size() > 4) {
    friction = Double.parseDouble(args[4]);
}

double mass = 1.0e-18;
if (args.size() > 5) {
    mass = Double.parseDouble(args[5]);
}

// Restart File
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn"); 
if (!dyn.exists()) {
    dyn = null;
} 

// Run 1000 molecular dynamics steps to build up a biasing potential.
int nSteps = 1000;
// Time step in femtoseconds.
double timeStep = 1.0;
// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;
// Choose to save coordinates only beyond the simulation length (avoid file cleanup).
double saveInterval = 10.0;
// Temperature in degrees Kelvin.
double temperature = 300.0;
// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;
// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Testing OSRW biasing potential for " + filename);

System.setProperty("lambdaterm", "true");

open(filename);

ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
int n = atoms.length;

// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
}

// Wrap the potential energy inside an OSRW instance.
OSRW osrw = new OSRW(energy, energy, active.getProperties(), atoms, temperature, timeStep);
osrw.setLambda(initialLambda);
osrw.setThetaFrication(friction);
osrw.setThetaMass(mass);

MolecularDynamics molDyn = new MolecularDynamics(active, osrw, active.getProperties(), null, thermostat);

molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);

/**
 * Stop propagating lambda to prevent adding new Gaussians 
 * to the biasing potential, which would introduce artifacts into the 
 * finite-difference derivatives.
 */
osrw.setPropagateLambda(false);

// Finite-difference step size.
double step = 1.0e-5;

double[] x = new double[n*3];
double[] analytic = new double[3*n];
double[] g = new double[3*n];
double[] numeric = new double[3];
energy.getCoordinates(x);

// Test Lambda gradients.
for (int j=0; j<3; j++) {
    double lambda = initialLambda - 0.001 + 0.001 * j;
   
    if (lambda - step < 0.0) {
        continue;
    }
    else if (lambda + step > 1.0) {
        continue;
    } else {
        osrw.setLambda(lambda);
    }

    // Calculate the energy and analytic dE/dX
    double eL = osrw.energyAndGradient(x, g);

    // Analytic dEdL
    double dEdLambda = osrw.getTotaldEdLambda();

    // Calculate the finite-difference dEdL
    osrw.setLambda(lambda + step);
    double lp = osrw.energyAndGradient(x, g);

    osrw.setLambda(lambda - step);
    double lm = osrw.energyAndGradient(x, g);

    double dedl = (lp - lm) / (2.0 * step);

    println(String.format(" Analytic dE/dL:   %15.8f", dEdLambda));
    println(String.format(" Numeric  dE/dL:   %15.8f\n", dedl));
    
    // Calculate analytic dE/dX/dL
    osrw.setLambda(lambda);
    double e = 0.0;
    double orig = 0.0;
    double gradientTolerance = 1.0e-3;
   
    // Calculate finite-difference coordinate gradient
    for (int i=ligandStart-1; i<ligandStop; i++) {
        Atom a0 = atoms[i];
        int i3 = i*3;
        int i0 = i3 + 0;
        int i1 = i3 + 1;
        int i2 = i3 + 2;

        // Calculate the analytic dE/dX
        osrw.energyAndGradient(x, analytic);
        
        // Find numeric dX
        orig = x[i0];
        x[i0] += step;
        e = osrw.energyAndGradient(x,g);
        x[i0] -= 2.0 * step;
        e -= osrw.energyAndGradient(x,g);
        x[i0] = orig;
        numeric[0] = e / (2.0 * step);

        // Find numeric dY
        orig = x[i1];
        x[i1] += step;
        e = osrw.energyAndGradient(x,g);
        x[i1] -= 2.0 * step;
        e -= osrw.energyAndGradient(x,g);
        x[i1] = orig;
        numeric[1] = e / (2.0 * step);

        // Find numeric dZ
        orig = x[i2];
        x[i2] += step;
        e = osrw.energyAndGradient(x,g);
        x[i2] -= 2.0 * step;
        e -= osrw.energyAndGradient(x,g);
        x[i2] = orig;
        numeric[2] = e / (2.0 * step);

        double dx = analytic[i0] - numeric[0];
        double dy = analytic[i1] - numeric[1];
        double dz = analytic[i2] - numeric[2];
        double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (len > gradientTolerance) {
            println(" " + a0.toShortString() + String.format(" failed: %10.6f.", len) + String.format(
                 "\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2]) 
                + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
            return;
        } else {
            println(" " + a0.toShortString() + String.format(" passed: %10.6f.", len) + String.format(
                 "\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2]) 
                + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
        }
    }
}
