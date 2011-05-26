// TEST LAMBDA GRADIENTS

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

if (args.size() < 3) {
   println("Usage: ffx testLambda filename ligandStart ligandStop");
   return;
}

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("\n Usage: ffxc testLambda filename ligandStart ligandStop lambda fiction mass ");
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

double initialFriction = 1.0e-18;
if (args.size() > 4) {
    initialFriction = Double.parseDouble(args[4]);
}

double mass = 1.0e-20;
if (args.size() > 5) {
    mass = Double.parseDouble(args[5]);
}

// Restart File
File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn"); 
if (!dyn.exists()) {
   dyn = null;
} 

// Run 500 molecular dynamics steps to build up a biasing potential.
int nSteps = 10000;
// Time step in femtoseconds.
double timeStep = 1.0;
// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.001;
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

println("\n Running molecular dynmaics on " + filename);

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
molDyn.setThetaFrication(initialFriction);
molDyn.setThetaMass(mass);
molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);

// Stop counting energy evaluations to prevent adding new Gaussians to the biasing potential
energy.doLambdaCounting(false);

// Compute the λ = 1.0 energy.
// double lambda = 1.0;
// energy.setLambda(lambda);
// double e1 = energy.energy(false, true);
// println(String.format(" E(λ=1.0) : %20.8f.", e1));

// Compute the λ = 0.0 energy.
// lambda = 0.0;
// energy.setLambda(lambda);
// double e0 = energy.energy(false, true);
// println(String.format(" E(λ=0.0) : %20.8f.", e0));
// println(String.format(" E(1)-E(0): %20.8f.", e1-e0));

// Test λ gradients.
for (int j=0; j<3; j++) {
   lambda = initialLambda - 0.001 + 0.001 * j;
   energy.setLambda(lambda);

   // Calculate the energy, dE/dX, dE/dλ, d2E/dλ2 and dE/dλ/dX
   double eL = energy.energy(true, false);

   // Analytic dEdλ and d2E/dλ2
   double dEdLambda = energy.getdEdLambda();
   double dE2dLambda2 = energy.getd2EdLambda2();

   // Analytic dEdλdX
   double[] dEdLdX = new double[n * 3];
   energy.getdEdLambdaGradient(dEdLdX);

   // Calculate the finite-difference dEdλ, d2Edλ2 and dEdλdX
   double step = 1.0e-5;
   double[][] dedx = new double[2][n * 3];

   energy.setLambda(lambda + step);
   double lp = energy.energy(true, false);
   double dedlp = energy.getdEdLambda();
   energy.getGradients(dedx[0]);

   energy.setLambda(lambda - step);
   double lm = energy.energy(true, false);
   double dedlm = energy.getdEdLambda()
   energy.getGradients(dedx[1]);

   double dedl = (lp - lm) / (2.0 * step);
   double d2edl2 = (dedlp - dedlm) / (2.0 * step);

   println(String.format(" Analytic dE/dλ:   %15.8f", dEdLambda));
   println(String.format(" Numeric  dE/dλ:   %15.8f\n", dedl));

   //println(String.format(" Analytic d2E/dλ2: %15.8f", dE2dLambda2));
   //println(String.format(" Numeric  d2E/dλ2: %15.8f\n", d2edl2));

   /*
   for (int i = ligandStart - 1; i < ligandStop; i++) {
       println(" dE/dX/dλ for Ligand Atom " + (i + 1));
       int index = i * 3;
       double dX = (dedx[0][index] - dedx[1][index]) / (2.0 * step);
       index++;
       double dY = (dedx[0][index] - dedx[1][index]) / (2.0 * step);
       index++;
       double dZ = (dedx[0][index] - dedx[1][index]) / (2.0 * step);
       println(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", 
            dEdLdX[i*3],dEdLdX[i*3+1],dEdLdX[i*3+2]));
       println(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX,dY,dZ));
   }
   printl("\n");
   */

   double[] x = new double[n*3];
   double[] analytic = new double[3*n];
   energy.getCoordinates(x);
   energy.energyAndGradient(x,analytic);
   double[] g = new double[3*n];

   double e = 0.0;
   double orig = 0.0;
   double gradientTolerance = 1.0e-3;
   double[] numeric = new double[3];

   for (int i=ligandStart-1; i<ligandStop; i++) {
      Atom a0 = atoms[i];
      int i3 = i*3;
      int i0 = i3 + 0;
      int i1 = i3 + 1;
      int i2 = i3 + 2;

      // Find numeric dX
      orig = x[i0];
      x[i0] += step;
      e = energy.energyAndGradient(x,g);
      x[i0] -= 2.0 * step;
      e -= energy.energyAndGradient(x,g);
      x[i0] = orig;
      numeric[0] = e / (2.0 * step);

      // Find numeric dY
      orig = x[i1];
      x[i1] += step;
      e = energy.energyAndGradient(x,g);
      x[i1] -= 2.0 * step;
      e -= energy.energyAndGradient(x,g);
      x[i1] = orig;
      numeric[1] = e / (2.0 * step);

      // Find numeric dZ
      orig = x[i2];
      x[i2] += step;
      e = energy.energyAndGradient(x,g);
      x[i2] -= 2.0 * step;
      e -= energy.energyAndGradient(x,g);
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
