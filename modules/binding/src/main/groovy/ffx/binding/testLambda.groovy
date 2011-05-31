// TEST LAMBDA GRADIENTS

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
int ligandStart = Integer.parseInt(args[1]);
int ligandStop = Integer.parseInt(args[2]);

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Testing the Lambda gradients of " + filename);
open(filename);

ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
int n = atoms.length;

// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
}

// Compute the Lambda = 1.0 energy.
double lambda = 1.0;
energy.setLambda(lambda);
double e1 = energy.energy(false, false);
println(String.format(" E(Lambda=1.0) : %20.8f.", e1));

// Compute the Lambda = 0.0 energy.
lambda = 0.0;
energy.setLambda(lambda);
double e0 = energy.energy(false, false);
println(String.format(" E(Lambda=0.0) : %20.8f.", e0));
println(String.format(" E(1)-E(0): %20.8f.", e1-e0));

// Set Lambda to be 0.5.
lambda = 0.5;
energy.setLambda(lambda);
// Turn on calculation of dE/dLambda, d2E/dLambda2 and dE/dLambda/dX.
energy.lambdaGradient(true);

// Calculate the energy, dE/dX, dE/dLambda, d2E/dLambda2 and dE/dLambda/dX
double e = energy.energy(true, false);

// Analytic dEdLambda and d2E/dLambda2
double dEdLambda = energy.getdEdLambda();
double dE2dLambda2 = energy.getd2EdLambda2();

// Analytic dEdLambdadX
double[] dEdLdX = new double[n * 3];
energy.getdEdLambdaGradient(dEdLdX);

// Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
double step = 1.0e-5;
double[][] dedx = new double[2][n * 3];

energy.setLambda(lambda + step);
double dedlp = energy.energy(true, false);
double d2edl2p = energy.getdEdLambda();
energy.getGradients(dedx[0]);

energy.setLambda(lambda - step);
double dedlm = energy.energy(true, false);
double d2edl2m = energy.getdEdLambda()
energy.getGradients(dedx[1]);

double dedl = (dedlp - dedlm) / (2.0 * step);
double d2edl2 = (d2edl2p - d2edl2m) / (2.0 * step);

println(String.format(" Analytic vs. Numeric dE/dLambda:   %15.8f %15.8f", dEdLambda, dedl));
println(String.format(" Analytic vs. Numeric d2E/dLambda2: %15.8f %15.8f", dE2dLambda2, d2edl2));

for (int i = ligandStart - 1; i < ligandStop; i++) {
    println(" dE/dX/dLambda for Ligand Atom " + (i + 1));
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
