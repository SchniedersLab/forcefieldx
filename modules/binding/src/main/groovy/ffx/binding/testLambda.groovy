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

println("\n Testing the λ gradients of " + filename);

// For now - only use the van der Waals potential.

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

// Compute the λ = 1.0 energy.
double lambda = 1.0;
energy.setSoftCoreLambda(lambda);
double e1 = energy.energy(false, false);
println(String.format(" E(λ=1.0) : %20.8f.", e1));

// Compute the λ = 0.0 energy.
lambda = 0.0;
energy.setSoftCoreLambda(lambda);
double e0 = energy.energy(false, false);
println(String.format(" E(λ=0.0) : %20.8f.", e0));
println(String.format(" E(1)-E(0): %20.8f.", e1-e0));

// Set λ to be 0.5.
lambda = 0.5;
energy.setSoftCoreLambda(lambda);
// Turn on calculation of dE/dλ and dE/dλ/dX.
energy.lambdaGradients(true);

// Calculate the energy, dE/dX, dE/dλ and dE/dλ/dX
energy.energy(true, false);

// Analytic dEdλ
double dEdLambda = energy.getdEdLambda();

// Analytic dEdλdX
double[] dEdLdX = new double[n * 3];
energy.getdEdLambdadX(dEdLdX);

// Turn off calculation of dE/dλ and dE/dλ/dX.
energy.lambdaGradients(false);

// Calculate the finite-difference dEdλ and dEdλdX
double step = 1.0e-5;
double dedl = 0.0;
double[][] dedx = new double[2][n * 3];

energy.setSoftCoreLambda(lambda + step);
dedl = energy.energy(true, false);
energy.getGradients(dedx[0]);

energy.setSoftCoreLambda(lambda - step);
dedl -= energy.energy(true, false);
energy.getGradients(dedx[1]);

dedl /= (2.0 * step);
println(String.format(" Analytic vs. Numeric dE/dλ: %15.8f %15.8f", dEdLambda, dedl));

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
