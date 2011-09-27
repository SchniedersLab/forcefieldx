// TEST LAMBDA GRADIENT

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

// First ligand atom.
int ligandStart = 1;

// Last ligand atom.
int ligandStop = ligandStart;

// Initial lambda value.
double initialLambda = 0.5;

// Print out the energy for each step.
boolean print = false;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testLambdaGradient [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName:'1', 'Starting ligand atom.');
cli.f(longOpt:'final', args:1, argName:'start', 'Final ligand atom.');
cli.l(longOpt:'lambda', args:1, argName:'0.5', 'Lambda value.');
cli.v(longOpt:'verbose', 'Print out the energy for each step');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line file.
String filename = arguments.get(0);

// Starting ligand atom.
if (options.s) {
    ligandStart = Integer.parseInt(options.s);
}

// Final ligand atom.
if (options.f) {
    ligandStop = Integer.parseInt(options.f);
}

// Starting lambda value.
if (options.l) {
    initialLambda = Double.parseDouble(options.l);
}

// Print the energy for each step.
if (options.v) {
    print = Boolean.parseBoolean(options.v);
}

logger.info("\n Testing lambda gradient for " + filename);

System.setProperty("lambdaterm","true");
System.setProperty("ligand-start", Integer.toString(ligandStart));
System.setProperty("ligand-stop", Integer.toString(ligandStop));

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
double e1 = energy.energy(false, print);
logger.info(String.format(" E(1):      %20.8f.", e1));

// Compute the Lambda = 0.0 energy.
lambda = 0.0;
energy.setLambda(lambda);
double e0 = energy.energy(false, print);
logger.info(String.format(" E(0):      %20.8f.", e0));
logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1-e0));

// Finite-difference step size.
double step = 1.0e-5;

// Error tolerence
double errTol = 1.0e-3;

// Test Lambda gradient in the neighborhood of the lambda variable.
for (int j=0; j<3; j++) {
    lambda = initialLambda - 0.01 + 0.01 * j;

    if (lambda - step < 0.0) {
        continue;
    }
    if (lambda + step > 1.0) {
        continue;
    }

    logger.info(String.format(" Current lambda value %6.4f", lambda));

    energy.setLambda(lambda);

    // Calculate the energy, dE/dX, dE/dLambda, d2E/dLambda2 and dE/dLambda/dX
    double e = energy.energy(true, print);

    // Analytic dEdLambda and d2E/dLambda2
    double dEdLambda = energy.getdEdL();
    double dE2dLambda2 = energy.getd2EdL2();

    // Analytic dEdLambdadX
    double[] dEdLdX = new double[n * 3];
    energy.getdEdXdL(dEdLdX);

    // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
    double[][] dedx = new double[2][n * 3];

    energy.setLambda(lambda + step);
    double lp = energy.energy(true, print);
    double dedlp = energy.getdEdL();
    energy.getGradients(dedx[0]);

    energy.setLambda(lambda - step);
    double lm = energy.energy(true, print);
    double dedlm = energy.getdEdL();
    energy.getGradients(dedx[1]);

    double dedl = (lp - lm) / (2.0 * step);
    double d2edl2 = (dedlp - dedlm) / (2.0 * step);

    double err = Math.abs(dedl - dEdLambda);
    if (err < errTol) {
        logger.info(String.format(" dE/dL passed:   %10.6f", err));
    } else {
        logger.info(String.format(" dE/dL failed: %10.6f", err));
    }
    logger.info(String.format(" Numeric:   %15.8f", dedl));
    logger.info(String.format(" Analytic:  %15.8f", dEdLambda));

    err = Math.abs(d2edl2 - dE2dLambda2);
    if (err < errTol) {
        logger.info(String.format(" d2E/dL2 passed: %10.6f", err));
    } else {
        logger.info(String.format(" d2E/dL2 failed: %10.6f", err));
    }
    logger.info(String.format(" Numeric:   %15.8f", d2edl2));
    logger.info(String.format(" Analytic:  %15.8f", dE2dLambda2));

    for (int i = ligandStart - 1; i < ligandStop; i++) {
        int ii = i * 3;
        double dX = (dedx[0][ii] - dedx[1][ii]) / (2.0 * step);
        double dXa = dEdLdX[ii];
        double eX = dX - dXa;
        ii++;
        double dY = (dedx[0][ii] - dedx[1][ii]) / (2.0 * step);
        double dYa = dEdLdX[ii];
        double eY = dY - dYa;
        ii++;
        double dZ = (dedx[0][ii] - dedx[1][ii]) / (2.0 * step);
        double dZa = dEdLdX[ii];
        double eZ = dZ - dZa;

        double error = Math.sqrt(eX * eX + eY * eY + eZ * eZ);
        if (error < errTol) {
            logger.info(String.format(" dE/dX/dL for Ligand Atom %d passed: %10.6f", i + 1, error));
        } else {
            logger.info(String.format(" dE/dX/dL for Ligand Atom %d failed: %10.6f", i + 1, error));
        }
        logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa,dYa,dZa));
        logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ));
    }

    logger.info("");
}
