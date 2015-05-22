/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */

// ORTHOGONAL SPACE RANDOM WALK

// Apache Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Paralle Java Imports
import edu.rit.pj.Comm;

// Force Field X Imports
import ffx.algorithms.Barostat;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// First atom of the ligand.
int ligandStart = 1;

// Last atom of the ligand.
int ligandStop = -1;

// First atom of the ligand for the 2nd topology.
int ligandStart2 = 1;

// Last atom of the ligand for the 2nd topology.
int ligandStop2 = -1;

// Initial lambda value (0 is ligand in vacuum; 1 is ligand in PBC).
double initiaLambda = 0.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testOSRWBias [options] <filename> [filename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName:'1', 'Starting ligand atom.');
cli.s2(longOpt:'start2', args:1, argName:'1', 'Starting ligand atom for the 2nd topology.');
cli.f(longOpt:'final', args:1, argName:'-1', 'Final ligand atom.');
cli.f2(longOpt:'final2', args:1, argName:'-1', 'Final ligand atom for the 2nd topology.');
cli.l(longOpt:'lambda', args:1, argName:'0.0', 'Initial lambda value (> 1.0 distributes lambda across walkers)');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() > 2) {
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

// No electrostatics on ligand 1.
if (options.e) {
    noElec = true;
}

// Starting ligand atom for the 2nd topology.
if (options.s2) {
    ligandStart2 = Integer.parseInt(options.s2);
}

// Final ligand atom for the 2nd topology.
if (options.f2) {
    ligandStop2 = Integer.parseInt(options.f2);
}

// No electrostatics on ligand 2.
if (options.e2) {
    noElec2 = true;
}

// Starting lambda value.
if (options.l) {
    initialLambda = Double.parseDouble(options.l);
}

println("\n Testing Orthogonal Space Random Walk on " + filename);

File structureFile = new File(FilenameUtils.normalize(filename));
structureFile = new File(structureFile.getAbsolutePath());
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File histogramRestart = new File(baseFilename + ".his");
File lambdaRestart = null;
File dyn = null;

if (!histogramRestart.exists()) {
    logger.severe("\n Histogram restart file does not exist.");
} else if (!histogramRestart.canRead()) {
    logger.severe("\n Histogram restart file can not be read.");
}

// For a single process job, try to get the restart files from the current directory.
lambdaRestart = new File(baseFilename + ".lam");
dyn = new File(baseFilename + ".dyn");

if (!dyn.exists()) {
    dyn = null;
}

// Turn on computation of lambda derivatives.
System.setProperty("lambdaterm", "true");

// Relative free energies via the DualTopologyEnergy class require different
// default OSRW parameters than absolute free energies.
if (arguments.size() == 2) {
    // Condensed phase polarization is evaluated over the entire range.
    System.setProperty("polarization-lambda-start","0.0");
    // Polarization energy is not scaled individually by lambda, but
    // along with the overall potential energy of a topology.
    System.setProperty("polarization-lambda-exponent","0.0");
    // Ligand vapor electrostatics are not calculated. This cancels when the
    // difference between protein and water environments is considered.
    System.setProperty("ligand-vapor-elec","false");
    // Condensed phase polarization, without the ligand present, is unecessary.
    System.setProperty("no-ligand-condensed-scf","false");
}

// Open the first system
open(filename);

// Get a reference to the first system's ForceFieldEnergy and atom array.
ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
    ai.print();
}

// Turn off checks for overlapping atoms, which is expected for lambda=0.
energy.getCrystal().setSpecialPositionCutoff(0.0);
// OSRW will be configured for either single or dual topology.
OSRW osrw = null;
// Save a reference to the first topology.
topology1 = active;

if (arguments.size() == 1) {
    // Wrap the single topology ForceFieldEnergy inside an OSRW instance.
    osrw = new OSRW(energy, energy, lambdaRestart, histogramRestart, active.getProperties(),
        298.15, 1.0, 1.0, 1.0, false, sh);
} else {
    // Open the 2nd topology.
    filename = arguments.get(1);
    open(filename);
    energy = active.getPotentialEnergy();
    atoms = active.getAtomArray();
    // Apply the ligand atom selection for the 2nd topology.
    for (int i = ligandStart2; i <= ligandStop2; i++) {
        Atom ai = atoms[i - 1];
        ai.setApplyLambda(true);
        ai.print();
    }
    // Save a reference to the second topology.
    topology2 = active;
    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff(0.0);
    // Create the DualTopology potential energy.
    DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topology1, active);
    // Wrap the DualTopology potential energy inside an OSRW instance.
    osrw = new OSRW(dualTopologyEnergy, dualTopologyEnergy, lambdaRestart, histogramRestart, active.getProperties(),
        298.15, 1.0, 1.0, 1.0, false, sh);
}

/**
 * Stop propagating lambda to prevent adding new Gaussians
 * to the biasing potential, which would introduce artifacts into the
 * finite-difference derivatives.
 */
osrw.setPropagateLambda(false);
osrw.setLambda(initialLambda);
n = osrw.getNumberOfVariables();

assert(n%3 == 0);
n = n / 3;

// Finite-difference step size.
double step = 1.0e-5;

double[] x = new double[3*n];
double[] analytic = new double[3*n];
double[] g = new double[3*n];
double[] numeric = new double[3];
osrw.getCoordinates(x);

// Test Lambda gradients.
for (int j=0; j<3; j++) {
    double lambda = initialLambda - 0.001 + 0.001 * j;

    if (lambda - step < 0.0) {
        continue;
    } else if (lambda + step > 1.0) {
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

    logger.info(String.format(" Analytic dE/dL:   %15.8f", dEdLambda));
    logger.info(String.format(" Numeric  dE/dL:   %15.8f\n", dedl));

    // Calculate analytic dE/dX/dL
    osrw.setLambda(lambda);
    double e = 0.0;
    double orig = 0.0;
    double gradientTolerance = 1.0e-3;

    // Calculate finite-difference coordinate gradient
    for (int i=0; i<n; i++) {
        //Atom a0 = atoms[i];
        int i3 = i*3;
        int i0 = i3 + 0;
        int i1 = i3 + 1;
        int i2 = i3 + 2;

        // Calculate the analytic dE/dX
        osrw.energyAndGradient(x, analytic);

        // Find numeric dX
        orig = x[i0];
        x[i0] = orig + step;
        e = osrw.energyAndGradient(x,g);
        x[i0] = orig - step;
        e = e - osrw.energyAndGradient(x,g);
        x[i0] = orig;
        numeric[0] = e / (2.0 * step);

        // Find numeric dY
        orig = x[i1];
        x[i1] = orig + step;
        e = osrw.energyAndGradient(x,g);
        x[i1] = orig - step;
        e = e - osrw.energyAndGradient(x,g);
        x[i1] = orig;
        numeric[1] = e / (2.0 * step);

        // Find numeric dZ
        orig = x[i2];
        x[i2] = orig + step;
        e = osrw.energyAndGradient(x,g);
        x[i2] = orig - step;
        e = e - osrw.energyAndGradient(x,g);
        x[i2] = orig;
        numeric[2] = e / (2.0 * step);

        double dx = analytic[i0] - numeric[0];
        double dy = analytic[i1] - numeric[1];
        double dz = analytic[i2] - numeric[2];
        double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (len > gradientTolerance) {
            logger.info(" Atom " + i + String.format(" failed: %10.6f.", len)
                + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2])
                + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
            return;
        } else {
            logger.info(" Atom " + i + String.format(" passed: %10.6f.", len)
                + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2])
                + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
        }
    }
}
