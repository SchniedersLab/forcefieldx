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

// TIMER

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData
import ffx.xray.DiffractionFile
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode

// The number of iterations.
int nEvals = 5;

// Compute the atomic coordinate gradient.
boolean gradient = true;

// Print the energy for each iteraction.
boolean print = true;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.timer [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.n(longOpt:'iterations', args:1, argName:'5', 'Number of iterations');
cli.c(longOpt:'threads', args:1, argName:'all', 'Number of SMP threads (ie. default uses all CPU cores)');
cli.g(longOpt:'gradient', args:1, argName:'true', 'Compute the atomic coordinats gradeint');
cli.v(longOpt:'verbose', args:1, argName:'true', 'Print out the energy for each step');
cli.d(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.r(longOpt:'mode', args:1, argName:'coordinates', 'type of refinement: [coordinates / bfactors / coordinates_and_bfactors / occupancies / bfactors_and_occupancies / coordinates_and_occupancies / coordinates_and_bfactors_and_occupancies]');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Load the number iterations.
if (options.n) {
    nEvals = Integer.parseInt(options.n);
}

// Load the number of threads.
if (options.c) {
    System.setProperty("pj.nt", options.c);
}

// Compute the gradient for each step.
if (options.g) {
    gradient = Boolean.parseBoolean(options.g);
}

// Print the energy for each step.
if (options.v) {
    print = Boolean.parseBoolean(options.v);
}

// Name of the file (PDB or XYZ).
String filename = arguments.get(0);

// set up diffraction data (can be multiple files)
List diffractionfiles = new ArrayList();
if (arguments.size() > 1) {
    DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false);
    diffractionfiles.add(diffractionfile);
}
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=3) {
	double wA = Double.parseDouble(options.ds[i+1]);
	boolean neutron = Boolean.parseBoolean(options.ds[i+2]);
	DiffractionFile diffractionfile = new DiffractionFile(options.ds[i], wA, neutron);
	diffractionfiles.add(diffractionfile);
    }
}

if (options.r) {
    try {
	refinementmode = RefinementMode.valueOf(options.r.toUpperCase());
    } catch (Exception e) {
	refinementmode = RefinementMode.COORDINATES;
    }
}

if (options.p) {
    System.setProperty("polarization", options.p);
}


logger.info("\n Timing energy and gradient for " + filename);

open(filename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(active, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(active, active.getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));
diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

RefinementEnergy refinementEnergy = RefinementEnergy.refinementEnergyFactory(diffractiondata, refinementmode);
int n = refinementEnergy.getNumberOfVariables();
double[] x = new double[n];
double[] g = new double[n];
refinementEnergy.getCoordinates(x);

Potential energy = refinementEnergy.getDataEnergy();

for (int i=0; i<nEvals; i++) {
    long time = -System.nanoTime();
    energy.energyAndGradient(x,g);
    time += System.nanoTime();
    logger.info(String.format(" Time %12.8f", time * 1.0e-9));
}