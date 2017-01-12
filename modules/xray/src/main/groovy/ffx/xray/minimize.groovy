/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

// XRAY MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.parsers.DiffractionFile;

// RMS gradient per atom convergence criteria
double eps = -1.0;

// RMS gradient convergence criteria for three stage refinement
double coordeps = -1.0;
double beps = -1.0;
double occeps = -1.0;

// maximum number of refinement cycles
int maxiter = 1000;

// type of refinement
RefinementMode refinementmode = RefinementMode.COORDINATES;

// do 3 stage refinement (coordinates, then B, then occupancies)?
boolean threestage = false;

// suffix to append to output data
String suffix = "_refine";


// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.minimize [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.e(longOpt:'eps', args:1, argName:'-1.0', 'RMS gradient convergence criteria (negative: automatically determine based on refinement type)');
cli.f(longOpt:'threeeps', args:3, valueSeparator:',', argName:'-1.0,-1.0,-1.0', 'RMS gradient convergence criteria for three stage refinement (negative: automatically determine for each stage)');
cli.m(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed refinement iterations');
cli.r(longOpt:'mode', args:1, argName:'coordinates', 'type of refinement: [coordinates / bfactors / coordinates_and_bfactors / occupancies / bfactors_and_occupancies / coordinates_and_occupancies / coordinates_and_bfactors_and_occupancies]');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');
cli.t(longOpt:'threestage', 'set to perform refinement in 3 stages: coordinates, bfactors, then occupancies - overrides mode setting if true');
cli.s(longOpt:'suffix', args:1, argName:'_refine', 'output suffix');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

List<String> arguments = options.arguments();
String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    modelfilename = arguments.get(0);
} else {
    return cli.usage();
}

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

if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.f) {
    coordeps = Double.parseDouble(options.fs[0]);
    beps = Double.parseDouble(options.fs[1]);
    occeps = Double.parseDouble(options.fs[2]);
}

if (options.m) {
    maxiter = Integer.parseInt(options.m);
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

if (options.t) {
    threestage = true;
}

if (options.s) {
    suffix = options.s;
}

logger.info("\n Running x-ray minimize on " + modelfilename);
systems = open(modelfilename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(),
    SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

if (threestage) {
    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES);
    if (coordeps < 0.0) {
	coordeps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + coordeps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(coordeps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();
    energy();

    refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.BFACTORS);
    if (beps < 0.0) {
	beps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + beps + "\n Maximum number of iterations: " + maxiter + "\n");
    refinementMinimize.minimize(beps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();

    if (diffractiondata.getAltResidues().size() > 0
	|| diffractiondata.getAltMolecules().size() > 0){
	refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.OCCUPANCIES);
	if (occeps < 0.0){
	    occeps = refinementMinimize.getEps();
	}
	logger.info("\n RMS gradient convergence criteria: " + occeps + " max number of iterations: " + maxiter);
	refinementMinimize.minimize(occeps, maxiter);
	diffractiondata.scaleBulkFit();
	diffractiondata.printStats();
    } else {
	logger.info("Occupancy refinement not necessary, skipping");
    }
} else {
    RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode);
    if (eps < 0.0) {
	eps = refinementMinimize.getEps();
    }
    logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
    refinementMinimize.minimize(eps, maxiter);
    diffractiondata.scaleBulkFit();
    diffractiondata.printStats();
}

energy();

diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb");
diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + suffix + ".mtz");
