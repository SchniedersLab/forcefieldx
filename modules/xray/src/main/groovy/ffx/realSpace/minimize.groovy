/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */

// REALSPACE MINIMIZE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// RMS gradient per atom convergence criteria
double eps = 1.0;

// maximum number of refinement cycles
int maxiter = 1000;

// suffix to append to output data
String suffix = "_rsrefine";


// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc realspace.minimize [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:2, valueSeparator:',', argName:'data.map,1.0', 'specify input data filename (or simply provide the datafilename argument after the PDB file) and weight to apply to the data (wA)');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.m(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed refinement iterations');
cli.s(longOpt:'suffix', args:1, argName:'_rsrefine', 'output suffix');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Name of the file (PDB or XYZ).
String modelFilename = arguments.get(0);

// set up real space map data (can be multiple files)
List mapFiles = new ArrayList();
if (arguments.size() > 1) {
    RealSpaceFile realSpaceFile = new RealSpaceFile(arguments.get(1), 1.0);
    mapFiles.add(realSpaceFile);
}
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=2) {
	double wA = Double.parseDouble(options.ds[i+1]);
	RealSpaceFile realSpaceFile = new RealSpaceFile(options.ds[i], wA);
	mapFiles.add(realSpaceFile);
    }
}

if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.m) {
    maxiter = Integer.parseInt(options.m);
}

if (options.s) {
    suffix = options.s;
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

logger.info("\n Running real-space minimize on " + modelFilename);
systems = open(modelFilename);
if (mapFiles.size() == 0) {
    RealSpaceFile realspacefile = new RealSpaceFile(systems);
    mapFiles.add(realspacefile);
}
RealSpaceData realspacedata = new RealSpaceData(systems, systems[0].getProperties(), mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
energy();

RefinementMinimize refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES);

if (eps < 0.0) {
    eps = 1.0;
}
logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter);
refinementMinimize.minimize(eps, maxiter);

energy();
saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelFilename) + suffix + ".pdb"));