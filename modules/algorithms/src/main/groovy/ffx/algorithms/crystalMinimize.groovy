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

// MINIMIZE

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.CrystalMinimize;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;
import ffx.algorithms.CrystalMinimize
import ffx.potential.XtalEnergy

// Default convergence criteria.
double eps = 1.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc crystalMinimize [options] <filename1> [filename2]');
cli.h(longOpt:'help', 'Print this help message.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.p(longOpt:'polarization', args:1, 'polarization model: [none / direct / mutual]');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Read in the first command line file.
String filename = arguments.get(0);

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Open the first topology.
systems = open(filename);

logger.info("\n Running minimize on " + filename);
logger.info(" RMS gradient convergence criteria: " + eps);

ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
XtalEnergy xtalEnergy = new XtalEnergy(forceFieldEnergy, active);

// Do the minimization
CrystalMinimize crystalMinimize = new CrystalMinimize(active, xtalEnergy, sh);


e = crystalMinimize.minimize(eps);

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(systems, new File(filename + ".pdb"));
}