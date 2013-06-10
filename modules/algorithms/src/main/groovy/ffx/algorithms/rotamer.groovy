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

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// ENERGY
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.RotamerLibrary;
import ffx.potential.Rotamer;
import ffx.potential.ResidueEnumerations;
import ffx.potential.ForceFieldEnergy;

// Groovy Imports
import groovy.util.CliBuilder;
import ffx.algorithms.RotamerOptimization

// Things below this line normally do not need to be changed.
// ===============================================================================================
def startResID = -1;
def finalResID = -1;
def algorithm = 1;
def min = false;
def eps = 0.01;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rotamer [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'algorithm', args:1, argName:'1', 'Choices are independent residues (1) or all permuations (2).');
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting residue to perform the rotamer search on (-1 exits).');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Final residue to perform the rotamer search on (-1 exits).');
cli.m(longOpt:'minimize', args:1, argName:'0.01', 'Minimize the final structure to the given RMS gradient (Kcal/mole/A).');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Algorithm.
if (options.a) {
    algorithm = Integer.parseInt(options.a);
}

// Starting residue.
if (options.s) {
    startResID = Integer.parseInt(options.s);
}
// Load the number iterations.
if (options.f) {
    finalResID = Integer.parseInt(options.f);
}

if (options.m) {
    min = true;
    eps = Double.parseDouble(options.m);
}

if (finalResID < startResID || startResID < 0 || finalResID < 0) {
    return;
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Evaluating rotamers for residues " + startResID + " to " + finalResID);

open(filename);

logger.info(" Beginning Energy\n");
energy();

RotamerOptimization rotamerOptimization = new RotamerOptimization(active, sh);

if (algorithm == 1) {
    rotamerOptimization.optimize(startResID, finalResID, RotamerOptimization.Algorithm.INDEPENDENT);
} else if (algorithm == 2) {
    rotamerOptimization.optimize(startResID, finalResID, RotamerOptimization.Algorithm.GLOBAL);
} else if (algorithm == 3) {
    rotamerOptimization.optimize(startResID, finalResID, RotamerOptimization.Algorithm.SLIDING_WINDOW);
}

if (min) {
    minimize(eps);
}

logger.info(" Final Energy\n");
energy();

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);
if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(new File(filename + ".pdb"));
}
