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

package ffx.potential

import ffx.potential.parsers.PDBFileFilter
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils
import org.biojava.bio.structure.Structure
import org.biojava.bio.structure.io.PDBFileReader

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc biojLoader [options] <PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.p(longOpt:'print', args:1, argName:'true', 'Test by printing to file.');
cli.f(longOpt:'file', args:1, argName:'<output file>.pdb', 'Specifies file to print to.');
cli.m(longOpt:'minimize', args:1, argName:'false', 'Not implemented: tests minimization of structure.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criterion');
cli.r(longOpt:'rmsd', args:1, argName:'false', 'Not implemented: calculates RMSD between initial and final structure.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

String filename = arguments.get(0);
File inFile = new File(filename);
if (inFile == null || !inFile.exists() || inFile.isDirectory()) {
    logger.severe(" Invalid file name %s", filename);
}
if (!new PDBFileFilter().acceptDeep(inFile)) {
    logger.severe(" File %s is not a recognized PDB file", filename);
}

boolean print = true;
File outFile = inFile;
boolean minimize = false;
double eps = 1.0;
boolean rmsd = false;

if (options.f) {
    outFile = new File(options.f);
    if (outFile.exists()) {
        if (outFile.isDirectory()) {
            logger.warning(" Output file is a directory: using versioned input file for output.");
            outFile = inFile;
        } else {
            logger.warning(" Output file already exists");
        }
    }
}

PDBFileReader reader = new PDBFileReader();
Structure struct = reader.getStructure(inFile);
PotentialsFunctions functions;
try {
    functions = getPotentialsFunctions();
} catch (MissingMethodException ex) {
    functions = new PotentialsUtils();
}
// Compare to just opening locally.
MolecularAssembly[] assemblies = functions.open(filename);
MolecularAssembly assembly = assemblies[0];

double e1 = functions.returnEnergy(assembly);

MolecularAssembly[] biojAssemblies = functions.convertDataStructure(struct, outFile);
MolecularAssembly bioj = biojAssemblies[0];

double e2 = functions.returnEnergy(bioj);

logger.info(String.format(" Energy from original: %f", e1));
logger.info(String.format(" Energy from Biojava: %f", e2));

if (options.f && options.p) {
    String outFileName = options.f;
    outFile = new File(outFileName);
    int count = 1;
    while (outFile.exists() && count < 1000) {
        outFile = new File(outFileName + String.format("_%d", count++));
    }
    if (outFile.exists) {
        logger.severe(String.format(" Could not version file: %s and all files %s through %s already exist.", outFileName, outFileName + "_1", outFileName + "_999"));
    }
    logger.info(String.format(" Printing to file %s", outFileName + String.format("_%d", count)));
}

if (options.m) {
    // Minimize code
    // Also check eps.
}

if (options.r) {
    // RMSD code.
}