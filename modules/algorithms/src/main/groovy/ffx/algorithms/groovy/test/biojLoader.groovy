/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

package ffx.algorithms.groovy.test;

import org.apache.commons.io.FilenameUtils
import org.biojava.bio.structure.Structure
import org.biojava.bio.structure.io.PDBFileReader

import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.potential.MolecularAssembly
import ffx.potential.parsers.PDBFileFilter

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc test.biojLoader [options] <PDB> <Biojava PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.p(longOpt:'print', 'Test by printing to file.');
cli.m(longOpt:'minimize', 'Not implemented: tests minimization of structure.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criterion');
cli.l(longOpt:'local', 'Use the local Potentials package loader instead of the User Interfaces loader');
cli.d(longOpt:'biojDefault', 'Forces the use of default mechanisms to find the Biojava structure file');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 2) {
    return cli.usage();
}

String filename = arguments.get(0);
File inFile = new File(filename);
if (inFile == null || !inFile.exists() || inFile.isDirectory()) {
    logger.severe(String.format(" Invalid file name %s", filename));
}
if (!new PDBFileFilter().acceptDeep(inFile)) {
    logger.severe(String.format(" File %s is not a recognized PDB file", filename));
}

String bjFilename = arguments.get(1);
File biojFile = new File(bjFilename);
if (biojFile == null || !biojFile.exists() || biojFile.isDirectory()) {
    logger.severe(String.format(" Invalid file name %s", filename));
}
if (!new PDBFileFilter().acceptDeep(biojFile)) {
    logger.severe(String.format(" File %s is not a recognized PDB file", bjFilename));
}

boolean print = true;
boolean minimize = false;
double eps = 1.0;
boolean rmsd = false;

PDBFileReader reader = new PDBFileReader();
Structure struct = reader.getStructure(biojFile);
/*if (options.d) {
biojFile = null;
}*/

AlgorithmFunctions functions;
if (options.l) {
    functions = new AlgorithmUtils();
} else {
    try {
        functions = getAlgorithmUtils();
    } catch (MissingMethodException ex) {
        logger.severe(String.format("Could not get any upper-level AlgorithmFunctions implementation: %s", ex.toString()));
    }
}

MolecularAssembly[] assemblies = functions.open(filename);
MolecularAssembly assembly = assemblies[0];

double e1 = 0;
if (options.m) {
    if (options.e) {
        eps = Double.parseDouble(options.e);
    }
    e1 = functions.minimize(assembly, eps).getTotalEnergy();
} else {
    e1 = functions.returnEnergy(assembly);
}

if (options.p) {
    String outFileName = FilenameUtils.removeExtension(inFile.getName());
    File outFile = new File(outFileName + ".pdb");
    int count = 1;
    while (outFile.exists() && count < 1000) {
        outFile = new File(String.format("%s%s_%d", outFileName, ".pdb", count++));
    }
    if (outFile.exists()) {
        logger.severe(String.format(" Could not version file: %s and all files %s through %s already exist.", outFileName, outFileName + "_1", outFileName + "_999"));
    }
    logger.info(String.format(" Printing FFX-loaded structure to file %s", outFile.getName()));
    functions.saveAsPDB(assembly, outFile);
}

MolecularAssembly[] biojAssemblies;
if (options.d) {
    biojAssemblies = functions.convertDataStructure(struct);
    biojFile = biojAssemblies[0].getFile();
} else {
    biojAssemblies = functions.convertDataStructure(struct, biojFile);
}
MolecularAssembly bioj = biojAssemblies[0];

double e2 = 0;
if (options.m) {
    e2 = functions.minimize(bioj, eps).getTotalEnergy();
} else {
    e2 = functions.returnEnergy(bioj);
}

logger.info(String.format(" Energy from original: %f", e1));
logger.info(String.format(" Energy from Biojava: %f", e2));

if (options.p) {
    String outFileName = FilenameUtils.removeExtension(biojFile.getName());
    File outFile = new File(outFileName + ".pdb");
    int count = 1;
    while (outFile.exists() && count < 1000) {
        outFile = new File(String.format("%s%s_%d", outFileName, ".pdb", count++));
    }
    if (outFile.exists()) {
        logger.severe(String.format(" Could not version file: %s and all files %s through %s already exist.", outFileName, outFileName + "_1", outFileName + "_999"));
    }
    logger.info(String.format(" Printing Biojava-loaded structure to file %s", outFile.getName()));
    functions.saveAsPDB(bioj, outFile);
}
