
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

import org.apache.commons.io.FilenameUtils;
import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmUtils;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.parsers.PDBFileFilter;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc biojLoader [options] <PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.p(longOpt:'print', 'Test by printing to file.');
cli.f(longOpt:'file', args:1, argName:'<biojava file>.pdb', 'Specifies a separate file for Biojava to load from');
cli.m(longOpt:'minimize', 'Not implemented: tests minimization of structure.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criterion');
cli.l(longOpt:'local', 'Use the local Potentials package loader instead of the User Interfaces loader');
cli.d(longOpt:'biojDefault', 'Forces the use of default mechanisms to find the Biojava structure file');

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
File biojFile = inFile;
boolean minimize = false;
double eps = 1.0;
boolean rmsd = false;

if (options.f) {
    biojFile = new File(options.f);
    if (biojFile.isDirectory()) {
        logger.warning(" Biojava file is a directory: using other input file.");
        biojFile = inFile;
    }
}

if (options.d) {
    biojFile = null;
}

PDBFileReader reader = new PDBFileReader();
Structure struct = reader.getStructure(biojFile);

AlgorithmFunctions functions;
if (!options.l) {
    try {
        functions = getAlgorithmUtils();
    } catch (MissingMethodException ex) {
        logger.severe(String.format("Could not get any upper-level AlgorithmFunctions implementation: %s", ex.toString()));
    }
} else {
    functions = new AlgorithmUtils();
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
    File outFile = new File(String.format("%s_%d", outFileName, count));
    while (outFile.exists() && count < 1000) {
        outFile = new File(String.format("%s_%d", outFileName, count++));
    }
    if (outFile.exists()) {
        logger.severe(String.format(" Could not version file: %s and all files %s through %s already exist.", outFileName, outFileName + "_1", outFileName + "_999"));
    }
    logger.info(String.format(" Printing FFX structure to file %s", outFile.getName()));
    functions.saveAsPDB(bioj, outFile);
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
    e1 = functions.minimize(assembly, eps).getTotalEnergy();
} else {
    e1 = functions.returnEnergy(assembly);
}

logger.info(String.format(" Energy from original: %f", e1));
logger.info(String.format(" Energy from Biojava: %f", e2));

if (options.p) {
    String outFileName = FilenameUtils.removeExtension(biojFile.getName());
    File outFile = new File(outFileName);
    int count = 1;
    while (outFile.exists() && count < 1000) {
        outFile = new File(String.format("%s_%d", outFileName, count++));
    }
    if (outFile.exists()) {
        logger.severe(String.format(" Could not version file: %s and all files %s through %s already exist.", outFileName, outFileName + "_1", outFileName + "_999"));
    }
    logger.info(String.format(" Printing Biojava structure to file %s", outFile.getName()));
    functions.saveAsPDB(bioj, outFile);
}