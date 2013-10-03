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

package ffx.potentials;

import org.apache.commons.io.FilenameUtils;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.Rotamer;
import ffx.potential.RotamerLibrary;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc saveRotamers [options] <PDB|XYZ>');
cli.h(longOpt:'help', 'Print this help message.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name (default is \' \').');
cli.l(longOpt:'library', args:1, argName:'1', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.i(longOpt:'independent', args:1, argName: 'false', 'Independent draws nucleic acid rotamers independently of chain context.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

int resID = 1;
String chain = " ";
int library = 1;
boolean independent = false;

// Residue number.
if (options.r) {
    resID = Integer.parseInt(options.r);
}

// Chain Name.
if (options.c) {
    chain = options.c;
}

// Rotamer Library.
if (options.l) {
    library = Integer.parseInt(options.l);
}

// Rotamer independence
if (options.i) {
    independent = Boolean.parseBoolean(options.i);
}

logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".");

// Read in command line.
String filename = arguments.get(0);

// Open the file.
systems = open(filename);

MolecularAssembly molecularAssembly = (MolecularAssembly) active;
RotamerLibrary.loadPriorAtomicCoordinates(molecularAssembly);
Polymer polymer = molecularAssembly.getChain(chain);
if (polymer == null) {
    logger.info(" Polymer + " + chain + " does not exist.");
    return;
}
Residue residue = polymer.getResidue(resID);
if (residue == null) {
    logger.info(" Residue + " + resID + " does not exist.");
    return;
}

if (library == 1) {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
} else {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
}

Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
if (rotamers == null) {
    logger.info(" There are no rotamers for residue + " + residue.toString());
}

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

for (int i = 0; i < rotamers.length; i++){
    RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
    if (ext.toUpperCase().contains("XYZ")) {
        saveAsXYZ(new File(filename + ".xyz"));
    } else {
        saveAsPDB(systems, new File(filename + ".pdb"));
    }
}

