/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

package ffx.potentials;

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.GenerateRotamers
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc saveRotamers [options] <PDB|XYZ>');
cli.h(longOpt:'help', 'Print this help message.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name (default is \' \').');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.t(longOpt:'torsions', args:1, 'Torsions to apply as "chi1,chi2,chi3;chi1,chi2,chi3;...');
cli.n(longOpt:'nChi', args:1, 'Number of torsions (will fill rotamers with 0).');
cli.vf(longOpt:'videoFile', args:1, argName:'file_rots.pdb', 'File to print torsion snapshots to.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line and open the file.
String filename = arguments.get(0);
open(filename);

int resID = 1;
String chain = " ";
String[] torStrings;
int nChi = 0;
File outFile;

if (options.t) {
    torStrings = options.t.split(";");
} else {
    logger.severe(" No torsion list provided.");
}

// Residue number.
if (options.r) {
    resID = Integer.parseInt(options.r);
} else {
    logger.severe(" No residue provided.");
}

// Chain Name.
if (options.c) {
    chain = options.c;
} else {
    logger.severe(" No chain provided.");
}

// Chain Name.
if (options.n) {
    nChi = Integer.parseInt(options.n);
} else {
    logger.severe(" No nChi provided.");
}

String newName = FilenameUtils.getBaseName(filename);
if (options.vw) {
    videoFile = options.vw;
} else {
    videoFile = newName + "_rots.pdb";
}

outFile = new File(newName + ".rotout.tmp");
outFile.deleteOnExit();

logger.info("\n Saving torsions for residue number " + resID + " of chain " + chain + ".");

MolecularAssembly molecularAssembly = (MolecularAssembly) active;
RotamerLibrary.initializeDefaultAtomicCoordinates(molecularAssembly.getChains());
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

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

GenerateRotamers genr = new GenerateRotamers(active, active.getPotentialEnergy(), residue, outFile, nChi, sh);
genr.setVideo(videoFile);
genr.applyAndSaveTorsions(torStrings);

