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

package ffx.potentials;

import org.apache.commons.io.FilenameUtils;

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
cli.l(longOpt:'library', args:1, argName:'1', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.i(longOpt:'independent', args:1, argName: 'false', 'Independent draws nucleic acid rotamers independently of chain context.');
cli.s(longOpt:'start', args:1, argName: '-1', 'First rotamer to draw. Indexed from rotamer 0.');
cli.f(longOpt:'finish', args:1, argName: '-1', 'Last rotamer to draw. Indexed from rotamer 0.');
cli.x(longOpt:'all', args:1, argName: '0', 'Draw all rotamers beginning from the passed rotamer number (overrides other options). Indexed from rotamer 0.');
cli.u(longOpt:'upstreamPucker', args:1, argName: 'true', 'Adjusts the pucker of the 5\' residue to match the rotamer.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

int resID = 1;
String chain = " ";
int library = 1;
boolean independent = false;
int start = -1;
int finish = -1;
int allStart = 0;
// Will be checked for validity, and set false if invalid.
boolean upstreamPucker = true;
RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();

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

// Default behavior: draw all rotamers starting from 0 (as though set -x 0).
boolean saveAllRotamers = true;

// Start point
if (options.s) {
    start = Integer.parseInt(options.s);
    saveAllRotamers = false;
}

// Finish point
if (options.f) {
    finish = Integer.parseInt(options.f);
}

// Start from to all
if (options.x) {
    allStart = Integer.parseInt(options.x);
    saveAllRotamers = true;
    // Over-rides options.s
}

// Upstream pucker adjustment
if (options.u) {
    upstreamPucker = Boolean.parseBoolean(options.u);
}

logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".");

// Read in command line.
String filename = arguments.get(0);

// Open the file.
systems = open(filename);

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

if (library == 1) {
    rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
} else {
    rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
}

Rotamer[] rotamers = residue.getRotamers(rLib);
if (rotamers == null) {
    logger.severe(" There are no rotamers for residue + " + residue.toString());
}
int nrotamers = rotamers.length;

boolean isDeoxy = false; // Applies to prevResidue.
Residue prevResidue;
// If upstreamPucker is specified true, ensure that it is valid to apply it
// (residue is nucleic acid with an upstream partner), and set isDeoxy.
// If it is invalid, set upstreamPucker false.
if (upstreamPucker) {
    // Exception gets thrown if it's an amino acid, since "NA" is undefined.
    try {
        if (residue.getResidueType() == NA) {
            prevResidue = (Residue) residue.getPreviousResidue();
            // If no previous residue, set upstream pucker false.
            // The method used will ensure prevResidue is a nucleic acid.
            if (prevResidue == null) {
                upstreamPucker = false;
            } else {
                Atom HOs = (Atom) prevResidue.getAtomNode("HO\'");
                if (HOs == null) {
                    isDeoxy = true;
                }
            }
        }  else {
            upstreamPucker = false;
        }
    } catch (Exception e) {
        upstreamPucker = false;
    }
}

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);

if (saveAllRotamers) {
    if (allStart >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.");
    } else {
        for (int i = allStart; i < nrotamers; i++) {
            RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
            if (upstreamPucker) {
                double prevDelta = rotamers[i].chi1;
                if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                    // North pucker
                    RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true);
                } else {
                    RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true);
                }
            }
            if (ext.toUpperCase().contains("XYZ")) {
                saveAsXYZ(new File(filename + ".xyz"));
            } else {
                saveAsPDB(systems, new File(filename + ".pdb"));
            }
        }
    }
} else {
    if (start >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.");
    } else {
        if (finish >= nrotamers) {
            finish = nrotamers - 1;
        } else if (finish < start) {
            logger.info(" Specified finish point is before the start point; drawing only rotamer " + start);
            finish = start;
        }
        for (int i = start; i <= finish; i++) {
            RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
            if (upstreamPucker) {
                double prevDelta = rotamers[i].chi1;
                if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                    // North pucker
                    RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true);
                } else {
                    RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true);
                }
            }
            if (ext.toUpperCase().contains("XYZ")) {
                saveAsXYZ(new File(filename + ".xyz"));
            } else {
                saveAsPDB(systems, new File(filename + ".pdb"));
            }
        }
    }
}

