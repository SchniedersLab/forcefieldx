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

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// ENERGY
import ffx.potential.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3;
import ffx.potential.bonded.Residue.ResidueType;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.RotamerOptimization
import ffx.algorithms.RotamerOptimization.Direction;
import edu.rit.pj.Comm
import java.util.Scanner;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rotamer [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.l(longOpt:'library', args:1, argName:'2', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.t(longOpt:'threeBody', args:1, argName:'true', 'Include 3-Body interactions in elimination criteria.');
cli.lR(longOpt:'listResidues', args: 1, argName: '-1', 'Choose a list of individual residues to optimize (eg. A11,A24,B40).');
cli.x(longOpt:'all', args:1, argName:'1', 'Optimize all residues in the system beginning from the passed residue number (overrides other options); for box optimization, optimizes all boxes from the passed index.');

int algorithm = 2;
int library = 2;
boolean threeBodyTerm = false;

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

List<String> resList = new ArrayList<>();
if (options.lR) {
    def tok = (options.lR).tokenize('.');
    for (String t : tok) {
        logger.info(" Adding " + t);
        resList.add(t);
    }
}

// Rotamer Library.
if (options.l) {
    library = Integer.parseInt(options.l);
}

if (options.t) {
    threeBodyTerm = Boolean.parseBoolean(options.t);
}

if (options.x) {
    allStartResID = Integer.parseInt(options.x);
}

if (options.x && options.lR) {
    logger.info("\n Incompatible options: x, lR.")
}
if (!options.x && !options.lR) {
    logger.info("\n Missing a selection option: x, lR.")
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Rotamer optimizing... ");
open(filename);

RotamerOptimization rotamerOptimization = new RotamerOptimization(active, active.getPotentialEnergy(), sh);
rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
RotamerLibrary.setUseOrigCoordsRotamer(true);
if (library == 1) {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
} else {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
}

ArrayList<Residue> residueList = new ArrayList<>();
Polymer[] polymers = active.getChains();

if (options.x) {
    int counter = 1;
    int nPolymers = polymers.length;
    for (int p=0; p<nPolymers; p++) {
        Polymer polymer = polymers[p];
        ArrayList<Residue> residues = polymer.getResidues();
        int nResidues = residues.size();
        for (int i=0; i<nResidues; i++) {
            Residue residue = residues.get(i);
            Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
            if (rotamers != null) {
                int nrot = rotamers.length;
                if (nrot == 1) {
                    RotamerLibrary.applyRotamer(residue, rotamers[0]);
                } else if (nrot > 1) {
                    if (counter >= allStartResID) {
                        residueList.add(residue);
                    }
                }
            }
            counter++;
        }
    }
    rotamerOptimization.setResidues(residueList);
} else if (options.lR) {
    int n = 0;
    for (String s : resList) {
        Character chainID = s.charAt(0);
        int i = Integer.parseInt(s.substring(1));
        for (Polymer p : polymers) {
            if (p.getChainID() == chainID) {
                List<Residue> rs = p.getResidues();
                for (Residue r : rs) {
                    if (r.getResidueNumber() == i) {
                        residueList.add(r);
                        Rotamer[] rotamers = RotamerLibrary.getRotamers(r);
                        if (rotamers != null && rotamers.size() > 1) {
                            n++;
                        }
                    }
                }
            }
        }
    }
    rotamerOptimization.setResiduesIgnoreNull(residueList);
    if (n < 1) {
        return;
    }
}

residueList = rotamerOptimization.getResidues();
energy();
RotamerLibrary.measureRotamers(residueList, false);
rotamerOptimization.optimize(RotamerOptimization.Algorithm.GLOBAL_DEE);

logger.info(" Final Minimum Energy");
energy();
String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);
if (ext.toUpperCase().contains("XYZ")) {
    saveAsXYZ(new File(filename + ".xyz"));
} else {
    saveAsPDB(new File(filename + ".pdb"));
}
