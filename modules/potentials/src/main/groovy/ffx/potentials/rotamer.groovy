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

// Things below this line normally do not need to be changed.
// ===============================================================================================
def startResID = -1;
def finalResID = -1;
def algorithm = 1;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rotamer [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'algorithm', args:1, argName:'1', 'Choices are independent residues (1) or all permuations (2).');
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting residue to perform the rotamer search on (-1 exits).');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Final residue to perform the rotamer search on (-1 exits).');

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

if (finalResID < startResID || startResID < 0 || finalResID < 0) {
    return;
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Evaluating rotamers for residues " + startResID + " to " + finalResID);

open(filename);

energy();

ForceFieldEnergy potential = (ForceFieldEnergy) active.getPotentialEnergy();
Polymer[] polymers = active.getChains();

if (algorithm == 1) {
    for (def i=startResID; i<=finalResID; i++) {
        Residue residue = polymers[0].getResidue(i);
        print residue.toString();
        def name = ResidueEnumerations.AminoAcid3.valueOf(residue.getName());
        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
        if (rotamers == null) {
            continue;
        }
        e = potential.energy(false, true);
        bestRotamer = -1;
        for (j=0; j<rotamers.length;j++) {
            Rotamer rotamer = rotamers[j];
            RotamerLibrary.applyRotamer(name, residue, rotamer);
            newE = potential.energy(false, true);
            if (newE < e) {
                bestRotamer = j;
            }
        }
        if (bestRotamer > -1) {
            Rotamer rotamer = rotamers[bestRotamer];
            RotamerLibrary.applyRotamer(name, residue, rotamer);
        }
    }
} else {
    ArrayList<Residue> residues = new ArrayList<Residue>();
    int permutations = 1;
    for (def i=startResID; i<=finalResID; i++) {
        Residue residue = polymers[0].getResidue(i);
        residues.add(residue);
        def name = ResidueEnumerations.AminoAcid3.valueOf(residue.getName());
        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
        if (rotamers != null) {
            permutations *= rotamers.length;
        }
    }
    logger.info(" Number of permutations: " + permutations);
    ArrayList<Integer> optimum = new ArrayList<Integer>();
    double minEnergy = RotamerLibrary.rotamerOptimization(active, residues, Double.MAX_VALUE, optimum);
    for (def i=startResID; i<=finalResID; i++) {
        Residue residue = polymers[0].getResidue(i);
        def name = ResidueEnumerations.AminoAcid3.valueOf(residue.getName());
        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
        int j = optimum.remove(0);
        if (rotamers != null) {
            Rotamer rotamer = rotamers[j];
            RotamerLibrary.applyRotamer(name, residue, rotamer);
        }
    }
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
