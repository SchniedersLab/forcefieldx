
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

// Groovy Imports
import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.algorithms.NQFlipper
import ffx.algorithms.RotamerOptimization
import ffx.potential.bonded.Residue
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import groovy.util.CliBuilder;

// FFX Imports

// Things below this line normally do not need to be changed.
// ===============================================================================================
int algorithm = 1;
int library = 2;
int start = -1;
int end = -1;
int allStart = 1;
double radius = 4.0;
String chain = "A";

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rotamer [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'algorithm', args:1, argName:'1', 'Only try a simple flip (1), allow all rotamers (2), or repack around each NQ (3)');
cli.l(longOpt:'library', args:1, argName:'2', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.x(longOpt:'all', args:1, argName:'1', 'Optimize all NQ residues in the system beginning from the passed residue number (overrides other options)');
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting residue to perform the NQ search on (-1 exits).');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Final residue to perform the NQ search on (-1 exits).');
cli.r(longOpt:'radius', args:1, argName:'4.0', 'Radius around NQ residues to repack.');
cli.c(longOpt:'chain', args:1, argName:'A', 'Chain the residue selection belongs to.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

if (options.a) {
    algorithm = Integer.parseInt(options.a);
}
if (options.l) {
    library = Integer.parseInt(options.l);
}
if (options.r) {
    radius = Double.parseDouble(options.r);
}
if (options.s) {
    start = Integer.parseInt(options.s);
}
if (options.f) {
    end = Integer.parseInt(options.f);
}
if (options.x) {
    allStart = Integer.parseInt(options.x);
}
if (options.c) {
    chain = options.c;
}

if (!options.x) {
    logger.info("\n Flipping NQ for residues " + start + " to " + end);
} else {
    logger.info("\n Flipping NQ for all residues beginning at " + allStart);
}

AlgorithmFunctions functions;
try {
    functions = getAlgorithmUtils();
} catch (MissingMethodException e) {
    functions = new AlgorithmUtils();
}
MolecularAssembly[] systems = functions.open(arguments.get(0));
MolecularAssembly active = systems[0];

int counter = 1;
List<Residue> residueList = new ArrayList<>();
if (options.x) {
    Polymer[] polymers = active.getChains();
    int nPolymers = polymers.length;
    for (int p = 0; p < nPolymers; p++) {
        Polymer polymer = polymers[p];
        List<Residue> residues = polymer.getResidues();
        int nResidues = residues.size();
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues.get(i);
            Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
            if (rotamers != null) {
                int nrot = rotamers.length;
                if (nrot == 1) {
                    RotamerLibrary.applyRotamer(residue, rotamers[0]);
                } else if (nrot > 1) {
                    if (counter >= allStart) {
                        residueList.add(residue);
                    }
                }
            }
            counter++;
        }
    }
}

if (algorithm == 1) {
    NQFlipper flipper = new NQFlipper(active, active.getPotentialEnergy());
    if (options.x) {
        flipper.setResidues(residueList);
    } else if (options.c) {
        flipper.setResidues(chain, start, end);
    } else {
        flipper.setResidues(start, end);
    }
    flipper.flipNQs();
    functions.saveAsPDB(active, active.getFile());
}