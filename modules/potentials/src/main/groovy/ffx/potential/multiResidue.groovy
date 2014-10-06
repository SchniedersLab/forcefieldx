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

package ffx.potentials

// Groovy Imports
import groovy.util.CliBuilder;

import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.MolecularAssembly;

// Things below this line normally do not need to be changed.
// ===============================================================================================

int resID = 1;
Character chain = ' ';
MultiResidue multiResidue;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc multiResidue [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.r(longOpt:'resID', args:1, argName:'1', 'Residue number.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name (default is \' \').');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Residue number.
if (options.r) {
    resID = Integer.parseInt(options.r);
}

// Chain Name.
if (options.c) {
    chain = options.c.toCharacter();
}

// Read in command line.
String filename = arguments.get(0);

open(filename);

Residue residue;
Polymer[] polymers = active.getChains();
for (int i = 0; i < polymers.length; i++) {
    Polymer polymer = polymers[i];
    if (chain.equals(polymer.getChainID())) {
        residue = polymer.getResidue(resID);
        if (residue != null) {
            multiResidue = new MultiResidue(residue, active);
            polymer.addMultiResidue(multiResidue);
        }
    }
}

Residue newResidue = new Residue("ALA", residue.getResidueNumber(), residue.getResidueType());
multiResidue.addResidue(newResidue);
