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

import org.apache.commons.configuration.CompositeConfiguration

import ffx.potential.bonded.MolecularAssembly
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.ForceFieldFilter
import ffx.potential.parsers.PDBFilter
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Atom
import ffx.utilities.Keyword;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc mutatePDB [options] <PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.n(longOpt:'resname', args:1, argName:'GLY', 'New residue name.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name (default is \' \').');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

int resID = 1;
String resName = "ALA";
Character chain = ' ';

// Residue number.
if (options.r) {
    resID = Integer.parseInt(options.r);
}

// Residue Name.
if (options.n) {
    resName = options.n;
}

// Chain Name.
if (options.c) {
    chain = options.c.toCharacter();
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Mutating residue number " + resID + " of chain " + chain + " to " + resName);

File structure = new File(filename);
int index = filename.lastIndexOf(".");
String name = filename.substring(0, index);
MolecularAssembly molecularAssembly = new MolecularAssembly(name);
molecularAssembly.setFile(structure);

CompositeConfiguration properties = Keyword.loadProperties(structure);
ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
ForceField forceField = forceFieldFilter.parse();
molecularAssembly.setForceField(forceField);

PDBFilter pdbFilter = new PDBFilter(structure, molecularAssembly, forceField, properties);
pdbFilter.mutate(chain,resID,resName);
pdbFilter.readFile();
molecularAssembly.finalize(true);
pdbFilter.writeFile(structure, false);

