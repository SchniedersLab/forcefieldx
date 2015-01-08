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

package ffx.potentials;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Atom;
import ffx.utilities.Keyword;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc mutatePDB [options] <PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.n(longOpt:'resname', args:1, argName:'ALA', 'New residue name.');
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
molecularAssembly.finalize(true, forceField);
pdbFilter.writeFile(structure, false);

