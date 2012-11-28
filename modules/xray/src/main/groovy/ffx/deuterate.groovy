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

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.MSNode;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc deuterate [options] <pdbfilename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Name of the PDB with crystal header information
String modelfilename = arguments.get(0);

systems = open(modelfilename);

for (int i=0; i<systems.length; i++) {
    Atom[] atoms = systems[i].getAtomArray();
    for (Atom a : atoms) {
	if (a.getAtomicNumber() == 1) {
	    Atom b = a.getBonds().get(0).get1_2(a);

	    // criteria for converting H to D
	    if (b.getAtomicNumber() == 7
		|| b.getAtomicNumber() == 8) {
		String name = a.getName().replaceFirst("H", "D");
		a.setName(name);
	    }
	}
    }

    ArrayList<MSNode> waters = systems[i].getWaters();
    for (MSNode node : waters) {
    	Molecule water = (Molecule) node;
	water.setName("DOD");
    }
}

saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + "_deuterate.pdb"));
