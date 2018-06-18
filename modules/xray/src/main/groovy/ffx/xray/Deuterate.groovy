package xray

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Molecule
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * Deuterate changes exchangeable hydrogen atoms to deuterium atoms for a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Deuterate &lt;pdbfile1&gt;
 */
class Deuterate extends Script {

	/**
	 * Options for the X-ray Deuterate Script.
	 * <br>
	 * Usage:
	 * <br>
	 * ffxc xray.Deuterate &lt;pdbfile1&gt;
	 */
	class Options {
		/**
		 * -h or --help to print a help message
		 */
		@Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
		boolean help
		/**
		 * The final argument should be a PDB filename.
		 */
		@Unparsed(description = 'A PDB filename.')
		List<String> filename
	}

	/**
	 * Execute the script.
	 */
	def run() {
		def cli = new CliBuilder()
		cli.name = "ffxc xray.Deuterate"

		def options = new Options()
		cli.parseFromInstance(options, args)

		if (options.help == true) {
			return cli.usage()
		}

		List<String> arguments = options.filename

		// Name of the PDB with crystal header information
		String modelfilename
		if (arguments != null && arguments.size() > 0) {
			// Read in command line.
			modelfilename = arguments.get(0)
		}  else if (active == null) {
			return cli.usage()
		} else {
			modelfilename = active.getFile()
		}

		logger.info("\n Running xray.Deuterate on " + modelfilename)

		// This is an interface specifying the closure-like methods.
		PotentialsFunctions functions
		try {
			// Use a method closure to try to get an instance of UIUtils (the User Interfaces
			// implementation, which interfaces with the GUI, etc.).
			functions = getPotentialsUtils()
		} catch (MissingMethodException ex) {
			// If Groovy can't find the appropriate closure, catch the exception and build
			// an instance of the local implementation.
			functions = new PotentialsUtils()
		}

		// Use PotentialsFunctions methods instead of Groovy method closures to do work.
		MolecularAssembly[] systems = functions.open(modelfilename)

		for (int i=0; i<systems.length; i++) {
			Atom[] atoms = systems[i].getAtomArray()
			for (Atom a : atoms) {
				if (a.getAtomicNumber() == 1) {
					Atom b = a.getBonds().get(0).get1_2(a)

					// Criteria for converting H to D
					if (b.getAtomicNumber() == 7
							|| b.getAtomicNumber() == 8) {
						String name = a.getName().replaceFirst("H", "D")
						a.setName(name)
					}
				}
			}

			ArrayList<MSNode> waters = systems[i].getWaters()
			for (MSNode node : waters) {
				Molecule water = (Molecule) node
				water.setName("DOD")
			}
		}

		functions.saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + "_deuterate.pdb"))
	}
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
