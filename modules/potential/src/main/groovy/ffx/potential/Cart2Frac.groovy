package ffx.potential

// Groovy Imports
import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.util.CliBuilder

// Apache Imports
import org.apache.commons.io.FilenameUtils

// FFX Imports
import ffx.crystal.Crystal
import ffx.potential.bonded.Atom
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The Cart2Frac script converts Cartesian coordinates to Fractional.
 * <br>
 * Usage:
 * <br>
 * ffxc Cart2Frac &lt;filename&gt;
 */
class Cart2Frac extends Script {

    /**
     * Options for the Cart2Frac Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Cart2Frac &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', description='Print this help message.') boolean help
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc Cart2Frac <filename>', header:' Options:');
        def options = new Options()
        cli.parseFromInstance(options, args)
        if (options.help) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Converting from Cartesian to fractional coordinates for " + modelFilename + "\n");

        // This is an interface specifying the closure-like methods.
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsFunctions()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }
        // Use PotentialsFunctions methods instead of Groovy method closures to do work.
        MolecularAssembly[] assemblies = functions.open(modelFilename)

        // Loop over each system.
        for (int i=0; i<assemblies.length; i++) {
            def system = assemblies[i];
            Crystal crystal = system.getCrystal().getUnitCell();

            List<Atom> atoms = system.getAtomList();
            int nAtoms = atoms.size();
            double[] frac = new double[3];
            double[] cart = new double[3];

            for (Atom atom in atoms) {
                atom.getXYZ(cart);
                crystal.toFractionalCoordinates(cart, frac);
                atom.moveTo(frac);
            }
        }

        String ext = FilenameUtils.getExtension(filename);
        filename = FilenameUtils.removeExtension(filename);
        if (ext.toUpperCase().contains("XYZ")) {
            functions.saveAsXYZ(assemblies[0], new File(filename + ".xyz"))
        } else {
            functions.saveAsPDB(assemblies, new File(filename + ".pdb"))
        }

    }
}

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