
package ffx.potential

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.potential.bonded.Loop
import ffx.potential.parsers.PDBFilter
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The LoopClosure script tests a loop closure algorithm for placing loops.
 * <br>
 * Usage:
 * <br>
 * ffxc LoopClosure [options] &lt;filename&gt;
 */
class LoopClosure extends Script {

    /**
     * Options for the LoopClosure Script.
     * <br>
     * Usage:
     * <br>
     * ffxc LoopClosure [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', description='Print this help message.') boolean help
        /**
         * -f or --first to set the first residue to draw (indexed from residue 0).
         */
        @Option(shortName='f', defaultValue='2', description='First residue to draw. Indexed from residue 0.') boolean first
        /**
         * -l or --last to set the last residue to draw (indexed from residue 0).
         */
        @Option(shortName='l', defaultValue='4', description='Last residue to draw. Indexed from residue 0.') boolean last
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
        def cli = new CliBuilder(usage:' ffxc FracToCart <filename>', header:' Options:');
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

        logger.info("\n Running LoopClosure on " + modelFilename + "\n");

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
        MolecularAssembly[] assemblies = functions.open(modelFilename)

        int firstResidue = options.first
        int lastResidue = options.last

        if (lastResidue - firstResidue != 3)
        {
            logger.info("\n\nInvalid residue range. Residue range must consist of 3 residues.\n\n");
        }

        boolean writeFile = true;
        ForceFieldEnergy forceFieldEnergy = assemblies[0].getPotentialEnergy();
        forcefieldEnergy.setPrintOnFailure(false, false);
        Loop loop = new Loop(assemblies[0]);
        List<double[]> loopSolutions = loop.generateLoops(firstResidue, lastResidue);

        for(int i = 0; i < loopSolutions.size(); i++){
            //Test for using alternative coordinates with generateLoops method.
            // loopSolutions = loop.generateLoops(stt_res+1, end_res+1,loopSolutions.get(i));
            forceFieldEnergy.setCoordinates(loopSolutions.get(i));
            File modifiedFile = new File("loop_" + modelFilename + "_" + i);
            PDBFilter modFilter = new PDBFilter(modifiedFile, assemblies[0], null, null);
            if (writeFile) {
                modFilter.writeFile(modifiedFile, true);
            }
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