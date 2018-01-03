
package ffx.potential

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The Analysis script provides analysis of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Analysis [options] &lt;filename&gt;
 */
class Analysis extends Script {

    /**
     * Options for the Analysis Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Energy [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {

        def cli = new CliBuilder(usage:' ffxc Analysis [options] <filename>', header:' Options:')

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
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

        logger.info("\n Analysis of " + modelFilename);

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
        MolecularAssembly activeAssembly = assemblies[0]
        functions.analysis(activeAssembly)
        functions.energy(activeAssembly)
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
