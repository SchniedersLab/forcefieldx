
package ffx.realspace

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.potential.MolecularAssembly
import ffx.xray.RefinementMinimize
import ffx.xray.RefinementMinimize.RefinementMode

/**
 * The Real Space Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Minimize [options] &lt;filename [file2...]&gt;
 */
class Minimize extends Script {

    /**
     * Options for the Real Space Minimize Script.
     * <br>
     * Usage:
     * <br>
     * ffxc realspace.Minimize [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -e or --eps RMS gradient convergence criteria (default of -1 automatically determines eps based on refinement type).
         */
        @Option(shortName = 'e', longName = 'eps', defaultValue = '-1.0', description = 'RMS gradient convergence criteria (default of -1 automatically determines eps based on refinement type)')
        double eps
        /**
         * -i or --iterations Maximum number of optimization iterations.
         */
        @Option(shortName = 'i', longName = 'iterations ', defaultValue = '1000', description = ' Maximum number of optimization iterations.')
        int iterations
        /**
         * -D or --data Specify input data filename and weight applied to the data (wA).
         */
        @Option(shortName = 'D', longName = 'data', defaultValue = '', numberOfArguments = 2, valueSeparator = ',',
                description = 'Specify input data filename and weight applied to the data (wA).')
        String[] data
        /**
         * The final arguments should be a PDB filename and data filename (CIF or MTZ).
         */
        @Unparsed(description = "PDB file and a Real Space Map file.")
        List<String> filenames
    }
    
    def run() {

        def cli = new CliBuilder()
        cli.name = "ffxc realspace.Minimize"

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        AlgorithmFunctions aFuncts
        try {
            // getAlgorithmUtils is a magic variable/closure passed in from ModelingShell
            aFuncts = getAlgorithmUtils()
        } catch (MissingMethodException ex) {
            // This is the fallback, which does everything necessary without magic names
            aFuncts = new AlgorithmUtils()
        }

        List<String> arguments = options.filenames

        // RMS gradient per atom convergence criteria.
        double eps = options.eps
        
        // Maximum number of refinement cycles.
        int maxiter = options.iterations
        
        // Suffix to append to output data
        String suffix = "_refine"

        String modelfilename
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelfilename = active.getFile()
        }

        logger.info("\n Running simulated annealing on " + modelfilename)

        MolecularAssembly[] systems = aFuncts.open(modelfilename)

        // set up real space map data (can be multiple files)
        List mapfiles = new ArrayList()
        if (arguments.size() > 1) {
            ffx.realspace.parsers.RealSpaceFile realspacefile = new ffx.realspace.parsers.RealSpaceFile(arguments.get(1), 1.0)
            mapfiles.add(realspacefile)
        }
        if (options.data) {
            for (int i=0; i<options.data.size(); i+=2) {
                double wA = Double.parseDouble(options.data[i+1])
                ffx.realspace.parsers.RealSpaceFile realspacefile = new ffx.realspace.parsers.RealSpaceFile(options.data[i], wA)
                mapfiles.add(realspacefile)
            }
        }

        if (mapfiles.size() == 0) {
            ffx.realspace.parsers.RealSpaceFile realspacefile = new ffx.realspace.parsers.RealSpaceFile(systems[0], 1.0)
            mapfiles.add(realspacefile)
        }

        RealSpaceData realspacedata = new RealSpaceData(systems[0], systems[0].getProperties(), systems[0].getParallelTeam(),
                mapfiles.toArray(new ffx.realspace.parsers.RealSpaceFile[mapfiles.size()]))

        aFuncts.energy(systems[0])
        
        RefinementMinimize refinementMinimize = new RefinementMinimize(realspacedata, RefinementMode.COORDINATES)

        if (eps < 0.0) {
            eps = 1.0
        }
        
        logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter)
        refinementMinimize.minimize(eps, maxiter)

        aFuncts.energy(systems[0])
        aFuncts.saveAsPDB(systems, new File(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb"))
        
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
