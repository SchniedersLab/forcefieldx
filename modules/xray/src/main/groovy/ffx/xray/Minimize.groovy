
package ffx.xray

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.potential.MolecularAssembly
import ffx.xray.CrystalReciprocalSpace.SolventModel
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.parsers.DiffractionFile

/**
 * The X-ray Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Minimize [options] &lt;filename&gt;
 */
class Minimize extends Script {

    /**
     * Options for the Anneal Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.Minimize [options] &lt;filename [file2...]&gt;
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
        @Option(shortName='e', longName='eps', defaultValue='-1.0', description='RMS gradient convergence criteria (default of -1 automatically determines eps based on refinement type)') double eps;
        /**
         * -i or --iterations Maximum number of optimization iterations.
         */
        @Option(shortName='i', longName='iterations ', defaultValue='1000', description=' Maximum number of optimization iterations.') int iterations;
        /**
         * -r or --mode sets the desired refinement mode
         * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
         */
        @Option(shortName = 'r', longName = 'mode', convert = { s -> return RefinementMinimize.parseMode(s) }, defaultValue = 'COORDINATES',
                description = 'Refinement mode: coordinates, bfactors, occupancies.')
        RefinementMode mode
        /**
         * -t or --threeStage Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true).
         */
        @Option(shortName='t', longName='threeStage', defaultValue='false',
                description='Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true)') boolean threeStage;
        /**
         * -E or --eps3 RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determine eps for each stage).
         */
        @Option(shortName = 'E', longName = 'eps3', defaultValue = '-1.0,-1.0,-1.0', numberOfArguments = 3, valueSeparator = ',',
                description = 'RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determine eps for each stage).')
        double[] eps3
        /**
         * -D or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
         */
        @Option(shortName = 'D', longName = 'data', defaultValue = '', numberOfArguments = 3, valueSeparator = ',',
                description = 'Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.')
        String[] data
        /**
         * The final arguments should be a PDB filename and data filename (CIF or MTZ).
         */
        @Unparsed(description = "PDB file and a CIF or MTZ file.")
        List<String> filenames
    }

    def run() {

        def cli = new CliBuilder()
        cli.name = "ffxc xray.Minimize"

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

        // Do 3 stage refinement (coordinates, then B, then occupancies).
        boolean threestage = options.threeStage

        // RMS gradient convergence criteria for three stage refinement
        double coordeps = options.eps3[0];
        double beps = options.eps3[1];
        double occeps = options.eps3[2];

        // Maximum number of refinement cycles.
        int maxiter = options.iterations

        // Type of refinement.
        RefinementMode refinementmode = options.mode

        // Suffix to append to output data
        String suffix = "_refine"

        String modelfilename = null
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

        // Set up diffraction data (can be multiple files)
        List diffractionfiles = new ArrayList()
        if (arguments.size() > 1) {
            DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        if (options.data) {
            for (int i = 0; i < options.data.size(); i += 3) {
                double wA = Double.parseDouble(options.data[i + 1])
                boolean neutron = Boolean.parseBoolean(options.data[i + 2])
                DiffractionFile diffractionfile = new DiffractionFile(options.data[i], wA, neutron)
                diffractionfiles.add(diffractionfile)
            }
        }

        if (diffractionfiles.size() == 0) {
            DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(),
                SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

        diffractiondata.scaleBulkFit()
        diffractiondata.printStats()

        aFuncts.energy(systems[0]);

        if (threestage) {
            RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES)
            if (coordeps < 0.0) {
                coordeps = refinementMinimize.getEps();
            }
            logger.info("\n RMS gradient convergence criteria: " + coordeps + " max number of iterations: " + maxiter)
            refinementMinimize.minimize(coordeps, maxiter)
            diffractiondata.scaleBulkFit()
            diffractiondata.printStats()

            aFuncts.energy(systems[0])

            refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.BFACTORS);
            if (beps < 0.0) {
                beps = refinementMinimize.getEps()
            }
            logger.info("\n RMS gradient convergence criteria: " + beps + "\n Maximum number of iterations: " + maxiter + "\n")
            refinementMinimize.minimize(beps, maxiter)
            diffractiondata.scaleBulkFit()
            diffractiondata.printStats()

            if (diffractiondata.getAltResidues().size() > 0
                    || diffractiondata.getAltMolecules().size() > 0) {
                refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.OCCUPANCIES)
                if (occeps < 0.0) {
                    occeps = refinementMinimize.getEps()
                }
                logger.info("\n RMS gradient convergence criteria: " + occeps + " max number of iterations: " + maxiter)
                refinementMinimize.minimize(occeps, maxiter)
                diffractiondata.scaleBulkFit()
                diffractiondata.printStats()
            } else {
                logger.info("Occupancy refinement not necessary, skipping")
            }
        } else {
            RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, refinementmode)
            if (eps < 0.0) {
                eps = refinementMinimize.getEps()
            }
            logger.info("\n RMS gradient convergence criteria: " + eps + " max number of iterations: " + maxiter)
            refinementMinimize.minimize(eps, maxiter)
            diffractiondata.scaleBulkFit()
            diffractiondata.printStats()
        }

        aFuncts.energy(systems[0])

        diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb")
        diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + suffix + ".mtz")
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