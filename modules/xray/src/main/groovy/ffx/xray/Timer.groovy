package ffx.xray

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.CrystalReciprocalSpace.SolventModel
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.parsers.DiffractionFile

/**
 * The X-ray Timer script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Timer [options] &lt;filename&gt;
 */
class Timer extends Script {

    /**
     * Options for the Timer Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.Timer [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -n or --iterations to set the number of iterations
         */
        @Option(longName='iterations', shortName='n', defaultValue='5', description='Number of iterations.')
        int iterations
        /**
         * -c or --threads to set the number of SMP threads (the default of 0 specifies use of all CPU cores)
         */
        @Option(longName='threads', shortName='c', defaultValue='0', description='Number of SMP threads (the default of 0 specifies use of all CPU cores)')
        int threads
        /**
         * -g or --gradient to ignore computation of the atomic coordinates gradient
         */
        @Option(longName='gradient', shortName='g', defaultValue='false', description='Ignore computation of the atomic coordinates gradient')
        boolean gradient
        /**
         * -r or --mode sets the desired refinement mode
         * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
         */
        @Option(shortName = 'r', longName = 'mode', convert = { s -> return RefinementMinimize.parseMode(s) }, defaultValue = 'COORDINATES',
                description = 'Refinement mode: coordinates, bfactors, occupancies.')
        RefinementMode mode
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
        cli.name = "ffxc xray.Timer"

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

        String modelfilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelfilename = active.getFile()
        }

        logger.info("\n Running xray.Timer on " + modelfilename)

        // Set the number of threads (needs to be done before opening the files).
        if (options.threads > 0) {
            System.setProperty("pj.nt", Integer.toString(options.threads))
        }

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

        aFuncts.energy(systems[0])

        // Type of refinement.
        RefinementMode refinementmode = options.mode

        RefinementEnergy refinementEnergy = new RefinementEnergy(diffractiondata, refinementmode)
        int n = refinementEnergy.getNumberOfVariables()
        double[] x = new double[n]
        double[] g = new double[n]
        refinementEnergy.getCoordinates(x)

        Potential energy = refinementEnergy.getDataEnergy()

        int nEvals = options.iterations

        for (int i=0; i<nEvals; i++) {
            long time = -System.nanoTime()

            if (options.gradient) {
                energy.energyAndGradient(x,g)
            } else {
                energy.energy(x)
            }

            time += System.nanoTime()

            logger.info(String.format(" Time %12.8f", time * 1.0e-9))
        }
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