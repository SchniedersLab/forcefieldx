
package ffx.xray

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.algorithms.MolecularDynamics
import ffx.algorithms.integrators.Integrator
import ffx.algorithms.integrators.IntegratorEnum
import ffx.algorithms.thermostats.Thermostat
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.potential.MolecularAssembly
import ffx.xray.CrystalReciprocalSpace.SolventModel
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.parsers.DiffractionFile

/**
 * The X-ray Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Dynamics [options] &lt;filename&gt;
 */
class Dynamics extends Script {

    /**
     * Options for the Anneal Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.Dynamics [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -n or --steps sets the number of molecular dynamics steps (default is 1 nsec).
         */
        @Option(shortName = 'n', defaultValue = '1000000', description = 'Number of molecular dynamics steps')
        int steps
        /**
         * -d or --dt Time step in femtosceonds (1.0).
         */
        @Option(shortName = 'd', longName = 'dt', defaultValue = '1.0', description = 'Time step (fsec).')
        double d
        /**
         * -l or --log sets the thermodynamics logging frequency in picoseconds (0.1 psec default).
         */
        @Option(shortName = 'l', longName = 'log', defaultValue = '0.25', description = 'Interval to report thermodynamics (psec).')
        double log;
        /**
         * -w or --write sets snapshot save frequency in picoseconds (1.0 psec default).
         */
        @Option(shortName = 'w', longName = 'write', defaultValue = '10.0', description = 'Interval to write out coordinates (psec).')
        double write
        /**
         * -t or --temperature sets the simulation temperature (Kelvin).
         */
        @Option(shortName = 't', longName = 'temperature', defaultValue = '298.15', description = 'Temperature (Kelvin)')
        double temperature
        /**
         * -b or --thermostat sets the desired thermostat [Adiabatic, Berendsen, Bussi].
         */
        @Option(shortName = 't', longName = 'thermostat', convert = { s -> return Thermostat.parseThermostat(s) }, defaultValue = 'Berendsen',
                description = 'Thermostat: Adiabatic, Berendsen or Bussi.')
        ThermostatEnum thermostat
        /**
         * -i or --integrator sets the desired integrator [Beeman, RESPA, Stochastic].
         */
        @Option(shortName = 'i', longName = 'integrator', convert = { s -> return Integrator.parseIntegrator(s) }, defaultValue = 'Beeman',
                description = 'Integrator: Beeman, RESPA or Stochastic.')
        IntegratorEnum integrator
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
        cli.name = "ffxc xray.Dynamics"

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

        String modelfilename
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelfilename = active.getFile()
        }

        logger.info("\n Running dynamics on " + modelfilename)

        MolecularAssembly[] systems = aFuncts.open(modelfilename)

        // Set up diffraction data (can be multiple files)
        List diffractionfiles = new ArrayList()
        if (arguments.size() > 1) {
            DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        if (options.data) {
            for (int i=0; i<options.data.size(); i+=3) {
                double wA = Double.parseDouble(options.data[i+1])
                boolean neutron = Boolean.parseBoolean(options.data[i+2])
                DiffractionFile diffractionfile = new DiffractionFile(options.data[i], wA, neutron)
                diffractionfiles.add(diffractionfile)
            }
        }

        RefinementMode refinementmode = options.mode

        if (diffractionfiles.size() == 0) {
            DiffractionFile diffractionfile = new DiffractionFile(systems[0], 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        DiffractionData diffractiondata = new DiffractionData(systems[0], systems[0].getProperties(),
                SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

        diffractiondata.scaleBulkFit()
        diffractiondata.printStats()

        RefinementEnergy refinementEnergy = new RefinementEnergy(diffractiondata, refinementmode)

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn")
        if (!dyn.exists()) {
            dyn = null
        }
        MolecularDynamics molDyn = new MolecularDynamics(systems[0], refinementEnergy, systems[0].getProperties(), refinementEnergy,
                options.thermostat, options.integrator)
        refinementEnergy.setThermostat(molDyn.getThermostat())

        boolean initVelocities = true;
        molDyn.dynamic(options.steps, options.d, options.log, options.write, options.temperature, initVelocities, dyn);
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
