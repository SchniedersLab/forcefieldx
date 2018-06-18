
package ffx.realspace

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
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode

/**
 * The Real Space Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Dynamics [options] &lt;filename&gt;
 */
class Dynamics extends Script {

    /**
     * Options for the Real Space Dynamics Script.
     * <br>
     * Usage:
     * <br>
     * ffxc realspace.Dynamics [options] &lt;filename [file2...]&gt;
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
        @Option(shortName = 'b', longName = 'thermostat', convert = { s -> return Thermostat.parseThermostat(s) }, defaultValue = 'Berendsen',
                description = 'Thermostat: Adiabatic, Berendsen or Bussi.')
        ThermostatEnum thermostat
        /**
         * -i or --integrator sets the desired integrator [Beeman, RESPA, Stochastic].
         */
        @Option(shortName = 'i', longName = 'integrator', convert = { s -> return Integrator.parseIntegrator(s) }, defaultValue = 'Beeman',
                description = 'Integrator: Beeman, RESPA or Stochastic.')
        IntegratorEnum integrator
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
        cli.name = "ffxc realspace.Dynamics"

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

        logger.info("\n Running Real Space Dynamics on " + modelfilename)

        MolecularAssembly[] systems = aFuncts.open(modelfilename)

        // set up real space map data (can be multiple files)
        List mapfiles = new ArrayList();
        if (arguments.size() > 1) {
            RealSpaceFile realspacefile = new RealSpaceFile(arguments.get(1), 1.0);
            mapfiles.add(realspacefile);
        }
        if (options.data) {
            for (int i=0; i<options.data.size(); i+=2) {
                double wA = Double.parseDouble(options.data[i+1]);
                RealSpaceFile realspacefile = new RealSpaceFile(options.data[i], wA);
                mapfiles.add(realspacefile);
            }
        }

        if (mapfiles.size() == 0) {
            RealSpaceFile realspacefile = new RealSpaceFile(systems[0], 1.0);
            mapfiles.add(realspacefile);
        }

        RealSpaceData realspacedata = new RealSpaceData(systems[0], systems[0].getProperties(), systems[0].getParallelTeam(),
                mapfiles.toArray(new RealSpaceFile[mapfiles.size()]));

        aFuncts.energy(systems[0]);

        RefinementEnergy refinementEnergy = new RefinementEnergy(realspacedata, RefinementMode.COORDINATES);

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
        if (!dyn.exists()) {
            dyn = null;
        }
        MolecularDynamics molDyn = new MolecularDynamics(systems[0], refinementEnergy,
                systems[0].getProperties(), refinementEnergy, options.thermostat, options.integrator, null);
        refinementEnergy.setThermostat(molDyn.getThermostat());

        // Reset velocities (ignored if a restart file is given)
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