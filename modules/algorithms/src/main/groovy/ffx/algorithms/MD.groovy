
package ffx.algorithms;

// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly
import ffx.potential.OpenMMForceFieldEnergy
import ffx.potential.parameters.ForceField;

/**
 * The MD script implements stochastic and molecular dynamics algorithms.
 * <br>
 * Usage:
 * <br>
 * ffxc MD [options] &lt;filename&gt; [file2...]
 */
class MD extends Script {

    /**
     * Options for the MD Script.
     * <br>
     * Usage:
     * <br>
     * ffxc MD [options] &lt;filename&gt; [file2...]
     */
    class Options {

        private static parseThermo(String str) {
            try {
                return Thermostats.valueOf(str.toUpperCase());
            } catch (Exception e) {
                System.err.println(String.format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
                return Thermostats.BERENDSEN;
            }
        }

        private static parseIntegrator(String str) {
            try {
                return Integrators.valueOf(str.toUpperCase());
            } catch (Exception e) {
                System.err.println(String.format(" Could not parse %s as an integrator; defaulting to Beeman.", str));
                return Integrators.BEEMAN;
            }
        }

        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help;

        /**
         * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
         */
        @Option(shortName = 'b', longName = 'thermostat', convert = { s -> return parseThermo(s); },
                defaultValue = 'Berendsen', description = 'Thermostat: [Adiabatic / Berendsen / Bussi].')
        Thermostats tstat;
        /**
         * -i or --integrator sets the desired integrator: current choices are Beeman, RESPA, or Stochastic (AKA Langevin dynamics).
         */
        @Option(shortName = 'i', longName = 'integrator', convert = { s -> return parseIntegrator(s); },
                defaultValue = 'Beeman', description = 'Integrator: [Beeman / Respa / Stochastic]')
        Integrators integrator;
        /**
         * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of 2.0 is possible for the RESPA integrator.
         */
        @Option(shortName = 'd', defaultValue = '1.0', description = 'Time discretization step in femtoseconds.')
        double dt;
        /**
         * -r or --report sets the thermodynamics reporting frequency in picoseconds (0.1 psec default).
         */
        @Option(shortName = 'r', longName = 'report', defaultValue = '0.25', description = 'Interval to report thermodynamics (psec).')
        double report;
        /**
         * -w or --write sets snapshot save frequency in picoseconds (1.0 psec default).
         */
        @Option(shortName = 'w', longName = 'write', defaultValue = '10.0', description = 'Interval to write out coordinates (psec).')
        double write;
        /**
         * -t or --temperature sets the simulation temperature (Kelvin).
         */
        @Option(shortName = 't', longName = 'temperature', defaultValue = '298.15', description = 'Temperature (Kelvin)')
        double temp;
        /**
         * -n or --steps sets the number of molecular dynamics steps (default is 1 nsec).
         */
        @Option(shortName = 'n', defaultValue = '1000000', description = 'Number of molecular dynamics steps')
        int steps;
        /**
         * -ld or --minDensity sets a tin box constraint on the barostat, preventing over-expansion of the box (particularly in vapor phase), permitting an analytic correction.
         */
        @Option(shortName = 'ld', longName = 'minDensity', defaultValue = '0.75',
                description = 'Minimum density allowed by the barostat')
        double minDensity;
        /**
         * -hd or --maxDensity sets a maximum density on the barostat, preventing under-expansion of the box.
         */
        @Option(shortName = 'hd', longName = 'maxDensity', defaultValue = '1.6',
                description = 'Maximum density allowed by the barostat')
        double maxDensity;
        /**
         * -sm or --maxSideMove sets the width of proposed crystal side length moves (rectangularly distributed) in Angstroms.
         */
        @Option(shortName = 'sm', longName = 'maxSideMove', defaultValue = '0.25',
                description = 'Maximum side move allowed by the barostat in Angstroms')
        double maxSideMove;
        /**
         * -am or --maxAngleMove sets the width of proposed crystal angle moves (rectangularly distributed) in degrees.
         */
        @Option(shortName = 'am', longName = 'maxAngleMove', defaultValue = '0.5',
                description = 'Maximum angle move allowed by the barostat in degrees')
        double maxAngleMove;
        /**
         * -mi or --meanInterval sets the mean number of MD steps (Poisson distribution) between barostat move proposals.
         */
        @Option(shortName = 'mi', longName = 'meanInterval', defaultValue = '10',
                description = 'Mean number of MD steps between barostat move proposals.')
        int meanInterval;
        /**
         * -p or --npt Specify use of a MC Barostat at the given pressure (default 1.0 atm).
         */
        @Option(shortName = 'p', longName = 'npt', defaultValue = '0',
                description = 'Specify use of a MC Barostat at the given pressure (default of 0 = disabled)')
        double pressure;
        /**
         * -f or --file Choose the file type to write [PDB/XYZ].
         */
        @Option(shortName = 'f', longName = 'file', defaultValue = 'PDB',
                description = 'Choose file type to write [PDB/XYZ].')
        String fileType;

        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed
        List<String> filenames;
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc MD [options] <filename> [file2...]', header: ' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }

        // Load the number of molecular dynamics steps.
        nSteps = options.steps

        // Write dyn interval in picoseconds
        restartFrequency = options.write
        saveInterval = options.write

        fileType = options.fileType;

        // Load the time steps in femtoseconds.
        timeStep = options.dt

        // Report interval in picoseconds.
        printInterval = options.report

        // Temperature in degrees Kelvin.
        temperature = options.temp

        // Pressure in atm.
        pressure = options.pressure

        // Minimum density
        minDensity = options.minDensity

        // Maximum density
        maxDensity = options.maxDensity

        // Max side move
        maxSideMove = options.maxSideMove

        // Max angle move
        maxAngleMove = options.maxAngleMove

        // Mean interval to apply the Barostat
        meanInterval = options.meanInterval

        // Thermostat.
        thermostat = options.tstat

        // Integrator.
        integrator = options.integrator

        List<String> arguments = options.filenames;

        String modelfilename = null;
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0);
            open(modelfilename);
        } else if (active == null) {
            return cli.usage();
        } else {
            modelfilename = active.getFile();
        }

        Potential potential = active.getPotentialEnergy();
        logger.info(" Starting energy (before .dyn restart loaded):");
        boolean updatesDisabled = active.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);

        potential.energy(x, true);

        if (pressure > 0) {
            if (potential instanceof OpenMMForceFieldEnergy) {
                logger.warning(" NPT with OpenMM acceleration is still experimental and may not function correctly.");
            }
            logger.info(String.format(" Running NPT dynamics at pressure %7.4g", pressure));
            MolecularAssembly molecAssem = active; // Why this is needed is utterly inexplicable to me.
            CrystalPotential cpot = (CrystalPotential) potential;
            Barostat barostat = new Barostat(molecAssem, cpot);
            barostat.setPressure(pressure);
            barostat.setMaxDensity(maxDensity);
            barostat.setMinDensity(minDensity);
            barostat.setMaxSideMove(maxSideMove);
            barostat.setMaxAngleMove(maxAngleMove);
            barostat.setMeanBarostatInterval(meanInterval);
            potential = barostat;
        }

        logger.info("\n Running molecular dynamics on " + modelfilename);

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
        if (!dyn.exists()) {
            dyn = null;
        }

        MolecularDynamics molDyn = new MolecularDynamics(active, potential, active.getProperties(), sh, thermostat, integrator);
        molDyn.setFileType(fileType);
        molDyn.setRestartFrequency(restartFrequency);
        molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
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
