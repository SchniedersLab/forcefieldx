package ffx.algorithms;

// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.potential.ForceFieldEnergy;

/**
 * The Histogram script prints out a (TT-)OSRW histogram from a *.his file.
 * <br>
 * Usage:
 * <br>
 * ffxc Histogram [options] &lt;filename&gt;
 */
class Histogram extends Script {

    /**
     * Options for the Histogram Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Histogram [options] &lt;filename&gt;
     */
    class Options {

        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -p or --pmf Print out potential of mean force information.
         */
        @Option(shortName = 'p', longName = 'pmf', defaultValue = 'false', description = 'Print out potential of mean force information.')
        boolean pmf
        /**
         * -u or --untempered Histogram for untempered OSRW.
         */
        @Option(shortName = 'u', longName = 'untempered', defaultValue = 'false', description = 'Histogram for untempered OSRW.')
        boolean untempered

        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed
        List<String> filenames;
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc Histogram [options] <filename>', header: ' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }

        List<String> arguments = options.filenames;

        // Read in command line file.
        String filename = arguments.get(0);

        // Print out PMF information.
        boolean pmf = options.pmf;

        // The default Histogram is tempered (i.e. the counts are floating point values).
        boolean untempered = options.untempered

        println("\n Evaluating Histogram for " + filename);

        File structureFile = new File(FilenameUtils.normalize(filename));
        structureFile = new File(structureFile.getAbsolutePath());
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
        File histogramRestart = new File(baseFilename + ".his");
        File lambdaRestart = null;

        open(filename);

        // Get a reference to the active system's ForceFieldEnergy and atom array.
        ForceFieldEnergy energy = active.getPotentialEnergy();

        // Print the current energy
        energy.energy(true, true);

        // Asychronous communication between walkers.
        boolean asynchronous = false;

        // Time step in femtoseconds.
        double timeStep = 1.0;

        // Frequency to log thermodynamics information in picoseconds.
        double printInterval = 1.0;

        // Frequency to write out coordinates in picoseconds.
        double saveInterval = 100.0;

        // Temperture in degrees Kelvin.
        double temperature = 298.15;

        if (!untempered) {
            TransitionTemperedOSRW ttosrw = new TransitionTemperedOSRW(energy, energy, lambdaRestart, histogramRestart, active.getProperties(),
                    temperature, timeStep, printInterval, saveInterval, asynchronous, sh);
            if (pmf) {
                ttosrw.evaluatePMF();
            }
        } else {
            // Wrap the potential energy inside an OSRW instance.
            OSRW osrw = new OSRW(energy, energy, lambdaRestart, histogramRestart, active.getProperties(),
                    temperature, timeStep, printInterval, saveInterval, asynchronous, sh);
            if (pmf) {
                osrw.evaluatePMF();
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