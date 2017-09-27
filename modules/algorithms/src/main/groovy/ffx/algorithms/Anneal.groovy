
package ffx.algorithms

// FFX Imports
import ffx.algorithms.Integrator.Integrators
import ffx.algorithms.Thermostat.Thermostats

// Groovy Imports
import groovy.cli.Option
import groovy.cli.Unparsed

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils

/**
 * The Anneal script.
 * <br>
 * Usage:
 * <br>
 * ffxc Anneal [options] &lt;filename&gt;
 */
class Anneal extends Script {

    private static parseThermo(String str) {
        try {
            return Thermostats.valueOf(str.toUpperCase())
        } catch (Exception e) {
            System.err.println(String.format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
            return Thermostats.BERENDSEN;
        }
    }

    private static parseIntegrator(String str) {
        try {
            return Integrators.valueOf(str.toUpperCase())
        } catch (Exception e) {
            System.err.println(String.format(" Could not parse %s as an integrator; defaulting to Beeman.", str));
            return Integrators.BEEMAN;
        }
    }

    /**
     * Options for the Anneal Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Anneal [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.') boolean help
        /**
         * -n or --steps Number of molecular dynamics steps per annealing window (1000).
         */
        @Option(shortName='n', longName='steps', defaultValue='1000', description='Number of MD steps per annealing window.') int n
        /**
         * -d or --dt Time step in femtosceonds (1.0).
         */
        @Option(shortName='d', longName='dt', defaultValue='1.0', description='Time step (fsec).') double d
        /**
         * -w or --windows Number of annealing windows (10).
         */
        @Option(shortName='w', longName='windows', defaultValue='10', description='Number of annealing windows.') int w
        /**
         * -l or --low Low temperature limit in degrees Kelvin (10.0).
         */
        @Option(shortName='l', longName='low', defaultValue='10.0', description='Low temperature limit (Kelvin).') double l
        /**
         * -h or --high High temperature limit in degrees Kelvin (1000.0).
         */
        @Option(shortName='h', longName='high', defaultValue='1000.0', description='High temperature limit (Kelvin).') double h
        /**
         * -b or --thermostat sets the desired thermostat [Adiabatic, Berendsen, Bussi].
         */
        @Option(shortName='t', longName='thermostat', convert = {s -> return parseThermo(s);}, defaultValue='Berendsen',
                description='Thermostat: Adiabatic, Berendsen or Bussi.') Thermostats thermo
        /**
         * -i or --integrator sets the desired integrator [Beeman, RESPA, Stochastic].
         */
        @Option(shortName='i', longName='integrator', convert = {s -> return parseIntegrator(s);}, defaultValue='Beeman',
                description='Integrator: Beeman, RESPA or Stochastic.') Integrators integrate
        /**
         * The final argument should be a filename.
         */
        @Unparsed List<String> filenames
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc Anneal [options] <filename>', header: ' Options:')

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

        // Load the number of molecular dynamics steps at each temperature.
        int steps = options.n

        // Load the number of annealing steps.
        int windows = options.w

        // Load the low temperature end point.
        double low = options.l

        // Load the high temperature end point.
        double high = options.h

        // Load the time steps in femtoseconds.
        double timeStep = options.d

        // Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
        Thermostats thermostat = options.thermo;

        // Integrators [ BEEMAN, RESPA, STOCHASTIC]
        Integrators integrator = options.integrate;

        String filename = null;
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0);
            open(modelfilename);
        } else if (active == null) {
            return cli.usage();
        } else {
            modelfilename = active.getFile();
        }

        logger.info("\n Running simulated annealing on " + modelfilename)
        SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, active.getPotentialEnergy(),
                active.getProperties(), null, thermostat, integrator)
        simulatedAnnealing.anneal(high, low, windows, steps, timeStep)

        String ext = FilenameUtils.getExtension(modelfilename)
        modelfilename = FilenameUtils.removeExtension(modelfilename)

        if (ext.toUpperCase().contains("XYZ")) {
            saveAsXYZ(new File(modelfilename + ".xyz"));
        } else {
            saveAsPDB(new File(modelfilename + ".pdb"));
        }
    }
}