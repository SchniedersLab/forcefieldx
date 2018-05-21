
package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.algorithms.integrators.Integrator
import ffx.algorithms.integrators.IntegratorEnum
import ffx.algorithms.thermostats.Thermostat
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.potential.ForceFieldEnergy
import ffx.potential.parameters.ForceField

class DynamicsOpenMM extends Script{

    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.') boolean help;
        /**
         * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
         */
        @Option(shortName = 'b', longName = 'thermostat', convert = { s -> return Thermostat.parseThermostat(s); },
            defaultValue = 'Berendsen', description = 'Thermostat: [Berendsen, Anderson]')
        ThermostatEnum tstat;
        /**
         * -i or --integrator sets the desired integrator: current choices are Beeman, RESPA, Velocity Verlet, or Stochastic (AKA Langevin dynamics).
         */
        @Option(shortName = 'i', longName = 'integrator', convert = { s -> return Integrator.parseIntegrator(s); },
            defaultValue = 'VelocityVerlet', description = 'Integrator: [Beeman / Respa / Stochastic / VelocityVerlet ]')
        IntegratorEnum integrator;
        /**
         * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of 2.0 is possible for the RESPA integrator.
         */
        @Option(shortName = 'd', defaultValue = '1.0', description = 'Time discretization step in femtoseconds.')
        double dt;
        /**
         * -z or -intervalSteps sets the length of the MD trajectory run on the GPU in femtoseconds(defaul is 100 femtoseconds)
         */
        @Option(shortName = 'z', longName = 'trajSteps', defaultValue = '100', description = 'Number of steps for each MD Trajectory in femtoseconds')
        int trajSteps;
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
        @Option(shortName = 'n', longName = 'steps', defaultValue = '1000000', description = 'Number of molecular dynamics steps')
        int steps;
        /**
         * -f or --file Choose the file type to write [PDB/XYZ].
         */
        @Option(shortName = 'f', longName = 'file', defaultValue = 'XYZ',
            description = 'Choose file type to write [PDB/XYZ].')
        String fileType;
        /**
         * -cf or -coeffOfFriction Specifies what the coefficient of friction is to be used with Langevin and Brownian integrators
         */ 
        @Option(shortName = 'cf', longName = 'coeffOfFriction', defaultValue = '0.01',
            description = 'Coefficient of friction to be used with the Langevin and Brownian integrators')
        double coeffOfFriction;
        /**
         * -q or -collisionFreq Specifies the frequency for particle collision to be used with the Anderson thermostat
         */
        @Option(shortName = 'q', longName = 'collisionFreq', defaultValue = '91.0',
            description = 'Collision frequency to be set when Anderson Thermostat is created: Can be used with Verlet integrator')
        double collisionFreq;

        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed
        List<String> filenames;
    }
    
    def run(){
        logger.info("in the run method");
        def cli = new CliBuilder(usage: ' ffxc DynamicsOpenMM [options] <filename> [file2...]', header: ' Options:');

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
        
        // Number of MD steps per trajectory (offloaded to GPU)
        
        intervalSteps = options.trajSteps

        // Report interval in picoseconds.
        printInterval = options.report

        // Temperature in degrees Kelvin.
        temperature = options.temp

        // Thermostat.
        thermostat = options.tstat

        // Integrator.
        integrator = options.integrator
        
        // Coefficient of Friction
        frictionCoeff = options.coeffOfFriction
        
        // Collision Frequency
        collisionFreq = options.collisionFreq

        boolean initVelocities = true;

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

        ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
        switch (forceFieldEnergy.getPlatform()) {
            case ForceFieldEnergy.Platform.OMM:
            case ForceFieldEnergy.Platform.OMM_CUDA:
            case ForceFieldEnergy.Platform.OMM_OPENCL:
            case ForceFieldEnergy.Platform.OMM_OPTCPU:
            case ForceFieldEnergy.Platform.OMM_REF:
                logger.fine(" Platform is appropriate for OpenMM Dynamics.");
                break;
            case ForceFieldEnergy.Platform.FFX:
            default:
                logger.severe(String.format(" Platform %s is inappropriate for OpenMM dynamics. Please explicitly specify an OpenMM platform.", forceFieldEnergy.getPlatform()));
                break;
        }


        logger.info(" Starting energy (before .dyn restart loaded):");
        boolean updatesDisabled = active.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x, true);

        logger.info("\n Running molecular dynmaics on " + modelfilename);

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
        if (!dyn.exists()) {
            dyn = null;
        }

        MolecularDynamics moldyn = MolecularDynamics.dynamicsFactory(active, forceFieldEnergy, active.getProperties(), sh, thermostat, integrator)
        if (moldyn instanceof MolecularDynamicsOpenMM){
            moldyn.setRestartFrequency(restartFrequency);
            moldyn.setFileType(fileType);
            moldyn.setIntervalSteps(intervalSteps);
            moldyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
        }
        else{
            logger.severe(" Could not start OpenMM molecular dynamics.");
        }

    }
}

