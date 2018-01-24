/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.algorithms.integrators.Integrator.Integrators
import ffx.algorithms.thermostats.Thermostat.Thermostats
import ffx.crystal.CrystalPotential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.nonbonded.CoordRestraint
import ffx.potential.parameters.ForceField

class DynamicsOpenMM extends Script{
    /*
    // Number of molecular dynamics steps
    int nSteps = 1000000;

    // Time step in femtoseconds.
    double timeStep = 1.0;

    int intervalSteps = 100;

    // Frequency to print out thermodynamics information in picoseconds.
    double printInterval = 100.0;

    // Frequency to save out coordinates in picoseconds.
    double saveInterval = 2.0;

    // Temperature in degrees Kelvin.
    double temperature = 298.15;

    // Thermostats [ ANDERSON ]
    Thermostats thermostat = null;

    // Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET] for regular MolecularDynamics
    Integrators integrator = null;

    // Integrator string, sets integrator for MolecularDynamicsOpenMM
    //String stringIntegrator = "VERLET";

    // Reset velocities (ignored if a restart file is given)
    boolean initVelocities = true;

    // Interval to write out restart file (psec)
    double restartFrequency = 2.0;

    // Coefficient of Friction for Langevin itegrator
    double frictionCoeff = 91.0;

    // Collision Frequency for Anderson Thermostat used in conjucntion with Verlet Integrator
    double collisionFreq = 0.01;

    List<CoordRestraint> restraints = null;

    // Pressure in atm; invalid (such as negative) values mean do not run NPT.
    double pressure = -1.0;

    // Set min and max density and set the maximum MC move size
    int meanInterval = 10;
    double minDensity = 0.5;
    double maxDensity = 1.5;
    double maxSideMove = 0.25;
    double maxAngleMove = 1.0;

    int numThreads = 1;

    // File type of snapshots.
    String fileType = "PDB";
    */
    class Options {
        // Create the command line parser.
        //def cli = new CliBuilder(usage:' ffxc md [options] <filename>');
        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.') boolean help;
        /**
         * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
         */
        @Option(shortName = 'b', longName = 'thermostat', convert = { s -> return ffx.algorithms.thermostats.Thermostat.parseThermostat(s); },
            defaultValue = 'Berendsen', description = 'Thermostat: [Berendsen, Anderson]')
        Thermostats tstat;
        /**
         * -i or --integrator sets the desired integrator: current choices are Beeman, RESPA, Velocity Verlet, or Stochastic (AKA Langevin dynamics).
         */
        @Option(shortName = 'i', longName = 'integrator', convert = { s -> return ffx.algorithms.integrators.Integrator.parseIntegrator(s); },
            defaultValue = 'VelocityVerlet', description = 'Integrator: [Beeman / Respa / Stochastic / VelocityVerlet ]')
        Integrators integrator;
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
         * -nt or --numThreads Specifies the number of threads to be used per process
         */ 
        @Option(shortName = 'nt', longName = 'numThreads', defaultValue = '1',
            description = 'Number of threads to use per process')
        int numThreads;
        
        //cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
        //cli.z(longOpt: 'intervalSteps', args: 1, argName:'100', 'Number of steps for each MD trajectory (fsec)');
        //cli.si(longOpt:'integrate', args:1, argName:'Verlet', 'Integrator: [Langevin / Verlet]');
        //cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
        //cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
        //cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
        //cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
        //cli.a(longOpt:'pressure', args:1, argName:'Disabled', 'If numeric argument, then sets NPT pressure, else do not run NPT.')
        //cli.ld(longOpt:'minDensity', args:1, argName:'0.5', 'Minimum density allowed by the barostat.');
        //cli.hd(longOpt:'maxDensity', args:1, argName:'1.5', 'Maximum density allowed by the barostat.');
        //cli.sm(longOpt:'maxSideMove', args:1, argName:'0.25', 'Maximum side move allowed by the barostat.');
        //cli.am(longOpt:'maxAngleMove', args:1, argName:'1.0', 'Maximum angle move allowed by the barostat.');
        //cli.mi(longOpt:'meanInterval', args:1, argName:'10', 'Mean number of MD steps between applications of the barostat.');
        //cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
        //cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
        //cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
        //cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
        //cli.cf(longOpt:'Coefficient of friction', args:1, argName:'1.0', 'Coefficient of friction to be used with Langevin and Brownian integrators')
        //cli.q(longOpt:'Collision frequency in 1/ps', args:1, argName:'10.0', 'Collision frequency to be set when Anderson Thermostat is created: Verlet integrator requires this')
        //def options = cli.parse(args);
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
        
        // Coefficent of Friction
        
        frictionCoeff = options.coeffOfFriction
        
        // Collision Frequency
        collisionFreq = options.collisionFreq
        
        // Number of Threads
        
        numThreads = options.numThreads

        boolean initVelocities = true;
        
        logger.info("done parsing arguments and storing values");
        
        /*
        if (options.h) {
            return cli.usage();
        }

        // Thermostat.
        if (options.b) {
            try {
                thermostat = Thermostats.valueOf(options.b.toUpperCase());
            } catch (Exception e) {
                thermostat = null;
            }
        }

        // Load the time steps in femtoseconds.
        if (options.d) {
            timeStep = Double.parseDouble(options.d);
        }

        if (options.z) {
            intervalSteps = Integer.parseInt(options.z);
        }

        if (options.i) {
            try {
                integrator = Integrators.valueOf(options.i.toUpperCase());
            } catch (Exception e) {
                logger.warning(String.format(" Exception in parsing integrator %s: %s", options.i, e.toString()));
                integrator = null;
            }
        }

        // Report interval in picoseconds.
        if (options.l) {
            printInterval = Double.parseDouble(options.l);
        }

        // Load the number of molecular dynamics steps.
        if (options.n) {
            nSteps = Integer.parseInt(options.n);
        }

        if (options.p) {
            System.setProperty("polarization", options.p);
        }

        // Pressure in atm.
        if (options.a) {
            pressure = Double.parseDouble(options.a);
        }

        // Minimum density
        if (options.ld) {
            minDensity = Double.parseDouble(options.ld);
        }
        // Maximum density
        if (options.hd) {
            maxDensity = Double.parseDouble(options.hd);
        }
        // Max side move
        if (options.sm) {
            maxSideMove = Double.parseDouble(options.sm);
        }
        // Max angle move
        if (options.am) {
            maxAngleMove = Double.parseDouble(options.am);
        }
        // Max angle move
        if (options.mi) {
            meanInterval = Integer.parseInt(options.mi);
        }

        // Temperature in degrees Kelvin.
        if (options.t) {
            temperature = Double.parseDouble(options.t);
        }

        // Write snapshot interval in picoseconds.
        if (options.w) {
            saveInterval = Double.parseDouble(options.w);
        }

        // Write dyn interval in picoseconds
        if (options.s) {
            restartFrequency = Double.parseDouble(options.s);
        }

        //
        if (options.f) {
            fileType = options.f.toUpperCase();
        }

        // Coefficient of friction to be used with Langevin and Brownian integrators
        if (options.cf) {
            frictionCoeff = Double.parseDouble(options.cf);
        }

        //Collision frequency (in 1/ps) set when using Anderson Thermostat: Verlet integrator requires this
        if (options.q) {
            collisionFreq = Double.parseDouble(options.q);
        }
        */


        // List<String> arguments = options.arguments();
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
        logger.info(" Starting energy (before .dyn restart loaded):");
        boolean updatesDisabled = active.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x, true);

        if (pressure > 0) {
            if (potential instanceof ffx.potential.ForceFieldEnergyOpenMM) {
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

        logger.info("\n Running molecular dynmaics on " + modelfilename);

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
        if (!dyn.exists()) {
            dyn = null;
        }

        logger.info(" about to create dynamics object");

        //if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
        MolecularDynamics moldyn = MolecularDynamics.dynamicsFactory(active, forceFieldEnergy, active.getProperties(), sh, thermostat, integrator)
        if (moldyn instanceof MolecularDynamicsOpenMM){
            moldyn.setRestartFrequency(restartFrequency);
            //moldyn.setIntegratorString(stringIntegrator);
            //moldyn.setFrictionCoefficient(frictionCoeff);
            //moldyn.setCollisionFrequency(collisionFreq);
            moldyn.setIntervalSteps(intervalSteps);
            moldyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
        }
        else{
            logger.severe(" Could not start OpenMM molecular dynamics");
        }
    }
}
//}
