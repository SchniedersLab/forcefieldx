
package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.integrators.Integrator
import ffx.algorithms.integrators.IntegratorEnum
import ffx.algorithms.thermostats.Thermostat
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.parameters.ForceField

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

@Command(description = " Run dynamics on a system using OpenMM.", name = "ffxc DynamicsOpenMM")
class DynamicsOpenMM extends AlgorithmsScript{

    
    @Mixin
    DynamicsOptions dynamics;
    
    /**
     * -z or --trajSteps sets the length of the MD trajectory run on the GPU in femtoseconds(defaul is 100 femtoseconds)
     */
    @Option(names = ['-z', '--trajSteps'], paramLabel = '100', 
        description = 'Number of steps for each MD Trajectory in femtoseconds')
    int trajSteps;
    /**
     * --cf or --coeffOfFriction specifies what the coefficient of friction is to be used with Langevin and Brownian integrators
     */ 
    @Option(names = ['--cf', '--coeffOfFriction'], paramLabel = '0.01',
        description = 'Coefficient of friction to be used with the Langevin and Brownian integrators')
    double coeffOfFriction;
    /**
     * -q or --collisionFreq specifies the frequency for particle collision to be used with the Anderson thermostat
     */
    @Option(names = ['-q', '--collisionFreq'], paramLabel = '91.0',
        description = 'Collision frequency to be set when Anderson Thermostat is created: Can be used with Verlet integrator')
    double collisionFreq;

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
        description = "XYZ or PDB input files.")
    private List<String> filenames

    
    def run(){

        if (!init()) {
            return
        }
        
        dynamics.init()

        boolean initVelocities = true;

        /**
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
        }*/
        
        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
        }
        
        File structureFile = new File(FilenameUtils.normalize(modelfilename))
        structureFile = new File(structureFile.getAbsolutePath())
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName())

        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy();
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
        boolean updatesDisabled = activeAssembly.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x, true);

        logger.info("\n Running molecular dynamics on " + modelfilename);

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
        if (!dyn.exists()) {
            dyn = null;
        }

        MolecularDynamics moldyn = MolecularDynamics.dynamicsFactory(activeAssembly, forceFieldEnergy, activeAssembly.getProperties(), sh, dynamics.thermostat, dynamics.integrator)
        //MolecularDynamics molDyn = dynamics.getDynamics(potential, activeAssembly, sh)
        if (moldyn instanceof MolecularDynamicsOpenMM){
            moldyn.setRestartFrequency(dynamics.write);
            moldyn.setFileType(dynamics.fileType);
            moldyn.setIntervalSteps(trajSteps);
            moldyn.dynamic(dynamics.steps, dynamics.dt, dynamics.report, dynamics.write, dynamics.temp, initVelocities, dyn);
        }
        else{
            logger.severe(" Could not start OpenMM molecular dynamics.");
        }

    }
}

