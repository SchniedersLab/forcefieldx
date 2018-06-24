package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.AnnealOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.potential.MolecularAssembly

import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters
import picocli.CommandLine.Command

/**
 * The Anneal script.
 * <br>
 * Usage:
 * <br>
 * ffxc Anneal [options] &lt;filename&gt;
 */
@Command(description = " Run simulated annealing on a system.", name = "ffxc Anneal")
class Anneal extends AlgorithmsScript {

    @Mixin
    DynamicsOptions dynamics

    @Mixin
    AnnealOptions anneal

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files", description = "XYZ or PDB input files.")
    private List<String> filenames

    def run() {

        if (!init()) {
            return
        }

        dynamics.init()

        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath();
        }

        logger.info("\n Running simulated annealing on " + modelfilename + "\n")

        SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(activeAssembly,
                activeAssembly.getPotentialEnergy(), activeAssembly.getProperties(), null,
                dynamics.thermostat, dynamics.integrator)

        simulatedAnnealing.anneal(anneal.upper, anneal.low, anneal.windows, dynamics.steps, dynamics.dt)

        String ext = FilenameUtils.getExtension(modelfilename)
        modelfilename = FilenameUtils.removeExtension(modelfilename)

        if (ext.toUpperCase().contains("XYZ")) {
            algorithmFunctions.saveAsXYZ(new File(modelfilename + ".xyz"));
        } else {
            algorithmFunctions.saveAsPDB(new File(modelfilename + ".pdb"));
        }
    }
}
