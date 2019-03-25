package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.AnnealOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.optimize.SimulatedAnnealing
import ffx.potential.MolecularAssembly

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

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
    @Parameters(arity = "1", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    private SimulatedAnnealing simulatedAnnealing = null;

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    @Override
    Anneal run() {

        if (!init()) {
            return this
        }

        dynamics.init()

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath();
        }

        logger.info("\n Running simulated annealing on " + modelFilename + "\n")

        simulatedAnnealing = new SimulatedAnnealing(activeAssembly,
                activeAssembly.getPotentialEnergy(), activeAssembly.getProperties(),
                algorithmListener, dynamics.thermostat, dynamics.integrator)

        simulatedAnnealing.anneal(anneal.upper, anneal.low, anneal.windows, dynamics.steps, dynamics.dt)

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        String fileName = FilenameUtils.getName(modelFilename)
        String ext = FilenameUtils.getExtension(fileName)
        fileName = FilenameUtils.removeExtension(fileName)

        String dirName = FilenameUtils.getFullPath(saveDir.getAbsolutePath())

        if (ext.toUpperCase().contains("XYZ")) {
            algorithmFunctions.saveAsXYZ(activeAssembly, new File(dirName + fileName + ".xyz"))
        } else {
            algorithmFunctions.saveAsPDB(activeAssembly, new File(dirName + fileName + ".pdb"))
        }
        return this
    }

    public SimulatedAnnealing getAnnealing() {
        return simulatedAnnealing;
    }
}
