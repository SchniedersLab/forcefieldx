//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.groovy

import edu.rit.pj.ParallelTeam
import ffx.crystal.CrystalPotential
import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.EstimateBootstrapper
import ffx.numerics.estimator.MBARFilter
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.SystemFilter
import ffx.potential.bonded.LambdaInterface
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option
import picocli.CommandLine.Mixin
import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.estimator.MultistateBennettAcceptanceRatio
import ffx.numerics.estimator.MultistateBennettAcceptanceRatio.*
import static java.lang.String.format
import org.apache.commons.io.FilenameUtils

/**
 * Perform MBAR calculation or necessary energy evaluations. See PhEnergy for CpHMD energy evaluations.
 * <br>
 * Usage:
 * <br>
 * ffxc test.MBAR [options] &lt;path&gt
 */
@Command(description = " Evaluates a free energy change with the Multistate Bennett Acceptance Ratio algorithm.",
        name = "MBAR")
class MBAR extends AlgorithmsScript {

    @Mixin
    private AlchemicalOptions alchemical

    @Mixin
    private TopologyOptions topology

    @Option(names = ["--bar"], paramLabel = "true",
            description = "Run BAR calculation as well using a subset of the MBAR data.")
    boolean bar = true

    @Option(names = ["--convergence"], paramLabel = "false",
            description = "Run MBAR multiple times across different time periods of the data to examine the change in FE over time.")
    boolean convergence = false

    @Option(names = ["--numBootstrap", "--nb"], paramLabel = "0",
            description = "Number of bootstrap samples to use.")
    int numBootstrap = 0

    @Option(names = ["--numLambda", "--nL", "--nw"], paramLabel = "-1",
            description = "Required for lambda energy evaluations. Ensure numLambda is consistent with the trajectory lambdas, i.e. gaps between traj can be filled easily. nL >> nTraj is recommended.")
    int numLambda = -1

    @Option(names = ["--seed"], paramLabel = "BAR",
            description = "Seed MBAR calculation with this: ZEROS, ZWANZIG, BAR. Fallback to ZEROS if input is does not or is unlikely to converge.")
    String seedWith = "BAR"

    /**
     * The path to MBAR/BAR files.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'Path to MBAR/BAR files to analyze or an PDB/XYZ in a directory with archive(s).')
    List<String> fileList = null

    public MultistateBennettAcceptanceRatio mbar = null
    int numTopologies

    /**
     * MBAR Constructor.
     */
    MBAR() {
        this(new Binding())
    }

    /**
     * MBAR Constructor.
     * @param binding The Groovy Binding to use.
     */
    MBAR(Binding binding) {
        super(binding)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    MBAR run() {
        if (!init()) {
            return this
        }

        /**
         * Load user supplied fileNames into an array.
         */
        if (fileList == null) {
            logger.severe("No path to MBAR/BAR or trajectory(s) file names specified.")
            return this
        }
        int nFiles = fileList.size()
        String[] fileNames = new String[nFiles]
        File[] files = new File[nFiles]
        for (int i = 0; i < nFiles; i++) {
            fileNames[i] = fileList.get(i)
            files[i] = (new File(fileNames[i])).getAbsoluteFile()
            if (!files[i].exists()) {
                logger.severe("File does not exist: " + fileNames[i])
                return this
            }
        }
        boolean isArc = !files[0].isDirectory()

        // Write MBAR file
        if(isArc){
            if (numLambda == -1){
                logger.severe("numLambda must be specified for lambda energy evaluations.")
                return this
            }
            // Get list of fileNames & check validity
            File parent = files[0].getParentFile() // Run directory
            int window = Integer.parseInt(parent.getName()) // Run name should be int
            File outputDir = new File(parent.getParentFile(), "mbarFiles") // Make mbarFiles
            if(!outputDir.exists()) {
                outputDir.mkdir()
            }
            // Write MBAR file with window number, although this will be reassigned by the file filter based on
            // placement relative to other fileNames with energy values.
            File outputFile = new File(outputDir, "energy_" + window + ".mbar")
            //TODO: Fix atrocious setting of temperatures
            double[][] energies = getEnergyForLambdas(files, numLambda) // Long step!
            MultistateBennettAcceptanceRatio.writeFile(energies, outputFile, 298) // Assume 298 K
            return this
        }

        // Run MBAR calculation
        File path = new File(fileList.get(0))
        if (!path.exists()) {
            logger.severe("Path to MBAR/BAR fileNames does not exist: " + path)
            return this
        }
        if (!path.isDirectory() && !(path.isFile() && path.canRead())) {
            logger.severe("Path to MBAR/BAR fileNames is not accessible: " + path)
            return this
        }
        MBARFilter filter = new MBARFilter(path);
        seedWith = seedWith.toUpperCase()
        SeedType seed = SeedType.valueOf(seedWith) as SeedType
        if (seed == null) {
            logger.severe("Invalid seed type: " + seedWith)
            return this
        }

        mbar = filter.getMBAR(seed as MultistateBennettAcceptanceRatio.SeedType)
        this.mbar = mbar
        if (mbar == null) {
            logger.severe("Could not create MBAR object.")
            return this
        }
        logger.info("\n MBAR Results:")
        logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol", mbar.getFreeEnergy(),
                mbar.getUncertainty()))
        double[] dGs = mbar.getBinEnergies()
        double[] uncertainties = mbar.getBinUncertainties()
        double[][] uncertaintyMatrix = mbar.getDiffMatrix()
        for (int i = 0; i < dGs.length; i++) {
            logger.info(format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
        }
        logger.info("\n MBAR uncertainty between all i & j: ")
        for(int i = 0; i < uncertaintyMatrix.length; i++) {
            StringBuilder sb = new StringBuilder()
            sb.append("    [")
            for(int j = 0; j < uncertaintyMatrix[i].length; j++) {
                sb.append(format(" %6.5f ", uncertaintyMatrix[i][j]))
            }
            sb.append("]")
            logger.info(sb.toString())
        }
        logger.info("\n")
        if(bar){
            try {
                logger.info("\n BAR Results:")
                BennettAcceptanceRatio bar = mbar.getBAR()
                logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol", bar.getFreeEnergy(),
                        bar.getUncertainty()))
                dGs = bar.getBinEnergies()
                uncertainties = bar.getBinUncertainties()
                for (int i = 0; i < dGs.length; i++) {
                    logger.info(format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                }
            } catch (Exception ignored) {
                logger.warning(" BAR calculation failed to converge.")
            }
        }
        logger.info("\n")
        if(numBootstrap != 0) {
            EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar)
            bootstrapper.bootstrap(numBootstrap)
            logger.info("\n MBAR Bootstrap Results from " + numBootstrap + " Samples:")
            logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol", bootstrapper.getTotalFE(),
                    bootstrapper.getTotalUncertainty()))
            dGs = bootstrapper.getFE()
            uncertainties = bootstrapper.getUncertainty()
            for (int i = 0; i < dGs.length; i++) {
                logger.info(format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
            }
            logger.info("\n")
            if (bar) {
                try {
                    logger.info("\n BAR Bootstrap Results:")
                    bootstrapper = new EstimateBootstrapper(mbar.getBAR())
                    bootstrapper.bootstrap(numBootstrap)
                    logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol", bootstrapper.getTotalFE(),
                            bootstrapper.getTotalUncertainty()))
                    dGs = bootstrapper.getFE()
                    uncertainties = bootstrapper.getUncertainty()
                    for (int i = 0; i < dGs.length; i++) {
                        logger.info(format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                    }
                } catch (Exception ignored) {
                    logger.warning(" BAR calculation failed to converge.")
                }
            }
        }

        if (convergence){
            MultistateBennettAcceptanceRatio.FORCE_ZEROS_SEED = true;
            MultistateBennettAcceptanceRatio[] mbarTimeConvergence =
                    filter.getTimeConvergenceMBAR(seed as MultistateBennettAcceptanceRatio.SeedType, 1e-7)
            double[][] dGTime = new double[mbarTimeConvergence.length][]
            for (int i = 0; i < mbarTimeConvergence.length; i++) {
                dGTime[i] = mbarTimeConvergence[i].getBinEnergies()
            }
            MultistateBennettAcceptanceRatio[] mbarPeriodComparison =
                    filter.getPeriodComparisonMBAR(seed as MultistateBennettAcceptanceRatio.SeedType, 1e-7)
            double[][] dGPeriod = new double[mbarPeriodComparison.length][]
            for (int i = 0; i < mbarPeriodComparison.length; i++) {
                dGPeriod[i] = mbarPeriodComparison[i].getBinEnergies()
            }
            logger.info("\n MBAR Time Convergence Results:")
            logger.info(format("     %10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%% ",
                    10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
            for(int i = 0; i < dGTime[0].length; i++) {
                StringBuilder sb = new StringBuilder()
                sb.append(" dG_").append(i).append(": ")
                for(int j = 0; j < dGTime.length; j++) {
                    sb.append(format("%10.4f ", dGTime[j][i]))
                }
                logger.info(sb.toString())
            }
            logger.info("\n")
            logger.info("\n MBAR Period Comparison Results:")
            logger.info(format("     %10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%% ",
                    10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
            for(int i = 0; i < dGTime[0].length; i++) {
                StringBuilder sb = new StringBuilder()
                sb.append(" dG_").append(i).append(": ")
                for(int j = 0; j < dGTime.length; j++) {
                    sb.append(format("%10.4f ", dGTime[j][i]))
                }
                logger.info(sb.toString())
            }
        }
        return this
    }

    private double[][] getEnergyForLambdas(File[] files, int nLambda) {
        // Stuff copied from BAR.groovy
        /*
            Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
            Checking numTopologies == 1 should only be done for scripts that imply
            some sort of lambda scaling.
            The Minimize script, for example, may be running on a single, unscaled physical topology.
        */
        boolean lambdaTerm = (numTopologies == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())
        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }
        // Relative free energies via the DualTopologyEnergy class require different
        // default free energy parameters than absolute free energies.
        if (numTopologies == 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }
        if (numTopologies == 2) {
            logger.info(format(" Initializing two topologies for each window."))
        } else {
            logger.info(format(" Initializing a single topology for each window."))
        }
        // Read in files and meta-details
        int numTopologies = files.length
        int threadsAvail = ParallelTeam.getDefaultThreadCount()
        int numParallel = topology.getNumParallel(threadsAvail, numTopologies)
        int threadsPer = (int) (threadsAvail / numParallel)
        MolecularAssembly[] molecularAssemblies = new MolecularAssembly[files.length]
        SystemFilter[] openers = new SystemFilter[files.length]
        for (int i = 0; i < numTopologies; i++) {
            MolecularAssembly ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, files[i].getName(), i)
            molecularAssemblies[i] = ma
            openers[i] = algorithmFunctions.getFilter()
        }
        StringBuilder sb = new StringBuilder(format(
                "\n Performing FEP evaluations for: %s\n ", files))
        CrystalPotential potential = (CrystalPotential) topology.assemblePotential(molecularAssemblies, threadsAvail, sb)
        String[] arcFileName = new String[files.length]
        for (int j = 0; j < numTopologies; j++) {
            // Only use the first arc file (even if restarts happened)
            arcFileName[j] = FilenameUtils.removeExtension(files[j].getAbsolutePath()) + ".arc"
            File archiveFile = new File(arcFileName[j])
            openers[j].setFile(archiveFile)
            molecularAssemblies[j].setFile(archiveFile)
        }

        // Slightly modified from BAR.groovy
        int nSnapshots = openers[0].countNumModels()
        double[] x = new double[potential.getNumberOfVariables()]
        double[] lambdaValues = new double[nLambda]
        double[][] energy = new double[nLambda][nSnapshots]
        for (int k = 0; k < lambdaValues.length; k++) {
            lambdaValues[k] = k.toDouble() / (nLambda - 1)
            energy[k] = new double[nSnapshots]
        }
        LambdaInterface linter1 = (LambdaInterface) potential
        logger.info(format("\n\n Performing energy evaluations for %d snapshots.", nSnapshots))
        logger.info(format(" Using %d lambda values.", nLambda))
        logger.info(format(" Using %d topologies.", numTopologies))
        logger.info(" Lambda values: " + lambdaValues)
        logger.info("")
        for (int i = 0; i < nSnapshots; i++) {
            // Read coords
            boolean resetPosition = (i == 0)
            for (int n = 0; n < openers.length; n++) {
                openers[n].readNext(resetPosition, false)
            }
            x = potential.getCoordinates(x)

            // Compute energies for each lambda
            StringBuilder sb2 = new StringBuilder().append("Snapshot ").append(i).append(" Energy Evaluations: ")
            for (int k = 0; k < lambdaValues.length; k++) {
                double lambda = lambdaValues[k]
                linter1.setLambda(lambda)
                energy[k][i] = potential.energy(x, false)
                sb2.append(" ").append(energy[k][i])
            }
            logger.info(sb2.append("\n").toString())
        }

        return energy
    }
}
