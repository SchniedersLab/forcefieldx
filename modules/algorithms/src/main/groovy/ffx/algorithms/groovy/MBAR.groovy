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

    @Option(names = ["--lambdaDerivative", "--lD"], paramLabel = "false",
            description = "Calculate lambda derivatives for each snapshot.")
    boolean lambdaDerivative = false

    @Option(names = ["--continuousLambda", "--cL"], paramLabel = "false",
            description = "Data comes from continuous lambda source and only contains mbar file.")
    boolean continuousLambda = false

    @Option(names = ["--outputDir", "--oD"], paramLabel = "",
            description = "Where to place MBAR files. Default is ../mbarFiles/energy_(window#).mbar. Will write out a file called energy_0.mbar.")
    String outputDirectory = ""

    @Option(names = ["--seed"], paramLabel = "BAR",
            description = "Seed MBAR calculation with this: ZEROS, ZWANZIG, BAR. Fallback to ZEROS if input is does not or is unlikely to converge.")
    String seedWith = "BAR"

    @Option(names = ["--tol", "--tolerance"], paramLabel = "1e-7",
            description = "Iteration change tolerance.")
    double tol = 1e-7

    @Option(names = ["--ss", "--startSnapshot"], paramLabel = "-1",
            description = "Start at this snapshot when reading in tinker BAR files.")
    int startingSnapshot = -1

    @Option(names = ["--es", "--endSnapshot"], paramLabel = "-1",
            description = "End at this snapshot when reading in tinker BAR files.")
    private int endingSnapshot = -1

    @Option(names = ["--verbose"], paramLabel = "false",
            description = "Log weight matrices, iterations, and other details.")
    boolean verbose = false


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

        // Write MBAR file if option is set
        if(isArc){
            if (numLambda == -1){
                logger.severe("numLambda must be specified for lambda energy evaluations.")
                return this
            }
            // Get list of fileNames & check validity
            File parent = files[0].getParentFile() // Run directory
            int window;
            File outputDir;
            if (outputDirectory.isEmpty()) {
                window = Integer.parseInt(parent.getName()) // Run name should be int
                outputDir = new File(parent.getParentFile(), "mbarFiles") // Make mbarFiles
                if(!outputDir.exists()) {
                    outputDir.mkdir()
                }
            } else {
                outputDir = new File(outputDirectory)
                if(!outputDir.exists()) {
                    outputDir.mkdir()
                }
                window = 0
            }
            // Write MBAR file with window number, although this will be reassigned by the file filter based on
            // placement relative to other fileNames with energy values.
            File outputFile = new File(outputDir, "energy_" + window + ".mbar")
            //TODO: Fix atrocious setting of temperatures
            double[][][] energiesAndDerivatives = getEnergyForLambdas(files, numLambda)
            double[][] energies =  energiesAndDerivatives[0] // Long step!
            MultistateBennettAcceptanceRatio.writeFile(energies, outputFile, 298) // Assume 298 K
            if (lambdaDerivative) {
                double[][] lambdaDerivatives = energiesAndDerivatives[1]
                File outputDerivFile = new File(outputDir, "derivatives_" + window + ".mbar")
                MultistateBennettAcceptanceRatio.writeFile(lambdaDerivatives, outputDerivFile, 298)
            }
            return this
        }

        // Run MBAR calculation if file write-out is not requested & files are correct
        File path = new File(fileList.get(0))
        if (!path.exists()) {
            logger.severe("Path to MBAR/BAR fileNames does not exist: " + path)
            return this
        }
        if (!path.isDirectory() && !(path.isFile() && path.canRead())) {
            logger.severe("Path to MBAR/BAR fileNames is not accessible: " + path)
            return this
        }
        MBARFilter filter = new MBARFilter(path, continuousLambda);
        if (startingSnapshot >= 0){
            filter.setStartSnapshot(startingSnapshot)
            logger.info("Starting with snapshot index: " + startingSnapshot)
        }
        if (endingSnapshot >= 0){
            filter.setEndSnapshot(endingSnapshot)
            logger.info("Ending with snapshot index: " + endingSnapshot)
        }
        seedWith = seedWith.toUpperCase()
        SeedType seed = SeedType.valueOf(seedWith) as SeedType
        if (seed == null) {
            logger.severe("Invalid seed type: " + seedWith)
            return this
        }
        MultistateBennettAcceptanceRatio.VERBOSE = verbose

        // Runs calculation on class creation
        mbar = filter.getMBAR(seed as MultistateBennettAcceptanceRatio.SeedType, tol)
        this.mbar = mbar
        if (mbar == null) {
            logger.severe("Could not create MBAR object.")
            return this
        }

        // Print out results
        logger.info("\n MBAR Results:")
        double[] dGs = mbar.getBinEnergies()
        double[] uncertainties = mbar.getBinUncertainties()
        logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol\n", mbar.getFreeEnergy(),
                mbar.getUncertainty()))
        for (int i = 0; i < dGs.length; i++) {
            logger.info(format("   dG %3d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
        }

        logger.info("\n MBAR Enthalpy & Entropy Results:")
        double[] enthalpies = mbar.getBinEnthalpies()
        double[] entropies = mbar.getBinEntropies()
        double totalEnthalpy = sum(enthalpies)
        double totalEntropy = sum(entropies)
        logger.info(format(" Total dG = %10.4f (dH) - %10.4f (TdS) kcal/mol\n", totalEnthalpy, totalEntropy))
        for (int i = 0; i < enthalpies.length; i++) {
            logger.info(format("   dG %3d = %10.4f (dH) - %10.4f (TdS) kcal/mol", i, enthalpies[i], entropies[i]))
        }

        logger.info("\n MBAR uncertainty between all i & j: ")
        double[][] uncertaintyMatrix = mbar.getUncertaintyMatrix()
        for(int i = 0; i < uncertaintyMatrix.length; i++) {
            StringBuilder sb = new StringBuilder()
            sb.append("    [")
            for(int j = 0; j < uncertaintyMatrix[i].length; j++) {
                sb.append(format(" %6.5f ", uncertaintyMatrix[i][j]))
            }
            sb.append("]")
            logger.info(sb.toString())
        }

        // Read in and compute expectations for observable data
        boolean observableData = filter.readObservableData(true, false, true);
        if (observableData){
            logger.info("\n Observable data read in.")
        }
        boolean biasData = filter.readObservableData(true, true, false);
        if(biasData){
            logger.info(" Bias data read in.")
        }
        if(observableData) {
            logger.info("\n MBAR Observable Data: ")
            double[] observableValues = mbar.getObservationEnsembleAverages()
            for (int i = 0; i < observableValues.length; i++) {
                logger.info(format("     %3d = %10.4f ", i, observableValues[i]))
            }
            logger.info(" Integral:    " + mbar.getTIIntegral())
        }
        logger.info("\n")
        // BAR to compare (negligible difference if properly converged and doesn't take long at all)
        if(bar){
            try {
                logger.info("\n BAR Results:")
                BennettAcceptanceRatio bar = mbar.getBAR()
                logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol\n", bar.getFreeEnergy(),
                        bar.getUncertainty()))
                dGs = bar.getBinEnergies()
                uncertainties = bar.getBinUncertainties()
                for (int i = 0; i < dGs.length; i++) {
                    logger.info(format("   dG %3d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                }
                enthalpies = bar.getBinEnthalpies()
                totalEnthalpy = sum(enthalpies)
                logger.info(format("\n Total dH = %10.4f kcal/mol\n", totalEnthalpy))
                for (int i = 0; i < enthalpies.length; i++) {
                    logger.info(format("   dH %3d = %10.4f kcal/mol", i, enthalpies[i]))
                }
            } catch (Exception ignored) {
                logger.warning(" BAR calculation failed to converge.")
            }
        }
        logger.info("\n")
        // Bootstrapping MBAR & BAR
        if(numBootstrap != 0) {
            EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar)
            bootstrapper.bootstrap(numBootstrap)
            logger.info("\n MBAR Bootstrap Results from " + numBootstrap + " Samples:")
            logger.info(format(" Total dG = %10.4f +/- %10.4f kcal/mol", bootstrapper.getTotalFE(),
                    bootstrapper.getTotalUncertainty()))
            dGs = bootstrapper.getFE()
            uncertainties = bootstrapper.getUncertainty()
            for (int i = 0; i < dGs.length; i++) {
                logger.info(format("    dG %3d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
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
                        logger.info(format("    dG %3d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                    }
                } catch (Exception ignored) {
                    logger.warning(" BAR calculation failed to converge.")
                }
            }
        }
        // Compare MBAR calculation across different tenths of the data
        if (convergence){
            MultistateBennettAcceptanceRatio.FORCE_ZEROS_SEED = true;
            MultistateBennettAcceptanceRatio[] mbarPeriodComparison =
                    filter.getPeriodComparisonMBAR(seed as MultistateBennettAcceptanceRatio.SeedType, 1e-7)
            double[][] dGPeriod = new double[mbarPeriodComparison.length][]
            for (int i = 0; i < mbarPeriodComparison.length; i++) {
                dGPeriod[i] = mbarPeriodComparison[i].getBinEnergies()
            }
            logger.info("\n MBAR Period Comparison Results:")
            logger.info(format("     %10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%%%10d%% ",
                    10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
            double[] totals = new double[dGPeriod[0].length]
            for(int i = 0; i < dGPeriod[0].length; i++) {
                StringBuilder sb = new StringBuilder()
                sb.append(" dG_").append(i).append(": ")
                for(int j = 0; j < dGPeriod.length; j++) {
                    sb.append(format("%10.4f ", dGPeriod[j][i]))
                    totals[j] += dGPeriod[j][i]
                }
                logger.info(sb.toString())
            }
            StringBuilder totalsSB = new StringBuilder()
            for(int i = 0; i < totals.length; i++) {
                totalsSB.append(format("%10.4f ", totals[i]))
            }
            logger.info("")
            logger.info("  Tot: " + totalsSB.toString())
        }
        return this
    }

    private double sum(double[] values) {
        double sum = 0
        for (double value : values) {
            sum += value
        }
        return sum
    }

    private double[][][] getEnergyForLambdas(File[] files, int nLambda) {
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
        double[][] lambdaDerivatives = new double[nLambda][nSnapshots]
        for (int k = 0; k < lambdaValues.length; k++) {
            lambdaValues[k] = k.toDouble() / (nLambda - 1)
            energy[k] = new double[nSnapshots]
            lambdaDerivatives[k] = new double[nSnapshots]
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
            StringBuilder sb3 = new StringBuilder().append("Snapshot ").append(i).append(" Lambda Derivatives: ")
            for (int k = 0; k < lambdaValues.length; k++) {
                double lambda = lambdaValues[k]
                if (lambda <= 1E-6){
                    lambda += .00275
                }
                if (lambda - 1.0 < 1E-6){
                    lambda -= .00275
                }
                linter1.setLambda(lambda)
                energy[k][i] = potential.energyAndGradient(x, new double[x.length * 3])
                if (lambdaDerivative) {
                    lambdaDerivatives[k][i] = linter1.getdEdL()
                    sb3.append(" ").append(lambdaDerivatives[k][i])
                }
                sb2.append(" ").append(energy[k][i])
            }
            logger.info(sb2.append("\n").toString())
            if (lambdaDerivative) {
                logger.info(sb3.append("\n").toString())
            }
        }

        return new double[][][]{energy, lambdaDerivatives}
    }
}
