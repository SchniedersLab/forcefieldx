//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.Crystal
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.EstimateBootstrapper
import ffx.numerics.estimator.SequentialEstimator
import ffx.numerics.estimator.Zwanzig
import ffx.numerics.math.BootStrapStatistics
import ffx.numerics.math.SummaryStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.SystemFilter
import ffx.utilities.Constants
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils //look for ways to adjust file names
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.min

/**
 * The BAR script find the free energy difference across a lambda window. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt &lt;structures2&gt;
 */
@Command(description = " Evaluates a free energy change with the Bennett Acceptance Ratio algorithm using pregenerated snapshots.", name = "ffxc BAR")
class BAR extends AlgorithmsScript {

    @Mixin
    private AlchemicalOptions alchemical

    @Mixin
    private TopologyOptions topology

    @Option(names = ["--l2", "--lambdaTwo"], paramLabel = "1.0",
            description = "Lambda value for the upper edge of the window")
    private double lambda2 = 1.0

    @Option(names = ["--t1", "--temperature1"], paramLabel = "298.15",
            description = "Temperature for system 1")
    private double temp1 = 298.15

    @Option(names = ["--t2", "--temperature2"], paramLabel = "298.15",
            description = "Temperature for system 2")
    private double temp2 = 298.15

    @Option(names = ["--a1", "--archive1"], paramLabel = "filename.arc",
            description = "Archive for ensemble 1")
    private String arc1 = null

    @Option(names = ["--a2", "--archive2"], paramLabel = "filename.arc",
            description = "Archive for ensemble 2")
    private String arc2 = null

    @Option(names = ["--dV", "--volume"], paramLabel = "false",
            description = "Write out snapshot volumes to the Tinker BAR file.")
    private boolean includeVolume = false

    @Option(names = ["--tb", "--tinkerBAR"], paramLabel = "false",
            description = "Write out a Tinker BAR file.")
    private boolean tinkerBAR = false

    @Option(names = ["--nw", "--nWindows"], paramLabel = "-1",
            description = "If set, auto-determine lambda values and subdirectories (overrides other flags).")
    private int nWindows = -1

    /**
     * The final argument(s) should be filenames for lambda windows in order..
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
    List<String> filenames = null

    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail
    MolecularAssembly[] topologies
    SystemFilter[] openers

    CrystalPotential potential

    private Configuration additionalProperties
    private String fileNameAsString


    /**
     * Sets an optional Configuration with additional properties.
     * @param additionalProps
     */
    void setProperties(Configuration additionalProps) {
        this.additionalProperties = additionalProps
    }

    /**
     * BAR Constructor.
     */
    BAR() {
        this(new Binding())
    }

    /**
     * BAR Constructor.
     * @param binding The Groovy Binding to use.
     */
    BAR(Binding binding) {
        super(binding)
    }
    /* */

    List<String> windowFiles = new ArrayList<>()
    /**
     * {@inheritDoc}
     */
    @Override
    BAR run() {
        // Begin boilerplate code.
        if (!init()) {
            return this
        }

        double[] lambdaValues
        fileNameAsString = filenames.toString()
        fileNameAsString = fileNameAsString.replace("[", "")
        fileNameAsString = fileNameAsString.replace("]", "")

        logger.info(fileNameAsString)
        //make work for given files or nWindows
        if (nWindows != -1) {
            for (int i = 0; i < nWindows; i++) {
                String fullPathToFile = FilenameUtils.getFullPath(fileNameAsString)
                String directoryFullPath = fullPathToFile.replace(fileNameAsString,"") + i;
                windowFiles.add(directoryFullPath + "/" + i)
            }
            lambdaValues = new double[nWindows]
            for (int i = 0; i < nWindows; i++) {
                lambdaValues[i] = alchemical.getInitialLambda(nWindows, i, false);
            }
        } else {
            lambdaValues = new double[2]
            lambdaValues[0] = alchemical.getInitialLambda()
            lambdaValues[1] = lambda2
        }
        if (filenames == null) {
            return this
        }

        int nFiles = filenames.size()
        //number of files needs to equal 1 (single topology) or 2 (dual topology)
        if (nFiles != 1 && nFiles != 2) {
            return this
        }
        String directoryPath = FilenameUtils.getFullPath(fileNameAsString)
        logger.info(directoryPath)
        logger.info("the full path is" + FilenameUtils.getFullPath(fileNameAsString))
        String archiveName = FilenameUtils.getBaseName(fileNameAsString) + ".arc"
        String[] archiveFullPaths = new String[nWindows]
        for(int i=0; i<nWindows; i++){
            archiveFullPaths[i] = directoryPath + i + "/" + archiveName
            logger.info(archiveFullPaths[i])
        }


//walk through the different file names instead of listing all the files */

        for (int i = 0; i < nWindows; i++) {
            boolean directoryExist = FilenameUtils.directoryContains(directoryPath, i.toString())
            boolean fileExist = FilenameUtils.isExtension(archiveFullPaths[i], "arc")
            if (!directoryExist || !fileExist) {
                logger.info("Cannot find directories")
                return this
            }
        }


        topologies = new MolecularAssembly[nFiles]
        openers = new SystemFilter[nFiles]

        int numParallel = topology.getNumParallel(threadsAvail, nFiles)
        threadsPer = (int) (threadsAvail / numParallel)

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
    The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nFiles == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OST parameters than absolute free energies.
        if (nFiles >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        logger.info(format(" Initializing %d topologies for each end", nFiles))
        for (int i = 0; i < nFiles; i++) {
            MolecularAssembly ma =
                    alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i)

            topologies[i] = ma
            openers[i] = algorithmFunctions.getFilter()
        }
        /*
      openers of 0 setFile to the current archive
      openers.readNext(true)
      subsequent snapshots
      openers.readNext(false)
    */





           /* StringBuilder sb = new StringBuilder(format(
                    "\n Using BAR to analyze a free energy change between L=%.5f and L=%.5f for\n ", lambdaValues[k],
                    lambdaValues[k + 1]))
            potential = (CrystalPotential) topology.assemblePotential(topologies, threadsAvail, sb)
            sb.append(" and ")

            logger.info(sb.toString())

            // Check for periodic boundary conditions


            //3d array 1st column ensemble, second column window/lambda, third column energy values
            //method input (Molecular Assembly, lambdaValues[]) return 2D array with lambdaValue and energy

            String lamString1 = format("%.3f", lambdaValues[k])
            String lamString2 = format("%.3f", lambdaValues[k + 1])

            logger.info(format("\n Ensemble 1 collected at L=%s", lamString1))
            if (isPBC) {
                logger.info(
                        format(" Snapshot     E(L=%s)     E(L=%s)             dE         Volume", lamString1,
                                lamString2))
            } else {
                logger.info(format(" Snapshot     E(L=%s)     E(L=%s)             dE", lamString1, lamString2))
            }

            logger.info(format("\n Ensemble 2 collected at L=%s", lamString2))
            if (isPBC) {
                logger.info(
                        format(" Snapshot     E(L=%s)     E(L=%s)             dE         Volume", lamString1,
                                lamString2))
            } else {
                logger.info(format(" Snapshot     E(L=%s)     E(L=%s)             dE", lamString1, lamString2))
            }




            SummaryStatistics eDiff1Stats = new SummaryStatistics(eDiff1)
            logger.info(format(" Ensemble 1 Statistics (%d snapshots)", eDiff1Stats.count))
            logger.info(format("  Energy difference: %s", eDiff1Stats.describe()))
            if (isPBC) {
                SummaryStatistics vol1Stats = new SummaryStatistics(vol1)
                logger.info(format("  Volume:            %s", vol1Stats.describe()))
            }

            SummaryStatistics eDiff2Stats = new SummaryStatistics(eDiff2)
            logger.info(format(" Ensemble 2 Statistics (%d snapshots)", eDiff2Stats.count))
            logger.info(format("  Energy difference: %s", eDiff2Stats.describe()))
            if (isPBC) {
                SummaryStatistics vol2Stats = new SummaryStatistics(vol2)
                logger.info(format("  Volume:            %s", vol2Stats.describe()))
            }

            if (tinkerBAR) {
                String barFileName = FilenameUtils.removeExtension(filenames.get(0)) + ".bar"
                logger.info(format("\n Writing Tinker-compatible BAR file to %s.", barFileName))

                new File(barFileName).withWriter { bw ->
                    StringBuilder fnames1 = new StringBuilder()
                    for (int i = 0; i < nFiles; i++) {
                        fnames1.append("  ").append(filenames.get(i))
                    }
                    bw.write(format("%8d %9.3f%s\n", nSnapshots1, temp1, fnames1.toString()))
                    for (int i = 0; i < nSnapshots1; i++) {
                        if (isPBC) {
                            bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e1L1[i], e1L2[i], vol1[i]))
                        } else {
                            bw.write(format("%8d %20.10f %20.10f\n", i + 1, e1L1[i], e1L2[i]))
                        }
                    }

                    StringBuilder fnames2 = new StringBuilder()
                    for (int i = nFiles; i < nFiles; i++) {
                        fnames2.append("  ").append(filenames.get(i))
                    }
                    bw.write(format("%8d %9.3f  %s\n", nSnapshots2, temp2, fnames2.toString()))
                    for (int i = 0; i < nSnapshots2; i++) {
                        if (isPBC) {
                            bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e2L1[i], e2L2[i], vol2[i]))
                        } else {
                            bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e2L1[i], e2L2[i]))
                        }
                    }
                }
            }*/

            /*// Load energy values for FEP and BAR
            double[][] energyLow = new double[2][]
            double[][] energyAt = new double[2][]
            double[][] energyHigh = new double[2][]
            // First Lambda Value.
            energyLow[0] = new double[nSnapshots1]
            energyAt[0] = e1L1
            energyHigh[0] = e1L2
            // Second Lambda Value.
            energyLow[1] = e2L1
            energyAt[1] = e2L2
            energyHigh[1] = new double[nSnapshots2]

            if(k== 0 || k== nWindows-1){

            } else {
                currentLambdas = new double[3]
            }*/
        StringBuilder sb = new StringBuilder(format(
                "\n Using BAR to analyze a free energy change for %s\n ", filenames))
        potential = (CrystalPotential) topology.assemblePotential(topologies, threadsAvail, sb)
        Crystal unitCell1 = potential.getCrystal().getUnitCell()
        Crystal unitCell2 = potential.getCrystal().getUnitCell()
        boolean isPBC = includeVolume && !unitCell1.aperiodic()
        isPBC = isPBC && !unitCell2.aperiodic()
        isPBC = isPBC && (unitCell1.getNumSymOps() == unitCell2.getNumSymOps())
        int nSymm = 0
        if (isPBC) {
            nSymm = unitCell1.getNumSymOps()
        }
        double[] currentLambdas
        double[][] energyLow = new double[nWindows][]
        double[][] energyAt = new double [nWindows][]
        double[][] energyHigh = new double [nWindows][]
        double[] volume
        double[][] energy
        double[] energyMean = new double[nWindows]
        double[] energySD = new double[nWindows]
        double[] energyVar = new double[nWindows]
            for(int w=0; w<nWindows; w++){

                if (w==0){
                    currentLambdas = new double[2]
                    currentLambdas[0] = lambdaValues[w]
                    currentLambdas[1] = lambdaValues[w+1]
                } else if(w == nWindows -1){
                    currentLambdas = new double[2]
                    currentLambdas[0] = lambdaValues[w-1]
                    currentLambdas[1] = lambdaValues[w]
                } else{
                    currentLambdas = new double[3]
                    currentLambdas[0] = lambdaValues[w-1]
                    currentLambdas[1] = lambdaValues[w]
                    currentLambdas[2] = lambdaValues[w+1]
                }
                energy = new double[currentLambdas.length][]
                //opener pointing to the right archive
                //cant cast string to file look at documentation and adjust this line
                File archiveFile = new File(archiveFullPaths[w])
                openers[0].setFile(archiveFile)
                volume = getEnergyForLambdas(topologies, currentLambdas,
                        archiveFullPaths[w], energy, isPBC, nSymm)
                logger.info("volume works")
                if(w == 0){
                        energyLow[w] = null
                        energyAt[w] = energy[0]
                        energyHigh[w] = energy[1]

                } else if(w == nWindows-1){
                        energyLow[w] = energy[0]
                        energyAt[w] = energy[1]
                        energyHigh[w] = null
                } else {
                    energyLow[w] = energy[0]
                    energyAt[w] = energy[1]
                    energyHigh[w] = energy[2]
                }
                logger.info("energy sent")
                BootStrapStatistics energyStats = new BootStrapStatistics(energy[w])
                energyMean[w] = energyStats.mean
                energySD[w] = energyStats.sd
                energyVar[w] = energyStats.var

            }


        logger.info("Successfully parsed energy")
            double[] temperature = new double[nWindows]
            Arrays.fill(temperature, temp1)
            //reference the 2D array for energies
            SequentialEstimator bar = new BennettAcceptanceRatio(lambdaValues, energyLow, energyAt, energyHigh,
                    temperature)
            SequentialEstimator forwards = bar.getInitialForwardsGuess()
            SequentialEstimator backwards = bar.getInitialBackwardsGuess()

            EstimateBootstrapper barBS = new EstimateBootstrapper(bar)
            EstimateBootstrapper forBS = new EstimateBootstrapper(forwards)
            EstimateBootstrapper backBS = new EstimateBootstrapper(backwards)

            long MAX_BOOTSTRAP_TRIALS = 100000L
            long bootstrap = min(MAX_BOOTSTRAP_TRIALS, min(volume.length, volume.length))

            //make logging information look nicer (evenly spaced, new lines)
        //free energy estimates for BAR and bootstrapping do every window
            logger.info("\n Free Energy Difference via FEP Method\n")
            long time = -System.nanoTime()
            forBS.bootstrap(bootstrap)
            time += System.nanoTime()
            logger.fine(format(" Forward FEP Bootstrap Complete:      %7.4f sec", time * Constants.NS2SEC))
            double sumForeFE = forBS.getTotalFE()
            double sumEnthalpyFore = forBS.getTotalEnthalpy()
            double varForeFE = forBS.getTotalUncertainty()
            double varEnthalpyFore = forBS.getTotalEnthalpyUncertainty()
            logger.info(format(" Free energy via Forwards FEP:   %12.4f +/- %6.4f kcal/mol.", sumForeFE, varForeFE))


            time = -System.nanoTime()
            backBS.bootstrap(bootstrap)
            time += System.nanoTime()
            logger.fine(format(" Backward FEP Bootstrap Complete:     %7.4f sec", time * Constants.NS2SEC))

            double sumBackFE = backBS.getTotalFE()
            double sumEnthalpyBack = backBS.getTotalEnthalpy()
            double varBackFE = backBS.getTotalUncertainty()
            double varEnthalpyBack = backBS.getTotalEnthalpyUncertainty()
            logger.info(format(" Free energy via Backwards FEP:  %12.4f +/- %6.4f kcal/mol.", sumBackFE, varBackFE))


            logger.info("\n Free Energy Difference via BAR Method\n")
            logger.info(
                    format(" Free energy via BAR Iteration:  %12.4f +/- %6.4f kcal/mol.", bar.getFreeEnergy(),
                            bar.getUncertainty()))
            time = -System.nanoTime()
            barBS.bootstrap(bootstrap)
            time += System.nanoTime()
            logger.fine(format(" BAR Bootstrap Complete:              %7.4f sec", time * Constants.NS2SEC))

            double sumBARFE = barBS.getTotalFE()
            double varBARFE = barBS.getTotalUncertainty()
            double barEnthalpy = barBS.getTotalEnthalpy()
            double varEnthalpy = barBS.getTotalEnthalpyUncertainty()
            logger.info(format(" Free energy via BAR Bootstrap:  %12.4f +/- %6.4f kcal/mol.", sumBARFE, varBARFE))

            //not set up generally yet
        //adjust to work for windows
        //add statistics to loop then pass into this section
            logger.info(" Enthalpy from Potential Energy Averages")

        for(int t=0; t<nWindows; t++){
            logger.info(format(" Average Energy for State %d:          %12.4f +/- %6.4f kcal/mol.",
                   t, energyMean[t], energySD[t]))
            if(t+1>nWindows){
                logger.info("end")
            } else {
                double enthalpyDiff = energyMean[t+1] - energyMean[t]
                double enthalpyDiffSD = Math.sqrt(energyVar[nWindows-1] + energyVar[0])
                logger.info(format(" Enthalpy via Direct Estimate between Energy State %d and %d:  " +
                        " %12.4f +/- %6.4f kcal/mol.", t+1, t, enthalpyDiff, enthalpyDiffSD))
            }

        }
        double enthalpyDiff = energyMean[nWindows-1] - energyMean[0]
        double enthalpyDiffSD = Math.sqrt(energyVar[nWindows-1] + energyVar[0])
        logger.info(format(" Enthalpy via Direct Estimate:   %12.4f +/- %6.4f kcal/mol.",
                enthalpyDiff, enthalpyDiffSD))

            logger.info("Enthalpy and Entropy via FEP")

            double forwardsEntropy = (sumEnthalpyFore - sumForeFE) / temp1
            double backwardsEntropy = (sumEnthalpyBack - sumBackFE) / temp1

            logger.info(format(" Enthalpy via Forward FEP       %12.4f +/- %6.4f kcal/mol.", sumEnthalpyFore, varEnthalpyFore))
            logger.info(format(" Entropy via Forward FEP        %12.4f kcal/mol/K.", forwardsEntropy))
            logger.info(format(" Forward FEP -T*ds Value         %12.4f kcal/mol.", -(forwardsEntropy * temp1)))

            logger.info(format(" Enthalpy via Backward FEP     %12.4f +/- %6.4f kcal/mol.", sumEnthalpyBack, varEnthalpyBack))
            logger.info(format(" Entropy via Backward FEP     %12.4f kcal/mol/K.", backwardsEntropy))
            logger.info(format(" Backward FEP -T*ds Value      %12.4f kcal/mol.", -(backwardsEntropy * temp1)))

            double tsBar = barEnthalpy - sumBARFE
            double barEntropy = tsBar / (temp1)
            logger.info(format(" Enthalpy via BAR     %12.4f +/- %6.4f kcal/mol.", barEnthalpy, varEnthalpy))
            logger.info(format(" Entropy via BAR     %12.4f kcal/mol/K.", barEntropy))
            logger.info(format(" BAR Estimate of -T*ds      %12.4f kcal/mol.", -(tsBar)))



        return this
    }

    //logic to set up each lambda array
    //send energy array to a 3D array with windows
    //look at logging for statistics of energyAt, low, high
    private double[] getEnergyForLambdas(MolecularAssembly[] topologies, double[] lambdaValues,
                                          String arcFileName, double[][] energy, boolean isPBC, int nSymm) {
        int nSnapshots = openers[0].countNumModels()
        logger.info("Number of models:" + nSnapshots.toString())
        double[] x = new double[potential.getNumberOfVariables()]
        double[] vol = new double[nSnapshots]
        energy = new double[nSnapshots]
        LambdaInterface linter1 = (LambdaInterface) topologies[0].getPotentialEnergy()
        StringBuilder sb = new StringBuilder(format(
                "\n Evaluating energies for %s\n ",  arcFileName))

        logger.info(sb as String)

        for(int i=0; i<nSnapshots; i++){

            for (int j = 0; j < 1; j++) {
                openers[j].readNext(false, false)
            }
            for (int k = 0; k < lambdaValues.length; k++) {

                logger.info("energy snapshots")
                double lambda = lambdaValues[k]
                x = potential.getCoordinates(x)
                linter1.setLambda(lambda)
                energy[k][i] = potential.energy(x, false)



            if (isPBC) {
                Crystal unitCell = potential.getCrystal().getUnitCell()
                vol[i] = unitCell.volume / nSymm
                logger.info(format(" %8d %14.4f %14.4f %14.4",
                        i + 1,lambdaValues[k], energy[k][i], vol[i]))
            } else {
                logger.info(format(" %8d %14.4f %14.4f",
                        i + 1,lambdaValues[k], energy[k][i]))
            }
        }


        }

        return vol

    }
    /**
     * {@inheritDoc}
     */
    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (potential == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = new ArrayList<>()
            if (potential != null) {
                potentials.add(potential)
            }
        }
        return potentials
    }
}










