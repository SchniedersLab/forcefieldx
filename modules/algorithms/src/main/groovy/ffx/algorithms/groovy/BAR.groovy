//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import ffx.numerics.math.BootStrapStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.BARFilter
import ffx.potential.parsers.SystemFilter
import ffx.utilities.Constants
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
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
@Command(description = " Evaluates a free energy change with the Bennett Acceptance Ratio algorithm using pregenerated snapshots.", name = "BAR")
class BAR extends AlgorithmsScript {

  @Mixin
  private AlchemicalOptions alchemical

  @Mixin
  private TopologyOptions topology

  @Option(names = ["--l2", "--lambdaTwo"], paramLabel = "1.0",
      description = "Lambda value for the upper edge of the window")
  private double lambda2 = 1.0

  @Option(names = ["-t", "--temperature"], paramLabel = "298.15",
      description = "Temperature for system")
  private double temperature = 298.15

  @Option(names = ["--dV", "--volume"], paramLabel = "false",
      description = "Write out snapshot volumes to the Tinker BAR file.")
  private boolean includeVolume = false

  @Option(names = ["--tb", "--tinkerBAR"], paramLabel = "false",
      description = "Write out a Tinker BAR file.")
  private boolean tinkerBAR = false

  @Option(names = ["--nw", "--nWindows"], paramLabel = "-1",
      description = "If set, auto-determine lambda values and subdirectories (overrides other flags).")
  private int nWindows = -1

  @Option(names = ["--utb", "--useTinker"], paramLabel = "false",
      description = "If set, use tinker BAR files for energy snapshots.")
  private boolean useTinkerBAR = false

  @Option(names = ["--sa", "--sortedArc"], paramLabel = "false",
      description = "If set, use sorted archive values.")
  private boolean sortedArc = false

  @Option(names = ["--ss", "--startSnapshot"], paramLabel = "0",
      description = "Start at this snapshot when reading in tinker BAR files.")
  private int startingSnapshot = 0

  @Option(names = ["--es", "--endSnapshot"], paramLabel = "0",
      description = "End at this snapshot when reading in tinker BAR files.")
  private int endingSnapshot = 0

  /**
   * The final argument(s) should be filenames for lambda windows in order.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = 'A single PDB/XYZ when windows are auto-determined (or two for dual topology). Two trajectory files for BAR between two ensembles (or four for dual topology).')
  List<String> filenames = null

  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  private int threadsPer = threadsAvail
  MolecularAssembly[] molecularAssemblies
  CrystalPotential potential
  SystemFilter[] openers
  BARFilter[] barOpeners
  BARFilter[] barWriters
  List<String> windowDirectories = new ArrayList<>()
  private Configuration additionalProperties
  private String[] files

  /**
   * Free energy difference from forward FEP.
   */
  private double forwardFEP
  /**
   * Enthalpy difference from forward FEP.
   */
  private double forwardEnthalpy
  /**
   * Entropy difference from forward FEP.
   */
  private double forwardEntropy
  /**
   * Free energy difference from backward FEP.
   */
  private double backwardFEP
  /**
   * Enthalpy difference from backward FEP.
   */
  private double backwardEnthalpy
  /**
   * Entropy difference from backward FEP.
   */
  private double backwardEntropy
  /**
   * Free energy difference from BAR.
   */
  private double barEnergy
  /**
   * Free energy difference from BAR using boot-strapping.
   */
  private double barEnergyBS
  /**
   * Enthalpy difference from BAR.
   */
  private double barEnthalpy
  /**
   * Entropy difference from BAR.
   */
  private double barEntropy
  /**
   * Number of Topologies
   */
  private int numTopologies

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

  /**
   * {@inheritDoc}
   */
  @Override
  BAR run() {
    // Begin boilerplate code.
    if (!init()) {
      return this
    }

    /**
     * Load user supplied files into an array.
     */
    int nFiles = filenames.size()
    files = new String[nFiles]
    for (int i = 0; i < nFiles; i++) {
      files[i] = filenames.get(i)
    }

    numTopologies = 1
    if (nFiles <= 1 && nWindows < 2) {
      logger.info(' At least two ensembles must be specified')
      return this
    } else if (nFiles == 1 && nWindows >= 2) {
      logger.
          info(format(' Auto-detecting %d windows for single topology:\n %s.',
              nWindows, files[0]))
    } else if (nFiles == 2 && nWindows >= 2) {
      logger.info(format(' Auto-detecting %d windows for dual topology:\n %s\n %s.',
          nWindows, files[0], files[1]))
      numTopologies = 2
    } else if (nFiles == 2 && nWindows < 2) {
      logger.info(format(' Applying BAR between two single topology ensembles:\n %s\n %s.',
          files[0], files[1]))
    } else if (nFiles == 4 && nWindows < 2) {
      logger.info(format(' Applying BAR between two dual topology ensembles:\n %s %s\n %s %s.',
          files[0], files[1], files[2], files[3]))
      numTopologies = 2
    } else {
      logger.info(' Inconsistent input of files and/or windows.')
      return this
    }

    boolean autodetect = false
    double[] lambdaValues
    if (nWindows > 1) {
      autodetect = true
      // Auto-determine subdirectories and their lambda values.
      for (int i = 0; i < nWindows; i++) {
        for (int j = 0; j < nFiles; j++) {
          String fullPathToFile = FilenameUtils.getFullPath(files[j])
          String directoryFullPath = fullPathToFile.replace(files[j], "") + i;
          windowDirectories.add(directoryFullPath + File.separator + i)
        }
      }
      lambdaValues = new double[nWindows]
      for (int i = 0; i < nWindows; i++) {
        lambdaValues[i] = alchemical.getInitialLambda(nWindows, i, false);
      }
    } else {
      // Otherwise we assume two ensembles at then given lambda values.
      lambdaValues = new double[2]
      lambdaValues[0] = alchemical.getInitialLambda()
      lambdaValues[1] = lambda2
      nWindows = 2
    }

    logger.info(" Lambda values for each window: ")
    for (int i = 0; i < lambdaValues.length; i++) {
      double l = lambdaValues[i]
      logger.info(format(" %d: %6.4f", i, l))
    }

    // If reading Tinker BAR file(s), allocate storage for file paths and openers.
    String[][] fullFilePaths
    if (useTinkerBAR) {
      fullFilePaths = new String[nWindows][1]
      barOpeners = new BARFilter[nWindows]
    } else {
      fullFilePaths = new String[nWindows][nFiles]
    }

    // If saving Tinker BAR file(s), allocate storage for writer(s).
    if (tinkerBAR) {
      barWriters = new BARFilter[nWindows]
    }

    File file = new File(files[0])
    String directoryPath = file.getAbsoluteFile().getParent() + File.separator

    // Loop over user supplied files.
    for (int j = 0; j < nFiles; j++) {
      // Loop over ensembles.
      for (int i = 0; i < nWindows; i++) {
        String archiveName
        if (sortedArc) {
          archiveName = FilenameUtils.getBaseName(files[j]) + "_E" + i.toString() + ".arc"
        } else {
          archiveName = FilenameUtils.getBaseName(files[j]) + ".arc"
        }
        if (useTinkerBAR && i != nWindows-1) {
          // Path to Tinker BAR files.
          fullFilePaths[i][0] = directoryPath + "barFiles" + File.separator + "energy_" + i + ".bar"
        } else if (!autodetect) {
          // Path to a file in the same directory as supplied archives.
          fullFilePaths[i][j] = directoryPath + File.separator + archiveName
        } else {
          // Paths to auto-detected subdirectories.
          fullFilePaths[i][j] = directoryPath + i + File.separator + archiveName
        }
      }
    }

    // Allocate space for each topology
    molecularAssemblies = new MolecularAssembly[numTopologies]
    openers = new SystemFilter[numTopologies]

    int numParallel = topology.getNumParallel(threadsAvail, numTopologies)
    threadsPer = (int) (threadsAvail / numParallel)

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

    for (int i = 0; i < numTopologies; i++) {
      MolecularAssembly ma =
          alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i)
      molecularAssemblies[i] = ma
      openers[i] = algorithmFunctions.getFilter()
    }

    StringBuilder sb = new StringBuilder(format(
        "\n Using BAR to analyze a free energy change for %s\n ", filenames))
    potential = (CrystalPotential) topology.assemblePotential(molecularAssemblies, threadsAvail, sb)
    Crystal unitCell = potential.getCrystal().getUnitCell()

    boolean isPBC = includeVolume && !unitCell.aperiodic()
    isPBC = isPBC && !unitCell.aperiodic()

    int nSymm = 0
    if (isPBC) {
      nSymm = unitCell.getNumSymOps()
    }

    double[] currentLambdas
    double[][] energyLow = new double[nWindows][]
    double[][] energyAt = new double[nWindows][]
    double[][] energyHigh = new double[nWindows][]
    double[][] volume = new double[nWindows][]
    double[][] energy
    double[] energyMean = new double[nWindows]
    double[] energySD = new double[nWindows]
    double[] energyVar = new double[nWindows]

    for (int w = 0; w < nWindows; w++) {

      if (w == 0) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w]
        currentLambdas[1] = lambdaValues[w + 1]
      } else if (w == nWindows - 1) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
      } else {
        currentLambdas = new double[3]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
        currentLambdas[2] = lambdaValues[w + 1]
      }
      energy = new double[currentLambdas.length][]

      if (useTinkerBAR) {
        File barFile = new File(fullFilePaths[w][0])
        barOpeners[w] = new BARFilter(barFile, startingSnapshot, endingSnapshot)
        if (w != nWindows-1){
          barOpeners[w].readFile()
        }
        if (w == 0) {
          energyLow[w] = new double[barOpeners[w].getSnaps()]
          energyAt[w] = barOpeners[w].getE1l1()
          energyHigh[w] = barOpeners[w].getE1l2()
        } else if (w == nWindows - 1) {
          energyLow[w] = barOpeners[w-1].getE2l1()
          energyAt[w] = barOpeners[w-1].getE2l2()
          energyHigh[w] = new double[barOpeners[w].getSnaps()]
        } else if (w > 0 && w < nWindows - 1) {
          energyLow[w] = barOpeners[w - 1].getE2l1()
          energyAt[w] = barOpeners[w].getE1l1()
          energyHigh[w] = barOpeners[w].getE1l2()
        }

        if (isPBC) {
          volume[w] = barOpeners[w].getVolume1()
        }

      } else {
        volume[w] = getEnergyForLambdas(molecularAssemblies, currentLambdas,
            fullFilePaths[w], energy, isPBC, nSymm)

        if (w == 0) {
          energyLow[w] = new double[energy[0].length]
          energyAt[w] = energy[0]
          energyHigh[w] = energy[1]
        }else if (w == nWindows - 1) {
          energyLow[w] = energy[0]
          energyAt[w] = energy[1]
          energyHigh[w] = new double[energy[0].length]
        } else if (w > 0 && w < nWindows - 1) {
          energyLow[w] = energy[0]
          energyAt[w] = energy[1]
          energyHigh[w] = energy[2]
        }
      }
      BootStrapStatistics energyStats = new BootStrapStatistics(energyAt[w])
      energyMean[w] = energyStats.mean
      energySD[w] = energyStats.sd
      energyVar[w] = energyStats.var
    }

    String tinkerFilePath = ""
    if (tinkerBAR) {
      String tinkerDirectoryPath = directoryPath + File.separator + "barFiles"
      File directory = new File(tinkerDirectoryPath)
      tinkerFilePath = tinkerDirectoryPath + File.separator
      directory.mkdir()
    }

    File xyzFile = new File(filenames.get(0))
    for (int w = 0; w < nWindows; w++) {
      if (tinkerBAR) {
        if (w == 0) {
          barWriters[w] = new BARFilter(xyzFile, energyAt[w], energyHigh[w], energyLow[w + 1],
              energyAt[w + 1], volume[w], volume[w + 1],
              this.temperature)
        } else if (w != nWindows - 1) {
          barWriters[w] = new BARFilter(xyzFile, energyAt[w], energyHigh[w], energyLow[w + 1],
              energyAt[w + 1], volume[w], volume[w + 1],
              this.temperature)
        }
        if (w != nWindows-1){
          String barFileName = tinkerFilePath + "energy_" + w.toString() + ".bar"
          barWriters[w].writeFile(barFileName, isPBC)
        }
      }
    }

    for (int w = 0; w < nWindows + 1; w++) {
      if (w == nWindows) {
        logger.info("\n\n Evaluating Overall:")
      } else {
        logger.info(format("\n\n Evaluating Window %d:", w))
      }

      if (w == 0) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w]
        currentLambdas[1] = lambdaValues[w + 1]
      } else if (w == nWindows - 1) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
      } else if (w == nWindows) {
        currentLambdas = lambdaValues
      } else {
        currentLambdas = new double[3]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
        currentLambdas[2] = lambdaValues[w + 1]
      }
      double[][] energyWindowLow
      double[][] energyWindowAt
      double[][] energyWindowHigh

      energyWindowLow = new double[currentLambdas.length][]
      energyWindowAt = new double[currentLambdas.length][]
      energyWindowHigh = new double[currentLambdas.length][]

      if (w == 0) {
        energyWindowLow[0] = energyLow[w]
        energyWindowLow[1] = energyLow[w + 1]

        energyWindowAt[0] = energyAt[w]
        energyWindowAt[1] = energyAt[w + 1]

        energyWindowHigh[0] = energyHigh[w]
        energyWindowHigh[1] = energyHigh[w + 1]
      } else if (w == nWindows - 1) {
        energyWindowLow[0] = energyLow[w - 1]
        energyWindowLow[1] = energyLow[w]

        energyWindowAt[0] = energyAt[w - 1]
        energyWindowAt[1] = energyAt[w]

        energyWindowHigh[0] = energyHigh[w - 1]
        energyWindowHigh[1] = energyHigh[w]
      } else if (w == nWindows) {
        energyWindowLow = energyLow
        energyWindowAt = energyAt
        energyWindowHigh = energyHigh
      } else {
        energyWindowLow[0] = energyLow[w - 1]
        energyWindowLow[1] = energyLow[w]
        energyWindowLow[2] = energyLow[w + 1]

        energyWindowAt[0] = energyAt[w - 1]
        energyWindowAt[1] = energyAt[w]
        energyWindowAt[2] = energyAt[w + 1]

        energyWindowHigh[0] = energyHigh[w - 1]
        energyWindowHigh[1] = energyHigh[w]
        energyWindowHigh[2] = energyHigh[w + 1]
      }

      SequentialEstimator bar = new BennettAcceptanceRatio(currentLambdas, energyWindowLow,
          energyWindowAt, energyWindowHigh,
          temperature)
      SequentialEstimator forwards = bar.getInitialForwardsGuess()
      SequentialEstimator backwards = bar.getInitialBackwardsGuess()

      EstimateBootstrapper barBS = new EstimateBootstrapper(bar)
      EstimateBootstrapper forBS = new EstimateBootstrapper(forwards)
      EstimateBootstrapper backBS = new EstimateBootstrapper(backwards)

      long MAX_BOOTSTRAP_TRIALS = 100000L
      long bootstrap = min(MAX_BOOTSTRAP_TRIALS, min(volume.length, volume.length))
      if (w == nWindows) {
        logger.info("\n Free Energy Difference:\n")
      } else {
        logger.info(format("\n Free Energy Difference for Window %d\n", w))
      }

      long time = -System.nanoTime()
      forBS.bootstrap(bootstrap)
      time += System.nanoTime()
      logger.fine(format(" Forward FEP Bootstrap Complete:      %7.4f sec", time * Constants.NS2SEC))
      forwardFEP = forBS.getTotalFE()
      forwardEnthalpy = forBS.getTotalEnthalpy()
      double varForeFE = forBS.getTotalUncertainty()
      double varEnthalpyFore = forBS.getTotalEnthalpyUncertainty()
      logger.info(format(" Free energy via Forwards FEP:   %12.4f +/- %6.4f kcal/mol.", forwardFEP,
          varForeFE))

      time = -System.nanoTime()
      backBS.bootstrap(bootstrap)
      time += System.nanoTime()
      logger.fine(format(" Backward FEP Bootstrap Complete:     %7.4f sec", time * Constants.NS2SEC))

      backwardFEP = backBS.getTotalFE()
      backwardEnthalpy = backBS.getTotalEnthalpy()
      double varBackFE = backBS.getTotalUncertainty()
      double varEnthalpyBack = backBS.getTotalEnthalpyUncertainty()
      logger.info(
          format(" Free energy via Backwards FEP:  %12.4f +/- %6.4f kcal/mol.", backwardFEP,
              varBackFE))
      barEnergy = bar.getFreeEnergy()

      logger.info(format(" Free energy via BAR Iteration:  %12.4f +/- %6.4f kcal/mol.", barEnergy,
          bar.getUncertainty()))
      time = -System.nanoTime()
      barBS.bootstrap(bootstrap)
      time += System.nanoTime()
      logger.fine(format(" BAR Bootstrap Complete:              %7.4f sec", time * Constants.NS2SEC))

      barEnergyBS = barBS.getTotalFE()
      double varBARFE = barBS.getTotalUncertainty()
      barEnthalpy = barBS.getTotalEnthalpy()
      double varEnthalpy = barBS.getTotalEnthalpyUncertainty()
      logger.info(
          format(" Free energy via BAR Bootstrap:  %12.4f +/- %6.4f kcal/mol.", barEnergyBS,
              varBARFE))

      if (w == nWindows) {
        logger.info("\n Enthalpy from Potential Energy Averages:\n")

        for (int n = 0; n < nWindows; n++) {
          logger.info(format(" Average Energy for State %d:       %12.4f +/- %6.4f kcal/mol.",
              n, energyMean[n], energySD[n]))

        }
        double enthalpyDiff = energyMean[nWindows - 1] - energyMean[0]
        double enthalpyDiffSD = Math.sqrt(energyVar[nWindows - 1] + energyVar[0])
        logger.info(format(" Enthalpy via Direct Estimate:     %12.4f +/- %6.4f kcal/mol.",
            enthalpyDiff, enthalpyDiffSD))

        logger.info("\n Enthalpy and Entropy:\n")
      } else {
        logger.info(format("\n Enthalpy and Entropy for Window %d\n", w))
      }

      forwardEntropy = (forwardEnthalpy - forwardFEP) / this.temperature
      backwardEntropy = (backwardEnthalpy - backwardFEP) / this.temperature

      logger.info(
          format(" Enthalpy via Forward FEP:       %12.4f +/- %6.4f kcal/mol.", forwardEnthalpy,
              varEnthalpyFore))
      logger.info(format(" Entropy via Forward FEP:        %12.4f kcal/mol/K.", forwardEntropy))
      logger.info(format(" Forward FEP -T*ds Value:        %12.4f kcal/mol.", -(forwardEntropy *
          this.temperature)))

      logger.info(
          format("\n Enthalpy via Backward FEP:      %12.4f +/- %6.4f kcal/mol.", backwardEnthalpy,
              varEnthalpyBack))
      logger.info(format(" Entropy via Backward FEP:       %12.4f kcal/mol/K.", backwardEntropy))
      logger.info(format(" Backward FEP -T*ds Value:       %12.4f kcal/mol.", -(backwardEntropy *
          this.temperature)))

      double tsBar = barEnthalpy - barEnergyBS
      double sBAR = tsBar / (this.temperature)
      logger.info(
          format("\n Enthalpy via BAR:               %12.4f +/- %6.4f kcal/mol.", barEnthalpy,
              varEnthalpy))
      logger.info(format(" Entropy via BAR:                %12.4f kcal/mol/K.", sBAR))
      logger.info(format(" BAR Estimate of -T*ds:          %12.4f kcal/mol.", -(tsBar)))

    }

    return this

  }

  private double[] getEnergyForLambdas(MolecularAssembly[] topologies, double[] lambdaValues,
      String[] arcFileName, double[][] energy, boolean isPBC, int nSymm) {
    for (int j = 0; j < numTopologies; j++) {
      File archiveFile = new File(arcFileName[j])
      openers[j].setFile(archiveFile)
      topologies[j].setFile(archiveFile)
      StringBuilder sb = new StringBuilder(format(
          "\n Evaluating energies for %s\n ", arcFileName[j]))
      logger.info(sb as String)
    }
    int nSnapshots = openers[0].countNumModels()

    double[] x = new double[potential.getNumberOfVariables()]
    double[] vol = new double[nSnapshots]
    for (int k = 0; k < lambdaValues.length; k++) {
      energy[k] = new double[nSnapshots]
    }

    LambdaInterface linter1 = (LambdaInterface) potential

    int endWindow = nWindows - 1
    String endWindows = endWindow + File.separator

    if (arcFileName[0].contains(endWindows)) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
          "Energy Low", "Lambda At", "Energy At"))
    } else if (arcFileName[0].contains("0/")) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda At",
          "Energy At", "Lambda High", "Energy High"))
    } else {
      logger.info(format(" %s     %s   %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
          "Energy Low", "Lambda At", "Energy At", "Lambda High", "Energy High"))
    }

    for (int i = 0; i < nSnapshots; i++) {
      boolean resetPosition = (i == 0)
      for (int n = 0; n < openers.length; n++) {
        openers[n].readNext(resetPosition, false)
      }

      x = potential.getCoordinates(x)
      for (int k = 0; k < lambdaValues.length; k++) {
        double lambda = lambdaValues[k]
        linter1.setLambda(lambda)
        energy[k][i] = potential.energy(x, false)
      }

      if (lambdaValues.length == 2) {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1, lambdaValues[0],
            energy[0][i], lambdaValues[1], energy[1][i]))
      } else {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1,
            lambdaValues[0],
            energy[0][i], lambdaValues[1], energy[1][i], lambdaValues[2], energy[2][i]))
      }

      if (isPBC) {
        Crystal unitCell = potential.getCrystal().getUnitCell()
        vol[i] = unitCell.volume / nSymm
        logger.info(format(" %8d %14.4f",
            i + 1, vol[i]))
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


  double getBarEnergy() {
    return barEnergy
  }

  double getBarEnergyBS() {
    return barEnergyBS
  }

  double getFepFor() {
    return forwardFEP
  }

  double getFepBack() {
    return backwardFEP
  }

  double gethFor() {
    return forwardEnthalpy
  }

  double gethBack() {
    return backwardEnthalpy
  }

  double gethBAR() {
    return barEnthalpy
  }

  double getsFor() {
    return forwardEntropy
  }

  double getsBack() {
    return backwardEntropy
  }

  double getsBAR() {
    return barEntropy
  }

}










