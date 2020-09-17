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
import ffx.numerics.math.SummaryStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
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

  @Option(names = ["--dV", "--volume"], paramLabel = "false",
      description = "Write out snapshot volumes to the Tinker BAR file.")
  private boolean includeVolume = false

  @Option(names = ["--tb", "--tinkerBAR"], paramLabel = "false",
      description = "Write out a Tinker BAR file.")
  private boolean tinkerBAR = false

  /**
   * The final argument(s) should be filenames for lambda windows in order..
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
  List<String> filenames = null

  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  private int threadsPer = threadsAvail
  MolecularAssembly[] topologies1
  MolecularAssembly[] topologies2
  SystemFilter[] openers1
  SystemFilter[] openers2

  CrystalPotential potential1
  CrystalPotential potential2

  private Configuration additionalProperties1
  private Configuration additionalProperties2

  /**
   * Sets an optional Configuration with additional properties.
   * @param additionalProps
   */
  void setProperties(Configuration additionalProps1, Configuration additionalProps2) {
    this.additionalProperties1 = additionalProps1
    this.additionalProperties2 = additionalProps2
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

    if (filenames == null || filenames.size() % 2 != 0) {
      return this
    }

    int nFiles = filenames.size()
    int nPer = (int) (nFiles / 2)

    topologies1 = new MolecularAssembly[nPer]
    topologies2 = new MolecularAssembly[nPer]
    openers1 = new SystemFilter[nPer]
    openers2 = new SystemFilter[nPer]

    int numParallel = topology.getNumParallel(threadsAvail, nPer)
    threadsPer = (int) (threadsAvail / numParallel)

    // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
    /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
    The Minimize script, for example, may be running on a single, unscaled physical topology. */
    boolean lambdaTerm = (nPer == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

    if (lambdaTerm) {
      System.setProperty("lambdaterm", "true")
    }

    // Relative free energies via the DualTopologyEnergy class require different
    // default OST parameters than absolute free energies.
    if (nPer >= 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    logger.info(format(" Initializing %d topologies for each end", nPer))
    for (int i = 0; i < nPer; i++) {
      MolecularAssembly ma =
          alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i)
      topologies1[i] = ma
      openers1[i] = algorithmFunctions.getFilter()
      int iSecond = i + nPer
      ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[iSecond], i)
      topologies2[i] = ma
      openers2[i] = algorithmFunctions.getFilter()
    }

    double lambda1 = alchemical.initialLambda

    StringBuilder sb = new StringBuilder(format(
        "\n Using BAR to analyze a free energy change between L=%.5f and L=%.5f for\n ", lambda1,
        lambda2))
    potential1 = (CrystalPotential) topology.assemblePotential(topologies1, threadsAvail, sb)
    sb.append(" and ")
    potential2 = (CrystalPotential) topology.assemblePotential(topologies2, threadsAvail, sb)

    LambdaInterface linter1 = (LambdaInterface) potential1
    LambdaInterface linter2 = (LambdaInterface) potential2
    logger.info(sb.toString())

    // Check for periodic boundary conditions
    Crystal unitCell1 = potential1.getCrystal().getUnitCell()
    Crystal unitCell2 = potential2.getCrystal().getUnitCell()
    boolean isPBC = includeVolume && !unitCell1.aperiodic()
    isPBC = isPBC && !unitCell2.aperiodic()
    isPBC = isPBC && (unitCell1.getNumSymOps() == unitCell2.getNumSymOps())
    int nSymm = 0
    if (isPBC) {
      nSymm = unitCell1.getNumSymOps()
    }

    int nSnapshots1 = openers1[0].countNumModels()
    int nSnapshots2 = openers2[0].countNumModels()
    double[] e1L1 = new double[nSnapshots1]
    double[] e1L2 = new double[nSnapshots1]
    double[] e2L1 = new double[nSnapshots2]
    double[] e2L2 = new double[nSnapshots2]
    double[] eDiff1 = new double[nSnapshots1]
    double[] eDiff2 = new double[nSnapshots2]
    double[] vol1 = null
    double[] vol2 = null
    if (isPBC) {
      vol1 = new double[nSnapshots1]
      vol2 = new double[nSnapshots2]
    }

    logger.info("\n Preliminary energy evaluation for first end of the window.")
    double[] x1 = new double[potential1.getNumberOfVariables()]
    x1 = potential1.getCoordinates(x1)
    linter1.setLambda(lambda1)
    e1L1[0] = potential1.energy(x1, true)
    linter1.setLambda(lambda2)
    e1L2[0] = potential1.energy(x1, false)
    eDiff1[0] = e1L2[0] - e1L1[0]
    if (isPBC) {
      Crystal unitCell = potential1.getCrystal().getUnitCell()
      vol1[0] = unitCell.volume / nSymm
    }

    logger.info("\n Preliminary energy evaluation for second end of the window.")
    double[] x2 = new double[potential2.getNumberOfVariables()]
    potential2.getCoordinates(x2)
    linter2.setLambda(lambda2)
    e2L2[0] = potential2.energy(x2, true)
    e2L1[0] = potential2.energy(x2, false)
    eDiff2[0] = e2L2[0] - e2L1[0]

    if (isPBC) {
      Crystal unitCell = potential2.getCrystal().getUnitCell()
      vol2[0] = unitCell.volume / nSymm
    }

    String lamString1 = format("%.3f", lambda1)
    String lamString2 = format("%.3f", lambda2)

    logger.info(format("\n Ensemble 1 collected at L=%s", lamString1))
    if (isPBC) {
      logger.info(
          format(" Snapshot     E(L=%s)     E(L=%s)             dE         Volume", lamString1,
              lamString2))
    } else {
      logger.info(format(" Snapshot     E(L=%s)     E(L=%s)             dE", lamString1, lamString2))
    }

    for (int i = 1; i < nSnapshots1; i++) {
      linter1.setLambda(lambda1)
      for (int j = 0; j < nPer; j++) {
        openers1[j].readNext(false, false)
      }
      x1 = potential1.getCoordinates(x1)
      e1L1[i] = potential1.energy(x1, false)
      linter1.setLambda(lambda2)
      e1L2[i] = potential1.energy(x1, false)
      eDiff1[i] = e1L2[i] - e1L1[i]
      if (isPBC) {
        Crystal unitCell = potential1.getCrystal().getUnitCell()
        vol1[i] = unitCell.volume / nSymm
        logger.info(format(" %8d %14.4f %14.4f %14.4f %14.4f",
            i + 1, e1L1[i], e1L2[i], eDiff1[i], vol1[i]))
      } else {
        logger.info(format(" %8d %14.4f %14.4f %14.4f",
            i + 1, e1L1[i], e1L2[i], eDiff1[i]))
      }
    }

    logger.info(format("\n Ensemble 2 collected at L=%s", lamString2))
    if (isPBC) {
      logger.info(
          format(" Snapshot     E(L=%s)     E(L=%s)             dE         Volume", lamString1,
              lamString2))
    } else {
      logger.info(format(" Snapshot     E(L=%s)     E(L=%s)             dE", lamString1, lamString2))
    }

    for (int i = 1; i < nSnapshots2; i++) {
      linter2.setLambda(lambda1)
      for (int j = 0; j < nPer; j++) {
        openers2[j].readNext(false, false)
      }
      x2 = potential2.getCoordinates(x2)
      e2L1[i] = potential2.energy(x2, false)
      linter2.setLambda(lambda2)
      e2L2[i] = potential2.energy(x2, false)
      eDiff2[i] = e2L2[i] - e2L1[i]
      if (isPBC) {
        Crystal unitCell = potential2.getCrystal().getUnitCell()
        vol2[i] = unitCell.volume / nSymm
        logger.info(format(" %8d %14.4f %14.4f %14.4f %14.4f",
            i + 1, e2L1[i], e2L2[i], eDiff2[i], vol2[i]))
      } else {
        logger.info(format(" %8d %14.4f %14.4f %14.4f",
            i + 1, e2L1[i], e2L2[i], eDiff2[i]))
      }
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

      new File(barFileName).withWriter {bw ->
        StringBuilder fnames1 = new StringBuilder()
        for (int i = 0; i < nPer; i++) {
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
        for (int i = nPer; i < nFiles; i++) {
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
    }

    // Load energy values for FEP and BAR
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

    double[] lambda = new double[] {lambda1, lambda2}
    SequentialEstimator bar = new BennettAcceptanceRatio(lambda, energyLow, energyAt, energyHigh,
        new double[] {temp1, temp2})
    SequentialEstimator forwards = bar.getInitialForwardsGuess()
    SequentialEstimator backwards = bar.getInitialBackwardsGuess()

    EstimateBootstrapper barBS = new EstimateBootstrapper(bar)
    EstimateBootstrapper forBS = new EstimateBootstrapper(forwards)
    EstimateBootstrapper backBS = new EstimateBootstrapper(backwards)

    long MAX_BOOTSTRAP_TRIALS = 100000L
    long bootstrap = min(MAX_BOOTSTRAP_TRIALS, min(nSnapshots1, nSnapshots2))

    logger.info("\n Free Energy Difference via FEP Method\n")
    long time = -System.nanoTime()
    forBS.bootstrap(bootstrap)
    time += System.nanoTime()
    logger.fine(format(" Forward FEP Bootstrap Complete:      %7.4f sec", time * Constants.NS2SEC))
    double sumFE = forBS.getTotalFE()
    double varFE = forBS.getTotalUncertainty()
    logger.info(format(" Free energy via Forwards FEP:   %12.4f +/- %6.4f kcal/mol.", sumFE, varFE))

    time = -System.nanoTime()
    backBS.bootstrap(bootstrap)
    time += System.nanoTime()
    logger.fine(format(" Backward FEP Bootstrap Complete:     %7.4f sec", time * Constants.NS2SEC))
    sumFE = backBS.getTotalFE()
    varFE = backBS.getTotalUncertainty()
    logger.info(format(" Free energy via Backwards FEP:  %12.4f +/- %6.4f kcal/mol.", sumFE, varFE))

    logger.info("\n Free Energy Difference via BAR Method\n")
    logger.info(
        format(" Free energy via BAR Iteration:  %12.4f +/- %6.4f kcal/mol.", bar.getFreeEnergy(),
            bar.getUncertainty()))
    time = -System.nanoTime()
    barBS.bootstrap(bootstrap)
    time += System.nanoTime()
    logger.fine(format(" BAR Bootstrap Complete:              %7.4f sec", time * Constants.NS2SEC))
    sumFE = barBS.getTotalFE()
    varFE = barBS.getTotalUncertainty()
    logger.info(format(" Free energy via BAR Bootstrap:  %12.4f +/- %6.4f kcal/mol.", sumFE, varFE))

    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential1 == null && potential2 == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = new ArrayList<>()
      if (potential1 != null) {
        potentials.add(potential1)
      }
      if (potential2) {
        potentials.add(potential2)
      }
    }
    return potentials
  }
}
