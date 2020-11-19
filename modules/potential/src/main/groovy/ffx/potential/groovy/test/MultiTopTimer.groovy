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
package ffx.potential.groovy.test

import edu.rit.pj.ParallelTeam
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.TimerOptions
import ffx.potential.cli.TopologyOptions
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The MultiTopTimer extends the functionality of the Timer script to handle
 * multi-topology energy functions.
 * <br>
 * Usage:
 * <br>
 * ffxc MultiTopTimer [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Time energy and gradients for 1-4 topologies with optional lambda scaling.", name = "ffxc MultiTopTimer")
class MultiTopTimer extends PotentialScript {

  @Mixin
  private AlchemicalOptions alchemical

  @Mixin
  private TopologyOptions topology

  @Mixin
  private TimerOptions timerOptions

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  private int threadsPer = threadsAvail
  MolecularAssembly[] topologies
  private Potential potential

  /**
   * MultiTopTimer Constructor.
   */
  MultiTopTimer() {
    this(new Binding())
  }

  /**
   * MultiTopTimer Constructor.
   * @param binding Groovy Binding to use.
   */
  MultiTopTimer(Binding binding) {
    super(binding)
  }

  @Override
  MultiTopTimer run() {

    if (!init()) {
      return this
    }

    List<String> arguments = filenames
    // Check nArgs should either be number of arguments (min 1), else 1.
    int nArgs = arguments ? arguments.size() : 1
    nArgs = (nArgs < 1) ? 1 : nArgs

    topologies = new MolecularAssembly[nArgs]

    int numParallel = topology.getNumParallel(threadsAvail, nArgs)
    threadsPer = (int) (threadsAvail / numParallel)

    // Turn on computation of lambda derivatives if softcore atoms exist.
    boolean lambdaTerm = alchemical.hasSoftcore() || topology.hasSoftcore()

    if (lambdaTerm) {
      System.setProperty("lambdaterm", "true")
    }

    double lambda = alchemical.getInitialLambda()

    // Relative free energies via the DualTopologyEnergy class require different
    // default OST parameters than absolute free energies.
    if (nArgs >= 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    List<MolecularAssembly> topologyList = new ArrayList<>(4)

    // Read in files.
    if (!arguments || arguments.isEmpty()) {
      MolecularAssembly molecularAssembly = potentialFunctions.getActiveAssembly()
      if (molecularAssembly == null) {
        logger.info(helpString())
        return this
      }
      arguments = new ArrayList<>()
      arguments.add(molecularAssembly.getFile().getName())
      topologyList.add(alchemical.processFile(topology, molecularAssembly, 0))
    } else {
      logger.info(format(" Initializing %d topologies...", nArgs))
      for (int i = 0; i < nArgs; i++) {
        topologyList.add(alchemical.openFile(potentialFunctions, topology, threadsPer,
            arguments.get(i), i))
      }
    }

    MolecularAssembly[] topologies =
        topologyList.toArray(new MolecularAssembly[topologyList.size()])

    // Configure the potential to test.
    StringBuilder sb = new StringBuilder("\n Testing energies ")
    boolean gradient = !timerOptions.getNoGradient()
    if (gradient) {
      sb.append("and gradients ")
    }
    sb.append("for ")
    potential = topology.assemblePotential(topologies, threadsAvail, sb)

    logger.info(sb.toString())

    LambdaInterface linter = (potential instanceof LambdaInterface) ? (LambdaInterface) potential :
        null
    linter?.setLambda(lambda)

    // End boilerplate opening code.

    long minTime = Long.MAX_VALUE
    int iterations = timerOptions.getIterations()
    int rmsIterations = (int) (((iterations - 1) / 2) + 1)
    double rmsTime = 0
    boolean print = timerOptions.getVerbose()

    int nVars = potential.getNumberOfVariables()
    double[] x = new double[nVars]
    potential.getCoordinates(x)
    double[] g = gradient ? new double[nVars] : null

    String baseString = gradient ? " Energy and gradient " : " Energy "

    for (int i = 0; i < iterations; i++) {
      long time = -System.nanoTime()
      if (gradient) {
        if (potential instanceof ForceFieldEnergyOpenMM) {
          ((ForceFieldEnergyOpenMM) potential).setLambda(lambda);
          ((ForceFieldEnergyOpenMM) potential).updateParameters();
        }
        potential.energyAndGradient(x, g, print)
      } else {
        potential.energy(x, print)
      }
      time += System.nanoTime()
      minTime = time < minTime ? time : minTime
      int rmsIndex = i + (rmsIterations - iterations)
      double thisTime = time * 1.0E-9
      logger.info(format("%s%d in %14.5g seconds.", baseString, (i + 1), thisTime))

      if (rmsIndex >= 0) {
        thisTime *= thisTime
        rmsTime += thisTime
      }
    }
    rmsTime /= ((double) rmsIterations)
    rmsTime = Math.sqrt(rmsTime)

    logger.info(format(" Minimum time: %14.5g (sec)", minTime * 1.0E-9))
    logger.info(format(" RMS time (latter half): %14.5g (sec)", rmsTime))
    for (int i = 0; i < topologies.size(); i++) {
      logger.info(format(" Number of threads for topology %d: %d", i, threadsPer))
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(potential)
    }
    return potentials
  }

}
