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
package ffx.algorithms.groovy.test

import edu.rit.pj.Comm
import ffx.algorithms.cli.RepexOSTOptions
import ffx.algorithms.cli.ThermodynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.groovy.Thermodynamics
import ffx.algorithms.thermodynamics.MonteCarloOST
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.algorithms.thermodynamics.RepExOST
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin

import java.util.stream.Collectors

/**
 * The RepexThermo command uses the Orthogonal Space Tempering with histogram replica exchange to estimate a free energy difference.
 * <br>
 * Usage:
 * <br>
 * ffxc test.RepexThermo [options] &lt;filename&gt [file2...];
 */
@Command(description = " Use Orthogonal Space Tempering with histogram replica exchange to estimate a free energy difference.", name = "test.RepexThermo")
class RepexThermo extends Thermodynamics {

  @Mixin
  RepexOSTOptions repex

  private RepExOST repExOST
  private CrystalPotential finalPotential

  /**
   * RepexThermo Constructor.
   */
  RepexThermo() {
    this(new Binding())
  }

  /**
   * RepexThermo Constructor.
   * @param binding The Groovy Binding to use.
   */
  RepexThermo(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  RepexThermo run() {

    // Begin boilerplate "make a topology" code.
    if (!init()) {
      return this
    }

    boolean fromActive
    List<String> arguments
    if (filenames) {
      arguments = filenames
      fromActive = false
    } else {
      logger.warning(" Untested: use of active assembly instead of provided filenames!")
      MolecularAssembly mola = algorithmFunctions.getActiveAssembly()
      if (mola == null) {
        logger.info(helpString())
        return this
      }
      arguments = Collections.singletonList(mola.getFile().getName())
      fromActive = true
    }

    int nArgs = arguments.size()

    topologies = new MolecularAssembly[nArgs]

    int numParallel = topologyOptions.getNumParallel(threadsAvail, nArgs)
    threadsPer = (int) (threadsAvail / numParallel)

    // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
    /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
    The Minimize script, for example, may be running on a single, unscaled physical topology. */
    boolean lambdaTerm = (
        nArgs == 1 || alchemicalOptions.hasSoftcore() || topologyOptions.hasSoftcore())

    if (lambdaTerm) {
      System.setProperty("lambdaterm", "true")
    }

    // Relative free energies via the DualTopologyEnergy class require different
    // default OST parameters than absolute free energies.
    if (nArgs >= 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    List<MolecularAssembly> topologyList = new ArrayList<>(nArgs)

    Comm world = Comm.world()
    int size = world.size()
    if (size < 2) {
      logger.severe(" RepexThermo requires multiple processes, found only one!")
    }
    int rank = (size > 1) ? world.rank() : 0

    double initLambda = alchemicalOptions.getInitialLambda(size, rank)

    // Segment of code for MultiDynamics and OST.
    List<File> structureFiles = arguments.stream().
        map {fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath())
        }.
        collect(Collectors.toList())

    File firstStructure = structureFiles.get(0)
    String filePathNoExtension = firstStructure.getAbsolutePath().replaceFirst(~/\.[^.]+$/, "")

    // SEGMENT DIFFERS FROM THERMODYNAMICS.

    String filepath = FilenameUtils.getFullPath(filePathNoExtension)
    String fileBase = FilenameUtils.getBaseName(FilenameUtils.getName(filePathNoExtension))
    String rankDirName = String.format("%s%d", filepath, rank)
    File rankDirectory = new File(rankDirName)
    if (!rankDirectory.exists()) {
      rankDirectory.mkdir()
    }
    // @formatter:off
    rankDirName = "${rankDirName}${File.separator}"
    String withRankName = "${rankDirName}${fileBase}"
    File lambdaRestart = new File("${withRankName}.lam")
    // @formatter:on

    boolean lamExists = lambdaRestart.exists()

    // Read in files.
    logger.info(String.format(" Initializing %d topologies...", nArgs))
    if (fromActive) {
      topologyList.add(alchemicalOptions.processFile(topologyOptions, activeAssembly, 0))
    } else {
      for (int i = 0; i < nArgs; i++) {
        topologyList.add(multiDynamicsOptions.openFile(algorithmFunctions, topologyOptions,
            threadsPer, arguments.get(i), i, alchemicalOptions, structureFiles.get(i), rank))
      }
    }

    MolecularAssembly[] topologies =
        topologyList.toArray(new MolecularAssembly[topologyList.size()])

    StringBuilder sb = new StringBuilder("\n Running ")
    switch (thermodynamicsOptions.getAlgorithm()) {
    // Labeled case blocks needed because Groovy (can't tell the difference between a closure and an anonymous code block).
      case ThermodynamicsOptions.ThermodynamicsAlgorithm.OST:
        ostAlg:
        {
          sb.append("Orthogonal Space Tempering")
        }
        break
      default:
        defAlg:
        {
          throw new IllegalArgumentException(
              " RepexThermo currently does not support fixed-lambda alchemy!")
        }
        break
    }
    sb.append(" with histogram replica exchange for ")

    potential = (CrystalPotential) topologyOptions.assemblePotential(topologies, threadsAvail, sb)

    LambdaInterface linter = (LambdaInterface) potential
    logger.info(sb.toString())

    double[] x = new double[potential.getNumberOfVariables()]
    potential.getCoordinates(x)
    linter.setLambda(initLambda)
    potential.energy(x, true)

    if (nArgs == 1) {
      randomSymopOptions.randomize(topologies[0])
    }

    multiDynamicsOptions.distribute(topologies, potential, algorithmFunctions, rank, size)

    boolean isMC = ostOptions.isMonteCarlo()
    boolean twoStep = ostOptions.isTwoStep()
    MonteCarloOST mcOST = null
    MolecularDynamics md

    if (thermodynamicsOptions.getAlgorithm() == ThermodynamicsOptions.ThermodynamicsAlgorithm.OST) {
      // @formatter:off
      File firstHisto = new File("${filepath}0${File.separator}${fileBase}.his")
      // @formatter:on

      orthogonalSpaceTempering =
          ostOptions.constructOST(potential, lambdaRestart, firstHisto, topologies[0],
              additionalProperties, dynamicsOptions, thermodynamicsOptions, lambdaParticleOptions,
              algorithmListener,
              false)
      finalPotential = ostOptions.applyAllOSTOptions(orthogonalSpaceTempering, topologies[0],
          dynamicsOptions, barostatOptions)

      if (isMC) {
        mcOST = ostOptions.
            setupMCOST(orthogonalSpaceTempering, topologies, dynamicsOptions, thermodynamicsOptions,
                verbose,
                algorithmListener)
        md = mcOST.getMD()
      } else {
        md = ostOptions.
            assembleMolecularDynamics(topologies, finalPotential, dynamicsOptions, algorithmListener)
      }
      if (!lamExists) {
        if (finalPotential instanceof LambdaInterface) {
          ((LambdaInterface) finalPotential).setLambda(initLambda)
        } else {
          orthogonalSpaceTempering.setLambda(initLambda)
        }
      }

      CompositeConfiguration allProperties = new CompositeConfiguration(
          topologies[0].getProperties())
      if (additionalProperties != null) {
        allProperties.addConfiguration(additionalProperties)
      }

      for (int i = 1; i < size; i++) {
        // @formatter:off
        File rankIHisto = new File("${filepath}${i}${File.separator}${fileBase}.his")
        // @formatter:on
        orthogonalSpaceTempering.addHistogram(ostOptions.generateHistogramSettings(rankIHisto,
            lambdaRestart.toString(), allProperties, i, dynamicsOptions, lambdaParticleOptions, true,
            false))
      }

      if (isMC) {
        repExOST = RepExOST.repexMC(orthogonalSpaceTempering, mcOST, dynamicsOptions, ostOptions,
            topologies[0].getProperties(), writeoutOptions.getFileType(), twoStep,
            repex.getRepexFrequency())
      } else {
        repExOST = RepExOST.repexMD(orthogonalSpaceTempering, md, dynamicsOptions, ostOptions,
            topologies[0].getProperties(), writeoutOptions.getFileType(), repex.getRepexFrequency())
      }

      long eSteps = thermodynamicsOptions.getEquilSteps()
      if (eSteps > 0) {
        repExOST.mainLoop(eSteps, true)
      }
      repExOST.mainLoop(dynamicsOptions.getNumSteps(), false)
    } else {
      logger.severe(" RepexThermo currently does not support fixed-lambda alchemy!")
    }

    // @formatter:off
    logger.info(" ${thermodynamicsOptions.getAlgorithm()} with Histogram Replica Exchange Done.")
    // @formatter:on

    return this
  }


  @Override
  OrthogonalSpaceTempering getOST() {
    return repExOST == null ? null : repExOST.getOST()
  }

  @Override
  CrystalPotential getPotential() {
    return (repExOST == null) ? potential : repExOST.getOST()
  }

  @Override
  List<Potential> getPotentials() {
    if (repExOST == null) {
      return
      potential == null ? Collections.emptyList() : Collections.singletonList((Potential) potential)
    }
    return Collections.singletonList((Potential) repExOST.getOST())
  }
}
