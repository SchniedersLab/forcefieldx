//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.cli.LambdaParticleOptions
import ffx.algorithms.cli.MultiDynamicsOptions
import ffx.algorithms.cli.OSTOptions
import ffx.algorithms.cli.ThermodynamicsOptions
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Histogram script prints out a Orthogonal Space histogram from a *.his file.
 * <br>
 * Usage:
 * <br>
 * ffxc Histogram [options] &lt;filename&gt;
 */
@Command(description = " Evaluate the Orthogonal Space Histogram.", name = "Histogram")
class Histogram extends AlgorithmsScript {

  /**
   * -s or --save Save the histogram, PMF and 2D bias to files.
   */
  @Option(names = ['-s', '--save'], paramLabel = 'false',
      description = 'Save the bias histogram to histogram.txt, the total PMF to pmf.txt, and 2D PMF to pmf.2D.txt')
  boolean save = false

  /**
   * -b or --bias By default, the PMF is saved. This flag flips the sign to give the OST bias.
   */
  @Option(names = ['-b', '--bias'], paramLabel = 'false',
      description = 'By default, the PMF is saved. This flag flips the sign to give the OST bias.')
  boolean bias = false

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private String filename

  private OrthogonalSpaceTempering orthogonalSpaceTempering
  private File saveDir = null

  /**
   * Histogram Constructor.
   */
  Histogram() {
    this(new Binding())
  }

  /**
   * Histogram Constructor.
   * @param binding The Groovy Binding to use.
   */
  Histogram(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  Histogram run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    println("\n Evaluating Histogram for " + filename)

    File structureFile = new File(FilenameUtils.normalize(filename))
    structureFile = new File(structureFile.getAbsolutePath())
    String baseFilename = FilenameUtils.removeExtension(structureFile.getName())
    File histogramRestart = new File(baseFilename + ".his")
    File lambdaRestart = null

    // Get a reference to the active system's ForceFieldEnergy and atom array.
    ForceFieldEnergy energy = activeAssembly.getPotentialEnergy()

    // Print the current energy
    energy.energy(true, true)
    logger.info("")

    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    // Construct some options with defaults.
    DynamicsOptions dynamicsOptions = new DynamicsOptions()
    LambdaParticleOptions lambdaParticleOptions = new LambdaParticleOptions()
    MultiDynamicsOptions multiDynamicsOptions = new MultiDynamicsOptions()
    ThermodynamicsOptions thermodynamicsOptions = new ThermodynamicsOptions()
    OSTOptions ostOptions = new OSTOptions()

    // Construct the Thermodynamics instance.
    orthogonalSpaceTempering = ostOptions.constructOST(energy,
        lambdaRestart, histogramRestart, activeAssembly, null,
        dynamicsOptions, thermodynamicsOptions, lambdaParticleOptions,
        algorithmListener, !multiDynamicsOptions.isSynchronous())

    if (save) {
      orthogonalSpaceTempering.setMolecularAssembly(activeAssembly)
      OrthogonalSpaceTempering.Histogram histogram = orthogonalSpaceTempering.getHistogram()
      histogram.updateFreeEnergyDifference(false, true)
      StringBuffer sb = histogram.evaluateTotalOSTBias(bias)

      String dirName = saveDir.toString() + File.separator
      String file = dirName + "pmf.txt"
      logger.info(" Writing " + file)
      FileWriter fileWriter = new FileWriter(file)
      fileWriter.write(sb.toString())
      fileWriter.close()

      sb = histogram.evaluate2DOSTBias(bias)
      file = dirName + "pmf.2D.txt"
      logger.info(" Writing " + file)
      fileWriter = new FileWriter(file)
      fileWriter.write(sb.toString())
      fileWriter.close()
    }
    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (orthogonalSpaceTempering == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList((Potential) orthogonalSpaceTempering)
    }
    return potentials
  }

}
