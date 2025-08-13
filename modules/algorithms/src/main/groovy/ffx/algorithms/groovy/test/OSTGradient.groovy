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
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.cli.LambdaParticleOptions
import ffx.algorithms.cli.MultiDynamicsOptions
import ffx.algorithms.cli.OSTOptions
import ffx.algorithms.cli.ThermodynamicsOptions
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.GradientOptions
import ffx.potential.utils.GradientUtils
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Option
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * OSTGradient tests the Orthogonal Space Tempering Potential.
 * <br>
 * Usage:
 * <br>
 * ffxc test.OSTGradient [options] &lt;filename&gt;
 */
@Command(description = " OSTGradient script tests the Orthogonal Space Tempering Potential.", name = "test.OSTGradient")
class OSTGradient extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  @Mixin
  GradientOptions gradientOptions

  @Mixin
  AlchemicalOptions alchemicalOptions

  /**
   * --meta or --metaDynamics Use a 1D metadynamics bias.
   */
  @Option(names = ["--meta", "--metaDynamics"],
      defaultValue = "false", description = "Use a 1D metadynamics style bias.")
  private boolean metaDynamics

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file", description = "XYZ or PDB input file.")
  private String filename

  private OrthogonalSpaceTempering orthogonalSpaceTempering
  private File saveDir = null

  public double dUdLError = 0.0
  public int nFailures = 0

  /**
   * CrystalSearch Constructor.
   */
  OSTGradient() {
    this(new groovy.lang.Binding())
  }

  /**
   * CrystalSearch Constructor.
   * @param binding The Groovy Binding to use.
   */
  OSTGradient(Binding binding) {
    super(binding)
  }

  @Override
  OSTGradient run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Turn on computation of lambda derivatives.
    System.setProperty("lambdaterm", "true")

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    File structureFile = activeAssembly.getFile()
    filename = structureFile.getAbsolutePath()

    logger.info("\n Evaluating OST Gradient for " + filename)

    String baseFilename = FilenameUtils.removeExtension(filename)
    File histogramRestart = new File(baseFilename + ".his")
    File lambdaRestart = null

    if (!histogramRestart.exists()) {
      logger.warning("\n Histogram restart file does not exist: " + histogramRestart.toString())
    } else if (!histogramRestart.canRead()) {
      logger.warning("\n Histogram restart file can not be read." + histogramRestart.toString())
    }

    // Get a reference to the active system's ForceFieldEnergy and atom array.
    ForceFieldEnergy energy = activeAssembly.getPotentialEnergy()

    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    // Construct some options with defaults.
    DynamicsOptions dynamicsOptions = new DynamicsOptions()
    LambdaParticleOptions lambdaParticleOptions = new LambdaParticleOptions()
    MultiDynamicsOptions multiDynamicsOptions = new MultiDynamicsOptions()
    ThermodynamicsOptions thermodynamicsOptions = new ThermodynamicsOptions()
    OSTOptions ostOptions = new OSTOptions()
    // Apply the metaDynamics flag.
    ostOptions.setMetaDynamics(metaDynamics)

    // Construct the Thermodynamics instance.
    orthogonalSpaceTempering = ostOptions.constructOST(energy, lambdaRestart, histogramRestart, activeAssembly, null,
        dynamicsOptions, thermodynamicsOptions, lambdaParticleOptions,
        null, !multiDynamicsOptions.isSynchronous())

    // Get the lambda value to test.
    double lambda = alchemicalOptions.getInitialLambda(false)
    logger.info("\n Lambda value: " + lambda)

    // Set the alchemical atoms.
    alchemicalOptions.setFirstSystemAlchemistry(activeAssembly)
    // Set the uncharged atoms.
    alchemicalOptions.setFirstSystemUnchargedAtoms(activeAssembly)
    // Set the active atoms.
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    /*
     * Stop propagating lambda to prevent adding new Gaussian potentials
     * to the bias, which would introduce artifacts into the
     * finite-difference derivatives.
     */
    orthogonalSpaceTempering.setPropagateLambda(false)
    orthogonalSpaceTempering.setLambda(lambda)
    int n = orthogonalSpaceTempering.getNumberOfVariables()

    assert (n % 3 == 0)
    n = (int) (n / 3)

    // Finite-difference step size.
    double[] x = new double[3 * n]
    orthogonalSpaceTempering.getCoordinates(x)
    // Calculate the energy and analytic dE/dX
    double[] g = new double[3 * n]
    double e = orthogonalSpaceTempering.energyAndGradient(x, g)
    double dEdLambda = orthogonalSpaceTempering.getTotaldEdLambda()

    // Calculate the finite-difference dEdL
    double step = gradientOptions.getDx()
    orthogonalSpaceTempering.setLambda(lambda + step)
    double lp = orthogonalSpaceTempering.energy(x)
    orthogonalSpaceTempering.setLambda(lambda - step)
    double lm = orthogonalSpaceTempering.energy(x)
    double dedl = (lp - lm) / (2.0 * step)

    // Report the results.
    logger.info(format(" Analytic dE/dL:   %15.8f", dEdLambda))
    logger.info(format(" Numeric  dE/dL:   %15.8f\n", dedl))
    dUdLError = Math.abs(dEdLambda - dedl)

    // Check the gradient.
    orthogonalSpaceTempering.setLambda(lambda)
    GradientUtils gradientUtils = new GradientUtils(orthogonalSpaceTempering)
    nFailures = gradientUtils.testGradient(gradientOptions)

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
