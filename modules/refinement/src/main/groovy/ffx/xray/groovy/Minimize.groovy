//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.xray.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.DiffractionData
import ffx.xray.RefinementMinimize
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.cli.XrayOptions
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The X-ray Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Minimize [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Refine an X-ray/Neutron target.", name = "xray.Minimize")
class Minimize extends AlgorithmsScript {

  @Mixin
  MinimizeOptions minimizeOptions

  @Mixin
  XrayOptions xrayOptions

  /**
   * -t or --threeStage Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true).
   */
  @Option(names = ['-t', '--threeStage'], paramLabel = 'false',
      description = 'Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true)')
  boolean threeStage = false

  /**
   * -E or --eps3 RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determine eps for each stage).
   */
  @Option(names = ['-E', '--eps3'], paramLabel = '-1.0', arity = '3',
      description = 'RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determines eps for each stage).')
  double[] eps3 = [-1.0, -1.0, -1.0]

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames
  private MolecularAssembly[] molecularAssemblies
  private DiffractionData diffractionData

  /**
   * Minimize constructor.
   */
  Minimize() {
    this(new Binding())
  }

  /**
   * Minimize constructor.
   * @param binding The Groovy Binding to use.
   */
  Minimize(Binding binding) {
    super(binding)
  }

  @Override
  Minimize run() {

    if (!init()) {
      return this
    }

    xrayOptions.init()

    String filename
    if (filenames != null && filenames.size() > 0) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = molecularAssemblies[0]
      filename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      molecularAssemblies = [activeAssembly]
      filename = activeAssembly.getFile().getAbsolutePath()
    }

    logger.info("\n Running xray.Minimize on " + filename)

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties()
    xrayOptions.setProperties(parseResult, properties)

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties)
    diffractionData.scaleBulkFit()
    diffractionData.printStats()
    algorithmFunctions.energy(molecularAssemblies)

    // RMS gradient convergence criteria for three stage refinement
    double coordinateEPS = eps3[0]
    double bfactorEPS = eps3[1]
    double occupancyEPS = eps3[2]

    // The number of corrections used in the BFGS update.
    int nBFGS = minimizeOptions.getNBFGS()

    // Maximum number of refinement cycles.
    int maxIterations = minimizeOptions.getIterations()

    if (threeStage) {
      // Coordinates
      RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData,
          RefinementMode.COORDINATES)
      if (coordinateEPS < 0.0) {
        coordinateEPS = refinementMinimize.getEps()
      }
      if (maxIterations < Integer.MAX_VALUE) {
        logger.info(format(
            "\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", coordinateEPS,
            maxIterations))
      } else {
        logger.info(format("\n RMS gradient convergence criteria: %8.5f", coordinateEPS))
      }
      refinementMinimize.minimize(nBFGS, coordinateEPS, maxIterations)
      diffractionData.scaleBulkFit()
      diffractionData.printStats()
      algorithmFunctions.energy(molecularAssemblies)

      // B-factors
      refinementMinimize = new RefinementMinimize(diffractionData, RefinementMode.BFACTORS)
      if (bfactorEPS < 0.0) {
        bfactorEPS = refinementMinimize.getEps()
      }
      if (maxIterations < Integer.MAX_VALUE) {
        logger.info(
            format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", bfactorEPS,
                maxIterations))
      } else {
        logger.info(format("\n RMS gradient convergence criteria: %8.5f", bfactorEPS))
      }
      refinementMinimize.minimize(nBFGS, bfactorEPS, maxIterations)
      diffractionData.scaleBulkFit()
      diffractionData.printStats()

      // Occupancies
      if (
          diffractionData.getAltResidues().size() > 0 || diffractionData.getAltMolecules().size() >
              0) {
        refinementMinimize = new RefinementMinimize(diffractionData, RefinementMode.OCCUPANCIES)
        if (occupancyEPS < 0.0) {
          occupancyEPS = refinementMinimize.getEps()
        }
        if (maxIterations < Integer.MAX_VALUE) {
          logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d",
              occupancyEPS, maxIterations))
        } else {
          logger.info(format("\n RMS gradient convergence criteria: %8.5f", occupancyEPS))
        }
        refinementMinimize.minimize(occupancyEPS, maxIterations)
        diffractionData.scaleBulkFit()
        diffractionData.printStats()
      } else {
        logger.info("Occupancy refinement not necessary, skipping")
      }

    } else {
      // Type of refinement.
      RefinementMode refinementMode = xrayOptions.refinementMode
      RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData, refinementMode)
      double eps = minimizeOptions.eps
      if (eps < 0.0) {
        eps = refinementMinimize.getEps()
      }

      if (maxIterations < Integer.MAX_VALUE) {
        logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps,
            maxIterations))
      } else {
        logger.info(format("\n RMS gradient convergence criteria: %8.5f", eps))
      }
      refinementMinimize.minimize(eps, maxIterations)
      diffractionData.scaleBulkFit()
      diffractionData.printStats()
    }

    // Print the final energy of each conformer.
    algorithmFunctions.energy(molecularAssemblies)

    logger.info(" ")
    diffractionData.writeModel(removeExtension(filename) + ".pdb")
    diffractionData.writeData(removeExtension(filename) + ".mtz")

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (molecularAssemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(molecularAssemblies).
          filter {a -> a != null
          }.map {a -> a.getPotentialEnergy()
      }.filter {e -> e != null
      }.collect(Collectors.toList())
    }
  }

  @Override
  boolean destroyPotentials() {
    return diffractionData == null ? true : diffractionData.destroy()
  }
}
