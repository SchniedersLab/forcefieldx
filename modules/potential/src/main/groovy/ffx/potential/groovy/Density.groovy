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

package ffx.potential.groovy

import ffx.crystal.Crystal
import ffx.numerics.math.SummaryStatistics
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.WriteoutOptions
import ffx.potential.parsers.SystemFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Density script calculates system density for single topology systems.
 * <br>
 * Usage:
 * <br>
 * ffxc Density [options] &lt;filename&gt;
 */
@Command(description = " Calculates system density.", name = "ffxc Density")
class Density extends PotentialScript {

  @Mixin
  private WriteoutOptions writeout

  /**
   * -s or --start First frame to evaluate (1-indexed).
   */
  @Option(names = ['-s', '--start'], paramLabel = "1", defaultValue = "1",
      description = 'First frame to evaluate (1-indexed).')
  private int start = 1

  /**
   * -f or --final Last frame to evaluate (1-indexed); values less than 1 evaluate to end of trajectory.
   */
  @Option(names = ['-f', '--final'], paramLabel = "all frames", defaultValue = "0",
      description = 'Last frame to evaluate (1-indexed); values less than 1 evaluate to end of trajectory.')
  private int finish = 0

  /**
   * --st or --stride Stride: evaluate density every N frames. Must be positive.
   */
  @Option(names = ['--st', '--stride'], paramLabel = "1", defaultValue = "1",
      description = "Stride: evaluate density every N frames. Must be positive.")
  private int stride = 1

  /**
   * -p or --printout writes out a file with density adjusted to match mean calculated density.
   */
  @Option(names = ['-p', '--printout'], defaultValue = "false",
      description = "Print out a file with density adjusted to match mean calculated density.")
  private boolean doPrint = false

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'An atomic coordinate file in PDB or XYZ format.')
  String filename = null

  /**
   * Density constructor.
   */
  Density() {
    this(new Binding())
  }

  /**
   * Density constructor.
   * @param binding The Groovy Binding to use.
   */
  Density(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Density run() {

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

    SystemFilter openFilter = potentialFunctions.getFilter()
    double totMass = activeAssembly.getMass()
    Crystal crystal = activeAssembly.getCrystal().getUnitCell()

    if (crystal.aperiodic()) {
      logger.info(format(" System %s appears aperiodic: total mass %16.7g g/mol", filename, totMass))
    } else {
      double volume = crystal.volume / crystal.numSymOps
      double density = crystal.getDensity(totMass)
      int nFrames = openFilter.countNumModels()
      double[] densities = new double[nFrames]
      int lastFrame = (finish < 1) ? nFrames : finish

      densities[0] = density
      logger.info(
          format(" Evaluating density for system %s with mass %16.7g (g/mol).", filename, totMass))
      logger.info(
          format(" Density at frame %9d is %16.7g (g/mL) from a volume of %16.7g (A^3)", 1, density,
              volume))

      int ctr = 1
      // TODO: Optimize by skipping frames by stride.
      while (openFilter.readNext(false, false)) {
        volume = crystal.volume / crystal.numSymOps
        density = crystal.getDensity(totMass)
        logger.info(
            format(" Density at frame %9d is %16.7g (g/mL) from a volume of %16.7g (A^3)", ++ctr,
                density, volume))
        int i = ctr - 1
        densities[i] = density
      }

      SummaryStatistics densStats = new SummaryStatistics(densities, start - 1, lastFrame, stride)
      logger.info(" Summary statistics for density:")
      logger.info(densStats.toString())

      if (doPrint) {
        crystal.setDensity(densStats.mean, totMass)
        String outFileName = FilenameUtils.removeExtension(filename)
        writeout.saveFile(outFileName, potentialFunctions, activeAssembly)
      }
    }

    return this
  }
}
