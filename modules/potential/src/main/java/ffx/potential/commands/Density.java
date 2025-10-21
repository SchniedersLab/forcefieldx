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
package ffx.potential.commands;

import ffx.crystal.Crystal;
import ffx.numerics.math.SummaryStatistics;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.cli.WriteoutOptions;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import static java.lang.String.format;

/**
 * Calculate the mean system density for single topology systems.
 * <p>
 * Usage:
 * ffxc Density [options] &lt;filename&gt;
 */
@Command(name = "Density", description = " Calculate the mean system density for single topology systems.")
public class Density extends PotentialCommand {

  @Mixin
  private WriteoutOptions writeoutOptions = new WriteoutOptions();

  /**
   * First frame to evaluate (1-indexed).
   */
  @Option(names = {"-s", "--start"}, paramLabel = "1", defaultValue = "1",
      description = "First frame to evaluate (1-indexed).")
  private int start = 1;

  /**
   * Last frame to evaluate (1-indexed); values less than 1 evaluate to end of trajectory.
   */
  @Option(names = {"-f", "--final"}, paramLabel = "all frames", defaultValue = "0",
      description = "Last frame to evaluate (1-indexed); values less than 1 evaluate to end of trajectory.")
  private int finish = 0;

  /**
   * Stride: evaluate density every N frames. Must be positive.
   */
  @Option(names = {"--st", "--stride"}, paramLabel = "1", defaultValue = "1",
      description = "Stride: evaluate density every N frames. Must be positive.")
  private int stride = 1;

  /**
   * Print out a file with density adjusted to match mean calculated density.
   */
  @Option(names = {"-p", "--printout"}, defaultValue = "false",
      description = "Print out a file with density adjusted to match mean calculated density.")
  private boolean doPrint = false;

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "An atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  public Density() {
    super();
  }

  public Density(FFXBinding binding) {
    super(binding);
  }

  public Density(String[] args) {
    super(args);
  }

  @Override
  public Density run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    SystemFilter openFilter = potentialFunctions.getFilter();
    double totMass = activeAssembly.getMass();
    Crystal crystal = activeAssembly.getCrystal().getUnitCell();

    if (crystal.aperiodic()) {
      logger.info(format(" System %s is aperiodic: total mass %16.7g g/mol", filename, totMass));
    } else {
      double volume = crystal.volume / crystal.getNumSymOps();
      double density = crystal.getDensity(totMass);
      int nFrames = openFilter.countNumModels();
      double[] densities = new double[nFrames];
      int lastFrame = (finish < 1) ? nFrames : finish;

      densities[0] = density;
      logger.info(format("\n Evaluating density for system %s with mass %16.7g (g/mol).", filename, totMass));
      logger.info(format(" Density at frame %9d is %16.7g (g/mL) from a volume of %16.7g (A^3)", 1, density, volume));

      int ctr = 1;
      while (openFilter.readNext(false, false)) {
        volume = crystal.volume / crystal.getNumSymOps();
        density = crystal.getDensity(totMass);
        logger.info(format(" Density at frame %9d is %16.7g (g/mL) from a volume of %16.7g (A^3)", ++ctr, density, volume));
        int i = ctr - 1;
        densities[i] = density;
      }

      SummaryStatistics densStats = new SummaryStatistics(densities, start - 1, lastFrame, stride);
      logger.info(" Summary statistics for density:");
      logger.info(densStats.toString());

      if (doPrint) {
        activeAssembly.setFractionalMode(MolecularAssembly.FractionalMode.MOLECULE);
        activeAssembly.computeFractionalCoordinates();
        crystal.setDensity(densStats.mean, totMass);
        activeAssembly.moveToFractionalCoordinates();
        saveByExtension(activeAssembly, filename, writeoutOptions.fileType);
        // Update filename to the written file.
        filename = activeAssembly.getFile().getAbsolutePath();
        logger.info("\n Saved density adjusted coordinates: " + filename);
        logger.info("\n Molecule based fractional coordinates were maintained.");
      }
    }

    return this;
  }
}
