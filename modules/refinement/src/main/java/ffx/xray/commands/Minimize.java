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
package ffx.xray.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.MinimizeOptions;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementMinimize;
import ffx.xray.cli.XrayOptions;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * The X-ray Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Minimize [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Refine an X-ray/Neutron target.", name = "xray.Minimize")
public class Minimize extends AlgorithmsCommand {

  @Mixin
  private XrayOptions xrayOptions;

  @Mixin
  AtomSelectionOptions atomSelectionOptions;

  @Mixin
  private MinimizeOptions minimizeOptions;

  /**
   * -t or --threeStage Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true).
   */
  @Option(names = {"-t", "--threeStage"}, paramLabel = "false",
      description = "Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true)")
  private boolean threeStage = false;

  /**
   * -E or --eps3 RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determine eps for each stage).
   */
  @Option(names = {"-E", "--eps3"}, paramLabel = "-1.0", arity = "3",
      description = "RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determines eps for each stage).")
  private double[] eps3 = {-1.0, -1.0, -1.0};

  /**
   * -c or --cycles Number of refinement cycles.
   */
  @Option(names = {"-c", "--cycles"}, paramLabel = "1", defaultValue = "1",
      description = "Number of refinement cycles.")
  private int ncycles = 1;

  /**
   * --mtz Write out MTZ containing structure factor coefficients.
   */
  @Option(names = {"--mtz"}, paramLabel = "false",
      description = "Write out an MTZ containing structure factor coefficients.")
  private boolean mtz = false;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames;

  private MolecularAssembly[] molecularAssemblies;
  private DiffractionData diffractionData;

  /**
   * Minimize constructor.
   */
  public Minimize() {
    super();
  }

  /**
   * Minimize constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Minimize(String[] args) {
    super(args);
  }

  /**
   * Minimize constructor.
   *
   * @param binding The Binding to use.
   */
  public Minimize(FFXBinding binding) {
    super(binding);
  }

  @Override
  public Minimize run() {

    if (!init()) {
      return this;
    }

    xrayOptions.init();

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.getFirst());
      activeAssembly = molecularAssemblies[0];
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
    }

    // Update the active filename
    filename = activeAssembly.getFile().getAbsolutePath();

    // Apply active atom flags.
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      atomSelectionOptions.setActiveAtoms(molecularAssembly);
    }

    logger.info("\n Running xray.Minimize on " + filename);

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties();
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);

    // Log the energy of each MolecularAssembly
    algorithmFunctions.energy(molecularAssemblies);

    // RMS gradient convergence criteria for three stage refinement
    double coordinateEPS = eps3[0];
    double bfactorEPS = eps3[1];
    double occupancyEPS = eps3[2];

    // The number of corrections used in the BFGS update.
    int nBFGS = minimizeOptions.getNBFGS();

    // Maximum number of refinement cycles.
    int maxIterations = minimizeOptions.getIterations();

    for (int cycle = 0; cycle < ncycles; cycle++) {
      if (ncycles > 1) {
        logger.info(format("\n Refinement Cycle: %d\n", cycle + 1));
      }
      diffractionData.scaleBulkFit();
      diffractionData.printStats();

      if (threeStage) {
        // Coordinates
        RefinementModel refinementModel = diffractionData.getRefinementModel();
        RefinementMode refinementMode = RefinementMode.COORDINATES;
        refinementModel.setRefinementMode(refinementMode);
        if (refinementModel.getNumCoordParameters() > 0) {
          RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData);
          if (coordinateEPS < 0.0) {
            coordinateEPS = refinementMode.getDefaultEps();
          }
          if (maxIterations < Integer.MAX_VALUE) {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d",
                coordinateEPS, maxIterations));
          } else {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f", coordinateEPS));
          }
          refinementMinimize.minimize(nBFGS, coordinateEPS, maxIterations);
          algorithmFunctions.energy(molecularAssemblies);
        }

        // B-factors
        refinementMode = RefinementMode.BFACTORS;
        refinementModel.setRefinementMode(refinementMode);
        if (refinementModel.getNumBFactorParameters() > 0) {
          RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData);
          if (bfactorEPS < 0.0) {
            boolean hasAnisous = refinementModel.getNumANISOU() > 0;
            bfactorEPS = refinementMode.getDefaultEps(hasAnisous);
          }
          if (maxIterations < Integer.MAX_VALUE) {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d",
                bfactorEPS, maxIterations));
          } else {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f", bfactorEPS));
          }
          refinementMinimize.minimize(nBFGS, bfactorEPS, maxIterations);
        }

        // Occupancies
        refinementMode = RefinementMode.OCCUPANCIES;
        refinementModel.setRefinementMode(refinementMode);
        if (refinementModel.getNumOccupancyParameters() > 0) {
          RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData);
          if (occupancyEPS < 0.0) {
            occupancyEPS = refinementMode.getDefaultEps();
          }
          if (maxIterations < Integer.MAX_VALUE) {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d",
                occupancyEPS, maxIterations));
          } else {
            logger.info(format("\n RMS gradient convergence criteria: %8.5f", occupancyEPS));
          }
          refinementMinimize.minimize(occupancyEPS, maxIterations);
        } else {
          logger.info("Occupancy refinement not necessary, skipping");
        }
      } else {
        RefinementMinimize refinementMinimize = new RefinementMinimize(diffractionData);
        RefinementModel refinementModel = diffractionData.getRefinementModel();
        RefinementMode refinementMode = refinementModel.getRefinementMode();
        double eps = minimizeOptions.getEps();
        if (eps < 0.0) {
          boolean hasAnisous = refinementModel.getNumANISOU() > 0;
          eps = refinementMode.getDefaultEps(hasAnisous);
        }
        if (maxIterations < Integer.MAX_VALUE) {
          logger.info(format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps,
              maxIterations));
        } else {
          logger.info(format("\n RMS gradient convergence criteria: %8.5f", eps));
        }
        refinementMinimize.minimize(eps, maxIterations);
      }
    }

    diffractionData.scaleBulkFit();
    diffractionData.printStats();

    // Print the final energy of each conformer.
    algorithmFunctions.energy(molecularAssemblies);

    logger.info(" ");
    diffractionData.writeModel(removeExtension(filename) + ".pdb");
    if (mtz) {
      diffractionData.writeData(removeExtension(filename) + ".mtz");
    }

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }

  @Override
  public boolean destroyPotentials() {
    return diffractionData == null ? true : diffractionData.destroy();
  }
}
