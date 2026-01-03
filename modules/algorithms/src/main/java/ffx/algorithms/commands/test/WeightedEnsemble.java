//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.WeightedEnsembleManager;
import ffx.algorithms.dynamics.WeightedEnsembleManager.OneDimMetric;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;

/**
 * WeightedEnsemble
 * <br>
 * Usage: Accelerated Sampling with Weighted Ensemble
 * <br>
 * ffxc test.WeightedEnsemble [options] &lt;filename&gt; [file2...];
 */
@Command(description = " Runs parallel simulations with intermittent resampling.", name = "test.WeightedEnsemble")
public class WeightedEnsemble extends AlgorithmsCommand {

  @Mixin
  private AtomSelectionOptions atomSelectionOptions;

  @Mixin
  private DynamicsOptions dynamicsOptions;

  @Mixin
  private WriteoutOptions writeOutOptions;

  @Option(names = {"--stepsPer"}, paramLabel = "10000", defaultValue = "10000",
      description = "Number of steps to take between resampling cycles.")
  private int stepsPer;

  @Option(names = {"--initDynamics"}, paramLabel = "10000", defaultValue = "10000",
      description = "Number of initialization steps to take before windows start. This is good for getting diverse starting structures.")
  private int initDynamics;

  @Option(names = {"--numPerBin"}, paramLabel = "2", defaultValue = "2",
      description = "Number of walkers per bin.")
  private int numPerBin;

  @Option(names = {"--oneDimensionalMetric"}, paramLabel = "RMSD", defaultValue = "RMSD",
      description = "Bin across this metric. Options: RMSD, POTENTIAL, RESIDUE_DISTANCE, ATOM_DISTANCE, COM_DISTANCE, RADIUS_OF_GYRATION")
  private String oneDimensionalMetric;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private String filename;

  /**
   * Constructor.
   */
  public WeightedEnsemble() {
    super();
  }

  /**
   * Constructor.
   * @param binding The Binding to use.
   */
  public WeightedEnsemble(FFXBinding binding) {
    super(binding);
  }

  /**
   * WeightedEnsemble constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public WeightedEnsemble(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public WeightedEnsemble run() {
    if (!init()) {
      return this;
    }
    dynamicsOptions.init();

    // Set up FFX integrator
    MolecularAssembly assembly = getActiveAssembly(filename);
    File file = assembly.getFile();
    if (file == null) {
      logger.severe(" No file found for assembly: " + assembly);
    } else {
      logger.info(" Running Weighted Ensemble on " + file);
    }
    Potential potential = assembly.getPotentialEnergy();
    double[] x = new double[potential.getNumberOfVariables()];
    potential.getCoordinates(x);
    potential.energy(x, true);
    MolecularDynamics md = dynamicsOptions.getDynamics(writeOutOptions, potential, assembly, algorithmListener);

    // Set up & run Weighted Ensemble
    OneDimMetric metric = null;
    try {
      metric = OneDimMetric.valueOf(oneDimensionalMetric);
    } catch (IllegalArgumentException e) {
      logger.severe(" Invalid oneDimensionalMetric: " + oneDimensionalMetric);
      return this;
    }
    WeightedEnsembleManager weightedEnsemble = new WeightedEnsembleManager(metric, numPerBin, md, file, true);
    weightedEnsemble.run(dynamicsOptions.getSteps(), stepsPer, dynamicsOptions.getTemperature(), dynamicsOptions.getDt());

    return this;
  }
}
