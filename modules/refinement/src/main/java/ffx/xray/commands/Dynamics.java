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
package ffx.xray.commands;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementEnergy;
import ffx.xray.cli.XrayOptions;
import ffx.xray.refine.RefinementMode;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * The X-ray Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Dynamics [options] &lt;filename&gt;
 */
@Command(description = " Run Dynamics on an X-ray target.", name = "xray.Dynamics")
public class Dynamics extends AlgorithmsCommand implements AlgorithmListener {

  @Mixin
  private XrayOptions xrayOptions;

  @Mixin
  private DynamicsOptions dynamicsOptions;

  @Mixin
  AtomSelectionOptions atomSelectionOptions;

  @Mixin
  private WriteoutOptions writeoutOptions;

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
  private RefinementEnergy refinementEnergy;

  /**
   * Dynamics constructor.
   */
  public Dynamics() {
    super();
  }

  /**
   * Dynamics constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Dynamics(String[] args) {
    super(args);
  }

  /**
   * Dynamics constructor.
   *
   * @param binding The Binding to use.
   */
  public Dynamics(FFXBinding binding) {
    super(binding);
  }

  @Override
  public Dynamics run() {

    if (!init()) {
      return this;
    }

    dynamicsOptions.init();
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

    logger.info("\n Running xray.Dynamics on " + filename);

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties();
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);
    refinementEnergy = xrayOptions.toXrayEnergy(diffractionData);

    // Log the energy of each MolecularAssembly
    algorithmFunctions.energy(molecularAssemblies);

    // Restart is currently only supported for COORDINATES mode.
    File dyn = null;
    if (xrayOptions.refinementMode == RefinementMode.COORDINATES) {
      dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn");
      if (!dyn.exists()) {
        dyn = null;
      }
    }

    MolecularDynamics molecularDynamics = dynamicsOptions.getDynamics(writeoutOptions,
        refinementEnergy, activeAssembly, this);
    refinementEnergy.setThermostat(molecularDynamics.getThermostat());
    boolean initVelocities = true;
    molecularDynamics.dynamic(dynamicsOptions.getSteps(), dynamicsOptions.getDt(), dynamicsOptions.getReport(),
        dynamicsOptions.getWrite(), dynamicsOptions.getTemperature(), initVelocities, dyn);

    // Print the final refinement statistics.
    diffractionData.scaleBulkFit();
    diffractionData.printStats();

    // Print the final energy of each conformer.
    algorithmFunctions.energy(molecularAssemblies);

    logger.info(" ");
    diffractionData.writeModel(FilenameUtils.removeExtension(filename) + ".pdb");

    if (mtz) {
      diffractionData.writeData(FilenameUtils.removeExtension(filename) + ".mtz");
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

  @Override
  public boolean algorithmUpdate(MolecularAssembly active) {
    logger.info(" R/Rfree " + diffractionData.printOptimizationUpdate());
    return true;
  }
}
