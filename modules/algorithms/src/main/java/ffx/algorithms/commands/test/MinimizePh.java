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
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.MinimizeOptions;
import ffx.algorithms.optimize.PhMinimize;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XPHFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import ffx.utilities.FileUtils;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;

/**
 * The MinimizePh command uses a limited-memory BFGS algorithm to minimize the energy of a CpHMD molecular system.
 * <br>
 * Usage:
 * <br>
 * ffxc test.MinimizePh [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Run L-BFGS minimization on a CpHMD system.", name = "test.MinimizePh")
public class MinimizePh extends AlgorithmsCommand {

  @Mixin
  private MinimizeOptions minimizeOptions;

  /**
   * --pH or --constantPH Constant pH value for the test.
   */
  @Option(names = {"--pH", "--constantPH"}, paramLabel = "7.4", defaultValue = "7.4",
      description = "pH value for the energy evaluation. (Only applies when esvTerm is true)")
  private double pH;

  /**
   * --coords
   */
  @Option(names = {"--coords"}, paramLabel = "false", defaultValue = "false",
      description = "Minimize spatial coordinates along with titration")
  private boolean coords;

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "Atomic coordinate files in PDB or XYZ format.")
  private String filename;

  private ForceFieldEnergy forceFieldEnergy;

  /**
   * MinimizePh Constructor.
   */
  public MinimizePh() {
    super();
  }

  /**
   * MinimizePh Constructor.
   *
   * @param binding The Binding to use.
   */
  public MinimizePh(FFXBinding binding) {
    super(binding);
  }

  /**
   * MinimizePh constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public MinimizePh(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public MinimizePh run() {

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

    logger.info("\n Running MinimizePh on " + filename);
    forceFieldEnergy = activeAssembly.getPotentialEnergy();

    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, pH, null);
    esvSystem.setConstantPh(pH);
    int numESVs = esvSystem.getExtendedResidueList().size();
    forceFieldEnergy.attachExtendedSystem(esvSystem);
    logger.info(format(" Attached extended system with %d residues.", numESVs));

    double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
    forceFieldEnergy.getCoordinates(x);
    forceFieldEnergy.energy(x, true);

    SystemFilter systemFilter = algorithmFunctions.getFilter();
    if (systemFilter instanceof XYZFilter) {
      XPHFilter xphFilter = new XPHFilter(activeAssembly.getFile(), activeAssembly,
          activeAssembly.getForceField(), activeAssembly.getProperties(), esvSystem);
      xphFilter.readFile();
      logger.info(" Reading ESV lambdas from XPH file");
      forceFieldEnergy.getCoordinates(x);
      forceFieldEnergy.energy(x, true);
    }
    PhMinimize minimize = new PhMinimize(activeAssembly, forceFieldEnergy, algorithmListener,
        esvSystem);

    double energy = minimize.getEnergy();
    double tolerance = 1.0e-10;

    if (coords) {
      while (true) {
        // Complete a round of coordinate optimization.
        minimize.minimizeCoordinates(minimizeOptions.getNBFGS(),
            minimizeOptions.getEps(), minimizeOptions.getIterations());
        double newEnergy = minimize.getEnergy();
        int status = minimize.getStatus();
        if (status != 0) {
          break;
        }
        if (abs(newEnergy - energy) <= tolerance) {
          break;
        }
        energy = newEnergy;

        // Complete a round of titration optimization.
        minimize.minimizeTitration(minimizeOptions.getNBFGS(),
            minimizeOptions.getEps(), minimizeOptions.getIterations());
        newEnergy = minimize.getEnergy();
        status = minimize.getStatus();
        if (status != 0) {
          break;
        }
        if (abs(newEnergy - energy) <= tolerance) {
          break;
        }
        energy = newEnergy;
      }
    } else {
      minimize.minimizeTitration(minimizeOptions.getNBFGS(),
          minimizeOptions.getEps(), minimizeOptions.getIterations());
    }

    forceFieldEnergy.getCoordinates(x);
    forceFieldEnergy.energy(x, true);
    esvSystem.writeRestart();

    // Handle Single Topology Cases.
    String modelFilename = activeAssembly.getFile().getAbsolutePath();
    if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
      baseDir = new File(FilenameUtils.getFullPath(modelFilename));
    }

    String dirName = baseDir.toString() + File.separator;
    String fileName = FilenameUtils.getName(modelFilename);
    String ext = FilenameUtils.getExtension(fileName);
    fileName = FilenameUtils.removeExtension(fileName);
    File saveFile = null;
    SystemFilter writeFilter = null;
    if (ext.toUpperCase().contains("XYZ")) {
      saveFile = new File(dirName + fileName + ".xyz");
      writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties());
      algorithmFunctions.saveAsXYZ(activeAssembly, saveFile);
    } else if (ext.toUpperCase().contains("ARC")) {
      saveFile = new File(dirName + fileName + ".arc");
      saveFile = algorithmFunctions.versionFile(saveFile);
      writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties());
      algorithmFunctions.saveAsXYZ(activeAssembly, saveFile);
    } else {
      saveFile = new File(dirName + fileName + ".pdb");
      saveFile = algorithmFunctions.versionFile(saveFile);

      PDBFilter pdbFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties());
      writeFilter = pdbFilter;
      int numModels = systemFilter.countNumModels();
      if (numModels > 1) {
        pdbFilter.setModelNumbering(0);
      }
      pdbFilter.writeFile(saveFile, true, false, false);
    }

    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
      while (systemFilter.readNext()) {
        Crystal crystal = activeAssembly.getCrystal();
        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy();
        forceFieldEnergy.setCrystal(crystal);
        if (systemFilter instanceof PDBFilter) {
          FileUtils.append(saveFile, "ENDMDL\n");
          minimize.minimizeTitration(minimizeOptions.getNBFGS(), minimizeOptions.getEps(), minimizeOptions.getIterations());
          PDBFilter pdbFilter = (PDBFilter) systemFilter;
          pdbFilter.writeFile(saveFile, true, false, false);
        } else {
          minimize.minimizeTitration(minimizeOptions.getNBFGS(), minimizeOptions.getEps(), minimizeOptions.getIterations());
          writeFilter.writeFile(saveFile, true);
        }
      }
      if (systemFilter instanceof PDBFilter) {
        FileUtils.append(saveFile, "END\n");
      }
    }

    return this;
  }
}
