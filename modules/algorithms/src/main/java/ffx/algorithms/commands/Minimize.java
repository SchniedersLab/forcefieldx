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
package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsScript;
import ffx.algorithms.cli.MinimizeOptions;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FileUtils;
import groovy.lang.Binding;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.String.format;

/**
 * The Minimize script uses a limited-memory BFGS algorithm to minimize the energy of a molecular system.
 * <br>
 * Usage:
 * <br>
 * ffxc Minimize [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run L-BFGS minimization on a system.", name = "Minimize")
public class Minimize extends AlgorithmsScript {

  @Mixin
  private MinimizeOptions minimizeOptions;

  @Mixin
  private AtomSelectionOptions atomSelectionOptions;

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Mixin
  private TopologyOptions topologyOptions;

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "Atomic coordinate files in PDB or XYZ format.")
  private List<String> filenames = null;

  private MolecularAssembly[] topologies;
  private Potential potential;

  /**
   * The minimization algorithm.
   */
  private ffx.algorithms.optimize.Minimize minimize;

  /**
   * Minimize Constructor.
   */
  public Minimize() {
    super();
  }

  /**
   * Minimize Constructor.
   * @param binding The Groovy Binding to use.
   */
  public Minimize(Binding binding) {
    super(binding);
  }

  /**
   * Minimize constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public Minimize(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Minimize run() {

    if (!init()) {
      return this;
    }

    // Determine the number of topologies to be read and allocate the array.
    int numTopologies = topologyOptions.getNumberOfTopologies(filenames);
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];

    // Turn on computation of lambda derivatives if softcore atoms exist.
    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    // Read in files.
    if (filenames == null || filenames.isEmpty()) {
      activeAssembly = getActiveAssembly(null);
      if (activeAssembly == null) {
        logger.info(helpString());
        return this;
      }
      filenames = new ArrayList<>();
      filenames.add(activeAssembly.getFile().getName());
      topologies[0] = alchemicalOptions.processFile(topologyOptions, activeAssembly, 0);
    } else {
      logger.info(format(" Initializing %d topologies...", numTopologies));
      for (int i = 0; i < numTopologies; i++) {
        topologies[i] = alchemicalOptions.openFile(algorithmFunctions,
            topologyOptions, threadsPerTopology, filenames.get(i), i);
      }
      activeAssembly = topologies[0];
    }

    if (topologies.length == 1) {
      atomSelectionOptions.setActiveAtoms(topologies[0]);
    }

    // Configure the potential to test.
    StringBuilder sb = new StringBuilder("\n Minimizing energy of ");
    potential = topologyOptions.assemblePotential(topologies, sb);
    logger.info(sb.toString());

    LambdaInterface linter = (potential instanceof LambdaInterface) ? (LambdaInterface) potential : null;
    if (linter != null) {
      linter.setLambda(alchemicalOptions.getInitialLambda());
    }

    SystemFilter systemFilter = algorithmFunctions.getFilter();

    double[] x = new double[potential.getNumberOfVariables()];
    potential.getCoordinates(x);
    potential.energy(x, true);

    minimize = new ffx.algorithms.optimize.Minimize(topologies[0], potential, algorithmListener);
    minimize.minimize(minimizeOptions.getNBFGS(), minimizeOptions.getEps(), minimizeOptions.getIterations());

    potential.getCoordinates(x);
    activeAssembly = systemFilter.getActiveMolecularSystem();
    updateTitle(potential.energy(x, true));
    if (topologies.length > 1) {
      // Handle Multiple Topology Cases.
      for (MolecularAssembly molecularAssembly : topologies) {
        String modelFilename = molecularAssembly.getFile().getAbsolutePath();

        if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
          baseDir = new File(FilenameUtils.getFullPath(modelFilename));
        }

        String dirName = baseDir.toString() + File.separator;
        String fileName = FilenameUtils.getName(modelFilename);
        String ext = FilenameUtils.getExtension(fileName);
        fileName = FilenameUtils.removeExtension(fileName);

        if (ext.toUpperCase().contains("XYZ")) {
          algorithmFunctions.saveAsXYZ(molecularAssembly, new File(dirName + fileName + ".xyz"));
        } else {
          algorithmFunctions.saveAsPDB(molecularAssembly, new File(dirName + fileName + ".pdb"));
        }
      }
    } else {
      // Handle Single Topology Cases.
      setActiveAssembly(topologies[0]);
      String modelFilename = activeAssembly.getFile().getAbsolutePath();
      if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
        baseDir = new File(FilenameUtils.getFullPath(modelFilename));
      }

      String dirName = baseDir.toString() + File.separator;
      String fileName = FilenameUtils.getName(modelFilename);
      String ext = FilenameUtils.getExtension(fileName);
      fileName = FilenameUtils.removeExtension(fileName);
      File saveFile;
      SystemFilter writeFilter;
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
            minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations());
            PDBFilter pdbFilter = (PDBFilter) systemFilter;
            pdbFilter.writeFile(saveFile, true, false, false);
          } else if (systemFilter instanceof XYZFilter) {
            minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations());
            writeFilter.writeFile(saveFile, true);
          }
        }
        if (systemFilter instanceof PDBFilter) {
          FileUtils.append(saveFile, "END\n");
        }
      }
    }

    return this;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    if (minimize != null) {
      sb.append(minimize.toString());
    } else {
      sb.append("No minimization has been performed.");
    }
    return sb.toString();
  }

  /**
   * Get the list of energies recorded during minimization.
   * @return A List of Double energies, or null if minimization was not run.
   */
  public List<Double> getEnergyList() {
    if (minimize != null) {
      return minimize.getEnergyList();
    }
    return null;
  }

  /**
   * Get the final RMS gradient after minimization.
   * @return The final RMS gradient or NaN if minimization was not run.
   */
  public double getRMSGradient() {
    if (minimize != null) {
      return minimize.getRMSGradient();
    }
    return Double.NaN;
  }

  /**
   * Get the final energy after minimization.
   * @return The final energy or NaN if minimization was not run.
   */
  public double getEnergy() {
    if (minimize != null) {
      return minimize.getEnergy();
    }
    return Double.NaN;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials;
    if (potential == null) {
      potentials = Collections.emptyList();
    } else {
      potentials = Collections.singletonList(potential);
    }
    return potentials;
  }

}
