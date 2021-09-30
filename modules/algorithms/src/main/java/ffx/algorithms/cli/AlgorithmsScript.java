// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.algorithms.cli;

import static java.lang.String.format;

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.AlgorithmUtils;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.utilities.FFXScript;
import groovy.lang.Binding;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Base class for scripts in the Algorithms package, providing some key functions.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AlgorithmsScript extends FFXScript {

  /** An instance of AlgorithmFunctions passed into the current context. */
  public AlgorithmFunctions algorithmFunctions;

  /**
   * An active MolecularAssembly passed into the current context or loaded by the Script from a file
   * argument.
   */
  public MolecularAssembly activeAssembly;

  /** An instance of the AlgorithmListener interface. */
  public AlgorithmListener algorithmListener;

  /** The directory in which to place output files. Mostly for tests. */
  protected File baseDir;

  public AlgorithmsScript() {
    this(new Binding());
  }

  public AlgorithmsScript(Binding binding) {
    super(binding);
  }

  /**
   * Reclaims resources associated with all Potential objects associated with this script.
   *
   * @return If all Potentials had resources reclaimed.
   */
  public boolean destroyPotentials() {
    boolean allSucceeded = true;
    for (Potential potent : getPotentials()) {
      logger.fine(format(" Potential %s is being destroyed. ", potent));
      allSucceeded = allSucceeded && potent.destroy();
    }
    return allSucceeded;
  }

  /**
   * Returns a List of all Potential objects associated with this script.
   *
   * @return All Potentials. Sometimes empty, never null.
   */
  public List<Potential> getPotentials() {
    List<Potential> plist = new ArrayList<>();
    if (activeAssembly != null && activeAssembly.getPotentialEnergy() != null) {
      plist.add(activeAssembly.getPotentialEnergy());
    }
    return plist;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Execute the BaseScript init method, then load algorithm functions.
   */
  @Override
  public boolean init() {
    if (!super.init()) {
      return false;
    }

    Binding binding = getBinding();

    if (binding.hasVariable("functions")) {
      algorithmFunctions = (AlgorithmFunctions) binding.getVariable("functions");
    } else {
      algorithmFunctions = new AlgorithmUtils();
      binding.setVariable("functions", algorithmFunctions);
    }

    activeAssembly = null;
    if (binding.hasVariable("active")) {
      activeAssembly = (MolecularAssembly) binding.getVariable("active");
    }

    algorithmListener = null;
    if (binding.hasVariable("listener")) {
      algorithmListener = (AlgorithmListener) binding.getVariable("listener");
    }

    if (binding.hasVariable("baseDir")) {
      baseDir = (File) binding.getVariable("baseDir");
    }

    return true;
  }

  /**
   * Sets the directory this script should save files to. Mostly used for tests.
   *
   * @param baseDir Directory to save output to.
   */
  public void setBaseDir(File baseDir) {
    this.baseDir = baseDir;
  }

  /**
   * Gets a File in the save directory with the same name as the input file. Can just be the original
   * file if saveDir was never set, which is the case for production runs.'
   *
   * @param file File to find a save location for.
   * @return File to save to
   */
  protected File saveDirFile(File file) {
    if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
      return file;
    } else {
      String baseName = file.getName();
      String newName = baseDir.getAbsolutePath() + File.separator + baseName;
      return new File(newName);
    }
  }

  /**
   * If a filename is supplied, open it and return the MolecularAssembly. Otherwise, the current
   * activeAssembly is returned (which may be null).
   *
   * @param filename Filename to open.
   * @return The active assembly.
   */
  public MolecularAssembly getActiveAssembly(String filename) {
    if (filename != null) {
      // Open the supplied file.
      MolecularAssembly[] assemblies = {algorithmFunctions.open(filename)};
      activeAssembly = assemblies[0];
    }
    return activeAssembly;
  }

  /**
   * If a filename is supplied, open it and return the MolecularAssemblies. Otherwise, the current
   * activeAssembly is returned (which may be null).
   *
   * @param filename Filename to open.
   * @return The active assemblies.
   */
  public MolecularAssembly[] getActiveAssemblies(String filename) {
    MolecularAssembly[] assemblies;
    if (filename != null) {
      // Open the supplied file.
      assemblies = algorithmFunctions.openAll(filename);
      activeAssembly = assemblies[0];
      return assemblies;
    } else {
      assemblies = new MolecularAssembly[] {activeAssembly};
    }
    return assemblies;
  }

  /**
   * Set the Active Assembly. This is a work-around for a strange Groovy static compilation bug where
   * direct assignment of activeAssembly in Groovy scripts that extend AlgorithmsScript fails (a NPE
   * results).
   *
   * @param molecularAssembly The MolecularAssembly that should be active.
   */
  public void setActiveAssembly(MolecularAssembly molecularAssembly) {
    activeAssembly = molecularAssembly;
  }

}
