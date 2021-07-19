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
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Minimize script uses OpenMM accelerated L-BFGS algorithm to minimize the 
 * energy of a molecular system.
 * <br>
 * Usage:
 * <br>
 * ffxc MinimizeOpenMM [options] &lt;filename&gt;
 */
@Command(description = " Run OpenMM Accelerated L-BFGS minimization on a system.", name = "ffxc MinimizeOpenMM")
class MinimizeOpenMM extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  @Mixin
  MinimizeOptions minimizeOptions


  /**
   * A PDB or XYZ filename.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private String filename

  private ForceFieldEnergy forceFieldEnergy

  /**
   * MinimizeOpenMM Constructor.
   */
  MinimizeOpenMM() {
    this(new Binding())
  }

  /**
   * MinimizeOpenMM Constructor.
   * @param binding The Groovy Binding to use.
   */
  MinimizeOpenMM(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  MinimizeOpenMM run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    if (System.getProperty("platform") != null && !System.getProperty("platform").isEmpty()) {
      System.setProperty("platform", System.getProperty("platform"))
    } else {
      System.setProperty("platform", "OMM")
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    atomSelectionOptions.setActiveAtoms(activeAssembly)

    forceFieldEnergy = activeAssembly.getPotentialEnergy()
    switch (forceFieldEnergy.getPlatform()) {
      case ForceFieldEnergy.Platform.OMM:
      case ForceFieldEnergy.Platform.OMM_CUDA:
      case ForceFieldEnergy.Platform.OMM_OPENCL:
      case ForceFieldEnergy.Platform.OMM_OPTCPU:
      case ForceFieldEnergy.Platform.OMM_REF:
        logger.fine(" Platform is appropriate for OpenMM Minimization.")
        break
      case ForceFieldEnergy.Platform.FFX:
      default:
        logger.severe(String.format(
            " Platform %s is inappropriate for OpenMM minimization. Please explicitly specify an OpenMM platform.",
            forceFieldEnergy.getPlatform()))
        break
    }

    if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
      ffx.algorithms.optimize.MinimizeOpenMM minimizeOpenMM = new ffx.algorithms.optimize.MinimizeOpenMM(
          activeAssembly)
      minimizeOpenMM.minimize(minimizeOptions.eps, minimizeOptions.iterations)

      if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
        baseDir = new File(FilenameUtils.getFullPath(filename))
      }

      String dirName = baseDir.toString() + File.separator
      String name = FilenameUtils.getName(filename)
      String ext = FilenameUtils.getExtension(name)
      name = FilenameUtils.removeExtension(name)

      File saveFile
      if (ext.toUpperCase().contains("XYZ")) {
        saveFile = new File(dirName + name + ".xyz")
        algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
      } else if (ext.toUpperCase().contains("ARC")) {
        saveFile = new File(dirName + name + ".arc")
        algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
      } else {
        saveFile = new File(dirName + name + ".pdb")
        algorithmFunctions.saveAsPDB(activeAssembly, saveFile)
      }

      SystemFilter systemFilter = algorithmFunctions.getFilter()
      saveFile = activeAssembly.getFile()

      if (systemFilter instanceof XYZFilter) {
        XYZFilter xyzFilter = (XYZFilter) systemFilter
        while (xyzFilter.readNext()) {
          minimizeOpenMM.minimize(minimizeOptions.eps, minimizeOptions.iterations)
          boolean append = true
          xyzFilter.writeFile(saveFile, append)
        }
      }
    } else {
      logger.severe(" Could not start OpenMM minimization.")
    }

    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (forceFieldEnergy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(forceFieldEnergy)
    }
    return potentials
  }

}
