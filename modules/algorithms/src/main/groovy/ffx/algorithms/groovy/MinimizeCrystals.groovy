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
import ffx.algorithms.optimize.CrystalMinimize
import ffx.algorithms.optimize.Minimize
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly.FractionalMode
import ffx.potential.XtalEnergy
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.math3.util.FastMath.abs

/**
 * The CrystalMin script uses a limited-memory BFGS algorithm to minimize the
 * energy of a crystal, including both coordinates and unit cell parameters.
 * <br>
 * Usage:
 * <br>
 * ffxc CrystalMin [options] &lt;filename&gt;
 */
@Command(description = " Minimize crystal unit cell parameters.", name = "ffxc CrystalMin")
class MinimizeCrystals extends AlgorithmsScript {

  @Mixin
  MinimizeOptions minimizeOptions

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * -c or --coords to cycle between lattice and coordinate optimization until both satisfy the convergence criteria.
   */
  @Option(names = ["-c", "--coords"], paramLabel = 'false', defaultValue = "false",
      description = 'Cycle between lattice and coordinate optimization until both satisfy the convergence criteria.')
  boolean coords
  /**
   * -f or --fractional to set the optimization to maintain fractional coordinates [ATOM/MOLECULE/OFF].
   */
  @Option(names = ["-f", "--fractional"], paramLabel = 'molecule', defaultValue = "MOLECULE",
      description = 'Maintain fractional coordinates during lattice optimization [OFF/MOLECULE/ATOM].')
  String fractional
  /**
   * -t or --tensor to print out partial derivatives of the energy with respect to unit cell parameters.
   */
  @Option(names = ["-t", "--tensor"], paramLabel = 'false', defaultValue = "false",
      description = 'Compute partial derivatives of the energy with respect to unit cell parameters.')
  boolean tensor

  /**
   * The final argument(s) should be an XYZ or PDB filename.
   */
  @Parameters(arity = "1", paramLabel = "file", description = 'Atomic coordinate files in PDB or XYZ format.')
  String filename = null

  private XtalEnergy xtalEnergy
  private CrystalMinimize crystalMinimize

  /**
   * CrystalMin Constructor.
   */
  MinimizeCrystals() {
    this(new Binding())
  }

  /**
   * CrystalMin Constructor.
   * @param binding The Groovy Binding to use.
   */
  MinimizeCrystals(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  MinimizeCrystals run() {

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

    // Set the filename
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Running CrystalMinimize on " + filename)
    logger.info("\n RMS gradient convergence criteria: " + minimizeOptions.eps)

    atomSelectionOptions.setActiveAtoms(activeAssembly)

    ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
    xtalEnergy = new XtalEnergy(forceFieldEnergy, activeAssembly)
    xtalEnergy.setFractionalCoordinateMode(FractionalMode.MOLECULE)

    SystemFilter systemFilter = algorithmFunctions.getFilter()

    // Apply fractional coordinate mode.
    if (fractional) {
      try {
        FractionalMode mode = FractionalMode.valueOf(fractional.toUpperCase())
        xtalEnergy.setFractionalCoordinateMode(mode)
      } catch (Exception e) {
        logger.info(" Unrecognized fractional coordinate mode: " + fractional)
        logger.info(" Fractional coordinate mode is set to MOLECULE.")
      }
    }

    runMinimize()

    if (tensor) {
      crystalMinimize.printTensor()
    }

    // Handle Single Topology Cases.
    if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
      baseDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = baseDir.toString() + File.separator
    String name = FilenameUtils.getName(filename)
    String ext = FilenameUtils.getExtension(name)
    name = FilenameUtils.removeExtension(name)
    File saveFile
    PDBFilter pdbFilter = null
    XYZFilter xyzFilter = null

    if (ext.toUpperCase().contains("XYZ")) {
      saveFile = new File(dirName + name + ".xyz")
      xyzFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
    } else if (ext.toUpperCase().contains("ARC")) {
      saveFile = new File(dirName + name + ".arc")
      saveFile = algorithmFunctions.versionFile(saveFile)
      xyzFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
    } else {
      saveFile = new File(dirName + name + ".pdb")
      saveFile = algorithmFunctions.versionFile(saveFile)
      pdbFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      int numModels = systemFilter.countNumModels()
      if (numModels > 1) {
        pdbFilter.setModelNumbering(0)
      }
      pdbFilter.writeFile(saveFile, true, false, false)
    }

    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
      while (systemFilter.readNext()) {
        if (systemFilter instanceof PDBFilter) {
          saveFile.append("ENDMDL\n")
          runMinimize()
          pdbFilter.writeFile(saveFile, true, false, false)
        } else if (systemFilter instanceof XYZFilter) {
          runMinimize()
          xyzFilter.writeFile(saveFile, true)
        }
      }

      if (systemFilter instanceof PDBFilter) {
        saveFile.append("END\n")
      }
    }

    return this
  }

  void runMinimize() {
    crystalMinimize = new CrystalMinimize(activeAssembly, xtalEnergy, algorithmListener)
    crystalMinimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
    double energy = crystalMinimize.getEnergy()

    double tolerance = 1.0e-10

    // Complete rounds of coordinate and lattice optimization.
    if (coords) {
      ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
      Minimize minimize = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      while (true) {
        // Complete a round of coordinate optimization.
        minimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
        double newEnergy = minimize.getEnergy()
        int status = minimize.getStatus()
        if (status != 0) {
          break
        }
        if (abs(newEnergy - energy) <= tolerance) {
          break
        }
        energy = newEnergy

        // Complete a round of lattice optimization.
        crystalMinimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
        newEnergy = crystalMinimize.getEnergy()
        status = crystalMinimize.getStatus()
        if (status != 0) {
          break
        }
        if (abs(newEnergy - energy) <= tolerance) {
          break
        }
        energy = newEnergy
      }
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (xtalEnergy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList((Potential) xtalEnergy)
    }
    return potentials
  }

}
