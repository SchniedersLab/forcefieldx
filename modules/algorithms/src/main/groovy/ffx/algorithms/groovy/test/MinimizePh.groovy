//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy.test

import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.algorithms.optimize.PhMinimize
import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Minimize script uses a limited-memory BFGS algorithm to minimize the energy of a molecular system.
 * <br>
 * Usage:
 * <br>
 * ffxc Minimize [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Run L-BFGS minimization on a system.", name = "ffxc Minimize")
class MinimizePh extends AlgorithmsScript {

  @Mixin
  MinimizeOptions minimizeOptions

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * --pH or --constantPH Constant pH value for the test.
   */
  @CommandLine.Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
          description = 'pH value for the energy evaluation. (Only applies when esvTerm is true)')
  double pH = 7.4

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate files in PDB or XYZ format.')
  List<String> filename = null

  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  private int threadsPer = threadsAvail
  MolecularAssembly[] topologies
  private ForceFieldEnergy forceFieldEnergy

  /**
   * Minimize Constructor.
   */
    MinimizePh() {
    this(new Binding())
  }

  /**
   * Minimize Constructor.
   * @param binding The Groovy Binding to use.
   */
    MinimizePh(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  MinimizePh run() {

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

    logger.info("\n Running Energy on " + filename)
    forceFieldEnergy = activeAssembly.getPotentialEnergy()

    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, null)
    esvSystem.setConstantPh(pH)
    int numESVs = esvSystem.extendedResidueList.size()
    forceFieldEnergy.attachExtendedSystem(esvSystem)
    logger.info(format(" Attached extended system with %d residues.", numESVs))

    SystemFilter systemFilter = algorithmFunctions.getFilter()

    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x, true)

    PhMinimize minimize = new PhMinimize(topologies[0], algorithmListener, esvSystem)
    minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())

    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x, true)

    if (topologies.length > 1) {
      // Handle Multiple Topology Cases.
      for (molecularAssembly in topologies) {
        String modelFilename = molecularAssembly.getFile().getAbsolutePath()

        if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
          baseDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        String dirName = baseDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        String ext = FilenameUtils.getExtension(fileName)
        fileName = FilenameUtils.removeExtension(fileName)

        if (ext.toUpperCase().contains("XYZ")) {
          algorithmFunctions.saveAsXYZ(molecularAssembly, new File(dirName + fileName + ".xyz"))
        } else {
          algorithmFunctions.saveAsPDB(molecularAssembly, new File(dirName + fileName + ".pdb"))
        }
      }
    } else {
      // Handle Single Topology Cases.
      setActiveAssembly(topologies[0])
      String modelFilename = activeAssembly.getFile().getAbsolutePath()
      if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
        baseDir = new File(FilenameUtils.getFullPath(modelFilename))
      }

      String dirName = baseDir.toString() + File.separator
      String fileName = FilenameUtils.getName(modelFilename)
      String ext = FilenameUtils.getExtension(fileName)
      fileName = FilenameUtils.removeExtension(fileName)
      File saveFile
      SystemFilter writeFilter
      if (ext.toUpperCase().contains("XYZ")) {
        saveFile = new File(dirName + fileName + ".xyz")
        writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
            activeAssembly.getProperties())
        algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
      } else if (ext.toUpperCase().contains("ARC")) {
        saveFile = new File(dirName + fileName + ".arc")
        saveFile = algorithmFunctions.versionFile(saveFile)
        writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
            activeAssembly.getProperties())
        algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
      } else {
        saveFile = new File(dirName + fileName + ".pdb")
        saveFile = algorithmFunctions.versionFile(saveFile)
        writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
            activeAssembly.getProperties())
        int numModels = systemFilter.countNumModels()
        if (numModels > 1) {
          writeFilter.setModelNumbering(0)
        }
        writeFilter.writeFile(saveFile, true, false, false)
      }

      if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
        while (systemFilter.readNext()) {
          Crystal crystal = activeAssembly.getCrystal()
          ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
          forceFieldEnergy.setCrystal(crystal)
          if (systemFilter instanceof PDBFilter) {
            saveFile.append("ENDMDL\n")
            minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())
            PDBFilter pdbFilter = (PDBFilter) systemFilter
            pdbFilter.writeFile(saveFile, true, false, false)
          } else if (systemFilter instanceof XYZFilter) {
            minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())
            writeFilter.writeFile(saveFile, true)
          }
        }
        if (systemFilter instanceof PDBFilter) {
          saveFile.append("END\n")
        }
      }
    }

    return this
  }
}
