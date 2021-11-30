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
package ffx.potential.groovy

import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The MoveIntoUnitCell script moves the center of mass of each molecule into the unit cell.
 * <br>
 * Usage:
 * <br>
 * ffxc MoveIntoUnitCell &lt;filename&gt;
 */
@Command(description = " Move all molecules into the unit cell.", name = "ffxc MoveIntoUnitCell")
class MoveIntoUnitCell extends PotentialScript {

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  MolecularAssembly[] molecularAssemblies = null
  public double[][] origCoordinates = null
  public double[][] unitCellCoordinates = null

  /**
   * MoveIntoUnitCell Constructor.
   */
  MoveIntoUnitCell() {
    this(new Binding())
  }

  /**
   * MoveIntoUnitCell Constructor.
   * @param binding Groovy Binding to use.
   */
  MoveIntoUnitCell(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  MoveIntoUnitCell run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load one or more MolecularAssembly instances.
    molecularAssemblies = getActiveAssemblies(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Moving molecular centers of mass into the unit cell for " + filename + "\n")

    // Loop over each system.
    for (int i = 0; i < molecularAssemblies.length; i++) {
      MolecularAssembly molecularAssembly = molecularAssemblies[i]

      Atom[] atoms = molecularAssembly.getAtomArray()
      int nAtoms = atoms.length
      origCoordinates = new double[nAtoms][3]
      unitCellCoordinates = new double[nAtoms][3]

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index]
        origCoordinates[index][0] = atom.getX()
        origCoordinates[index][1] = atom.getY()
        origCoordinates[index][2] = atom.getZ()
      }
      molecularAssembly.moveAllIntoUnitCell()

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index]
        unitCellCoordinates[index][0] = atom.getX()
        unitCellCoordinates[index][1] = atom.getY()
        unitCellCoordinates[index][2] = atom.getZ()
      }
    }

    // Configure the base directory if it has not been set.
    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = saveDir.toString() + File.separator
    String name = FilenameUtils.getName(filename)
    String ext = FilenameUtils.getExtension(name)
    name = FilenameUtils.removeExtension(name)

    if (ext.toUpperCase().contains("XYZ")) {
      potentialFunctions.saveAsXYZ(molecularAssemblies[0], new File(dirName + name + ".xyz"))
    } else {
      potentialFunctions.saveAsPDB(molecularAssemblies, new File(dirName + name + ".pdb"))
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (molecularAssemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(molecularAssemblies).
          filter {a -> a != null
          }.map {a -> a.getPotentialEnergy()
      }.filter {e -> e != null
      }.collect(Collectors.toList())
    }
  }
}
