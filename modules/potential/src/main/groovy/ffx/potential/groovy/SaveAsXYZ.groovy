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
package ffx.potential.groovy

import ffx.crystal.Crystal
import ffx.crystal.SymOp
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.logging.log4j.core.util.FileUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static ffx.crystal.SymOp.applyCartesianSymOp
import static org.apache.commons.io.FilenameUtils.concat
import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The SaveAsXYZ script saves a file as an XYZ file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsXYZ [options] &lt;filename&gt;
 */
@Command(description = " Save the system as an XYZ file.", name = "SaveAsXYZ")
class SaveAsXYZ extends PotentialScript {

  @Mixin
  SaveOptions saveOptions

  /**
   * -p or --pos-offset to set the positive atom type offset
   */
  @Option(names = ['-p', '--pos-offset'], paramLabel = "0",
          description = 'Positive offset of atom types in the new file')
  int posOffset = 0

  /**
   * -n or --neg-offset to set the negative atom type offset
   */
  @Option(names = ['-n', '--neg-offset'], paramLabel = "0",
          description = 'Negative offset of atom types in the new file.')
  int negOffset = 0

  /**
   * -r or --random to apply a random Cartesian symmetry operator with the specified translation range -X .. X (no default).
   */
  @Option(names = ['-r', '--random'], paramLabel = "X",
          description = 'Apply a random Cartesian SymOp with translation range -X .. X.')
  double scalar = -1

  /**
   * --fs or --firstSnapshot Provide the number of the first snapshot to be written.
   */
  @Option(names = ['--fs', '--firstSnapshot'], paramLabel = "-1", defaultValue = "-1",
          description = 'First snapshot to write out (indexed from 0).')
  private int firstSnapshot = 0

  /**
   * --ls or ---lastSnapshot Provide the number of the last snapshot to be written.
   */
  @Option(names = ['--ls', '--lastSnapshot'], paramLabel = "-1", defaultValue = "-1",
          description = 'Last snapshot to write out (indexed from 0).')
  private int lastSnapshot

  /**
   * --si or --snapshotIncrement Provide the number of the snapshot increment.
   */
  @Option(names = ['--si', '--snapshotIncrement'], paramLabel = "1", defaultValue = "1",
          description = 'Increment between written snapshots.')
  private int snapshotIncrement

  /**
   * --wd or --writeToDirectories Provide the number of the snapshot increment.
   */
  @Option(names = ['--wd', '--writeToDirectories'], paramLabel = "false", defaultValue = "false",
          description = 'Write snapshots to numbered subdirectories.')
  private boolean writeToDirectories = false

  /**
   * --alt or --alternateLocation Choose an alternate location for a PDB file.
   */
  @Option(names = ['--alt', '--alternateLocation'], paramLabel = "A", defaultValue = "A",
          description = 'Choose an alternate location for the PDB file (not supported for PDBs with multiple models.')
  private Character alternateLocation = 'A'

  /**
   * The final argument is a PDB coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
          description = 'The atomic coordinate file PDB format.')
  private String filename = null

  /**
   * SaveAsXYZ Constructor.
   */
  SaveAsXYZ() {
    this(new Binding())
  }

  /**
   * SaveAsXYZ Constructor.
   * @param binding Groovy Binding to use.
   */
  SaveAsXYZ(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveAsXYZ run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    MolecularAssembly[] molecularAssemblies = getActiveAssemblies(filename)
    if (molecularAssemblies == null) {
      logger.info(helpString())
      return this
    }

    if (molecularAssemblies.length > 1) {
      Character currentAltLoc = activeAssembly.getAlternateLocation()
      logger.info("\n Current alternate location: " + currentAltLoc)
      if (currentAltLoc != alternateLocation) {
        for (MolecularAssembly molecularAssembly : molecularAssemblies) {
          Character altLoc = molecularAssembly.getAlternateLocation()
          if (altLoc == alternateLocation) {
            activeAssembly = molecularAssembly
            logger.info(" Switching to alternate location: " + altLoc)
            break
          }
        }
      }
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    SystemFilter openFilter = potentialFunctions.getFilter()
    int numModels = openFilter.countNumModels()

    int offset = 0
    // Positive offset atom types.
    if (posOffset > 0) {
      offset = posOffset
    }

    // Negative offset atom types.
    if (negOffset > 0) {
      offset = negOffset
      offset = -offset
    }

    // Offset atom type numbers.
    if (offset != 0) {
      logger.info("\n Offset atom types by " + offset)
      ForceField forceField = activeAssembly.getForceField()
      forceField.renumberForceField(0, offset, 0)
    }

    if (scalar > 0.0) {
      SymOp symOp = SymOp.randomSymOpFactory(scalar)
      logger.info(format("\n Applying random Cartesian SymOp\n: %s", symOp.toString()))
      Atom[] atoms = activeAssembly.getAtomArray()
      double[] xyz = new double[3]
      for (int i = 0; i < atoms.length; i++) {
        atoms[i].getXYZ(xyz)
        applyCartesianSymOp(xyz, xyz, symOp)
        atoms[i].setXYZ(xyz)
      }
    }

    String dirString = getBaseDirString(filename)
    String name = getName(filename)

    // Choose a single snapshot to write out from an archive.
    if (firstSnapshot >= 0) {
      XYZFilter snapshotFilter = new XYZFilter(new File(dirString + name),
              activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
      openFilter.readNext(true)
      int counter = 0
      int snapshotCounter = 0
      logger.info(" Writing snapshots from " + firstSnapshot + " to " + lastSnapshot + " with increment " + snapshotIncrement)
      boolean reset = true

      while (openFilter.readNext(reset)) {
//        counter++
        // No more resets
        reset = false
        int counterOffset = counter - firstSnapshot
        // Write out the snapshot if it is within the range and the increment is met.
        if (counter >= firstSnapshot && counter <= lastSnapshot && counterOffset % snapshotIncrement == 0) {
          File snapshotFile
          if (writeToDirectories) {
            String subdirectory = concat(dirString, snapshotCounter.toString())
            FileUtils.mkdir(new File(subdirectory), true)
            snapshotFile = new File(concat(subdirectory, name))
          } else {
            snapshotFile = new File(concat(dirString, removeExtension(name) + "." + counter.toString() + ".xyz"))
          }

          potentialFunctions.versionFile(snapshotFile)
          saveOptions.preSaveOperations(activeAssembly)
//          snapshotFilter.setModelNumbering(snapshotCounter) // todo this isn't a method in XYZFilter - don't need it?
          logger.info("\n Writing out XYZ for " + snapshotFile.toString())
          snapshotFilter.writeFile(snapshotFile, true)
          // Increment the snapshot counter used to name the file or create a subdirectory.
          snapshotCounter++
        }
        // Increment the counter used to iterate through the snapshots.
        counter++
      }
      return this
    }

    logger.info("\n Writing out XYZ for " + filename)

    if (numModels <= 1) {
      // Just save a single snapshot.
      name = removeExtension(name) + ".xyz"
      File saveFile = new File(dirString + name)
      saveOptions.preSaveOperations(activeAssembly)
      potentialFunctions.save(activeAssembly, saveFile)
    } else {
      // Save to an arc file rather than an xyz file if more than one model exists.
      name = removeExtension(name) + ".arc"
      File saveFile = new File(dirString + name)
      saveFile = potentialFunctions.versionFile(saveFile)
      saveOptions.preSaveOperations(activeAssembly)
      potentialFunctions.save(activeAssembly, saveFile)

      XYZFilter saveFilter = (XYZFilter) potentialFunctions.getFilter()
      while (openFilter.readNext(false)) {
        saveOptions.preSaveOperations(activeAssembly)
        saveFilter.writeFile(saveFile, true)
      }
    }

    return this
  }
}
