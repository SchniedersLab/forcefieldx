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

import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import ffx.potential.parsers.XYZFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FileUtils.copyFile
import static org.apache.commons.io.FileUtils.createParentDirectories
import static org.apache.commons.io.FilenameUtils.*

/**
 * The SaveAsPDB script saves a file as a PDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsPDB [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "SaveAsPDB")
class SaveAsPDB extends PotentialScript {

  @Mixin
  SaveOptions saveOptions

  /**
   * --fs or ---firstSnapshot Provide the number of the first snapshot to be written.
   */
  @Option(names = ['--fs', '--firstSnapshot'], paramLabel = "-1", defaultValue = "-1",
      description = 'First snapshot to write out (indexed from 0).')
  private int firstSnapshot

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
   * --wp or --writeProperties Copy the property file to each subdirectory.
   */
  @Option(names = ['--cp', '--copyProperties'], paramLabel = "true", defaultValue = "true",
      description = 'Copy the property file to numbered subdirectories (ignored if not writing to subdirectories).')
  private boolean copyProperties = true

  /**
   * --esv Handle an extended system at the bottom of XYZ files using XPHFilter.
   */
  @Option(names = ['--esv'], paramLabel = "file", defaultValue = "",
      description = 'PDB file to build extended system from.')
  private String extended = ""

  /**
   * The final argument is an XYZ or ARC coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in XYZ or ARC format.')
  private String filename = null

  /**
   * SaveAsPDB Constructor.
   */
  SaveAsPDB() {
    this(new groovy.lang.Binding())
  }

  /**
   * SaveAsPDB Constructor.
   * @param binding Groovy Binding to use.
   */
  SaveAsPDB(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveAsPDB run() {

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
    SystemFilter openFilter = potentialFunctions.getFilter()
    ExtendedSystem esvSystem = null

    if (openFilter instanceof XYZFilter && extended != "") {
      logger.info("Building extended system from " + extended)
      activeAssembly = getActiveAssembly(extended) // Build from file with res info
      esvSystem = new ExtendedSystem(activeAssembly, 7.4, null)
      activeAssembly.setFile(new File(filename))
      openFilter = new XPHFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties(), esvSystem)
      openFilter.readFile()
      logger.info("Reading ESV lambdas from XPH file")
    }

    logger.info("\n Saving PDB for " + filename)

    // Use the current base directory, or update if necessary based on the given filename.
    String dirString = getBaseDirString(filename)

    String name = removeExtension(getName(filename)) + ".pdb"
    File saveFile = new File(dirString + name)

    if (firstSnapshot >= 0) {
      PDBFilter snapshotFilter = new PDBFilter(saveFile, activeAssembly,
          activeAssembly.getForceField(), activeAssembly.getProperties())
      openFilter.readNext(true)
      // Reset the filter to read the first snapshot.
      boolean resetPosition = true
      int counter = 0
      int snapshotCounter = 0
      logger.info(" Writing snapshots from " + firstSnapshot + " to " + lastSnapshot + " with increment " + snapshotIncrement)

      while (openFilter.readNext(resetPosition)) {
        // No more resets.
        resetPosition = false
        int offset = counter - firstSnapshot
        // Write out the snapshot if it is within the range and the increment is met.
        if (counter >= firstSnapshot && counter <= lastSnapshot && offset % snapshotIncrement == 0) {
          File snapshotFile
          if (writeToDirectories) {
            String subdirectory = concat(dirString, snapshotCounter.toString())
            snapshotFile = new File(concat(subdirectory, name))
            createParentDirectories(snapshotFile)
            if (copyProperties) {
              String propertyFile = activeAssembly.getProperties().getString("propertyFile")
              if (propertyFile != null) {
                String propertyFilename = getName(propertyFile)
                File copyOfPropFile = new File(concat(subdirectory, propertyFilename))
                logger.info("\n Copying properties to " + copyOfPropFile.toString())
                copyFile(new File(propertyFile), copyOfPropFile)
              }
            }
          } else {
            snapshotFile = new File(concat(dirString,
                removeExtension(name) + "." + counter.toString() + ".pdb"))
          }
          potentialFunctions.versionFile(snapshotFile)
          saveOptions.preSaveOperations(activeAssembly)
          logger.info(" Saving PDB to         " + snapshotFile.toString())
          snapshotFilter.writeFile(snapshotFile, true, false, false)
          snapshotFile.append("END\n")

          // Increment the snapshot counter used to name the file or create a subdirectory.
          snapshotCounter++
        }
        // Increment the counter used to iterate through the snapshots.
        counter++
      }
      return this
    }

    saveFile = potentialFunctions.versionFile(saveFile)

    int numModels = openFilter.countNumModels()

    if (numModels == 1) {
      // Write one model and return.
      saveOptions.preSaveOperations(activeAssembly)
      potentialFunctions.saveAsPDB(activeAssembly, saveFile)
      return this
    }

    // Write out the first model as "MODEL 1".
    saveFile.write("MODEL        1\n")
    saveOptions.preSaveOperations(activeAssembly)
    potentialFunctions.saveAsPDB(activeAssembly, saveFile, false, true)
    saveFile.append("ENDMDL\n")
    PDBFilter saveFilter = (PDBFilter) potentialFunctions.getFilter()
    saveFilter.setModelNumbering(1)

    // Iterate through the rest of the models in an arc or pdb.
    if (openFilter != null
        && (openFilter instanceof XYZFilter || openFilter instanceof PDBFilter || openFilter instanceof XPHFilter)) {
      try {
        while (openFilter.readNext(false)) {
          if (extended) {
            for (Atom atom : activeAssembly.getAtomList()) {
              int atomIndex = atom.getIndex() - 1
              atom.setOccupancy(esvSystem.getTitrationLambda(atomIndex))
              atom.setTempFactor(esvSystem.getTautomerLambda(atomIndex))
            }
          }

          saveOptions.preSaveOperations(activeAssembly)
          saveFilter.writeFile(saveFile, true, true, false)
        }
      } catch (Exception e) {
        // Do nothing.
      }
      // Add a final "END" record.
      saveFile.append("END\n")
    }
    return this
  }


}
