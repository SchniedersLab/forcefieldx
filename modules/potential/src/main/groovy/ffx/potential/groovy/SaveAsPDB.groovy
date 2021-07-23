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
package ffx.potential.groovy

import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The SaveAsPDB script saves a file as a PDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsPDB [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "ffxc SaveAsPDB")
class SaveAsPDB extends PotentialScript {

  @Mixin
  SaveOptions saveOptions

  /**
   * --wS or --writeSnapshot Write out a specific snapshot. Provide the number of the snapshot to be written.
   */
  @Option(names = ['--wS', '--writeSnapshot'], paramLabel = "0", defaultValue = "0",
      description = 'Write out a specific snapshot.')
  private int writeSnapshot = 0

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
    this(new Binding())
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

    logger.info("\n Saving PDB for " + filename)

    // Configure the base directory if it has not been set.
    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = saveDir.toString() + File.separator
    String fileName = FilenameUtils.getName(filename)
    fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
    File modelFile = new File(dirName + fileName)

    if (writeSnapshot >= 1) {
      PDBFilter snapshotFilter = new PDBFilter(modelFile, activeAssembly,
          activeAssembly.getForceField(), activeAssembly.getProperties())
      openFilter.readNext(true)
      int counter = 1
      if (counter == writeSnapshot) {
        File snapshotFile = new File(
            dirName + File.separator + "snapshot" + counter.toString() + ".pdb")
        potentialFunctions.versionFile(snapshotFile)
        saveOptions.preSaveOperations(activeAssembly)
        snapshotFilter.setModelNumbering(writeSnapshot - 1)
        logger.info("\n Writing out PDB for " + snapshotFile.toString())
        snapshotFilter.writeFile(snapshotFile, true)
      }
      while (openFilter.readNext(false)) {
        counter++
        if (counter == writeSnapshot) {
          File snapshotFile = new File(
              dirName + File.separator + "snapshot" + counter.toString() + ".pdb")
          potentialFunctions.versionFile(snapshotFile)
          saveOptions.preSaveOperations(activeAssembly)
          snapshotFilter.setModelNumbering(writeSnapshot - 1)
          logger.info("\n Writing out PDB for " + snapshotFile.toString())
          snapshotFilter.writeFile(snapshotFile, true, false, false)
          snapshotFile.append("END\n")
          break
        }
      }
      return this
    }

    File saveFile = potentialFunctions.versionFile(modelFile)

    int numModels = openFilter.countNumModels()
    if (numModels > 1) {
      new BufferedWriter(new FileWriter(saveFile)).withCloseable {bw ->
        bw.write("MODEL        1\n")
        bw.flush()
      }
      saveOptions.preSaveOperations(activeAssembly)
      potentialFunctions.saveAsPDB(activeAssembly, saveFile, false, true)
    } else {
      saveOptions.preSaveOperations(activeAssembly)
      potentialFunctions.saveAsPDB(activeAssembly, saveFile)
    }

    PDBFilter saveFilter
    try {
      saveFilter = (PDBFilter) potentialFunctions.getFilter()
    } catch (Throwable t) {
      logger.info(format(" The type of %s is %s", potentialFunctions.getFilter(),
          potentialFunctions.getFilter().class.getName()))
      throw t
    }

    //If SaveAsPDB is run on an arc file, iterate through the models in the arc file and save each as a pdb file.
    if (openFilter != null && (openFilter instanceof XYZFilter || openFilter instanceof PDBFilter)
        && numModels > 1) {
      saveFilter.setModelNumbering(1)
      try {
        while (openFilter.readNext(false)) {
          saveFile.append("ENDMDL\n")
          saveOptions.preSaveOperations(activeAssembly)
          saveFilter.writeFile(saveFile, true, true, false)
        }
      } catch (Exception e) {
        // Do nothing.
      }
      saveFile.append("END\n")
    }
    return this
  }
}
