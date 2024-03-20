//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFileFilter
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import ffx.potential.parsers.XYZFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension

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
   * --wS or --writeSnapshot Write out a specific snapshot. Provide the number of the snapshot to be written.
   */
  @Option(names = ['--wS', '--writeSnapshot'], paramLabel = "0", defaultValue = "0",
      description = 'Write out a specific snapshot.')
  private int writeSnapshot = 0

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
    ExtendedSystem esvSystem = null

    if(openFilter instanceof XYZFilter && extended != ""){
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

    if (writeSnapshot >= 1) {
      PDBFilter snapshotFilter = new PDBFilter(saveFile, activeAssembly,
          activeAssembly.getForceField(), activeAssembly.getProperties())
      openFilter.readNext(true)
      // Reset the filter to read the first snapshot.
      boolean resetPosition = true
      int counter = 0
      while (openFilter.readNext(resetPosition)) {
        // No more resets.
        resetPosition = false
        // Increment the snapshot counter.
        counter++
        if (counter == writeSnapshot) {
          File snapshotFile = new File(dirString + "snapshot" + counter.toString() + ".pdb")
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
    if (openFilter != null && (openFilter instanceof XYZFilter || openFilter instanceof PDBFilter || openFilter instanceof XPHFilter)) {
      try {
        while (openFilter.readNext(false)) {
          if(extended) {
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
