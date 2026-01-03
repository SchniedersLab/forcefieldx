//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.commands;

import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.cli.SaveOptions;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XPHFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

import static java.nio.file.StandardOpenOption.APPEND;
import static java.nio.file.StandardOpenOption.CREATE;
import static org.apache.commons.io.FilenameUtils.concat;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * Save the system as a PDB file.
 *
 * Usage:
 *   ffxc SaveAsPDB [options] &lt;filename&gt;
 */
@Command(name = "SaveAsPDB", description = " Save the system as a PDB file.")
public class SaveAsPDB extends PotentialCommand {

  @Mixin
  private SaveOptions saveOptions = new SaveOptions();

  /** --fs or --firstSnapshot Provide the number of the first snapshot to be written. */
  @Option(names = {"--fs", "--firstSnapshot"}, paramLabel = "-1", defaultValue = "-1",
      description = "First snapshot to write out (indexed from 0).")
  private int firstSnapshot = -1;

  /** --ls or --lastSnapshot Provide the number of the last snapshot to be written. */
  @Option(names = {"--ls", "--lastSnapshot"}, paramLabel = "-1", defaultValue = "-1",
      description = "Last snapshot to write out (indexed from 0).")
  private int lastSnapshot = -1;

  /** --si or --snapshotIncrement Provide the number of the snapshot increment. */
  @Option(names = {"--si", "--snapshotIncrement"}, paramLabel = "1", defaultValue = "1",
      description = "Increment between written snapshots.")
  private int snapshotIncrement = 1;

  /** --wd or --writeToDirectories Provide the number of the snapshot increment. */
  @Option(names = {"--wd", "--writeToDirectories"}, paramLabel = "false", defaultValue = "false",
      description = "Write snapshots to numbered subdirectories.")
  private boolean writeToDirectories = false;

  /** --cp or --copyProperties Copy the property file to each subdirectory. */
  @Option(names = {"--cp", "--copyProperties"}, paramLabel = "true", defaultValue = "true",
      description = "Copy the property file to numbered subdirectories (ignored if not writing to subdirectories).")
  private boolean copyProperties = true;

  /** --esv Handle an extended system at the bottom of XYZ files using XPHFilter. */
  @Option(names = {"--esv"}, paramLabel = "file", defaultValue = "",
      description = "PDB file to build extended system from.")
  private String extended = "";

  /** The final argument is an XYZ or ARC coordinate file. */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in XYZ or ARC format.")
  private String filename = null;

  public SaveAsPDB() { super(); }
  public SaveAsPDB(FFXBinding binding) { super(binding); }
  public SaveAsPDB(String[] args) { super(args); }

  @Override
  public SaveAsPDB run() {
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();
    SystemFilter openFilter = potentialFunctions.getFilter();
    ExtendedSystem esvSystem = null;

    if (openFilter instanceof XYZFilter && extended != null && !extended.isEmpty()) {
      logger.info(" Building extended system from " + extended);
      // Build extended system from file with residue info.
      activeAssembly = getActiveAssembly(extended);
      esvSystem = new ExtendedSystem(activeAssembly, 7.4, null);
      // Restore original filename while retaining extended data.
      activeAssembly.setFile(new File(filename));
      openFilter = new XPHFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties(), esvSystem);
      openFilter.readFile();
      logger.info(" Reading ESV lambdas from XPH file");
    }

    logger.info("\n Saving PDB for " + filename);

    // Use the current base directory, or update if necessary based on the given filename.
    String dirString = getBaseDirString(filename);

    String name = removeExtension(getName(filename)) + ".pdb";
    File saveFile = new File(dirString + name);

    if (firstSnapshot >= 0) {
      // Write selected snapshots to separate files.
      PDBFilter snapshotFilter = new PDBFilter(saveFile, activeAssembly,
          activeAssembly.getForceField(), activeAssembly.getProperties());
      openFilter.readNext(true);
      boolean resetPosition = true;
      int counter = 0;
      int snapshotCounter = 0;
      logger.info(" Writing snapshots from " + firstSnapshot + " to " + lastSnapshot + " with increment " + snapshotIncrement);

      while (openFilter.readNext(resetPosition)) {
        // No more resets after first read.
        resetPosition = false;
        int offset = counter - firstSnapshot;
        if (counter >= firstSnapshot && (lastSnapshot < 0 || counter <= lastSnapshot) && offset % snapshotIncrement == 0) {
          File snapshotFile;
          if (writeToDirectories) {
            String subdirectory = concat(dirString, Integer.toString(snapshotCounter));
            snapshotFile = new File(concat(subdirectory, name));
            // Ensure directory exists and optionally copy properties (handled by caller environment typically).
            //noinspection ResultOfMethodCallIgnored
            new File(subdirectory).mkdirs();
            if (copyProperties) {
              String propertyFile = activeAssembly.getProperties().getString("propertyFile");
              if (propertyFile != null) {
                File copyOfPropFile = new File(concat(subdirectory, getName(propertyFile)));
                try {
                  Files.createDirectories(copyOfPropFile.getParentFile().toPath());
                  Files.copy(new File(propertyFile).toPath(), copyOfPropFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException e) {
                  logger.info(" Could not copy properties file: " + e.getMessage());
                }
              }
            }
          } else {
            snapshotFile = new File(concat(dirString,
                removeExtension(name) + "." + counter + ".pdb"));
          }
          potentialFunctions.versionFile(snapshotFile);
          saveOptions.preSaveOperations(activeAssembly);
          logger.info(" Saving PDB to         " + snapshotFile);
          snapshotFilter.writeFile(snapshotFile, true, false, false);
          try {
            Files.writeString(snapshotFile.toPath(), "END\n", APPEND, CREATE);
          } catch (IOException e) { /* ignore */ }
          snapshotCounter++;
        }
        counter++;
      }
      return this;
    }

    // Version the save file for multi-model writing when needed.
    saveFile = potentialFunctions.versionFile(saveFile);

    int numModels = openFilter.countNumModels();
    if (numModels == 1) {
      // Write one model and return.
      saveOptions.preSaveOperations(activeAssembly);
      potentialFunctions.saveAsPDB(activeAssembly, saveFile);
      return this;
    }

    // Write out the first model as "MODEL 1".
    try {
      Files.writeString(saveFile.toPath(), "MODEL        1\n");
    } catch (IOException e) { /* ignore */ }
    saveOptions.preSaveOperations(activeAssembly);
    potentialFunctions.saveAsPDB(activeAssembly, saveFile, false, true);
    try {
      Files.writeString(saveFile.toPath(), "ENDMDL\n", APPEND, CREATE);
    } catch (IOException e) { /* ignore */ }

    PDBFilter saveFilter = (PDBFilter) potentialFunctions.getFilter();
    saveFilter.setModelNumbering(1);

    // Iterate through the rest of the models in an arc or pdb.
    if (openFilter instanceof XYZFilter || openFilter instanceof PDBFilter || openFilter instanceof XPHFilter) {
      try {
        while (openFilter.readNext(false)) {
          if (esvSystem != null) {
            for (Atom atom : activeAssembly.getAtomList()) {
              int atomIndex = atom.getIndex() - 1;
              atom.setOccupancy(esvSystem.getTitrationLambda(atomIndex));
              atom.setTempFactor(esvSystem.getTautomerLambda(atomIndex));
            }
          }
          saveOptions.preSaveOperations(activeAssembly);
          saveFilter.writeFile(saveFile, true, true, false);
        }
      } catch (Exception e) {
        // Do nothing.
      }
      // Add a final "END" record.
      try {
        Files.writeString(saveFile.toPath(), "END\n", APPEND, CREATE);
      } catch (IOException e) { /* ignore */ }
    }
    return this;
  }
}
