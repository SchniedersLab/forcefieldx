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

import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The SaveAsP1 script expands a specified file to P1
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsP1 [options] &lt;filename&gt;
 */
@Command(description = " Expand the system to P1 and then save it.", name = "ffxc SaveAsP1")
class SaveAsP1 extends PotentialScript {

  @CommandLine.Mixin
  SaveOptions saveOptions

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  /**
   * SaveAsP1 Constructor.
   */
  SaveAsP1() {
    this(new Binding())
  }

  /**
   * SaveAsP1 Constructor.
   * @param binding Groovy Binding to use.
   */
  SaveAsP1(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveAsP1 run() {
    if (!init()) {
      return this
    }

    MolecularAssembly[] assemblies
    if (filenames != null && filenames.size() > 0) {
      assemblies = potentialFunctions.openAll(filenames.get(0))
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    String modelFilename = activeAssembly.getFile().getAbsolutePath()
    logger.info("\n Expanding to P1 for " + modelFilename)

    // Configure the base directory if it has not been set.
    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(modelFilename))
    }

    String fileName = FilenameUtils.getName(modelFilename)
    String dirName = saveDir.getAbsolutePath()
    File saveLocation = new File(dirName + File.separator + fileName)

    logger.info(" Saving P1 file to: " + saveLocation)

    saveOptions.preSaveOperations(activeAssembly)
    potentialFunctions.saveAsP1(activeAssembly, saveLocation)

    return this
  }
}
