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

import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.PotentialsUtils
import picocli.CommandLine.Option
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.getExtension
import static org.apache.commons.io.FilenameUtils.removeExtension
import static org.apache.commons.io.FilenameUtils.getName

/**
 * The SaveAsP1 script expands a specified file to P1
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsP1 [options] &lt;filename&gt;
 */
@Command(description = " Expand the system to P1 and then save it.", name = "SaveAsP1")
class SaveAsP1 extends PotentialScript {

  @Mixin
  SaveOptions saveOptions

  /**
   * --lmn or --replicatesVector Number of unit cells in replicates crystal.
   */
  @Option(names = ['--lmn', '--replicatesVector'], paramLabel = "", defaultValue = "",
          description = "Dimension of replicates crystal (e.g., \"2,2,2\" for a 2 x 2 x 2).")
  private String lmn

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  private String filename = null


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

    logger.info("\n Expanding to P1 for " + filename)

    // Use the current base directory, or update if necessary based on the given filename.
    String dirString = getBaseDirString(filename)

    String name = getName(filename)
    String ext = getExtension(name)
    name = removeExtension(name)
    String[] tokens = lmn.trim().replaceAll(" +","").split(",")
    int numTokens = tokens.length
    int[] replicatesVector = new int[3]
    boolean noReplicate
    if(numTokens == 0 || numTokens == 1 && tokens[0].isEmpty()){
      noReplicate = true
    }else if(numTokens == 1){
      noReplicate = false
      replicatesVector[0] = Integer.parseInt(tokens[0])
      replicatesVector[1] =  Integer.parseInt(tokens[0])
      replicatesVector[2] = Integer.parseInt(tokens[0])
    }else if(numTokens == 3){
      noReplicate = false
      replicatesVector[0] = Integer.parseInt(tokens[0])
      replicatesVector[1] =  Integer.parseInt(tokens[1])
      replicatesVector[2] = Integer.parseInt(tokens[2])
    }else{
      logger.warning(" Replicates indices could not be parsed (should be a comma separated list of integers)." +
              " Saving as P1.")
      noReplicate = true;
    }
    if (ext.toUpperCase().contains("XYZ")) {
      File saveLocation = SystemFilter.version(new File(dirString + name + ".xyz"))
      logger.info(" Saving P1 file to: " + saveLocation)

      if (noReplicate) {
        potentialFunctions.saveAsXYZinP1(activeAssembly, saveLocation)
      } else {
        potentialFunctions.saveAsXYZasReplicates(activeAssembly, saveLocation, replicatesVector)
      }
    } else {
      File saveLocation = SystemFilter.version(new File(dirString + name + ".pdb"))
      logger.info(" Saving symmetry mates file to: " + saveLocation)
      if(noReplicate) {
        potentialFunctions.saveAsPDBinP1(activeAssembly, saveLocation)
      } else {
        if(potentialFunctions instanceof PotentialsUtils){
          potentialFunctions.saveAsPDBinP1(activeAssembly, saveLocation, replicatesVector) //Needs LMN support...
        }
      }
    }
    return this
  }
}
