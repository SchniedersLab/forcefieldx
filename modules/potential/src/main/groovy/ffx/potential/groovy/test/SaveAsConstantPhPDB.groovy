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
package ffx.potential.groovy.test

import ffx.potential.MolecularAssembly
import ffx.potential.parameters.TitrationUtils
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.ForceFieldFilter
import ffx.potential.parsers.PDBFilter
import ffx.utilities.Keyword
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Gradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Gradient [options] &lt;filename&gt;
 */
@Command(description = " Test the potential energy gradient.", name = "ffxc test.Gradient")
class SaveAsConstantPhPDB extends PotentialScript {

  @CommandLine.Option(names = ["--rt", "--rotamerTitration"], paramLabel = "false",
          description = "Prepare PDB for rotamer optimization with titration states.")
  boolean rotamerTitration = false
  /**
   * The final argument should be a PDB coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB format.')
  String filename = null

  /**
   * SaveAsConstantPhPDB constructor.
   */
  SaveAsConstantPhPDB() {
    this(new Binding())
  }

  /**
   * SaveAsConstantPhPDB constructor.
   * @param binding The Groovy Binding to use.
   */
  SaveAsConstantPhPDB(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveAsConstantPhPDB run() {

    if (!init()) {
      return this
    }

    if (rotamerTitration){
      logger.info("\n Adding rotamer optimization with titration protons to : " + filename + "\n")
    } else {
      logger.info("\n Adding constant pH protons to: " + filename + "\n")
    }


    // Read in command line.
    File structureFile = new File(filename)
    int index = filename.lastIndexOf(".")
    String name = filename.substring(0, index)
    activeAssembly = new MolecularAssembly(name)
    activeAssembly.setFile(structureFile)

    CompositeConfiguration properties = Keyword.loadProperties(structureFile)
    ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties)
    ForceField forceField = forceFieldFilter.parse()
    activeAssembly.setForceField(forceField)

    PDBFilter pdbFilter = new PDBFilter(structureFile, activeAssembly, forceField, properties)
    if (rotamerTitration){
      pdbFilter.setRotamerTitration(true)
    } else {
      pdbFilter.setConstantPH(true)
    }

    pdbFilter.readFile()
    pdbFilter.applyAtomProperties()
    activeAssembly.finalize(true, forceField)

    // Configure the base directory if it has not been set.
    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = saveDir.toString() + File.separator
    String fileName = FilenameUtils.getName(filename)
    File modelFile
    if (rotamerTitration) {
      fileName = FilenameUtils.removeExtension(fileName) + ".pdb_1"
      modelFile = new File(dirName + fileName)
    } else {
      fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
      modelFile = new File(dirName + fileName)
    }


    if (!pdbFilter.writeFile(modelFile, false, false, true)) {
      logger.info(format(" Save failed for %s", activeAssembly))
    }

    TitrationUtils constantPhUtils = new TitrationUtils(forceField)

    return this
  }

}
