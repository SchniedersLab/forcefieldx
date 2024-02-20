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

import ffx.numerics.Potential
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.CIFFilter
import ffx.potential.parsers.SystemFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.getExtension;
import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The ImportCIF script converts a CIF file to PDB or XYZ file(s) including atom types.
 * TODO: Move CIF parsing into a parsers.CIFFilter class.
 * <br>
 * Usage:
 * <br>
 * ffxc ImportCIF &lt;filename.cif&gt; &lt;filename.xyz&gt;
 */
@Command(description = " Convert a CIF file to PDB/XYZ format.", name = "ImportCIF")
class ImportCIF extends PotentialScript {

  /**
   * --zp or --zPrime Manually specify Z' (only affects writing CIF files)."
   */
  @Option(names = ['--zp', '--zPrime'], paramLabel = "-1", defaultValue = "-1",
          description = "Specify Z' when writing a CIF file.")
  private int zPrime

  /**
   * --sg or --spaceGroupNumber Override the CIF space group.
   */
  @Option(names = ['--sg', '--spaceGroupNumber'], paramLabel = "-1", defaultValue = "-1",
      description = 'Override the CIF space group.')
  private int sgNum

  /**
   * --name or --spaceGroupName Override the CIF space group.
   */
  @Option(names = ['--name', '--spaceGroupName'], paramLabel = "", defaultValue = "",
      description = 'Override the CIF space group.')
  private String sgName

  /**
   * --bt or --bondTolerance Tolerance added to covalent radius to bond atoms.
   */
  @Option(names = ['--bt', '--bondTolerance'], paramLabel = "0.2", defaultValue = "0.2",
      description = 'Tolerance added to covalent radius to determine if atoms should be bonded.')
  private double bondTolerance

  /**
   * --fl or --fixLattice Override CIF parameters to satisfy lattice conditions.
   */
  @Option(names = ['--fl', '--fixLattice'], paramLabel = "false", defaultValue = "false",
          description = 'Override CIF parameters to satisfy lattice conditions (Otherwise error).')
  private boolean fixLattice

  /**
   * --sc or --saveCIF Save file as a basic CIF.
   */
  @Option(names = ['--sc', '--saveCIF'], paramLabel = "false", defaultValue = "false",
          description = 'Attempt to save file in CIF format (input(s) in XYZ format).')
  private boolean saveCIF

  /**
   * --ca or --cifAppend Append cif files
   */
  @Option(names = ['--ca', '--cifAppend'], paramLabel = "false", defaultValue = "false",
          description = 'Append structures.')
  private boolean cifAppend

  /**
   * The final argument(s) should be a CIF file and a PDB or XYZ file that has been parameterized.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = "A CIF file and a PDB or XYZ file (already parameterized) containing one of each molecule from the CIF.")
  List<String> filenames = null

  /**
   * Array of strings containing files created.
   */
  public String[] createdFiles

  /**
   * ImportCIF Constructor.
   */
  ImportCIF() {
    this(new Binding())
  }

  /**
   * ImportCIF Constructor.
   * @param binding Groovy Binding to use.
   */
  ImportCIF(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  ImportCIF run() {

    // Turn off CDK logging.
    System.setProperty("cdk.logging.level", "fatal")

    // Turn off non-bonded interactions.
    System.setProperty("vdwterm", "false")

    if (!init()) {
      return this
    }
    if (filenames != null) {
      int fileInputs = filenames.size()
      System.clearProperty("mpoleterm")
      File saveFile
      String dirString = getBaseDirString(filenames.get(0))
      String name = removeExtension(getName(filenames.get(0)))
      if (saveCIF) {
        potentialFunctions.openAll(filenames.toArray() as String[])
        SystemFilter systemFilter = potentialFunctions.getFilter() as SystemFilter
        do{
          activeAssembly = systemFilter.getActiveMolecularSystem()
          saveFile = new File(dirString + name + ".cif")
          CIFFilter cifFilter = new CIFFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties(), saveCIF)
          cifFilter.setBondTolerance(bondTolerance)
          cifFilter.setFixLattice(fixLattice)
          cifFilter.setSgName(sgName)
          cifFilter.setSgNum(sgNum)
          cifFilter.setZprime(zPrime)
          if (cifFilter.writeFile(saveFile, cifAppend)) {
            createdFiles = cifFilter.getCreatedFileNames()
          }else{
            logger.warning(" Assembly " + activeAssembly.getName() + " was not successful...")
          }
        }while(systemFilter.readNext())
      } else if (fileInputs == 2) {
        getActiveAssembly(filenames[1])
        String ext = getExtension(filenames[1])
        saveFile = new File(dirString + name + "." + ext)
        CIFFilter cifFilter = new CIFFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties(), saveCIF)
        cifFilter.setBondTolerance(bondTolerance)
        cifFilter.setFixLattice(fixLattice)
        cifFilter.setSgName(sgName)
        cifFilter.setSgNum(sgNum)
        cifFilter.setZprime(zPrime)
        if (cifFilter.readFile()) {
          createdFiles = cifFilter.getCreatedFileNames()
        } else {
          logger.info(" Error occurred during conversion.")
        }
      }else{
        logger.info(helpString())
        logger.info(" Expected 2 files as input to convert CIF file(s).")
        return this
      }
    } else {
      logger.info(helpString())
      logger.info(" Expected 1 or 2 file(s) as input to ImportCIF.")
      return this
    }
    return this
  }

  @Override
  List<Potential> getPotentials() {
    return new ArrayList<Potential>()
  }
}