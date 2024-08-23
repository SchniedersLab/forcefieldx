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
import ffx.potential.bonded.Bond
import ffx.potential.bonded.Molecule
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.BioType
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The Biotype script prints out biotype properties.
 * <br>
 * Usage:
 * <br>
 * ffxc Biotype &lt;filename&gt;
 */
@Command(description = " Print out Biotype records for the atoms in an XYZ file.", name = "Biotype")
class Biotype extends PotentialScript {

  /**
   * --name or --moleculeName The molecule name to use for the Biotype records.
   */
  @Option(names = ['--name', '--moleculeName'], paramLabel = "MOL", defaultValue = "MOL",
      description = 'The molecule name to use for the Biotype records.')
  private String molName

  /**
   * -a or --useAtomNames Use the atom names in the XYZ file.
   */
  @Option(names = ['-a', '--useAtomNames'], paramLabel = "false", defaultValue = "false",
      description = 'Use the atom names in the XYZ file.')
  private boolean useAtomNames

  /**
   * -w or --writePDB Write out a PDB file with the updated atom and molecule names.
   */
  @Option(names = ['-c', '--writeCONECT'], paramLabel = "false", defaultValue = "false",
      description = 'Write out CONECT records to append to the PDB file.')
  private boolean writeCONNECT

  /**
   * -w or --writePDB Write out a PDB file with the updated atom and molecule names.
   */
  @Option(names = ['-w', '--writePDB'], paramLabel = "false", defaultValue = "false",
      description = 'Write out a PDB file with the updated atom and molecule names.')
  private boolean writePDB

  /**
   * The final argument is a single filename in XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = "An XYZ coordinate file.")
  String filename = null


  /**
   * The Script Binding contains the following List of BioType instances upon completion.
   */
  private List<BioType> bioTypes = null

  /**
   * Biotype Constructor.
   */
  Biotype() {
    this(new Binding())
  }

  /**
   * Biotype Constructor.
   * @param binding Groovy Binding to use.
   */
  Biotype(Binding binding) {
    super(binding)
  }

  /**
   * Get the list of created Biotype records.
   * @return
   */
  List<BioType> getBioTypes() {
    return bioTypes
  }

  @Override
  Biotype run() {

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

    logger.info("\n Running Biotype on " + filename)

    Molecule[] molecules = activeAssembly.getMoleculeArray()
    if (molecules.length > 1) {
      logger.info(" Biotype is intended for a system with one molecule.")
      return this
    }

    // Update the molecule name.
    molecules[0].setName(molName)
    Atom[] atoms = activeAssembly.getAtomArray()

    // Create a List of biotype String entries.
    bioTypes = new ArrayList<>()

    // Update atom names.
    if (!useAtomNames) {
      Atom.ElementSymbol[] symbols = Atom.ElementSymbol.values()
      int[] elementCounts = new int[symbols.length]
      for (Atom atom : atoms) {
        int element = atom.getAtomicNumber()
        int n = elementCounts[element]
        String name = symbols[element - 1].name().toUpperCase() + n
        atom.setName(name)
        elementCounts[element]++
      }
    }

    int index = 1
    for (Atom atom : atoms) {
      // Update the molecule name.
      atom.setResName(molName)

      // Collect the bond names.
      List bonds = atom.getBonds()
      String[] bondString = null
      if (bonds != null) {
        bondString = new String[bonds.size()]
        int i = 0
        for (Bond bond : bonds) {
          bondString[i++] = bond.get1_2(atom).getName()
        }
      }

      // Create a Biotype entry.
      BioType biotype = new BioType(index++, atom.getName(), molName,
          atom.getAtomType().type, bondString)

      bioTypes.add(biotype)
      logger.info(biotype.toString())
    }

    filename = activeAssembly.getFile().getAbsolutePath()
    String dirString = getBaseDirString(filename)
    String name = getName(filename)
    name = removeExtension(name)

    // Save out a PDB file with updated atom names.
    if (writePDB) {
      File pdbFile = new File(dirString + name + ".pdb")
      logger.info("\n Saving PDB file: " + pdbFile)
      potentialFunctions.saveAsPDB(activeAssembly, pdbFile)
    }

    // Return the bioTypes via the Binding.
    binding.setVariable("bioTypes", bioTypes)

    if (writeCONNECT) {
      for (Atom a : atoms) {
// =============================================================================
//  7 - 11        Integer        serial       Atom  serial number
// 12 - 16        Integer        serial       Serial number of bonded atom
// 17 - 21        Integer        serial       Serial number of bonded atom
// 22 - 26        Integer        serial       Serial number of bonded atom
// 27 - 31        Integer        serial       Serial number of bonded atom
        StringBuilder sb = new StringBuilder(String.format("CONECT%5s", Integer.toString(a.getXyzIndex())))
        Bond[] bonds = a.getBonds()
        for (Bond b : bonds) {
          sb.append(String.format("%5s", Integer.toString(b.get1_2(a).getXyzIndex())))
        }
        logger.info(sb.toString())
      }
    }

    return this
  }

}
