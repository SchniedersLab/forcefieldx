//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.BioType
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The Biotype script prints out biotype properties.
 * <br>
 * Usage:
 * <br>
 * ffxc Biotype &lt;filename&gt;
 */
@Command(description = " Print out Biotype records for the atoms in an XYZ file.", name = "ffxc Biotype")
class Biotype extends PotentialScript {

  /**
   * The final argument is a single filename in XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = "An XYZ coordinate file.")
  String filename = null

  // Create a List of bioptype String entries.
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

    Atom[] atoms = activeAssembly.getAtomArray()
    String mol = atoms[0].getAtomType().environment

    mol = mol.replaceAll("\"", "").trim()
    if (mol.length() > 3) {
      mol = mol.substring(0, 3)
    }

    // Create a List of bioptype String entries.
    bioTypes = new ArrayList<>()

    int index = 1
    for (Atom atom : atoms) {
      List bonds = atom.getBonds()
      String[] bondString = null
      if (bonds != null) {
        bondString = new String[bonds.size()]
        int i = 0
        for (Bond bond : bonds) {
          bondString[i++] = bond.get1_2(atom).getName()
        }
      }

      BioType biotype = new BioType(index++, atom.getName(), mol, atom.getAtomType().type,
          bondString)

      bioTypes.add(biotype)
      logger.info(biotype.toString())
    }

    // Return the bioTypes via the Binding.
    binding.setVariable("bioTypes", bioTypes)

    return this
  }

}
