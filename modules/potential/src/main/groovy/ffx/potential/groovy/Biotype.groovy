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
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.cli.PotentialScript
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

  @Parameters(arity = "1..*", paramLabel = "files",
      description = "An XYZ file.")
  List<String> xyzFile = null

  /**
   * Result of this script is an array of Biotype strings.
   */
  public List<String> biotypes

  @Override
  Biotype run() {
    if (!init()) {
      return null
    }

    if (xyzFile != null && xyzFile.size() > 0) {
      MolecularAssembly[] assemblies = [potentialFunctions.open(xyzFile.get(0))]
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return null
    }

    logger.info("\n Running Biotype on " + activeAssembly.toString())

    potentialFunctions.energy(activeAssembly)

    List atoms = activeAssembly.getAtomList()
    String mol = atoms.get(0).getAtomType().environment

    mol = mol.replaceAll("\"", "").trim()
    if (mol.length() > 3) {
      mol = mol.substring(0, 3)
    }

    // Create a List of bioptype String entries.
    biotypes = new ArrayList<>()

    int index = 1
    for (Atom atom : atoms) {
      StringBuilder sb = new StringBuilder()
      sb.append(String.format(" biotype %3d %4s \"%s\" %3d", index++, atom.getName(), mol,
          atom.getAtomType().type))
      List bonds = atom.getBonds()
      if (bonds != null) {
        for (Bond bond : bonds) {
          sb.append(String.format(" %4s", bond.get1_2(atom).getName()))
        }
      }
      biotypes.add(sb.toString())
      logger.info(sb.toString())
    }

    return this
  }

}
