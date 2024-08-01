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
package ffx.xray.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Molecule
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * Deuterate changes exchangeable hydrogen atoms to deuterium atoms for a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Deuterate &lt;pdbfile1&gt;
 */
@Command(description = " Deuterate exchangable hydrogen of the PDB model.", name = "xray.Deuterate")
class Deuterate extends AlgorithmsScript {

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB input file.")
  private List<String> filenames

  /**
   * Deuterate constructor.
   */
  Deuterate() {
    this(new Binding())
  }

  /**
   * Deuterate constructor.
   * @param binding The Groovy Binding to use.
   */
  Deuterate(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Deuterate run() {

    if (!init()) {
      return this
    }

    MolecularAssembly[] molecularAssemblies
    String filename
    if (filenames != null && filenames.size() > 0) {
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = molecularAssemblies[0]
      filename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      molecularAssemblies = [activeAssembly]
      filename = activeAssembly.getFile().getAbsolutePath()
    }

    logger.info("\n Running xray.Deuterate on " + filename)

    for (int i = 0; i < molecularAssemblies.length; i++) {
      Atom[] atoms = molecularAssemblies[i].getAtomArray()
      for (Atom a : atoms) {
        if (a.getAtomicNumber() == 1) {
          Atom b = a.getBonds().get(0).get1_2(a)

          // Criteria for converting H to D
          if (b.getAtomicNumber() == 7
              || b.getAtomicNumber() == 8) {
            String name = a.getName().replaceFirst("H", "D")
            a.setName(name)
          }
        }
      }

      List<MSNode> water = molecularAssemblies[i].getWater()
      for (MSNode node : water) {
        Molecule wat = (Molecule) node
        wat.setName("DOD")
      }
    }

    algorithmFunctions.saveAsPDB(molecularAssemblies, new File(removeExtension(filename) + "_deuterate.pdb"))

    return this
  }

  @Override
  List<Potential> getPotentials() {
    return Collections.emptyList()
  }
}
