//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import ffx.potential.bonded.*
import ffx.potential.bonded.RotamerLibrary.NucleicSugarPucker
import ffx.potential.cli.PotentialScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getExtension
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * WriteRestraints logs position restraints for a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc WriteRestraints [options] &lt;filename&gt;
 */
@Command(description = " Log position restraints for a PDB file.", name = "WriteRestraints")
class WriteRestraints extends PotentialScript {

  /**
   * -c or --chain to specify chain
   */
  @Option(names = ['--chain', '-c'], description = 'Single character chain name.')
  String chain = null

  /**
   * The final argument is an XYZ or PDB coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in XYZ or PDB format.')
  private String filename = null

  /**
   * SaveRotamers Constructor.
   */
  WriteRestraints() {
    this(new Binding())
  }

  /**
   * SaveRotamers Constructor.
   * @param binding Groovy Binding to use.
   */
  WriteRestraints(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  WriteRestraints run() {

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

    logger.info("\n Writing restraints for " + filename + "\n")


    Polymer[] polymers = activeAssembly.getChains()
    for (Polymer polymer : polymers) {
      if (chain != null && !chain.isEmpty()) {
        char requested = chain.charAt(0).toUpperCase()
        char current = polymer.getChainID().toUpperCase()
        if (current != requested) {
          logger.info(" Skipping chain " + current)
          continue
        } else {
          logger.info(" Restraints for chain " + current)
        }
      }
      Residue[] residues = polymer.getResidues()
      for (Residue residue : residues) {
        Atom CA = residue.getAtomByName("CA", true)
        if (CA != null) {
          double x = CA.getX()
          double y = CA.getY()
          double z = CA.getZ()
          logger.info(format("restrain-position %4d %18.15f %18.15f %18.15f", CA.getIndex(), x, y, z))
        }

        Atom P = residue.getAtomByName("P", true)
        if (P != null) {
          double x = P.getX()
          double y = P.getY()
          double z = P.getZ()
          logger.info(format("restrain-position %4d %18.15f %18.15f %18.15f", P.getIndex(), x, y, z))
        }
      }
    }

    return this
  }
}
