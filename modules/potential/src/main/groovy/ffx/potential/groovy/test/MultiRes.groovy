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
package ffx.potential.groovy.test

import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The MultiResidue script evaluates the energy of a MultiResidue system.
 * <br>
 * Usage:
 * <br>
 * ffxc test.MultiResidue [options] &lt;filename&gt;
 */
@Command(description = " Evaluates the energy of a MultiResidue system.", name = "ffxc test.MultiResidue")
class MultiRes extends PotentialScript {

  /**
   * -r or --resID to define the residue number (default is 1).
   */
  @Option(names = ['-r', '--resID'], paramLabel = "1", description = 'Define the residue number.')
  Integer resID

  /**
   * -c or --chain to set the single character chain name (default is \' \').'
   */
  @Option(names = ['-c', '--chain'], paramLabel = ' ', defaultValue = ' ',
      description = 'Single character chain name (default is \' \').')
  Character chain

  /**
   * -n or --name to set the name of residue to switch (default is ALA).
   */
  @Option(names = ['-n', '--resname'], paramLabel = 'ALA', defaultValue = 'ALA',
      description = 'Name of residue to switch to.')
  String name

  /**
   * -a or --all Add all amino acids at the given position.
   */
  @Option(names = ['-a', '--all'], paramLabel = 'false', defaultValue = 'false',
      description = 'Add all amino acids at the given position.')
  boolean all

  /**
   * The final argument is a single filename in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "A PDB or XYZ coordinate file.")
  String filename = null

  /**
   * MultiResidue Constructor.
   */
  MultiRes() {
    this(new Binding())
  }

  /**
   * MultiResidue Constructor.
   * @param binding Groovy Binding to use.
   */
  MultiRes(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  MultiRes run() {

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

    if (name.toUpperCase() == 'PRO') {
      logger.info(" Proline is not supported with a MultiResidue.")
      return this
    }

    logger.info("\n Running MultiResidue on " + filename + "\n")

    ForceField forceField = activeAssembly.getForceField()
    ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()

    MultiResidue multiResidue = null
    Residue residue = null
    AminoAcid3 startingResidueAA3 = null
    Polymer[] polymers = activeAssembly.getChains()
    for (int i = 0; i < polymers.length; i++) {
      Polymer polymer = polymers[i]
      if (chain == polymer.getChainID()) {
        residue = polymer.getResidue(resID)
        if (residue != null) {
          startingResidueAA3 = AminoAcid3.valueOf(residue.getName())
          if (!startingResidueAA3.useWithMultiResidue) {
            logger.info(format(" %s is not supported with a MultiResidue.", startingResidueAA3))
            return this
          }
          multiResidue = new MultiResidue(residue, forceField)
          polymer.addMultiResidue(multiResidue)
        }
      }
    }

    if (residue == null) {
      logger.info(" Chain " + chain + " residue " + resID + " was not found.")
      return this
    }

    if (all) {
      for (AminoAcid3 aa3 : AminoAcid3.getEnumConstants()) {
        if (aa3 != startingResidueAA3 && aa3.useWithMultiResidue) {
          logger.info(" Adding Residue: " + aa3.toString())
          multiResidue.addResidue(new Residue(aa3.toString(), resID, Residue.ResidueType.AA))
        }
      }
    } else {
      logger.info(" Adding Residue: " + name.toUpperCase())
      multiResidue.addResidue(new Residue(name.toUpperCase(), resID, Residue.ResidueType.AA))
    }

    int numResidues = multiResidue.getResidueCount()
    for (int i = 0; i < numResidues; i++) {
      multiResidue.setActiveResidue(i)
      logger.info("\n Active Residue: " + multiResidue.toString())
      forceFieldEnergy.reInit()
      forceFieldEnergy.energy(true, true)
    }

    return this
  }
}
