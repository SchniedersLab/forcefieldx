//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential.commands;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.cli.PotentialScript;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import groovy.lang.Binding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.List;

import static java.lang.String.format;

/**
 * WriteRestraints logs "restrain-position" properties for a PDB/XYZ file.
 *
 * Usage:
 *   ffxc WriteRestraints [options] &lt;filename&gt;
 */
@Command(name = "WriteRestraints", description = " Log position restraints for a PDB file.")
public class WriteRestraints extends PotentialScript {

  @Option(names = {"--chain", "-c"}, description = "Single character chain name.")
  private String chain = null;

  @Option(names = {"-k", "--forceConstant"}, defaultValue = "100.0", paramLabel = "100.0",
      description = "The force constant (kcal/mol/A^2).")
  private double forceConstant = 100.0;

  @Option(names = {"-d", "--flatBottom"}, defaultValue = "0.0", paramLabel = "0.0",
      description = "The flat bottom distance in Angstroms.")
  private double fbDistance = 0.0;

  @Option(names = {"--eh", "--excludeHydrogen"}, defaultValue = "false", paramLabel = "false",
      description = "Exclude writing restraints for hydrogen atoms.")
  private boolean excludeHydrogen = false;

  @Option(names = {"--ca", "--cAlphas"}, defaultValue = "false", paramLabel = "false",
      description = "Only write restraints for alpha carbons and/or phosphorus's.")
  private boolean onlyCalphas = false;

  @Option(names = {"-s", "--select"}, defaultValue = "1", paramLabel = "1",
      description = "Select every ith matching restraint and ignore the rest.")
  private int select = 1;

  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in XYZ or PDB format.")
  private String filename = null;

  public WriteRestraints() { super(); }
  public WriteRestraints(Binding binding) { super(binding); }
  public WriteRestraints(String[] args) { super(args); }

  @Override
  public WriteRestraints run() {
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    logger.info("\n Writing restraints for " + filename + "\n");

    SystemFilter systemFilter = this.potentialFunctions.getFilter();
    int count = 0;

    if (systemFilter instanceof PDBFilter) {
      Polymer[] polymers = activeAssembly.getChains();
      for (Polymer polymer : polymers) {
        if (chain != null && !chain.isEmpty()) {
          char requested = Character.toUpperCase(chain.charAt(0));
          Character chainID = polymer.getChainID();
          char current = (chainID == null) ? '?' : Character.toUpperCase(chainID);
          if (current != requested) {
            logger.info(" Skipping chain " + current);
            continue;
          } else {
            logger.info(" Restraints for chain " + current);
          }
        }
        java.util.List<Residue> residues = polymer.getResidues();
        for (Residue residue : residues) {
          if (onlyCalphas) {
            // Check for an amino acid C-alpha or a nucleic acid phosphate.
            Atom atom = residue.getAtomByName("CA", true);
            if (atom == null) {
              atom = residue.getAtomByName("P", true);
            }
            if (atom == null) {
              continue;
            }
            if (count % select == 0) {
              writeRestraints(atom);
            }
            count++;
          } else {
            List<Atom> atoms = residue.getAtomList();
            for (Atom atom : atoms) {
              if (excludeHydrogen && atom.isHydrogen()) {
                continue;
              }
              if (count % select == 0) {
                writeRestraints(atom);
              }
              count++;
            }
          }
        }
      }
    } else {
      Atom[] atoms = activeAssembly.getAtomArray();
      for (Atom atom : atoms) {
        if (excludeHydrogen && atom.isHydrogen()) {
          continue;
        }
        if (count % select == 0) {
          writeRestraints(atom);
        }
        count++;
      }
    }

    return this;
  }

  private void writeRestraints(Atom atom) {
    double x = atom.getX();
    double y = atom.getY();
    double z = atom.getZ();
    logger.info(format("restrain-position %4d %19.15f %19.15f %19.15f %12.8f %12.8f",
        atom.getIndex(), x, y, z, forceConstant, fbDistance));
  }
}
