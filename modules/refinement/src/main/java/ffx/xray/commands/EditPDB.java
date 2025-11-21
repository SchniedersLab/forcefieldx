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
package ffx.xray.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Option;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.List;

import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * EditPDB edits a PDB model for use in experimental refinement.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.EditPDB &lt;pdbfile1&gt;
 */
@Command(description = " Edit a PDB model for use in refinement.", name = "xray.EditPDB")
public class EditPDB extends AlgorithmsCommand {

  /**
   * Replace deuterium atoms with hydrogen.
   */
  @Option(names = {"--rd", "--repaceDeuterium"}, paramLabel = "false", defaultValue = "false",
      description = "Replace deuterium with hydrogen.")
  private boolean replaceDeuterium = false;

  /**
   * --altLoc Save a specific alternate location. All atoms are assigned the root alternate
   * location ' '.
   */
  @Option(names = {"--altLoc", "--saveAltLoc"}, paramLabel = " ", defaultValue = " ",
      description = "Save a specific alternate location")
  private Character altLoc = ' ';

  /**
   * If a single alternate conformation is being saved, set occupancy values to 1.0.
   */
  @Option(names = {"-o", "--occupancy"}, paramLabel = "false", defaultValue = "false",
      description = "If a single alternate conformation is being saved, set occupancy values to 1.0.")
  private boolean occupancy = false;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "file", description = "PDB input file.")
  private String filename;

  private MolecularAssembly[] molecularAssemblies;

  /**
   * Deuterate constructor.
   */
  public EditPDB() {
    super();
  }

  /**
   * Deuterate constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public EditPDB(String[] args) {
    super(args);
  }

  /**
   * Deuterate constructor.
   *
   * @param binding The Binding to use.
   */
  public EditPDB(FFXBinding binding) {
    super(binding);
  }

  /**
   * Execute the script.
   */
  @Override
  public EditPDB run() {

    if (!init()) {
      return this;
    }

    if (filename != null) {
      molecularAssemblies = algorithmFunctions.openAll(filename);
      activeAssembly = molecularAssemblies[0];
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
      filename = activeAssembly.getFile().getAbsolutePath();
    }

    logger.info("\n Running xray.EditPDB on " + filename);

    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      // Replace deuterium with hydrogen.
      if (replaceDeuterium) {
        List<Atom> atoms = molecularAssembly.getAtomList();
        for (Atom atom : atoms) {
          if (atom.isDeuterium()) {
            String name = atom.getName().toUpperCase();
            atom.setName(name.replaceFirst("D", "H"));
          }
        }
        List<MSNode> water = molecularAssembly.getWater();
        for (MSNode node : water) {
          Molecule wat = (Molecule) node;
          wat.setName("HOH");
        }
      }

      if (altLoc == ' ' || molecularAssembly.getAlternateLocation().equals(altLoc)) {
        if (altLoc != ' ') {
          // We will save this alternate location with the alternate location set to ' '.
          molecularAssembly.setAtomicAltLoc(' ');
          if (occupancy) {
            // Set all residues to have an occupancy of 1.0.
            molecularAssembly.setOccupancy(1.0);
          }
          algorithmFunctions.saveAsPDB(molecularAssembly, new File(filename));
        }
      }
    }

    // Only save all alternate locations if the altLoc flag is not set.
    if (altLoc == ' ') {
      algorithmFunctions.saveAsPDB(molecularAssemblies, new File(filename));
    }
    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }
}
