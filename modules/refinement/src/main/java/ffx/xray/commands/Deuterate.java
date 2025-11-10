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
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.List;

import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * Deuterate changes exchangeable hydrogen atoms to deuterium atoms for a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Deuterate &lt;pdbfile1&gt;
 */
@Command(description = " Deuterate exchangable hydrogen of the PDB model.", name = "xray.Deuterate")
public class Deuterate extends AlgorithmsCommand {

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB input file.")
  private List<String> filenames;

  private MolecularAssembly[] molecularAssemblies;

  /**
   * Deuterate constructor.
   */
  public Deuterate() {
    super();
  }

  /**
   * Deuterate constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Deuterate(String[] args) {
    super(args);
  }

  /**
   * Deuterate constructor.
   *
   * @param binding The Binding to use.
   */
  public Deuterate(FFXBinding binding) {
    super(binding);
  }

  /**
   * Execute the script.
   */
  @Override
  public Deuterate run() {

    if (!init()) {
      return this;
    }

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0));
      activeAssembly = molecularAssemblies[0];
      filename = filenames.get(0);
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
      filename = activeAssembly.getFile().getAbsolutePath();
    }

    logger.info("\n Running xray.Deuterate on " + filename);

    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      for (Atom a : atoms) {
        if (a.isHydrogen()) {
          Atom b = a.getBonds().getFirst().get1_2(a);

          // Criteria for converting H to D
          if (b.getAtomicNumber() == 7 || b.getAtomicNumber() == 8) {
            String name = a.getName().replaceFirst("H", "D");
            a.setName(name);
          }
        }
      }

      List<MSNode> water = molecularAssembly.getWater();
      for (MSNode node : water) {
        Molecule wat = (Molecule) node;
        wat.setName("DOD");
      }
    }

    algorithmFunctions.saveAsPDB(molecularAssemblies, new File(removeExtension(filename) + "_deuterate.pdb"));

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }
}
