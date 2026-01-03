//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.apache.commons.io.FilenameUtils.getExtension;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * Move the center of mass of each molecule into the unit cell.
 *
 * Usage:
 *   ffxc MoveIntoUnitCell &lt;filename&gt;
 */
@Command(name = "MoveIntoUnitCell", description = " Move all molecules into the unit cell.")
public class MoveIntoUnitCell extends PotentialCommand {

  /** The final argument is a PDB or XYZ coordinate file. */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  /** Save a reference to the MolecularAssembly instances to destroy their potentials. */
  protected MolecularAssembly[] molecularAssemblies = null;

  /** Original Cartesian coordinates before moving into the unit cell. */
  public double[][] origCoordinates = null;
  /** Cartesian coordinates after moving into the unit cell. */
  public double[][] unitCellCoordinates = null;

  public MoveIntoUnitCell() { super(); }
  public MoveIntoUnitCell(FFXBinding binding) { super(binding); }
  public MoveIntoUnitCell(String[] args) { super(args); }

  @Override
  public MoveIntoUnitCell run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Load one or more MolecularAssembly instances.
    molecularAssemblies = getActiveAssemblies(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    logger.info("\n Moving molecular centers of mass into the unit cell for " + filename + "\n");

    // Loop over each system.
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      int nAtoms = atoms.length;
      origCoordinates = new double[nAtoms][3];
      unitCellCoordinates = new double[nAtoms][3];

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index];
        origCoordinates[index][0] = atom.getX();
        origCoordinates[index][1] = atom.getY();
        origCoordinates[index][2] = atom.getZ();
      }

      molecularAssembly.moveAllIntoUnitCell();

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index];
        unitCellCoordinates[index][0] = atom.getX();
        unitCellCoordinates[index][1] = atom.getY();
        unitCellCoordinates[index][2] = atom.getZ();
      }
    }

    // Save the updated coordinates.
    saveByOriginalExtension(molecularAssemblies, filename);
    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }

}
