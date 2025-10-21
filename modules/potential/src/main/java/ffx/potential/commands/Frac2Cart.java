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

import ffx.crystal.Crystal;
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
 * Convert from fractional to Cartesian coordinates.
 *
 * Usage:
 *   ffxc Frac2Cart &lt;filename&gt;
 */
@Command(name = "Frac2Cart", description = " Convert from fractional to Cartesian coordinates.")
public class Frac2Cart extends PotentialCommand {

  /** The final argument should be a file in PDB or XYZ format. */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  /** Save a reference to the MolecularAssembly instances to destroy their potentials. */
  private MolecularAssembly[] molecularAssemblies;

  /** Cartesian coordinate output. */
  private double[][] cartCoordinates = null;
  /** Fractional coordinate input. */
  private double[][] fracCoordinates = null;

  public Frac2Cart() { super(); }
  public Frac2Cart(FFXBinding binding) { super(binding); }
  public Frac2Cart(String[] args) { super(args); }

  /** Return Cartesian Coordinate output. */
  public double[][] getCart() { return cartCoordinates; }
  /** Return Fractional Coordinate input. */
  public double[][] getFrac() { return fracCoordinates; }

  @Override
  public Frac2Cart run() {
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

    logger.info("\n Converting from fractional to Cartesian coordinates for " + filename);

    // Loop over each system. The original Groovy implementation only exposed the last system's
    // coordinates as 2D arrays; we preserve that API for test compatibility.
    for (MolecularAssembly ma : molecularAssemblies) {
      Crystal crystal = ma.getCrystal().getUnitCell();

      Atom[] atoms = ma.getAtomArray();
      int nAtoms = atoms.length;
      fracCoordinates = new double[nAtoms][3];
      cartCoordinates = new double[nAtoms][3];

      double[] frac = new double[3];
      double[] cart = new double[3];

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index];
        atom.getXYZ(frac);
        crystal.toCartesianCoordinates(frac, cart);

        // If the atom is at a special position, make sure it's active so the coordinates are updated.
        boolean active = atom.isActive();
        atom.setActive(true);
        atom.moveTo(cart);
        atom.setActive(active);

        cartCoordinates[index][0] = cart[0];
        cartCoordinates[index][1] = cart[1];
        cartCoordinates[index][2] = cart[2];

        fracCoordinates[index][0] = frac[0];
        fracCoordinates[index][1] = frac[1];
        fracCoordinates[index][2] = frac[2];
      }
    }

    saveByOriginalExtension(molecularAssemblies, filename);

    // Export results via Binding for compatibility (match Cart2Frac binding names).
    binding.setVariable("cart", cartCoordinates);
    binding.setVariable("frac", fracCoordinates);

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }
}
