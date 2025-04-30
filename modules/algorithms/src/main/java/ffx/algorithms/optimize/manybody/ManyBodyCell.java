// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.optimize.manybody;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.nonbonded.NeighborList;
import ffx.utilities.StringOutputStream;

import java.util.*;

import static ffx.crystal.SymOp.applyFracSymOp;
import static ffx.potential.bonded.RotamerLibrary.applyRotamer;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

/**
 * A cell used for optimization of a subdomain, its residues, its extent in fractional coordinates,
 * its overall (linear) index, and its indices along the a, b, and c axes.
 */
public class ManyBodyCell {

  /**
   * Cell extent in fractional coordinates.
   * fracCoords[0-2] contains min a, b and c
   * fracCoords[3-5] contains max a, b and c
   */
  private final double[] fracCoords = new double[6];
  /**
   * The a, b, and c indices of this cell.
   */
  private final int[] indices = new int[3];
  /**
   * The index of this cell.
   */
  private final int cellIndex;
  /**
   * The Residues contained within this cell.
   */
  private final ArrayList<Residue> residues = new ArrayList<>();

  /**
   * Constructs a ManyBodyCell instance, which takes up a set of fractional coordinates within the
   * Crystal, the Residues contained within, and the index of the cell along the crystal's a, b, and
   * c axes.
   *
   * @param fractionalCoordinates Extent in fractional coordinates (0-2 are min a,b,c and 3-5 are max).
   * @param indices               Index of the cell along a, b, and c axes.
   * @param cellIndex             Index of the cell in linear array.
   */
  public ManyBodyCell(double[] fractionalCoordinates, int[] indices, int cellIndex) {
    arraycopy(fractionalCoordinates, 0, fracCoords, 0, fracCoords.length);
    arraycopy(indices, 0, this.indices, 0, this.indices.length);
    this.cellIndex = cellIndex;
  }

  /**
   * Returns a string representation of this BoxOptCell.
   *
   * @return String representation.
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append(format(" Optimization Cell %2d (%2d, %2d, %2d)\n", cellIndex + 1, indices[0] + 1, indices[1] + 1, indices[2] + 1));
    int n = residues.size();
    if (n == 0) {
      sb.append("  Cell is Empty\n");
    } else if (n == 1) {
      sb.append(format("  Single Residue: %s\n", residues.getFirst()));
    } else {
      Residue firstResidue = residues.getFirst();
      Residue lastResidue = residues.getLast();
      sb.append(format("  %d Residues: %s ... %s\n", n, firstResidue.toString(), lastResidue.toString()));
    }
    sb.append("  Fractional Coordinates:\n");
    sb.append(format("   A-axis %5.3f to %5.3f\n", fracCoords[0], fracCoords[3]));
    sb.append(format("   B-axis %5.3f to %5.3f\n", fracCoords[1], fracCoords[4]));
    sb.append(format("   C-axis %5.3f to %5.3f", fracCoords[2], fracCoords[5]));
    return sb.toString();
  }

  /**
   * Add a residue to the box.
   *
   * @param residue Residue to be added.
   */
  public void addResidue(Residue residue) {
    residues.add(residue);
  }

  /**
   * Checks if any rotamer of a Residue is inside this BoxOptCell.
   *
   * @param residue      Residue to check.
   * @param crystal      A Crystal.
   * @param symOp        A symmetry operator to apply.
   * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
   * @return If contained inside this BoxOptCell.
   */
  public boolean anyRotamerInsideCell(Residue residue, Crystal crystal, SymOp symOp, boolean variableOnly) {
    ResidueState incomingState = residue.storeState();
    Rotamer[] rotamers = residue.getRotamers();
    boolean inside = Arrays.stream(rotamers).anyMatch((Rotamer r) -> {
      applyRotamer(residue, r);
      return residueInsideCell(residue, crystal, symOp, variableOnly);
    });
    residue.revertState(incomingState);
    return inside;
  }

  /**
   * Checks if an Atom would be contained inside this cell.
   *
   * @param atom    Atom to check.
   * @param crystal A Crystal.
   * @param symOp   A symmetry operator to apply.
   * @return If contained.
   */
  public boolean atomInsideCell(Atom atom, Crystal crystal, SymOp symOp) {
    double[] atXYZ = new double[3];
    atXYZ = atom.getXYZ(atXYZ);
    crystal.toFractionalCoordinates(atXYZ, atXYZ);
    applyFracSymOp(atXYZ, atXYZ, symOp);
    moveValuesBetweenZeroAndOne(atXYZ);
    return checkIfContained(atXYZ);
  }

  /**
   * Moves an array of doubles to be within 0.0 and 1.0 by addition or subtraction of a multiple of
   * 1.0. Typical use is moving an atom placed outside crystal boundaries from the symmetry mate back
   * into the crystal.
   *
   * @param valuesToMove Doubles to be moved between 0 and 1.
   */
  public static void moveValuesBetweenZeroAndOne(double[] valuesToMove) {
    for (int i = 0; i < valuesToMove.length; i++) {
      valuesToMove[i] = moveBetweenZeroAndOne(valuesToMove[i]);
    }
  }

  /**
   * Moves a double to be within 0.0 and 1.0 by addition or subtraction of a multiple of 1.0. Typical
   * use is moving an atom within unit cell boundaries.
   *
   * @param value Double to be moved between 0 and 1.
   * @return Shifted double.
   */
  private static double moveBetweenZeroAndOne(double value) {
    if (value < 0.0) {
      int belowZero = (int) (value);
      belowZero = 1 + (-1 * belowZero);
      value = value + belowZero;
    } else {
      value = value % 1.0;
    }
    return value;
  }

  /**
   * Returns the linear index of this Box.
   *
   * @return Linear index.
   */
  public int getCellIndex() {
    return cellIndex;
  }

  /**
   * Returns a copy of the ArrayList of residues.
   *
   * @return ArrayList of Residues in the cell.
   */
  public ArrayList<Residue> getResiduesAsList() {
    return new ArrayList<>(residues);
  }

  /**
   * Returns the a, b, and c axis indices of this cell.
   *
   * @return Cell indices.
   */
  public int[] getABCIndices() {
    int[] returnedIndices = new int[3];
    arraycopy(indices, 0, returnedIndices, 0, returnedIndices.length);
    return returnedIndices;
  }

  /**
   * Checks if a Residue is inside this BoxOptCell.
   *
   * @param residue      Residue to check.
   * @param crystal      A Crystal.
   * @param symOp        A symmetry operator to apply.
   * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
   * @return If contained inside this BoxOptCell.
   */
  public boolean residueInsideCell(Residue residue, Crystal crystal, SymOp symOp, boolean variableOnly) {
    List<Atom> atoms = variableOnly ? residue.getVariableAtoms() : residue.getAtomList();
    return atoms.stream().anyMatch(a -> atomInsideCell(a, crystal, symOp));
  }

  /**
   * Sorts residues in the box.
   */
  public void sortCellResidues() {
    Comparator<Residue> comparator = Comparator.comparing(Residue::getChainID)
        .thenComparingInt(Residue::getResidueNumber);
    residues.sort(comparator);
  }

  /**
   * Check if an atom's fractional coordinates would be contained by the box.
   *
   * @param atomFracCoords Atomic fractional coordinates
   * @return If contained.
   */
  private boolean checkIfContained(double[] atomFracCoords) {
    for (int i = 0; i < 3; i++) {
      if (atomFracCoords[i] < fracCoords[i] || atomFracCoords[i] > fracCoords[i + 3]) {
        // fracCoords[0-2] contain min a, b, and c-axes.
        // fracCoords[3-5] contain max a, b, and c-axes.
        return false;
      }
    }
    return true;
  }
}
