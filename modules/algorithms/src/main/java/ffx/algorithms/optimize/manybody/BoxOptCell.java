// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.optimize.manybody;

import static ffx.potential.bonded.RotamerLibrary.applyRotamer;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NeighborList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Contains a cell used for box optimization, its residues, the fractional coordinates within a
 * crystal it takes up, its overall (linear) index, and its indices along the a, b, and c crystal
 * axes.
 */
public class BoxOptCell {

  // fracCoords indexed by 1-3 min x,y,z, 4-6 are max x,y,z
  private final double[] fracCoords = new double[6];
  private final int[] indexXYZ = new int[3];
  private final int linearIndex;
  private final ArrayList<Residue> residues = new ArrayList<>();

  /**
   * Constructs a BoxOptCell object, which takes up a set of fractional coordinates within the
   * Crystal, the Residues contained within, and the index of the cell along the crystal's a, b, and
   * c axes.
   *
   * @param fractionalCoordinates Fractional coordinates contained, indexed by 1-3 min x,y,z, 4-6
   *     max x,y,z
   * @param indices Index of cell along a, b, and c (x, y, and z).
   * @param linearIndex Index of box in linear box array.
   */
  public BoxOptCell(double[] fractionalCoordinates, int[] indices, int linearIndex) {
    System.arraycopy(fractionalCoordinates, 0, fracCoords, 0, fracCoords.length);
    System.arraycopy(indices, 0, indexXYZ, 0, indexXYZ.length);
    this.linearIndex = linearIndex;
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
   * @param residue Residue to check.
   * @param crystal A Crystal.
   * @param symOp A symmetry operator to apply.
   * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
   * @return If contained inside this BoxOptCell.
   */
  public boolean anyRotamerInsideCell(
      Residue residue,
      Crystal crystal,
      SymOp symOp,
      boolean variableOnly) {
    ResidueState incomingState = residue.storeState();
    Rotamer[] rotamers = residue.getRotamers();
    boolean inside =
        Arrays.stream(rotamers)
            .anyMatch(
                (Rotamer r) -> {
                  applyRotamer(residue, r);
                  return residueInsideCell(residue, crystal, symOp, variableOnly);
                });
    residue.revertState(incomingState);
    return inside;
  }

  /**
   * Checks if an Atom would be contained inside this cell.
   *
   * @param atom Atom to check.
   * @param crystal A Crystal.
   * @param symOp A symmetry operator to apply.
   * @return If contained.
   */
  public boolean atomInsideCell(Atom atom, Crystal crystal, SymOp symOp) {
    double[] atXYZ = new double[3];
    atXYZ = atom.getXYZ(atXYZ);
    crystal.toFractionalCoordinates(atXYZ, atXYZ);
    NeighborList.moveValuesBetweenZeroAndOne(atXYZ);
    crystal.applyFracSymOp(atXYZ, atXYZ, symOp);
    return checkIfContained(atXYZ);
  }

  /**
   * Returns the linear index of this Box.
   *
   * @return Linear index.
   */
  public int getLinearIndex() {
    return linearIndex;
  }

  /**
   * Returns an array of the Residues contained within the cell.
   *
   * @return Array of Residues.
   */
  public Residue[] getResidues() {
    return (Residue[]) residues.toArray();
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
   * Returns the x, y, and z indices of this box.
   *
   * @return Box indices.
   */
  public int[] getXYZIndex() {
    int[] returnedIndices = new int[3];
    System.arraycopy(indexXYZ, 0, returnedIndices, 0, returnedIndices.length);
    return returnedIndices;
  }

  /**
   * Checks if a Residue is inside this BoxOptCell.
   *
   * @param residue Residue to check.
   * @param crystal A Crystal.
   * @param symOp A symmetry operator to apply.
   * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
   * @return If contained inside this BoxOptCell.
   */
  public boolean residueInsideCell(
      Residue residue, Crystal crystal, SymOp symOp, boolean variableOnly) {
    List<Atom> atoms = variableOnly ? residue.getVariableAtoms() : residue.getAtomList();
    return atoms.stream().anyMatch(a -> atomInsideCell(a, crystal, symOp));
  }

  /** Sorts residues in the box. */
  public void sortBoxResidues() {
    Comparator<Residue> comparator =
        Comparator.comparing(Residue::getChainID).thenComparingInt(Residue::getResidueNumber);
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
      if ((fracCoords[i] > atomFracCoords[i]) || (fracCoords[i + 3] < atomFracCoords[i])) {
        return false;
      }
    }
    return true;
  }
}
