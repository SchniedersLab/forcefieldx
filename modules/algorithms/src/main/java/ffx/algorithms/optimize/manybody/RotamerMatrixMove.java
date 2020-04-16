// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

import ffx.algorithms.mc.MCMove;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/** This implements single-rotamer changes in the framework of the rotamer energy matrices. */
public class RotamerMatrixMove implements MCMove {

  private final boolean useAllElims;
  /**
   * CurrentRots should point to the same array as being used in the overlying MetropolisMC
   * implementation.
   */
  private final int[] currentRots;

  private final int nRes;
  private final List<Integer> allowedRes;
  private final List<List<Integer>> allowedRots;
  private final int nAllowed;
  private RotamerOptimization rotamerOptimization;
  private EliminatedRotamers eR;
  private boolean monteCarloTesting;
  /** When we take a step, we need to remember which rotamer of which residue was changed. */
  private int changedRes;

  private int changedRot;

  /**
   * Constructs the RotamerMatrixMove set; at present, a new object must be made if rotamers or
   * residues are changed outside the scope of this class.
   *
   * @param useAllElims Use eliminated pair/triple info
   * @param rots Initial rotamer set
   */
  public RotamerMatrixMove(
      boolean useAllElims,
      int[] rots,
      Residue[] residues,
      RotamerLibrary library,
      RotamerOptimization rotamerOptimization,
      EliminatedRotamers eR,
      boolean monteCarloTesting) {
    this.useAllElims = useAllElims;
    this.rotamerOptimization = rotamerOptimization;
    this.eR = eR;
    this.monteCarloTesting = monteCarloTesting;

    nRes = rots.length;
    currentRots = rots;

    allowedRes = new ArrayList<>(nRes);
    allowedRots = new ArrayList<>(nRes);

    for (int i = 0; i < nRes; i++) {
      ArrayList<Integer> resAllowed = new ArrayList<>();

      int lenri = residues[i].getRotamers(library).length;
      for (int ri = 0; ri < lenri; ri++) {
        if (!eR.check(i, ri)) {
          resAllowed.add(ri);
        }
      }

      if (resAllowed.size() > 1) {
        resAllowed.trimToSize();
        allowedRes.add(i);
        allowedRots.add(resAllowed);
      }
    }

    ((ArrayList) allowedRes).trimToSize();
    nAllowed = allowedRes.size();
  }

  @Override
  public void move() {
    boolean validMove = !useAllElims;
    int indexI;
    int indexRI;
    do {
      // resi and roti correspond to their positions in allowedRes and
      // allowedRots. indexI and indexRI correspond to their numbers
      // in the rotamer matrix.

      Random rand = new Random();
      if (monteCarloTesting) {
        rand.setSeed(nAllowed);
      }
      int resi = rand.nextInt(nAllowed);
      indexI = allowedRes.get(resi);
      List<Integer> rotsi = allowedRots.get(resi);
      int lenri = rotsi.size();

      int roti = rand.nextInt(lenri);
      indexRI = rotsi.get(roti);
      if (useAllElims) {
        validMove = rotamerOptimization.checkValidMove(indexI, indexRI, currentRots);
      }
    } while (!validMove);

    changedRes = indexI;
    changedRot = currentRots[indexI];

    currentRots[indexI] = indexRI;
  }

  @Override
  public void revertMove() {
    currentRots[changedRes] = changedRot;
  }

  @Override
  public String toString() {
    return "Rotamer moves utilizing a rotamer energy matrix";
  }
}
