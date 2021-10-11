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

import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import java.util.logging.Level;
import java.util.logging.Logger;

/** Compute the minimum distance between each pair of residues for all rotamer permutations. */
public class DistanceRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(DistanceRegion.class.getName());
  private final DistanceLoop[] distanceLoops;
  private final int nResidues;
  private final Crystal crystal;
  private final int nSymm;
  private final int[][][] lists;
  private final IntegerSchedule pairwiseSchedule;
  /** AlgorithmListener who should receive updates as the optimization runs. */
  protected AlgorithmListener algorithmListener;

  private RotamerOptimization rotamerOptimization;
  private DistanceMatrix dM;
  /** MolecularAssembly to perform rotamer optimization on. */
  private MolecularAssembly molecularAssembly;
  /**
   * An array of all residues being optimized. Note that Box and Window optimizations operate on
   * subsets of this list.
   */
  private Residue[] allResiduesArray;
  /** RotamerLibrary instance. */
  private RotamerLibrary library;
  /**
   * The minimum distance between atoms of a residue pair, taking into account interactions with
   * symmetry mates.
   *
   * <p>[residue1][rotamer1][residue2][rotamer2]
   */
  private double[][][][] distanceMatrix;

  public DistanceRegion(
      int nt, int nResidues, Crystal crystal, int[][][] lists, IntegerSchedule schedule) {
    distanceLoops = new DistanceLoop[nt];
    this.nResidues = nResidues;
    this.crystal = crystal;
    this.nSymm = crystal.spaceGroup.getNumberOfSymOps();
    this.lists = lists;
    for (int i = 0; i < nt; i++) {
      distanceLoops[i] = new DistanceLoop();
    }
    pairwiseSchedule = schedule;
  }

  public void init(
      DistanceMatrix dM,
      RotamerOptimization rotamerOptimization,
      MolecularAssembly molecularAssembly,
      Residue[] allResiduesArray,
      RotamerLibrary library,
      AlgorithmListener algorithmListener,
      double[][][][] distanceMatrix) {
    this.dM = dM;
    this.rotamerOptimization = rotamerOptimization;
    this.molecularAssembly = molecularAssembly;
    this.allResiduesArray = allResiduesArray;
    this.library = library;
    this.algorithmListener = algorithmListener;
    this.distanceMatrix = distanceMatrix;
  }

  @Override
  public void run() {
    try {
      int threadID = getThreadIndex();
      execute(0, nResidues - 1, distanceLoops[threadID]);
    } catch (Exception e) {
      String message = " Exception computing residue-residue distances.";
      e.printStackTrace();
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Get the coordinates of a requested residue.
   *
   * @param i Residue index.
   * @param residues Array of residues.
   * @param rotamer Rotamer to apply.
   * @param forced True for a forced residue.
   * @return Returns the coordinates.
   */
  private double[][] getCoordinates(int i, Residue[] residues, Rotamer rotamer, boolean forced) {
    synchronized (residues[i]) {
      Residue residue = residues[i];
      if (!forced) {
        RotamerLibrary.applyRotamer(residue, rotamer);
        return residue.storeCoordinateArray();
      } else {
        return residue.storeCoordinateArray();
      }
    }
  }

  private class DistanceLoop extends IntegerForLoop {

    @Override
    public void run(int lb, int ub) {
      // Loop over symmetry operators.
      for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
        SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
        // Loop over residues.
        for (int i = lb; i <= ub; i++) {
          Residue residuei = allResiduesArray[i];
          Rotamer[] rotamersi = residuei.getRotamers();
          int lengthRi;
          boolean forcedResidueI = false;
          try {
            if (rotamerOptimization.checkIfForced(residuei)) {
              forcedResidueI = true;
              lengthRi = 1;
            } else {
              lengthRi = rotamersi.length;
            }
          } catch (IndexOutOfBoundsException ex) {
            logger.warning(format(" Exception in distance loop: %s", ex.toString()));
            continue;
          }

          int[] list = lists[iSymOp][i];
          int nList = list.length;

          // Loop over Residue i's rotamers
          for (int ri = 0; ri < lengthRi; ri++) {
            double[][] xi;
            if (forcedResidueI) {
              xi = getCoordinates(i, allResiduesArray, null, forcedResidueI);
            } else {
              xi = getCoordinates(i, allResiduesArray, rotamersi[ri], forcedResidueI);
            }

            // Loop over Residue i's neighbors.
            for (int j : list) {
              if (i == j) {
                continue;
              }

              Residue residuej = allResiduesArray[j];
              Rotamer[] rotamersj = residuej.getRotamers();
              int lengthRj;
              boolean forcedResidueJ = false;
              try {
                if (rotamerOptimization.checkIfForced(residuej)) {
                  forcedResidueJ = true;
                  lengthRj = 1;
                } else {
                  lengthRj = rotamersj.length;
                }
              } catch (IndexOutOfBoundsException ex) {
                continue;
              }

              // Loop over the neighbor's rotamers
              for (int rj = 0; rj < lengthRj; rj++) {
                double[][] xj;
                if (forcedResidueJ) {
                  xj = getCoordinates(j, allResiduesArray, null, forcedResidueJ);
                } else {
                  xj = getCoordinates(j, allResiduesArray, rotamersj[rj], forcedResidueJ);
                }
                if (getThreadIndex() == 0 && algorithmListener != null) {
                  algorithmListener.algorithmUpdate(molecularAssembly);
                }

                double r = dM.interResidueDistance(xi, xj, symOp);
                if (i < j) {
                  if (r < distanceMatrix[i][ri][j][rj]) {
                    distanceMatrix[i][ri][j][rj] = r;
                  }
                } else if (r < distanceMatrix[j][rj][i][ri]) {
                  distanceMatrix[j][rj][i][ri] = r;
                }
              }
            }
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return pairwiseSchedule;
    }
  }
}
