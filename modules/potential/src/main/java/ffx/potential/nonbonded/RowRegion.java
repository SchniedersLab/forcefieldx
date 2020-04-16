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
package ffx.potential.nonbonded;

import static ffx.potential.nonbonded.SpatialDensityRegion.logger;
import static java.util.Arrays.fill;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import java.nio.DoubleBuffer;
import java.util.logging.Level;

/**
 * The RowRegion class is used to parallelize placing onto a 3D grid
 *
 * <p>1) Multipoles using B-splines or
 *
 * <p>2) Diffraction form factors.
 *
 * <p>Each "row" of 3D grid (i.e. fixed values of the z and y-coordinates) is operated on by only a
 * single thread to logically enforce atomic updates of grid magnitudes.
 *
 * @author Armin Avdic
 */
public class RowRegion extends ParallelRegion {

  public int buff = 3;
  public boolean[][] select;
  protected RowLoop[] rowLoop;
  protected double[][][] coordinates;
  int nAtoms;
  int nSymm;
  private int gX, gY, gZ;
  private DoubleBuffer gridBuffer;
  private GridInitLoop[] gridInitLoop;
  private double initValue = 0.0;
  private int gridSize;
  private double[] grid;
  private boolean rebuildList;
  private int[][][] zyAtListBuild;

  /**
   * Constructor for RowRegion.
   *
   * @param gX a int.
   * @param gY a int.
   * @param gZ a int.
   * @param grid an array of {@link double} objects.
   * @param nSymm a int.
   * @param threadCount a int.
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param coordinates an array of {@link double} objects.
   */
  public RowRegion(
      int gX,
      int gY,
      int gZ,
      double[] grid,
      int nSymm,
      int threadCount,
      Atom[] atoms,
      double[][][] coordinates) {
    this.nAtoms = atoms.length;
    this.gX = gX;
    this.gY = gY;
    this.gZ = gZ;
    gridSize = gX * gY * gZ * 2;
    this.nSymm = nSymm;
    this.coordinates = coordinates;
    this.grid = grid;
    if (grid != null) {
      gridBuffer = DoubleBuffer.wrap(grid);
    }
    rowLoop = new RowLoop[threadCount];
    gridInitLoop = new GridInitLoop[threadCount];
    for (int i = 0; i < threadCount; i++) {
      gridInitLoop[i] = new GridInitLoop();
    }
    select = new boolean[nSymm][nAtoms];
    for (int i = 0; i < nSymm; i++) {
      fill(select[i], true);
    }
    zyAtListBuild = new int[nSymm][nAtoms][2];
    rebuildList = true;
  }

  /** {@inheritDoc} */
  @Override
  public void finish() {
    if (rebuildList) {
      rowLoop[0].saveZYValues(zyAtListBuild);
    }
    rebuildList = false;
  }

  /**
   * Getter for the field <code>grid</code>.
   *
   * @return an array of {@link double} objects.
   */
  public double[] getGrid() {
    return grid;
  }

  /**
   * getNatoms.
   *
   * @return a int.
   */
  public int getNatoms() {
    return nAtoms;
  }

  /**
   * getNsymm.
   *
   * @return a int.
   */
  public int getNsymm() {
    return nSymm;
  }

  /**
   * rowIndexForYZ.
   *
   * @param giy a int.
   * @param giz a int.
   * @return a int.
   */
  public int rowIndexForYZ(int giy, int giz) {
    return giy + gY * giz;
  }

  /** {@inheritDoc} */
  @Override
  public void run() throws Exception {
    int threadIndex = getThreadIndex();
    RowLoop loop = rowLoop[threadIndex];
    // This lets the same SpatialDensityLoops be used with different SpatialDensityRegions.
    loop.setNsymm(nSymm);
    try {
      execute(0, gridSize - 1, gridInitLoop[threadIndex]);
      execute(0, (gY * gZ) - 1, loop);
    } catch (Exception e) {
      String message = " Exception in RowRegion.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Select atoms that should be included. The default is to include all atoms, which is set up in
   * the constructor. This function should be over-ridden by subclasses that want finer control.
   */
  public void selectAtoms() {
    for (int i = 0; i < nSymm; i++) {
      fill(select[i], true);
    }
  }

  /**
   * Setter for the field <code>atoms</code>.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public void setAtoms(Atom[] atoms) {
    nAtoms = atoms.length;
    select = new boolean[nSymm][nAtoms];
    for (int i = 0; i < nSymm; i++) {
      fill(select[i], true);
    }
  }

  /**
   * Setter for the field <code>crystal</code>.
   *
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param gX a int.
   * @param gY a int.
   * @param gZ a int.
   */
  public final void setCrystal(Crystal crystal, int gX, int gY, int gZ) {
    // this.crystal = crystal.getUnitCell();
    this.gX = gX;
    this.gY = gY;
    this.gZ = gZ;
    gridSize = gX * gY * gZ * 2;
  }

  /**
   * setDensityLoop.
   *
   * @param loops an array of {@link ffx.potential.nonbonded.RowLoop} objects.
   */
  public void setDensityLoop(RowLoop[] loops) {
    rowLoop = loops;
  }

  /**
   * Setter for the field <code>initValue</code>.
   *
   * @param initValue a double.
   */
  public void setInitValue(double initValue) {
    this.initValue = initValue;
  }

  /** {@inheritDoc} */
  @Override
  public void start() {
    selectAtoms();
    rebuildList = (rebuildList || rowLoop[0].checkList(zyAtListBuild, buff));
  }

  /**
   * yFromRowIndex.
   *
   * @param i a int.
   * @return a int.
   */
  public int yFromRowIndex(int i) {
    return i % gY;
  }

  /**
   * zFromRowIndex.
   *
   * @param i a int.
   * @return a int.
   */
  public int zFromRowIndex(int i) {
    return i / gY;
  }

  /**
   * Setter for the field <code>gridBuffer</code>.
   *
   * @param grid a {@link java.nio.DoubleBuffer} object.
   */
  void setGridBuffer(DoubleBuffer grid) {
    gridBuffer = grid;
  }

  private class GridInitLoop extends IntegerForLoop {

    private final IntegerSchedule schedule = IntegerSchedule.fixed();
    // Extra padding to avert cache interference.
    long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    long pad8, pad9, pada, padb, padc, padd, pade, padf;

    @Override
    public void run(int lb, int ub) {
      if (gridBuffer != null) {
        // if (grid != null) {
        for (int i = lb; i <= ub; i++) {
          // grid[i] = initValue;
          gridBuffer.put(i, initValue);
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }
  }
}
