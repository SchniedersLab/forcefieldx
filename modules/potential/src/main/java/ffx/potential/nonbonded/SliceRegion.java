// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
 * The SliceLoop class is used to parallelize placing onto a 3D grid
 *
 * <p>1) Multipoles using B-splines or
 *
 * <p>2) Diffraction form factors.
 *
 * <p>Each "slice" of the grid (i.e. a fixed value of the z-coordinate) is operated on by only a
 * single thread to logically enforce atomic updates of grid magnitudes.
 *
 * @author Armin Avdic
 */
public class SliceRegion extends ParallelRegion {

  public int buff = 3;
  public boolean[][] select;
  protected SliceLoop[] sliceLoop;
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
  private int[][] zAtListBuild;

  /**
   * Constructor for SliceRegion.
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
  public SliceRegion(
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
    sliceLoop = new SliceLoop[threadCount];
    gridInitLoop = new GridInitLoop[threadCount];
    for (int i = 0; i < threadCount; i++) {
      gridInitLoop[i] = new GridInitLoop();
    }
    select = new boolean[nSymm][nAtoms];
    for (int i = 0; i < nSymm; i++) {
      fill(select[i], true);
    }
    zAtListBuild = new int[nSymm][nAtoms];
    rebuildList = true;
  }

  /** {@inheritDoc} */
  @Override
  public void finish() {
    if (rebuildList) {
      sliceLoop[0].saveZValues(zAtListBuild);
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

  /** {@inheritDoc} */
  @Override
  public void run() throws Exception {
    int threadIndex = getThreadIndex();
    SliceLoop loop = sliceLoop[threadIndex];
    // This lets the same SpatialDensityLoops be used with different SpatialDensityRegions.
    loop.setNsymm(nSymm);
    loop.setRebuildList(rebuildList);
    try {
      execute(0, gridSize - 1, gridInitLoop[threadIndex]);
      execute(0, gZ - 1, loop);
    } catch (Exception e) {
      String message = " Exception in SliceRegion.";
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
    zAtListBuild = new int[nSymm][nAtoms];
    for (int i = 0; i < nSymm; i++) {
      fill(select[i], true);
    }
    rebuildList = true;
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
   * @param loops an array of {@link ffx.potential.nonbonded.SliceLoop} objects.
   */
  public void setDensityLoop(SliceLoop[] loops) {
    sliceLoop = loops;
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
    rebuildList = (rebuildList || sliceLoop[0].checkList(zAtListBuild, buff));
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
