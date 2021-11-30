//******************************************************************************
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
//******************************************************************************
package ffx.potential.nonbonded.octree;

import java.util.ArrayList;

/**
 * OctreeCell: Object class for Octree method presented in the Fast Multipole Method (FMM) tutorial
 * from the Barba Group: https://github.com/barbagroup/FMM_tutorial
 */
public class OctreeCell {

  /**
   * Critical (maximum allowed) number of points allowed in any one cell: If a cell already contains
   * nCritical points, it needs to be split
   */
  private int nCritical;
  /** Number of leaves in a cell */
  private int numLeaves = 0;
  /** Array of leaf indices */
  private ArrayList<Integer> leaves = new ArrayList<>();
  /** Integer whose last 8 bits keep track of the empty child cells */
  private int nChild = 0;
  /** Array of child indices, length 8 */
  private int[] children = new int[] {0, 0, 0, 0, 0, 0, 0, 0};
  /** Parent cell index */
  private int parentIndex = 0;

  /** Coordinates for the center of the cell */
  private double x = 0.0;

  private double y = 0.0;
  private double z = 0.0;

  /** Radius of cell */
  private double r = 0.0;

  /** Multipole array, length 10 */
  private double[] multipole = new double[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  public OctreeCell(int nCritical) {
    setnCritical(nCritical);
  }

  public void addToMultipole(double[] calculatedMultipole) {
    for (int i = 0; i < 10; i++) {
      this.multipole[i] += calculatedMultipole[i];
    }
  }

  public int getChildAtIndex(int octant) {
    return children[octant];
  }

  public int getLeavesValueAtIndex(int index) {
    return this.leaves.get(index);
  }

  public double[] getMultipole() {
    return this.multipole;
  }

  public int getNumLeaves() {
    return this.numLeaves;
  }

  public void setNumLeaves(int num) {
    this.numLeaves = num;
  }

  public int getParentIndex() {
    return this.parentIndex;
  }

  public void setParentIndex(int p) {
    this.parentIndex = p;
  }

  /**
   * Returns cell radius
   *
   * @return Returns the cell radius.
   */
  public double getR() {
    return this.r;
  }

  /**
   * Sets cell radius
   *
   * @param r The radius.
   */
  public void setR(double r) {
    this.r = r;
  }

  /**
   * Gets x coordinate of center of cell
   *
   * @return Returns the x coordinate for center of the cell.
   */
  public double getX() {
    return this.x;
  }

  /**
   * Sets x coordinate of center of cell
   *
   * @param x Cell x-coordinate.
   */
  public void setX(double x) {
    this.x = x;
  }

  /**
   * Gets y coordinate of center of cell
   *
   * @return Returns the y coordinate for center of the cell.
   */
  public double getY() {
    return this.y;
  }

  /**
   * Sets y coordinate of center of cell
   *
   * @param y Cell y-coordinate.
   */
  public void setY(double y) {
    this.y = y;
  }

  /**
   * Gets z coordinate of center of cell
   *
   * @return Returns the z coordinate for center of the cell.
   */
  public double getZ() {
    return this.z;
  }

  /**
   * Sets z coordinate of center of cell
   *
   * @param z Cell z-coordinate.
   */
  public void setZ(double z) {
    this.z = z;
  }

  public int getnChild() {
    return this.nChild;
  }

  public void setnChild(int num) {
    this.nChild = num;
  }

  public void setChildren(int octant, int c) {
    this.children[octant] = c;
  }

  public void setLeaf(int index, int leaf) {
    this.leaves.set(index, leaf);
  }

  public void setnCritical(int nCrit) {
    this.nCritical = nCrit;
  }
}
