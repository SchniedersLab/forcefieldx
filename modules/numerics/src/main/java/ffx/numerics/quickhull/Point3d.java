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
package ffx.numerics.quickhull;

/**
 * A three-element spatial point. The only difference between a point and a
 * vector is in the way it is transformed by an affine transformation. Since
 * the transform method is not included in this reduced implementation for
 * QuickHull3D, the difference is purely academic.
 *
 * @author John E. Lloyd, Fall 2004
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Point3d extends Vector3d {

  /**
   * Creates a Point3d and initializes it to zero.
   */
  public Point3d() {
  }

  /**
   * Creates a Point3d by copying a vector
   *
   * @param v vector to be copied
   */
  public Point3d(Vector3d v) {
    set(v);
  }

  /**
   * Creates a Point3d with the supplied element values.
   *
   * @param x first element
   * @param y second element
   * @param z third element
   */
  public Point3d(double x, double y, double z) {
    set(x, y, z);
  }
}
