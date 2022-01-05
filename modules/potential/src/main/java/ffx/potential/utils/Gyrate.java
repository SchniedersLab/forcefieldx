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
package ffx.potential.utils;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.Atom;

/**
 * Gyrate computes the radius of gyration of a molecular system from its atomic coordinates.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Gyrate {

  /**
   * Compute the radius of gyration for all atoms in the supplied array.
   *
   * @param atoms Atom array.
   * @return The radius of gyration.
   */
  public static double radiusOfGyration(Atom[] atoms) {
    int nAtoms = atoms.length;
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    int index = 0;
    for (Atom atom : atoms) {
      x[index] = atom.getX();
      y[index] = atom.getY();
      z[index] = atom.getZ();
      index++;
    }

    return radiusOfGyration(x, y, z);
  }

  /**
   * Compute the radius of gyration.
   *
   * @param x Array of atomic X-coordinates.
   * @param y Array of atomic X-coordinates.
   * @param z Array of atomic X-coordinates.
   * @return The radius of gyration.
   */
  public static double radiusOfGyration(double[] x, double[] y, double[] z) {
    assert (x.length == y.length);
    assert (y.length == z.length);

    // Find the centroid of the atomic coordinates.
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    int nAtoms = x.length;
    for (int i = 0; i < nAtoms; i++) {
      xc += x[i];
      yc += y[i];
      zc += z[i];
    }
    xc /= nAtoms;
    yc /= nAtoms;
    zc /= nAtoms;

    // Compute the radius of gyration.
    double rg = 0.0;
    for (int i = 0; i < nAtoms; i++) {
      double dx = x[i] - xc;
      double dy = y[i] - yc;
      double dz = z[i] - zc;
      rg += dx * dx + dy * dy + dz * dz;
    }
    rg = sqrt(rg / nAtoms);

    return rg;
  }

  /**
   * Compute the radius of gyration.
   *
   * @param xyz Array of atomic coordinates (xyz = [X0, Y0, Z0, X1, Y1, Z1, ...].
   * @return The radius of gyration.
   */
  public static double radiusOfGyration(double[] xyz) {
    assert (xyz.length % 3 == 0);

    // Find the centroid of the atomic coordinates.
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    int nAtoms = xyz.length / 3;
    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      xc += xyz[index++];
      yc += xyz[index++];
      zc += xyz[index];
    }
    xc /= nAtoms;
    yc /= nAtoms;
    zc /= nAtoms;

    // Compute the radius of gyration
    double rg = 0.0;
    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      double dx = xyz[index++] - xc;
      double dy = xyz[index++] - yc;
      double dz = xyz[index] - zc;
      rg += dx * dx + dy * dy + dz * dz;
    }

    rg = sqrt(rg / nAtoms);

    return rg;
  }

}
