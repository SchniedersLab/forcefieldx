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

import static ffx.numerics.math.DoubleMath.dist2;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.sqrt;

import com.github.quickhull3d.QuickHull3D;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Constants;
import java.util.Arrays;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * This ConvexHullOps class uses the QuickHull3D package by John E. Lloyd to construct and operate
 * on 3D convex hulls: the minimal convex polyhedron that contains all points in a set of points.
 * This is especially useful for max-dist operations, as the most distant points in a set are
 * guaranteed to be part of the convex polyhedron.
 *
 * <p>The QuickHull3D package website is at quickhull3d.github.io/quickhull3d/ The algorithm it uses
 * is described in Barber, Dobkin, and Huhdanpaa, "The Quickhull Algorithm for Convex Hulls" (ACM
 * Transactions on Mathematical Software, Vol. 22, No. 4, December 1996), and the code is based on
 * the C package qhull.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0.0
 */
public class ConvexHullOps {
  private static final Logger logger = Logger.getLogger(ConvexHullOps.class.getName());

  /**
   * Constructs a convex hull from a MolecularAssembly.
   *
   * @param molecularAssembly MolecularAssembly to build a convex hull for.
   * @return A QuickHull3D implementation of convex hulls.
   */
  public static QuickHull3D constructHull(MolecularAssembly molecularAssembly) {
    Atom[] atoms = molecularAssembly.getAtomArray();
    return constructHull(atoms);
  }

  /**
   * Constructs a convex hull from a set of atoms.
   *
   * @param atoms Atoms to build a convex hull for.
   * @return A QuickHull3D implementation of convex hulls.
   */
  public static QuickHull3D constructHull(Atom[] atoms) {
    int nAts = atoms.length;
    if (nAts < 4) {
      throw new IllegalArgumentException(
          format(" 3D convex hull ill-defined for less than 4 points, found %d", nAts));
    }
    double[] xyz = new double[nAts * 3];
    for (int i = 0; i < nAts; i++) {
      Atom at = atoms[i];
      int i3 = 3 * i;
      xyz[i3] = at.getX();
      xyz[i3 + 1] = at.getY();
      xyz[i3 + 2] = at.getZ();
    }
    return new QuickHull3D(xyz);
  }

  /**
   * UNTESTED: Identifies atoms forming the convex hull.
   *
   * @param quickHull3D A QuickHull3D.
   * @param allAtoms Atoms used in building the QuickHull3D.
   * @return Atoms forming the convex hull.
   */
  public static Atom[] identifyHullAtoms(QuickHull3D quickHull3D, Atom[] allAtoms) {
    int[] indices = quickHull3D.getVertexPointIndices();
    return Arrays.stream(indices).mapToObj((int i) -> allAtoms[i]).toArray(Atom[]::new);
  }

  /**
   * Find the maximum pairwise distance between vertex points on a convex hull.
   *
   * @param quickHull3D A QuickHull3D object.
   * @return Maximum vertex-vertex distance.
   */
  public static double maxDist(QuickHull3D quickHull3D) {
    long time = -System.nanoTime();
    int nVerts = quickHull3D.getNumVertices();
    if (nVerts < 2) {
      return 0;
    }
    double[] vertPoints = new double[3 * nVerts];
    quickHull3D.getVertices(vertPoints);
    double maxDist =
        IntStream.range(0, nVerts)
            .parallel()
            .mapToDouble(
                (int i) -> {
                  double[] xyz = new double[3];
                  arraycopy(vertPoints, 3 * i, xyz, 0, 3);
                  double mij = 0;
                  for (int j = i + 1; j < nVerts; j++) {
                    double[] xyzJ = new double[3];
                    arraycopy(vertPoints, 3 * j, xyzJ, 0, 3);
                    double distIJ = dist2(xyz, xyzJ);
                    mij = max(mij, distIJ);
                  }
                  return mij;
                })
            .max()
            .getAsDouble();
    maxDist = sqrt(maxDist);
    time += System.nanoTime();
    if (time > 1E9) {
      logger.warning(
          format(
              " Required %12.6g sec to find max distance on a convex hull."
                  + " It may be time to further optimize this!",
              Constants.NS2SEC * time));
    }
    return maxDist;
  }

  /**
   * Maximum pairwise distance between atoms in an array. Uses either the convex hull method (more
   * than 10 atoms), or a brute-force loop (10 atoms or less).
   *
   * @param atoms Atoms to check max pairwise distance for.
   * @return Max pairwise distance in Angstroms, or 0 (0 or 1 atoms given).
   */
  public static double maxDist(Atom[] atoms) {
    int nAts = atoms.length;
    if (nAts < 2) {
      return 0;
    } else if (nAts > 10) {
      return maxDist(constructHull(atoms));
    } else {
      double maxDist = 0;
      for (int i = 0; i < nAts - 1; i++) {
        Atom atI = atoms[i];
        double[] xyzI = new double[3];
        xyzI = atI.getXYZ(xyzI);
        for (int j = i + 1; j < nAts; j++) {
          Atom atJ = atoms[j];
          double[] xyzJ = new double[3];
          xyzJ = atJ.getXYZ(xyzJ);
          double dist = dist2(xyzI, xyzJ);
          maxDist = max(dist, maxDist);
        }
      }
      return sqrt(maxDist);
    }
  }
}
