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

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.numerics.math.Double3;
import ffx.potential.bonded.Atom;
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.math3.geometry.euclidean.threed.Line;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

/**
 * Structure Metrics contains functionality to calculate characteristics of coordinate systems.
 * <p>
 * This includes: Gyrate computes the radius of gyration of a molecular system from its atomic
 * coordinates.
 * <p>
 * Inertia computes the principal moments of inertia for the system, and optionally translates the
 * center of mass to the origin and rotates the principal axes onto the global axes. Reference:
 * Herbert Goldstein, "Classical Mechanics, 2nd Edition", Addison-Wesley, Reading, MA, 1980; see the
 * Euler angle xyz convention in Appendix B
 *
 * @author Michael J. Schnieders
 * @author Aaron J. Nessler
 * @since 1.0
 */
public class StructureMetrics {

  private static final Logger logger = Logger.getLogger(StructureMetrics.class.getName());

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

  /**
   * Compute the components that make up the radius of gyration along three axes.
   *
   * @param atoms Atoms for calculation
   * @param pma Principal moment axes
   * @return radius of gyration along axes
   */
  public static double[][] radiusOfGyrationComponents(Atom[] atoms, boolean pma) {
    int nAtoms = atoms.length;
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    double[] mass = new double[nAtoms];

    int index = 0;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      mass[i] = atom.getMass();
      x[index] = atom.getX();
      y[index] = atom.getY();
      z[index] = atom.getZ();
      index++;
    }

    return radiusOfGyrationComponents(x, y, z, mass, pma);
  }

  /**
   * Compute the components that make up the radius of gyration along three axes.
   *
   * @param xyz Coordinates for calculation
   * @param mass Mass of each atom
   * @param pma Use principal moment axes
   * @return radius of gyration along axes
   */
  public static double[][] radiusOfGyrationComponents(double[] xyz, double[] mass, boolean pma) {
    assert (xyz.length % 3 == 0);

    int nAtoms = xyz.length / 3;
    // Find the centroid of the atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      x[i] = xyz[index++];
      y[i] = xyz[index++];
      z[i] = xyz[index];
    }

    return radiusOfGyrationComponents(x, y, z, mass, pma);
  }

  /**
   * Compute the components that make up the radius of gyration along three axes.
   *
   * @param x Coordinates for calculation
   * @param y Coordiantes for calculation
   * @param z Coordinates for calculation
   * @param mass Mass of each atom
   * @param pma Use principal moment axes
   * @return radius of gyration along axes
   */
  public static double[][] radiusOfGyrationComponents(double[] x, double[] y, double[] z,
      double[] mass, boolean pma) {
    assert (x.length <= mass.length);
    double massSum = Arrays.stream(mass).sum();

    double[][] inertia = momentsOfInertia(x, y, z, mass, false, false, pma);

    // Rg = sqrt(I/m)
    for (int i = 0; i < inertia.length; i++) {
      inertia[i][0] = sqrt(inertia[i][0] / massSum);
    }

    return inertia;
  }

  /**
   * Compute the components that make up the radius of gyration along three axes.
   *
   * @param xyz Coordinates for calculation
   * @param unitCell crystal used as system basis.
   * @return radius of gyration along axes
   */
  public static double[][] radiusOfGyrationComponents(double[] xyz, Crystal unitCell) {
    assert (xyz.length % 3 == 0);
    int nAtoms = xyz.length / 3;
    // Find the centroid of the atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      x[i] = xyz[index++];
      y[i] = xyz[index++];
      z[i] = xyz[index];
    }

    return radiusOfGyrationComponents(x, y, z, unitCell);
  }

  /**
   * Compute the components that make up the radius of gyration tensor about yz-, xz-, xy-planes.
   *
   * @param x Coordinates for calculation
   * @param y Coordinates for calculation
   * @param z Coordinates for calculation
   * @param unitCell Crystal to use as basis.
   * @return radius of gyration about planes.
   */
  public static double[][] radiusOfGyrationComponents(double[] x, double[] y, double[] z,
      Crystal unitCell) {
    assert (x.length == y.length);
    assert (y.length == z.length);

    // Define a Line for each unit cell vector.
//    Vector3D origin = new Vector3D(0.0, 0.0, 0.0);
//    Vector3D aAxis = new Vector3D(unitCell.Ai[0].clone());
//    Vector3D bAxis = new Vector3D(unitCell.Ai[1].clone());
//    Vector3D cAxis = new Vector3D(unitCell.Ai[2].clone());
//    Line[] axes = new Line[3];
//    axes[0] = new Line(origin, aAxis, 1.e0-8);
//    axes[1] = new Line(origin, bAxis, 1.e0-8);
//    axes[2] = new Line(origin, cAxis, 1.e0-8);

    // Find the centroid of the atomic coordinates.
    int nAtoms = x.length;
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    for (int i = 0; i < nAtoms; i++) {
      xc += x[i];
      yc += y[i];
      zc += z[i];
    }
    Vector3D centroid = new Vector3D(xc, yc, zc);
    centroid = centroid.scalarMultiply(1.0 / nAtoms);

    // Compute the radius of gyration about each unit cell vector.
    double[] radius = new double[3];
    for (int j = 0; j < nAtoms; j++) {
      Vector3D xyz = new Vector3D(x[j], y[j], z[j]);
      xyz = xyz.subtract(centroid);
      double[] coords = xyz.toArray();
      for (int i=0; i<3; i++) {
        // double mag = axes[i].distance(xyz);
        double mag = coords[i];
        radius[i] += mag * mag;
      }
    }

    for (int i = 0; i < 3; i++) {
      radius[i] = sqrt(radius[i] / nAtoms);
    }

    double[][] momentsAndVectors = new double[3][4];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 4; j++) {
        if (j == 0) {
          momentsAndVectors[i][j] = radius[i];
        } else {
          momentsAndVectors[i][j] = unitCell.Ai[i][j-1];
        }
      }
    }

    return momentsAndVectors;
  }

  /**
   * Compute the moments of inertia for all atoms in the supplied array.
   *
   * @param atoms Atom array.
   * @param moved Move coordinates
   * @param print Display values to user
   * @param pma Use principal moment axes.
   * @return The moments of inertia.
   */
  public static double[][] momentsOfInertia(Atom[] atoms, boolean moved, boolean print,
      boolean pma) {
    double[] mass = new double[atoms.length];
    int nAtoms = atoms.length;
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    int index = 0;
    for (Atom atom : atoms) {
      mass[index] = atom.getMass();
      x[index] = atom.getX();
      y[index] = atom.getY();
      z[index] = atom.getZ();
      index++;
    }

    return momentsOfInertia(x, y, z, mass, moved, print, pma);
  }

  /**
   * Compute the moments of inertia.
   *
   * @param xyz Array of atomic coordinates (xyz = [X0, Y0, Z0, X1, Y1, Z1, ...].
   * @param mass Mass of atoms
   * @param moved Move the from original coordinates to selection
   * @param print Display values to user
   * @param pma Use principal moment axes.
   * @return The radius of gyration.
   */
  public static double[][] momentsOfInertia(double[] xyz, double[] mass, boolean moved,
      boolean print, boolean pma) {
    assert (xyz.length % 3 == 0);
    int nAtoms = xyz.length / 3;
    // Find the centroid of the atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      x[i] = xyz[index++];
      y[i] = xyz[index++];
      z[i] = xyz[index];
    }

    return momentsOfInertia(x, y, z, mass, moved, print, pma);
  }

  /**
   * Compute the moments of inertia
   *
   * @param x Array of atomic X-coordinates.
   * @param y Array of atomic X-coordinates.
   * @param z Array of atomic X-coordinates.
   * @param mass mass of atoms
   * @param moved Move coordinates to principal axes
   * @param print Print out values to screen
   * @param pma Report moments of inertia to principal axes.
   * @return The moment of inertia.
   */
  public static double[][] momentsOfInertia(double[] x, double[] y, double[] z, double[] mass,
      boolean moved, boolean print, boolean pma) {
    assert (x.length == y.length);
    assert (y.length == z.length);

    // Find the centroid of the atomic coordinates.
    double total = 0.0;
    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    int nAtoms = x.length;
    for (int i = 0; i < nAtoms; i++) {
      double massValue = mass[i];
      total += massValue;
      xcm += x[i] * massValue;
      ycm += y[i] * massValue;
      zcm += z[i] * massValue;
    }
    xcm /= total;
    ycm /= total;
    zcm /= total;

    // Compute and then diagonalize the inertia tensor
    double xx = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yy = 0.0;
    double yz = 0.0;
    double zz = 0.0;

    double xterm;
    double yterm;
    double zterm;

    for (int i = 0; i < nAtoms; i++) {
      double massValue = mass[i];
      xterm = x[i] - xcm;
      yterm = y[i] - ycm;
      zterm = z[i] - zcm;
      xx += xterm * xterm * massValue;
      xy += xterm * yterm * massValue;
      xz += xterm * zterm * massValue;
      yy += yterm * yterm * massValue;
      yz += yterm * zterm * massValue;
      zz += zterm * zterm * massValue;
    }
    double[][] tensor = new double[3][3];
    tensor[0][0] = yy + zz;
    tensor[0][1] = -xy;
    tensor[0][2] = -xz;
    tensor[1][0] = -xy;
    tensor[1][1] = xx + zz;
    tensor[1][2] = -yz;
    tensor[2][0] = -xz;
    tensor[2][1] = -yz;
    tensor[2][2] = xx + yy;

    double[] moment;
    double[][] vec;
    if (pma) {
      // Diagonalize the matrix.
      Array2DRowRealMatrix cMatrix = new Array2DRowRealMatrix(tensor, false);
      EigenDecomposition eigenDecomposition = new EigenDecomposition(cMatrix);
      // Extract the quaternions.
      moment = eigenDecomposition.getRealEigenvalues();
      vec = eigenDecomposition.getV().getData();

      // Select the direction for each principal moment axis
      double dot = 0.0;
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < nAtoms; j++) {
          xterm = vec[i][0] * (x[j] - xcm);
          yterm = vec[i][1] * (y[j] - ycm);
          zterm = vec[i][2] * (z[j] - zcm);
          dot = xterm + yterm + zterm;
          if (dot < 0.0) {
            for (int k = 0; k < 3; k++) {
              vec[i][k] = -vec[i][k];
            }
          }
          if (dot != 0.0) {
            break;
          }
        }
        if (dot != 0.0) {
          break;
        }
      }

      // Moment axes must give a right-handed coordinate system.
      xterm = vec[0][0] * (vec[1][1] * vec[2][2] - vec[2][1] * vec[1][2]);
      yterm = vec[0][1] * (vec[2][0] * vec[1][2] - vec[1][0] * vec[2][2]);
      zterm = vec[0][2] * (vec[1][0] * vec[2][1] - vec[2][0] * vec[1][1]);
      dot = xterm + yterm + zterm;
      if (dot < 0.0) {
        for (int i = 0; i < 3; i++) {
          vec[2][i] = -vec[2][i];
        }
      }

      // Principal moment axes form rows of Euler rotation matrix.
      if (moved) {
        double[][] a = new double[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            a[j][i] = vec[j][i];
          }
        }
        // Translate to origin, then apply Euler rotation matrix.
        for (int i = 0; i < nAtoms; i++) {
          xterm = x[i] - xcm;
          yterm = y[i] - ycm;
          zterm = z[i] - zcm;
          x[i] = a[0][0] * xterm + a[1][0] * yterm + a[2][0] * zterm;
          y[i] = a[0][1] * xterm + a[1][1] * yterm + a[2][1] * zterm;
          z[i] = a[0][2] * xterm + a[1][2] * yterm + a[2][2] * zterm;
        }
      }
    } else {
      vec = new double[][] {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
      moment = new double[] {tensor[0][0], tensor[1][1], tensor[2][2]};
    }

    // print the center of mass and Euler angle values
    if (print) {
      logger.info(format("\n Center of Mass Coordinates: %8.4f %8.4f %8.4f", xcm, ycm, zcm));
      // invert vec
      double[] angles = new Rotation(vec, 1.0E-7).getAngles(RotationOrder.XYZ);
      double radian = 180 / PI;
      // Convert to degrees
      for (int i = 0; i < 3; i++) {
        angles[i] += radian;
      }
      logger.info(format(" Euler Angles (Phi/Theta/Psi): %8.3f %8.3f %8.3f", angles[0], angles[1],
          angles[2]));
      logger.info(
          " Moments of Inertia and Principle Axes:\n  Moments (amu Ang^2): \t X-, Y-, and Z-Components of Axes:");
      for (int i = 0; i < 3; i++) {
        logger.info(
            format("  %16.3f %12.6f %12.6f %12.6f", moment[i], vec[i][0], vec[i][1], vec[i][2]));
      }
    }

    double[][] momentsAndVectors = new double[3][4];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 4; j++) {
        if (j == 0) {
          momentsAndVectors[i][j] = moment[i];
        } else {
          momentsAndVectors[i][j] = vec[i][j - 1];
        }
      }
    }

    return momentsAndVectors;
  }

  /**
   * Compute the moments of inertia.
   *
   * @param xyz Array of atomic coordinates (xyz = [X0, Y0, Z0, X1, Y1, Z1, ...].
   * @param mass Mass of atoms
   * @param unitCell Crystal to use as system basis.
   * @return The radius of gyration.
   */
  public static double[][] momentsOfInertia(double[] xyz, double[] mass, Crystal unitCell) {
    assert (xyz.length % 3 == 0);
    int nAtoms = xyz.length / 3;
    // Find the centroid of the atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];

    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      x[i] = xyz[index++];
      y[i] = xyz[index++];
      z[i] = xyz[index];
    }

    return momentsOfInertia(x, y, z, mass, unitCell);
  }

  /**
   * Compute the moments of inertia
   *
   * @param x Array of atomic X-coordinates.
   * @param y Array of atomic X-coordinates.
   * @param z Array of atomic X-coordinates.
   * @param mass Mass of atoms
   * @param unitCell Crystal to use as system basis.
   * @return The moment of inertia.
   */
  public static double[][] momentsOfInertia(double[] x, double[] y, double[] z, double[] mass,
      Crystal unitCell) {
    assert (x.length == y.length);
    assert (y.length == z.length);
    // Crystal orientation
    double[][] radius = radiusOfGyrationComponents(x, y, z, unitCell);
    double massSum = Arrays.stream(mass).sum();

    double[][] momentsAndVectors = new double[3][4];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 4; j++) {
        if (j == 0) {
          momentsAndVectors[i][j] = massSum * radius[i][0] * radius[i][0];
        } else {
          momentsAndVectors[i][j] = radius[i][j - 1];
        }
      }
    }

    return momentsAndVectors;
  }

}