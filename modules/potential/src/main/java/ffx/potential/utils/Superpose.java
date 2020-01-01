//******************************************************************************
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
//******************************************************************************
package ffx.potential.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class Superpose {

    /**
     * Move the center of mass for both sets of atoms to the origin.
     *
     * @param x1    Cartesian coordinates of the first system.
     * @param mass1 The mass of each particle in the first system.
     * @param x2    Cartesian coordinates of the second system.
     * @param mass2 The mass of each particle in the second system.
     */
    public static void translate(double[] x1, double[] mass1, double[] x2, double[] mass2) {
        // Move the first set of atoms.
        translate(x1, mass1);
        // Move the second set of atoms.
        translate(x2, mass2);
    }

    /**
     * Move the center of mass for a set of atoms to the origin.
     *
     * @param x    Cartesian coordinates of the system; modified in-place.
     * @param mass The mass of each particle in the system.
     */
    private static void translate(double[] x, final double[] mass) {
        double[] translation = calculateTranslation(x, mass);
        applyTranslation(x, translation);
    }

    /**
     * Calculate a translation matrix [dx,dy,dz] to return a molecular system's
     * center of mass to the origin.
     *
     * @param x Coordinates of the system.
     * @param mass Mass of each atom.
     * @return Translation to the origin.
     */
    public static double[] calculateTranslation(double[] x, final double[] mass) {
        double xmid = 0.0;
        double ymid = 0.0;
        double zmid = 0.0;
        double norm = 0.0;
        int n = x.length / 3;

        for (int i = 0; i < n; i++) {
            int k = 3 * i;
            double weigh = mass[i];
            xmid = xmid + x[k] * weigh;
            ymid = ymid + x[k + 1] * weigh;
            zmid = zmid + x[k + 2] * weigh;
            norm = norm + weigh;
        }
        xmid = xmid / norm;
        ymid = ymid / norm;
        zmid = zmid / norm;

        return new double[]{xmid, ymid, zmid};
    }

    /**
     * Apply a translation matrix [dx,dy,dz] to a molecular system.
     *
     * @param x Coordinates to move; modified in-place.
     * @param translation Translation matrix.
     */
    public static void applyTranslation(double[] x, final double[] translation) {
        int n = x.length / 3;
        for (int i = 0; i < n; i++) {
            int i3 = 3 * i;
            for (int j = 0; j < 3; j++) {
                x[i3 + j] -= translation[j];
            }
        }
    }

    /**
     * Minimize the RMS distance between two sets of atoms using quaternions.
     *
     * @param x1   Cartesian coordinates of the first system. Unmodified.
     * @param x2   Cartesian coordinates of the second system. Modified in-place.
     * @param mass The mass of each particle in the system.
     */
    public static void rotate(double[] x1, double[] x2, double[] mass) {
        double[][] rotation = calculateRotation(x1, x2, mass);
        applyRotation(x2, rotation);
    }

    /**
     * Calculate a rotation to minimize RMS distance between two sets of atoms
     * using quaternions, overlapping x2 on x1.
     *
     * @param x1   Cartesian coordinates of the first system.
     * @param x2   Cartesian coordinates of the second system.
     * @param mass The mass of each particle in the system.
     * @return A rotation matrix.
     */
    public static double[][] calculateRotation(double[] x1, double[] x2, double[] mass) {

        // Build the upper triangle of the quadratic form matrix
        double xxyx = 0.0;
        double xxyy = 0.0;
        double xxyz = 0.0;
        double xyyx = 0.0;
        double xyyy = 0.0;
        double xyyz = 0.0;
        double xzyx = 0.0;
        double xzyy = 0.0;
        double xzyz = 0.0;

        int n = x1.length / 3;
        for (int i = 0; i < n; i++) {
            int k = i * 3;
            double weigh = mass[i];
            xxyx = xxyx + weigh * x1[k] * x2[k];
            xxyy = xxyy + weigh * x1[k + 1] * x2[k];
            xxyz = xxyz + weigh * x1[k + 2] * x2[k];
            xyyx = xyyx + weigh * x1[k] * x2[k + 1];
            xyyy = xyyy + weigh * x1[k + 1] * x2[k + 1];
            xyyz = xyyz + weigh * x1[k + 2] * x2[k + 1];
            xzyx = xzyx + weigh * x1[k] * x2[k + 2];
            xzyy = xzyy + weigh * x1[k + 1] * x2[k + 2];
            xzyz = xzyz + weigh * x1[k + 2] * x2[k + 2];
        }

        double[][] c = new double[4][4];
        c[0][0] = xxyx + xyyy + xzyz;
        c[0][1] = xzyy - xyyz;
        c[1][0] = c[0][1];
        c[1][1] = xxyx - xyyy - xzyz;
        c[0][2] = xxyz - xzyx;
        c[2][0] = c[0][2];
        c[1][2] = xxyy + xyyx;
        c[2][1] = c[1][2];
        c[2][2] = xyyy - xzyz - xxyx;
        c[0][3] = xyyx - xxyy;
        c[3][0] = c[0][3];
        c[1][3] = xzyx + xxyz;
        c[3][1] = c[1][3];
        c[2][3] = xyyz + xzyy;
        c[3][2] = c[2][3];
        c[3][3] = xzyz - xxyx - xyyy;

        // Diagonalize the quadratic form matrix
        Array2DRowRealMatrix cMatrix = new Array2DRowRealMatrix(c, false);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(cMatrix);

        // Extract the quaternions.
        double[] q = eigenDecomposition.getEigenvector(0).toArray();

        // Assemble rotation matrix that superimposes the molecules
        double[][] rot = new double[3][3];
        double q02 = q[0] * q[0];
        double q12 = q[1] * q[1];
        double q22 = q[2] * q[2];
        double q32 = q[3] * q[3];
        rot[0][0] = q02 + q12 - q22 - q32;
        rot[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
        rot[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
        rot[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
        rot[1][1] = q02 - q12 + q22 - q32;
        rot[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
        rot[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
        rot[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
        rot[2][2] = q02 - q12 - q22 + q32;

        return rot;
    }

    /**
     * Minimize the RMS distance between two sets of atoms using quaternions and
     * a pre-calculated rotation matrix; overlaps x2 onto x1.
     *
     * @param x2   Cartesian coordinates of the second system. Modified in-place.
     * @param rot  A pre-calculated rotation matrix.
     */
    public static void applyRotation(double[] x2, double[][] rot) {
        int n = x2.length / 3;

        // Rotate second molecule to best fit with first molecule
        for (int i = 0; i < n; i++) {
            int k = i * 3;
            double xrot = x2[k] * rot[0][0] + x2[k + 1] * rot[0][1] + x2[k + 2] * rot[0][2];
            double yrot = x2[k] * rot[1][0] + x2[k + 1] * rot[1][1] + x2[k + 2] * rot[1][2];
            double zrot = x2[k] * rot[2][0] + x2[k + 1] * rot[2][1] + x2[k + 2] * rot[2][2];
            x2[k] = xrot;
            x2[k + 1] = yrot;
            x2[k + 2] = zrot;
        }
    }

    /**
     * Compute the rms fit over superimposed atom pairs
     *
     * @param x1   Cartesian coordinates of the first system.
     * @param x2   Cartesian coordinates of the second system.
     * @param mass The mass of each particle in the system.
     * @return The RMSD.
     */
    public static double rmsd(double[] x1, double[] x2, double[] mass) {
        double rmsfit = 0.0;
        double norm = 0.0;
        int n = x1.length / 3;
        for (int i = 0; i < n; i++) {
            int k = 3 * i;
            double weigh = mass[i];
            double xr = x1[k] - x2[k];
            double yr = x1[k + 1] - x2[k + 1];
            double zr = x1[k + 2] - x2[k + 2];
            double dist2 = xr * xr + yr * yr + zr * zr;
            norm = norm + weigh;
            double rmsterm = dist2 * weigh;
            rmsfit = rmsfit + rmsterm;
        }
        return sqrt(rmsfit / norm);
    }

}
