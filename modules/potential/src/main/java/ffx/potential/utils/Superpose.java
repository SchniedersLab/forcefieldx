package ffx.potential.utils;

import java.util.logging.Logger;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.util.FastMath;

public class Superpose {

    private static final Logger logger = Logger.getLogger(Superpose.class.getName());

    /**
     * Move the center of mass for both sets of atoms to the origin.
     *
     * @param x1
     * @param mass1
     * @param x2
     * @param mass2
     */
    public static void center(double[] x1, double mass1[], double[] x2, double mass2[]) {
        // Move the first set of atoms.
        center(x1, mass1);
        // Move the second set of atoms.
        center(x2, mass2);
    }

    /**
     * Compute the rms fit over superimposed atom pairs
     *
     * @param x1
     * @param x2
     * @param mass
     * @return
     */
    public static double rmsd(double[] x1, double[] x2, double mass[]) {
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
        return Math.sqrt(rmsfit / norm);
    }

    /**
     * Move the center of mass for a set of atoms to the origin.
     *
     * @param x
     * @param mass
     */
    private static void center(double[] x, double mass[]) {
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
        for (int i = 0; i < n; i++) {
            int k = 3 * i;
            x[k] = x[k] - xmid;
            x[k + 1] = x[k + 1] - ymid;
            x[k + 2] = x[k + 2] - zmid;
        }
    }

    /**
     * Minimize the RMS distance between two sets of atoms using quaternions.
     *
     * @param x1
     * @param x2
     * @param mass
     */
    public static void quatfit(double[] x1, double[] x2, double mass[]) {
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
        c[1][1] = xxyx - xyyy - xzyz;
        c[0][2] = xxyz - xzyx;
        c[1][2] = xxyy + xyyx;
        c[2][2] = xyyy - xzyz - xxyx;
        c[0][3] = xyyx - xxyy;
        c[1][3] = xzyx + xxyz;
        c[2][3] = xyyz + xzyy;
        c[3][3] = xzyz - xxyx - xyyy;

//        c[0][0] = xxyx + xyyy + xzyz;
//        c[1][0] = xzyy - xyyz;
//        c[1][1] = xxyx - xyyy - xzyz;
//        c[2][0] = xxyz - xzyx;
//        c[2][1] = xxyy + xyyx;
//        c[2][2] = xyyy - xzyz - xxyx;
//        c[3][0] = xyyx - xxyy;
//        c[3][1] = xzyx + xxyz;
//        c[3][2] = xyyz + xzyy;
//        c[3][3] = xzyz - xxyx - xyyy;

        logger.info(String.format(" C Matrix (%8.3f %8.3f %8.3f %8.3f)",
                c[0][0], c[1][1], c[2][2], c[3][3]));

        // Diagonalize the quadratic form matrix
        Array2DRowRealMatrix cMatrix = new Array2DRowRealMatrix(c, false);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(cMatrix);

        // EigenVectors
        RealMatrix V = eigenDecomposition.getV();

        // EigenValues (on the diagonal)
        RealMatrix D = eigenDecomposition.getD();

        // Extract the desired quaternions.
        double[] q = new double[4];
        q[0] = V.getEntry(0,3);
        q[1] = V.getEntry(1,3);
        q[2] = V.getEntry(2,3);
        q[3] = V.getEntry(3,3);

        logger.info(String.format(" Quaternion (%8.3f %8.3f %8.3f %8.3f)",
                q[0], q[1], q[2], q[3]));

        double[] e = eigenDecomposition.getRealEigenvalues();

//        double[] e = new double[4];
//        e[0] = D.getEntry(0,0);
//        e[1] = D.getEntry(1,1);
//        e[2] = D.getEntry(2,2);
//        e[3] = D.getEntry(3,3);

        logger.info(String.format(" Eigenvalues (%8.3f %8.3f %8.3f %8.3f)",
                e[0], e[1], e[2], e[3]));

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

        // Rotate second molecule to best fit with first molecule
        for (int i = 0; i < n; i++) {
            int k = i * 3;
            double xrot = x2[k] * rot[0][0] + x2[k + 1] * rot[0][0] + x2[k + 2] * rot[0][2];
            double yrot = x2[k] * rot[1][0] + x2[k + 1] * rot[1][0] + x2[k + 2] * rot[1][2];
            double zrot = x2[k] * rot[2][0] + x2[k + 1] * rot[2][0] + x2[k + 2] * rot[2][2];
            x2[k] = xrot;
            x2[k + 1] = yrot;
            x2[k + 2] = zrot;
        }
    }


    /**
     * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
     *
     * @param householderMatrix Householder matrix of the transformation
     * to tridiagonal form.
     */
    private void findEigenVectors(final double[][] householderMatrix,
                                  ArrayRealVector[] eigenvectors, ArrayRealVector eigenValues) {

        final double[][] z = householderMatrix.clone();
        final byte MAX_ITER = 30;
        final long EXPONENT_OFFSET = 1023l;
        double EPSILON = Double.longBitsToDouble((EXPONENT_OFFSET - 53l) << 52);

        final int n = householderMatrix.length;

        /** Main diagonal of the tridiagonal matrix. */
        double[] main = new double[n];
        /** Secondary diagonal of the tridiagonal matrix. */
        double[] secondary = new double[n];

        double[] realEigenvalues = new double[n];
        double[] imagEigenvalues = new double[n];
        final double[] e = new double[n];
        for (int i = 0; i < n - 1; i++) {
            realEigenvalues[i] = main[i];
            e[i] = secondary[i];
        }
        realEigenvalues[n - 1] = main[n - 1];
        e[n - 1] = 0;

        // Determine the largest main and secondary value in absolute term.
        double maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(realEigenvalues[i]);
            }
            if (FastMath.abs(e[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(e[i]);
            }
        }
        // Make null any main and secondary value too small to be significant
        if (maxAbsoluteValue != 0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(realEigenvalues[i]) <= EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
                if (FastMath.abs(e[i]) <= EPSILON * maxAbsoluteValue) {
                    e[i]=0;
                }
            }
        }

        for (int j = 0; j < n; j++) {
            int its = 0;
            int m;
            do {
                for (m = j; m < n - 1; m++) {
                    double delta = FastMath.abs(realEigenvalues[m]) +
                            FastMath.abs(realEigenvalues[m + 1]);
                    if (FastMath.abs(e[m]) + delta == delta) {
                        break;
                    }
                }
                if (m != j) {
                    if (its == MAX_ITER) {
                        logger.severe(" Eigen decomposition failed.");
                    }
                    its++;
                    double q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j]);
                    double t = FastMath.sqrt(1 + q * q);
                    if (q < 0.0) {
                        q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t);
                    } else {
                        q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t);
                    }
                    double u = 0.0;
                    double s = 1.0;
                    double c = 1.0;
                    int i;
                    for (i = m - 1; i >= j; i--) {
                        double p = s * e[i];
                        double h = c * e[i];
                        if (FastMath.abs(p) >= FastMath.abs(q)) {
                            c = q / p;
                            t = FastMath.sqrt(c * c + 1.0);
                            e[i + 1] = p * t;
                            s = 1.0 / t;
                            c *= s;
                        } else {
                            s = p / q;
                            t = FastMath.sqrt(s * s + 1.0);
                            e[i + 1] = q * t;
                            c = 1.0 / t;
                            s *= c;
                        }
                        if (e[i + 1] == 0.0) {
                            realEigenvalues[i + 1] -= u;
                            e[m] = 0.0;
                            break;
                        }
                        q = realEigenvalues[i + 1] - u;
                        t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
                        u = s * t;
                        realEigenvalues[i + 1] = q + u;
                        q = c * t - h;
                        for (int ia = 0; ia < n; ia++) {
                            p = z[ia][i + 1];
                            z[ia][i + 1] = s * z[ia][i] + c * p;
                            z[ia][i] = c * z[ia][i] - s * p;
                        }
                    }
                    if (t == 0.0 && i >= j) {
                        continue;
                    }
                    realEigenvalues[j] -= u;
                    e[j] = q;
                    e[m] = 0.0;
                }
            } while (m != j);
        }

        //Sort the eigen values (and vectors) in increase order
        for (int i = 0; i < n; i++) {
            int k = i;
            double p = realEigenvalues[i];
            for (int j = i + 1; j < n; j++) {
                if (realEigenvalues[j] > p) {
                    k = j;
                    p = realEigenvalues[j];
                }
            }
            if (k != i) {
                realEigenvalues[k] = realEigenvalues[i];
                realEigenvalues[i] = p;
                for (int j = 0; j < n; j++) {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue=FastMath.abs(realEigenvalues[i]);
            }
        }
        // Make null any eigen value too small to be significant
        if (maxAbsoluteValue != 0.0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(realEigenvalues[i]) < EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
            }
        }
        eigenvectors = new ArrayRealVector[n];
        final double[] tmp = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp[j] = z[j][i];
            }
            eigenvectors[i] = new ArrayRealVector(tmp);
        }
    }

}
