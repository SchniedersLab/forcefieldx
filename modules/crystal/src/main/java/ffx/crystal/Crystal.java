/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.crystal;

import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;

import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.scalar;

import java.util.logging.Logger;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

/**
 * The Crystal class encapsulates the lattice parameters and space group that
 * describe the geometry and symmetry of a crystal. Methods are available to
 * apply the minimum image convention and space group symmetry operators.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Crystal {

    private static final Logger logger = Logger.getLogger(Crystal.class.getName());
    /**
     * Length of the cell edge in the direction of the <b>a</b> basis vector.
     */
    public final double a;
    /**
     * Length of the cell edge in the direction of the <b>b</b> basis vector.
     */
    public final double b;
    /**
     * Length of the cell edge in the direction of the <b>c</b> basis vector.
     */
    public final double c;
    /**
     * The interaxial lattice angle between <b>b</b> and <b>c</b>.
     */
    public final double alpha;
    /**
     * The interaxial lattice angle between <b>a</b> and <b>c</b>.
     */
    public final double beta;
    /**
     * The interaxial lattice angle between <b>a</b> and <b>b</b>.
     */
    public final double gamma;
    /**
     * The space group of the crystal.
     */
    public final SpaceGroup spaceGroup;
    /**
     * A refence the space group crystal system.
     */
    private final SpaceGroup.CrystalSystem crystalSystem;
    /**
     * The crystal unit cell volume.
     */
    public final double volume;
    /**
     * Matrix to convert from Cartesian to fractional coordinates.
     */
    public final double recip[][] = new double[3][3];
    /**
     * Entries in the "recip" array.
     */
    public final double r00, r01, r02, r10, r11, r12, r20, r21, r22;
    /**
     * Matric to convert from fractional to Cartesian coordinates.
     */
    public final double cart[][] = new double[3][3];
    /**
     * Entries in the "cart" array.
     */
    public final double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    /**
     * The direct space metric matrix.
     */
    public final double G[][] = new double[3][3];
    /**
     * the reciprocal space metric matrix.
     */
    public final double Gstar[][] = new double[3][3];
    /**
     * Direct space lattice vector <b>a</b>.
     */
    private final double ad[] = new double[3];
    /**
     * Direct space lattice vector <b>b</b>.
     */
    private final double bd[] = new double[3];
    /**
     * Direct space lattice vector <b>c</b>.
     */
    private final double cd[] = new double[3];
    /**
     * Reciprocal space lattice vector <b>a*</b>.
     */
    private final double as[] = new double[3];
    /**
     * Reciprocal space lattice vector <b>b*</b>.
     */
    private final double bs[] = new double[3];
    /**
     * Reciprocal space lattice vector <b>c*</b>.
     */
    private final double cs[] = new double[3];
    private final double half_a;
    private final double half_b;
    private final double half_c;
    private final double mhalf_a;
    private final double mhalf_b;
    private final double mhalf_c;
    private final double cos_alpha;
    private final double sin_beta;
    private final double cos_beta;
    private final double sin_gamma;
    private final double cos_gamma;
    private final double beta_term;
    private final double gamma_term;
    private final double tolerance = 1.0e-10;

    /**
     * The Crystal class encapsulates the lattice parameters and space group.
     * Methods are available to apply the minimum image convention and to apply
     * space group operators.
     *
     * @param a The a-axis length.
     * @param b The b-axis length.
     * @param c The c-axis lenght.
     * @param alpha The alpha angle.
     * @param beta The beta angle.
     * @param gamma The gamma angle.
     * @param sg The space group symbol.
     */
    public Crystal(double a, double b, double c, double alpha, double beta,
            double gamma, String sg) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        half_a = 0.5 * a;
        half_b = 0.5 * b;
        half_c = 0.5 * c;
        mhalf_a = -half_a;
        mhalf_b = -half_b;
        mhalf_c = -half_c;
        spaceGroup = SpaceGroup.spaceGroupFactory(sg);
        crystalSystem = spaceGroup.crystalSystem;

        if (!SpaceGroup.checkRestrictions(crystalSystem, a, b, c, alpha, beta, gamma)) {
            String message = "The lattice parameters do not satisfy the " + crystalSystem
                    + " crystal system restrictions:/n" + toString();
            logger.severe(message);
        }

        switch (crystalSystem) {
            case CUBIC:
            case ORTHORHOMBIC:
            case TETRAGONAL:
                cos_alpha = 0.0;
                sin_beta = 1.0;
                cos_beta = 0.0;
                sin_gamma = 1.0;
                cos_gamma = 0.0;
                beta_term = 0.0;
                gamma_term = 1.0;
                volume = a * b * c;
                break;
            case MONOCLINIC:
                cos_alpha = 0.0;
                sin_beta = sin(toRadians(beta));
                cos_beta = cos(toRadians(beta));
                sin_gamma = 1.0;
                cos_gamma = 0.0;
                beta_term = 0.0;
                gamma_term = sin_beta;
                volume = sin_beta * a * b * c;
                break;
            case HEXAGONAL:
            case TRICLINIC:
            case TRIGONAL:
            default:
                cos_alpha = cos(toRadians(alpha));
                sin_beta = sin(toRadians(beta));
                cos_beta = cos(toRadians(beta));
                sin_gamma = sin(toRadians(gamma));
                cos_gamma = cos(toRadians(gamma));
                beta_term = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
                gamma_term = sqrt(sin_beta * sin_beta - beta_term * beta_term);
                volume = sin_gamma * gamma_term * a * b * c;
                break;
        }

        // a basis vector.
        ad[0] = a;
        ad[1] = 0.0;
        ad[2] = 0.0;
        cart[0][0] = ad[0];
        cart[0][1] = ad[1];
        cart[0][2] = ad[2];
        // b basis vector.
        bd[0] = b * cos_gamma;
        bd[1] = b * sin_gamma;
        bd[2] = 0.0;
        cart[1][0] = bd[0];
        cart[1][1] = bd[1];
        cart[1][2] = bd[2];
        // c basis vector.
        cd[0] = c * cos_beta;
        cd[1] = c * beta_term;
        cd[2] = c * gamma_term;
        cart[2][0] = cd[0];
        cart[2][1] = cd[1];
        cart[2][2] = cd[2];
        c00 = cart[0][0];
        c10 = cart[1][0];
        c20 = cart[2][0];
        c01 = cart[0][1];
        c11 = cart[1][1];
        c21 = cart[2][1];
        c02 = cart[0][2];
        c12 = cart[1][2];
        c22 = cart[2][2];

        // a* basis vector = b x c / volume.
        cross(bd, cd, as);
        scalar(as, 1.0 / volume, as);
        recip[0][0] = as[0];
        recip[1][0] = as[1];
        recip[2][0] = as[2];
        // b* basis vector = c x a / volume.
        cross(cd, ad, bs);
        scalar(bs, 1.0 / volume, bs);
        recip[0][1] = bs[0];
        recip[1][1] = bs[1];
        recip[2][1] = bs[2];
        // c* basis vector = a x b / volume.
        cross(ad, bd, cs);
        scalar(cs, 1.0 / volume, cs);
        recip[0][2] = cs[0];
        recip[1][2] = cs[1];
        recip[2][2] = cs[2];
        r00 = recip[0][0];
        r01 = recip[0][1];
        r02 = recip[0][2];
        r10 = recip[1][0];
        r11 = recip[1][1];
        r12 = recip[1][2];
        r20 = recip[2][0];
        r21 = recip[2][1];
        r22 = recip[2][2];

        assert (abs(dot(ad, as) - 1.0) < tolerance);
        assert (abs(dot(bd, bs) - 1.0) < tolerance);
        assert (abs(dot(cd, cs) - 1.0) < tolerance);
        assert (abs(dot(ad, bs)) < tolerance);
        assert (abs(dot(ad, cs)) < tolerance);
        assert (abs(dot(bd, as)) < tolerance);
        assert (abs(dot(bd, cs)) < tolerance);
        assert (abs(dot(cd, as)) < tolerance);
        assert (abs(dot(cd, bs)) < tolerance);

        G[0][0] = a * a;
        G[0][1] = a * b * cos_gamma;
        G[0][2] = a * c * cos_beta;
        G[1][0] = G[0][1];
        G[1][1] = b * b;
        G[1][2] = b * c * cos_alpha;
        G[2][0] = G[0][2];
        G[2][1] = G[1][2];
        G[2][2] = c * c;

        // invert G to yield Gstar
        DenseMatrix A = new DenseMatrix(G);
        DenseMatrix I = Matrices.identity(3);
        DenseMatrix AI = I.copy();
        A.solve(I, AI);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Gstar[i][j] = AI.get(i, j);
            }
        }

    }

    /**
     * Apply the minimum image convention.
     *
     * @param xyz input distances that are over-written.
     * @return the output distance squared.
     */
    public double image(final double xyz[]) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];
        switch (crystalSystem) {
            case CUBIC:
            case ORTHORHOMBIC:
            case TETRAGONAL:
                while (x > half_a) {
                    x -= a;
                }
                while (x < mhalf_a) {
                    x += a;
                }
                while (y > half_b) {
                    y -= b;
                }
                while (y < mhalf_b) {
                    y += b;
                }
                while (z > half_c) {
                    z -= c;
                }
                while (z < mhalf_c) {
                    z += c;
                }
                break;
            case MONOCLINIC:
                double zm = z / sin_beta;
                double xm = x - zm * cos_beta;
                while (xm > half_a) {
                    xm -= a;
                }
                while (xm < mhalf_a) {
                    xm += a;
                }
                while (y > half_b) {
                    y -= b;
                }
                while (y < mhalf_b) {
                    y += b;
                }
                while (zm > half_c) {
                    zm -= c;
                }
                while (zm < mhalf_c) {
                    zm += c;
                }
                x = xm + zm * cos_beta;
                z = zm * sin_beta;
                break;
            case HEXAGONAL:
            case TRICLINIC:
            case TRIGONAL:
                double zt = z / gamma_term;
                double yt = (y - zt * beta_term) / sin_gamma;
                double xt = x - yt * cos_gamma - zt * cos_beta;
                while (xt > half_a) {
                    xt -= a;
                }
                while (xt < mhalf_a) {
                    xt += a;
                }
                while (yt > half_b) {
                    yt -= b;
                }
                while (yt < mhalf_b) {
                    yt += b;
                }
                while (zt > half_c) {
                    zt -= c;
                }
                while (zt < mhalf_c) {
                    zt += c;
                }
                x = xt + yt * cos_gamma + zt * cos_beta;
                y = yt * sin_gamma + zt * beta_term;
                z = zt * gamma_term;
        }
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        return x * x + y * y + z * z;
    }

    /**
     * Apply the minimum image convention.
     *
     * @param dx x-distance
     * @param dy y-distance
     * @param dz z-distance
     * @return the output distance squared.
     */
    public double image(double dx, double dy, double dz) {
        switch (crystalSystem) {
            case CUBIC:
            case ORTHORHOMBIC:
            case TETRAGONAL:
                while (dx > half_a) {
                    dx -= a;
                }
                while (dx < mhalf_a) {
                    dx += a;
                }
                while (dy > half_b) {
                    dy -= b;
                }
                while (dy < mhalf_b) {
                    dy += b;
                }
                while (dz > half_c) {
                    dz -= c;
                }
                while (dz < mhalf_c) {
                    dz += c;
                }
                break;
            case MONOCLINIC:
                double zm = dz / sin_beta;
                double xm = dx - zm * cos_beta;
                while (xm > half_a) {
                    xm -= a;
                }
                while (xm < mhalf_a) {
                    xm += a;
                }
                while (dy > half_b) {
                    dy -= b;
                }
                while (dy < mhalf_b) {
                    dy += b;
                }
                while (zm > half_c) {
                    zm -= c;
                }
                while (zm < mhalf_c) {
                    zm += c;
                }
                dx = xm + zm * cos_beta;
                dz = zm * sin_beta;
                break;
            case HEXAGONAL:
            case TRICLINIC:
            case TRIGONAL:
                double zt = dz / gamma_term;
                double yt = (dy - zt * beta_term) / sin_gamma;
                double xt = dx - yt * cos_gamma - zt * cos_beta;
                while (xt > half_a) {
                    xt -= a;
                }
                while (xt < mhalf_a) {
                    xt += a;
                }
                while (yt > half_b) {
                    yt -= b;
                }
                while (yt < mhalf_b) {
                    yt += b;
                }
                while (zt > half_c) {
                    zt -= c;
                }
                while (zt < mhalf_c) {
                    zt += c;
                }
                dx = xt + yt * cos_gamma + zt * cos_beta;
                dy = yt * sin_gamma + zt * beta_term;
                dz = zt * gamma_term;
        }
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * Apply a symmetry operator to an array of Cartesian coordinates. If the
     * arrays x, y or z are null or not of length n, the method returns
     * immediately. If xs, ys or zs are null or not of length n, new arrays are
     * allocated.
     *
     * @param n
     *            Number of atoms.
     * @param x
     *            Input x coordinates.
     * @param y
     *            Input y coordinates.
     * @param z
     *            Input z coordinates.
     * @param mateX
     *            Output x coordinates.
     * @param mateY
     *            Output y coordinates.
     * @param mateZ
     *            Output z coordinates.
     * @param symOp
     *            The symmetry operator.
     */
    public void applySymOp(int n, double x[], double y[], double z[],
            double mateX[], double mateY[], double mateZ[], SymOp symOp) {
        final double rot[][] = symOp.rot;
        final double trans[] = symOp.tr;
        final double rot00 = rot[0][0];
        final double rot10 = rot[1][0];
        final double rot20 = rot[2][0];
        final double rot01 = rot[0][1];
        final double rot11 = rot[1][1];
        final double rot21 = rot[2][1];
        final double rot02 = rot[0][2];
        final double rot12 = rot[1][2];
        final double rot22 = rot[2][2];
        final double t0 = trans[0];
        final double t1 = trans[1];
        final double t2 = trans[2];
        for (int i = 0; i < n; i++) {
            double xc = x[i];
            double yc = y[i];
            double zc = z[i];
            // Convert to fractional coordinates.
            double xi = xc * r00 + yc * r10 + zc * r20;
            double yi = xc * r01 + yc * r11 + zc * r21;
            double zi = xc * r02 + yc * r12 + zc * r22;
            // Apply Symmetry Operator.
            double fx = rot00 * xi + rot01 * yi + rot02 * zi + t0;
            double fy = rot10 * xi + rot11 * yi + rot12 * zi + t1;
            double fz = rot20 * xi + rot21 * yi + rot22 * zi + t2;
            // Convert back to Cartesian coordinates.
            mateX[i] = fx * c00 + fy * c10 + fz * c20;
            mateY[i] = fx * c01 + fy * c11 + fz * c21;
            mateZ[i] = fx * c02 + fy * c12 + fz * c22;
        }
    }

    /**
     * Apply a symmetry operator to one set of coordinates.
     *
     * @param xyz   Input coordinates.
     * @param mate  Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     */
    public void applySymOp(double xyz[], double mate[], SymOp symOp) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        // Convert to fractional coordinates.
        double xc = xyz[0];
        double yc = xyz[1];
        double zc = xyz[2];
        // Convert to fractional coordinates.
        double xi = xc * r00 + yc * r10 + zc * r20;
        double yi = xc * r01 + yc * r11 + zc * r21;
        double zi = xc * r02 + yc * r12 + zc * r22;
        // Apply Symmetry Operator.
        double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi + trans[0];
        double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi + trans[1];
        double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi + trans[2];
        // Convert back to Cartesian coordinates.
        mate[0] = fx * c00 + fy * c10 + fz * c20;
        mate[1] = fx * c01 + fy * c11 + fz * c21;
        mate[2] = fx * c02 + fy * c12 + fz * c22;
    }

    /**
     * Apply a symmetry operator to one HKL.
     *
     * @param hkl
     *            Input HKL.
     * @param mate
     *            Symmetry mate HKL.
     * @param symOp
     *            The symmetry operator.
     */
    public void applySymOp(HKL hkl, HKL mate, SymOp symOp) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        double h = hkl.h();
        double k = hkl.k();
        double l = hkl.l();
        // Apply Symmetry Operator.
        double hs = rot[0][0] * h + rot[0][1] * k + rot[0][2] * l + trans[0];
        double ks = rot[1][0] * h + rot[1][1] * k + rot[1][2] * l + trans[1];
        double ls = rot[2][0] * h + rot[2][1] * k + rot[2][2] * l + trans[2];
        // Convert back to HKL
        mate.h((int) hs);
        mate.k((int) ks);
        mate.l((int) ls);
    }

    /**
     * Apply a symmetry operator to one set of coordinates.
     *
     * @param xyz
     *            Input coordinates.
     * @param mate
     *            Symmetry mate coordinates.
     * @param symOp
     *            The symmetry operator.
     */
    public void applySymRot(double xyz[], double mate[], SymOp symOp) {
        double rot[][] = symOp.rot;
        // Convert to fractional coordinates.
        double xc = xyz[0];
        double yc = xyz[1];
        double zc = xyz[2];
        // Convert to fractional coordinates.
        double xi = xc * r00 + yc * r10 + zc * r20;
        double yi = xc * r01 + yc * r11 + zc * r21;
        double zi = xc * r02 + yc * r12 + zc * r22;
        // Apply Symmetry Operator.
        double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi;
        double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi;
        double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi;
        // Convert back to Cartesian coordinates.
        mate[0] = fx * c00 + fy * c10 + fz * c20;
        mate[1] = fx * c01 + fy * c11 + fz * c21;
        mate[2] = fx * c02 + fy * c12 + fz * c22;
    }

    /**
     * Apply a symmetry operator to one HKL.
     *
     * @param hkl
     *            Input HKL.
     * @param mate
     *            Symmetry mate HKL.
     * @param symOp
     *            The symmetry operator.
     */
    public void applySymRot(HKL hkl, HKL mate, SymOp symOp) {
        double rot[][] = symOp.rot;
        double h = hkl.h();
        double k = hkl.k();
        double l = hkl.l();
        // Apply Symmetry Operator.
        double hs = rot[0][0] * h + rot[0][1] * k + rot[0][2] * l;
        double ks = rot[1][0] * h + rot[1][1] * k + rot[1][2] * l;
        double ls = rot[2][0] * h + rot[2][1] * k + rot[2][2] * l;
        // Convert back to HKL
        mate.h((int) hs);
        mate.k((int) ks);
        mate.l((int) ls);
    }

    /**
     * Apply a symmetry operator to an array of Cartesian coordinates. If the
     * arrays x, y or z are null or not of length n, the method returns
     * immediately. If xs, ys or zs are null or not of length n, new arrays are
     * allocated.
     *
     * @param n
     *            Number of atoms.
     * @param x
     *            Input x coordinates.
     * @param y
     *            Input y coordinates.
     * @param z
     *            Input z coordinates.
     * @param mateX
     *            Output x coordinates.
     * @param mateY
     *            Output y coordinates.
     * @param mateZ
     *            Output z coordinates.
     * @param symOp
     *            The symmetry operator.
     */
    public void applySymRot(int n, double x[], double y[], double z[],
            double mateX[], double mateY[], double mateZ[], SymOp symOp) {
        double rot[][] = symOp.rot;
        final double rot00 = rot[0][0];
        final double rot10 = rot[1][0];
        final double rot20 = rot[2][0];
        final double rot01 = rot[0][1];
        final double rot11 = rot[1][1];
        final double rot21 = rot[2][1];
        final double rot02 = rot[0][2];
        final double rot12 = rot[1][2];
        final double rot22 = rot[2][2];
        for (int i = 0; i < n; i++) {
            // Convert to fractional coordinates.
            double xc = x[i];
            double yc = y[i];
            double zc = z[i];
            // Convert to fractional coordinates.
            double xi = xc * r00 + yc * r10 + zc * r20;
            double yi = xc * r01 + yc * r11 + zc * r21;
            double zi = xc * r02 + yc * r12 + zc * r22;
            // Apply Symmetry Operator.
            double fx = rot00 * xi + rot01 * yi + rot02 * zi;
            double fy = rot10 * xi + rot11 * yi + rot12 * zi;
            double fz = rot20 * xi + rot21 * yi + rot22 * zi;
            // Convert back to Cartesian coordinates.
            mateX[i] = fx * c00 + fy * c10 + fz * c20;
            mateY[i] = fx * c01 + fy * c11 + fz * c21;
            mateZ[i] = fx * c02 + fy * c12 + fz * c22;
        }
    }

    public void toFractionalCoordinates(int n, double x[], double y[],
            double z[], double xf[], double yf[], double zf[]) {
        for (int i = 0; i < n; i++) {
            double xc = x[i];
            double yc = y[i];
            double zc = z[i];
            xf[i] = xc * r00 + yc * r10 + zc * r20;
            yf[i] = xc * r01 + yc * r11 + zc * r21;
            zf[i] = xc * r02 + yc * r12 + zc * r22;
        }
    }

    public void toFractionalCoordinates(double x, double y,
            double z, double xf, double yf, double zf) {
        xf = x * r00 + y * r10 + z * r20;
        yf = x * r01 + y * r11 + z * r21;
        zf = x * r02 + y * r12 + z * r22;
    }

    public void toFractionalCoordinates(int n, double cart[], double frac[]) {
        int i3 = 0;
        for (int i = 0; i < n; i++) {
            // Convert to fractional coordinates.
            int iX = i3 + XX;
            int iY = i3 + YY;
            int iZ = i3 + ZZ;
            i3 += 3;
            double xc = cart[iX];
            double yc = cart[iY];
            double zc = cart[iZ];
            frac[iX] = xc * r00 + yc * r10 + zc * r20;
            frac[iY] = xc * r01 + yc * r11 + zc * r21;
            frac[iZ] = xc * r02 + yc * r12 + zc * r22;
        }
    }

    public void toFractionalCoordinates(double x[], double xf[]) {
        double xc = x[0];
        double yc = x[1];
        double zc = x[2];
        xf[0] = xc * r00 + yc * r10 + zc * r20;
        xf[1] = xc * r01 + yc * r11 + zc * r21;
        xf[2] = xc * r02 + yc * r12 + zc * r22;

    }

    public void toCartesianCoordinates(int n, double xf[], double yf[],
            double zf[], double x[], double y[], double z[]) {
        for (int i = 0; i < n; i++) {
            double xi = xf[i];
            double yi = yf[i];
            double zi = zf[i];
            x[i] = xi * c00 + yi * c10 + zi * c20;
            y[i] = xi * c01 + yi * c11 + zi * c21;
            z[i] = xi * c02 + yi * c12 + zi * c22;
        }
    }

    public void toCartesianCoordinates(double xf, double yf,
            double zf, double x, double y, double z) {
        x = xf * c00 + yf * c10 + zf * c20;
        y = xf * c01 + yf * c11 + zf * c21;
        z = xf * c02 + yf * c12 + zf * c22;
    }

    public void toCartesianCoordinates(int n, double frac[], double cart[]) {
        int i3 = 0;
        for (int i = 0; i < n; i++) {
            // Convert to cartesian coordinates.
            int iX = i3 + XX;
            int iY = i3 + YY;
            int iZ = i3 + ZZ;
            i3 += 3;
            double xf = frac[iX];
            double yf = frac[iY];
            double zf = frac[iZ];
            cart[iX] = xf * c00 + yf * c10 + zf * c20;
            cart[iY] = xf * c01 + yf * c11 + zf * c21;
            cart[iZ] = xf * c02 + yf * c12 + zf * c22;
        }
    }

    public void toCartesianCoordinates(double xf[], double x[]) {
        double fx = xf[0];
        double fy = xf[1];
        double fz = xf[2];
        x[0] = fx * c00 + fy * c10 + fz * c20;
        x[1] = fx * c01 + fy * c11 + fz * c21;
        x[2] = fx * c02 + fy * c12 + fz * c22;
    }

    public void toFractionalCoordinatesTINKER(double x[], double xf[]) {
        // Convert to fractional coordinates.
        xf[2] = (x[2] / gamma_term) / c;
        xf[1] = ((x[1] - xf[2] * c * beta_term) / sin_gamma) / b;
        xf[0] = (x[0] - xf[1] * b * cos_gamma - xf[2] * c * cos_beta) / a;
    }

    public static double quad_form(double v[], double mat[][]) {
        return (v[0] * (v[0] * mat[0][0] + 2 * (v[1] * mat[0][1] + v[2] * mat[0][2]))
                + v[1] * (v[1] * mat[1][1] + 2 * (v[2] * mat[1][2]))
                + v[2] * v[2] * mat[2][2]);
    }

    public static double quad_form(HKL hkl, double mat[][]) {
        return (hkl.h() * (hkl.h() * mat[0][0] + 2 * (hkl.k() * mat[0][1] + hkl.l() * mat[0][2]))
                + hkl.k() * (hkl.k() * mat[1][1] + 2 * (hkl.l() * mat[1][2]))
                + hkl.l() * hkl.l() * mat[2][2]);
    }

    public static double invressq(Crystal crystal, HKL hkl) {
        return quad_form(hkl, crystal.Gstar);
    }

    public static double res(Crystal crystal, HKL hkl) {
        return 1.0 / sqrt(quad_form(hkl, crystal.Gstar));
    }

    public static double sym_phase_shift(HKL hkl, SymOp symOp) {
        double trans[] = symOp.tr;
        // Apply translation
        return -2.0 * PI
                * (hkl.h() * trans[0] + hkl.k() * trans[1] + hkl.l() * trans[2]);
    }

    // this is here as its an atypical mod function used by xtal methods
    public static double mod(double a, double b) {
        double res = a % b;
        if (res < 0.0) {
            res += b;
        }
        return res;
    }

    @Override
    public String toString() {
        return String.format(
                "\n Unit cell: (%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f)\n", a,
                b, c, alpha, beta, gamma)
                + String.format(
                " Space group number %d (%s, Nsymm = %d)\n",
                spaceGroup.number,
                spaceGroup.shortName,
                spaceGroup.numSymEquiv);
    }
    private static final int XX = 0;
    private static final int YY = 1;
    private static final int ZZ = 2;
}
