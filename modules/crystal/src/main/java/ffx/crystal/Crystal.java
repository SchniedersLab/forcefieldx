/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import java.util.Vector;
import static java.lang.Math.*;

import ffx.utilities.HashCodeUtil;
import static ffx.numerics.VectorMath.mat3mat3;
import static ffx.numerics.VectorMath.mat3symvec6;
import static ffx.numerics.VectorMath.transpose3;

import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.util.MathUtils;

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
     * reference to the space group crystal system.
     */
    private final SpaceGroup.CrystalSystem crystalSystem;
    /**
     * Copy of symmetry operators in Cartesian coordinates.
     */
    private final Vector<SymOp> symOpsCartesian;
    /**
     * The crystal unit cell volume.
     */
    public final double volume;
    /**
     * Matrix to convert from fractional to Cartesian coordinates.
     */
    public final double Ai[][] = new double[3][3];
    /**
     * Entries in the Ai array.
     */
    public final double Ai00, Ai01, Ai02, Ai10, Ai11, Ai12, Ai20, Ai21, Ai22;
    /**
     * Matrix to convert from Cartesian to fractional coordinates.
     */
    public final double A[][];
    /**
     * Entries in the A array.
     */
    public final double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    /**
     * The direct space metric matrix.
     */
    public final double G[][] = new double[3][3];
    /**
     * the reciprocal space metric matrix.
     */
    public final double Gstar[][];
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
    private boolean aperiodic;
    public int scale_flag;
    public int scale_b[] = new int[6];
    public int scale_n;

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
        aperiodic = false;
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

        for (int i = 0; i < 6; i++) {
            scale_b[i] = -1;
        }
        SymOp symop;
        double rot[][];
        int index = 0;
        switch (crystalSystem) {
            case TRICLINIC:
                for (int i = 0; i < 6; i++) {
                    scale_b[i] = index++;
                }
                break;
            case MONOCLINIC:
                index = 0;
                scale_b[0] = index++;
                scale_b[1] = index++;
                scale_b[2] = index++;
                // determine unique axis
                symop = spaceGroup.symOps.get(1);
                rot = symop.rot;
                if (rot[0][0] > 0) {
                    scale_b[5] = index++;
                } else if (rot[1][1] > 0) {
                    scale_b[4] = index++;
                } else {
                    scale_b[3] = index++;
                }
                break;
            case ORTHORHOMBIC:
                index = 0;
                scale_b[0] = index++;
                scale_b[1] = index++;
                scale_b[2] = index++;
                break;
            case TETRAGONAL:
                index = 0;
                scale_b[0] = index++;
                scale_b[1] = scale_b[0];
                scale_b[2] = index++;
                break;
            case TRIGONAL:
            case HEXAGONAL:
                boolean hexagonal = false;
                for (int i = 1; i < spaceGroup.symOps.size(); i++) {
                    symop = spaceGroup.symOps.get(i);
                    rot = symop.rot;

                    index = 0;
                    if ((rot[1][1] * rot[1][2] == -1)
                            || (rot[2][1] * rot[2][2] == -1)
                            || (rot[1][1] * rot[1][2] == 1)
                            || (rot[2][1] * rot[2][2] == 1)) {
                        scale_b[0] = index++;
                        scale_b[1] = index++;
                        scale_b[2] = scale_b[1];
                        hexagonal = true;
                    } else if ((rot[0][0] * rot[0][2] == -1)
                            || (rot[2][0] * rot[2][2] == -1)
                            || (rot[0][0] * rot[0][2] == 1)
                            || (rot[2][0] * rot[2][2] == 1)) {
                        scale_b[0] = index++;
                        scale_b[1] = index++;
                        scale_b[2] = scale_b[0];
                        hexagonal = true;
                    } else if ((rot[0][0] * rot[0][1] == -1)
                            || (rot[1][0] * rot[1][1] == -1)
                            || (rot[0][0] * rot[0][1] == 1)
                            || (rot[1][0] * rot[1][1] == 1)) {
                        scale_b[0] = index++;
                        scale_b[1] = scale_b[0];
                        scale_b[2] = index++;
                        hexagonal = true;
                    }

                    if (hexagonal) {
                        break;
                    }
                }
                if (!hexagonal) {
                    // rhombohedral
                    index = 0;
                    scale_b[3] = index++;
                    scale_b[4] = scale_b[3];
                    scale_b[5] = scale_b[3];
                }
                break;
            case CUBIC:
                break;
        }
        scale_n = index;

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

        // a is the first row of A^(-1).
        Ai[0][0] = a;
        Ai[0][1] = 0.0;
        Ai[0][2] = 0.0;
        // b is the second row of A^(-1).
        Ai[1][0] = b * cos_gamma;
        Ai[1][1] = b * sin_gamma;
        Ai[1][2] = 0.0;
        // c is the third row of A^(-1).
        Ai[2][0] = c * cos_beta;
        Ai[2][1] = c * beta_term;
        Ai[2][2] = c * gamma_term;
        Ai00 = Ai[0][0];
        Ai10 = Ai[1][0];
        Ai20 = Ai[2][0];
        Ai01 = Ai[0][1];
        Ai11 = Ai[1][1];
        Ai21 = Ai[2][1];
        Ai02 = Ai[0][2];
        Ai12 = Ai[1][2];
        Ai22 = Ai[2][2];

        // Invert A^-1 to get A
        RealMatrix m = new Array2DRowRealMatrix(Ai, true);
        m = new LUDecompositionImpl(m).getSolver().getInverse();
        A = m.getData();
        // The columns of A are the reciprocal basis vectors
        A00 = A[0][0];
        A01 = A[0][1];
        A02 = A[0][2];
        A10 = A[1][0];
        A11 = A[1][1];
        A12 = A[1][2];
        A20 = A[2][0];
        A21 = A[2][1];
        A22 = A[2][2];

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
        m = new Array2DRowRealMatrix(G, true);
        m = new LUDecompositionImpl(m).getSolver().getInverse();
        Gstar = m.getData();

        symOpsCartesian = new Vector<SymOp>();
        Vector<SymOp> symOps = spaceGroup.symOps;
        int nSymm = symOps.size();
        RealMatrix toFrac = new Array2DRowRealMatrix(A);
        RealMatrix toCart = new Array2DRowRealMatrix(Ai);
        for (int i = 0; i < nSymm; i++) {
            SymOp symOp = symOps.get(i);
            m = new Array2DRowRealMatrix(symOp.rot);
            // rot_c = A^(-1).rot_f.A
            RealMatrix rotMat = m.preMultiply(toCart.transpose()).multiply(toFrac.transpose());
            // tr_c = tr_f.A^(-1)
            double tr[] = toCart.preMultiply(symOp.tr);
            symOpsCartesian.add(new SymOp(rotMat.getData(), tr));
        }

        /**
        double temp[] = {1.0/3.0, 1.0/7.0, 1.0/9.0};
        double ret[] = {0,0,0};
        for (int i=0; i < nSymm; i++) {
        SymOp symOp = symOps.get(i);
        SymOp symOpCart = symOpsCartesian.get(i);
        applySymOp(temp, ret, symOp);
        logger.info(String.format("Applied a fractional symmetry operator: %8.3f, %8.3f, %8.3f",
        ret[0], ret[1], ret[2]));
        applyCartesianSymOp(temp, ret, symOpCart);
        logger.info(String.format("Applied a Cartesian symmetry operator:  %8.3f, %8.3f, %8.3f",
        ret[0], ret[1], ret[2]));
        } */
    }

    public static Crystal checkProperties(CompositeConfiguration properties) {
        double a = properties.getDouble("a-axis", -1.0);
        double b = properties.getDouble("b-axis", -1.0);
        double c = properties.getDouble("c-axis", -1.0);
        double alpha = properties.getDouble("alpha", -1.0);
        double beta = properties.getDouble("beta", -1.0);
        double gamma = properties.getDouble("gamma", -1.0);
        String sg = properties.getString("spacegroup", null);

        sg = SpaceGroup.pdb2ShortName(sg);

        if (a < 0.0 || b < 0.0 || c < 0.0
                || alpha < 0.0 || beta < 0.0 || gamma < 0.0
                || sg == null) {
            return null;
        }

        // check the space group name is valid
        SpaceGroup spaceGroup = SpaceGroup.spaceGroupFactory(sg);
        if (spaceGroup == null) {
            sg = sg.replaceAll(" ", "");
            spaceGroup = SpaceGroup.spaceGroupFactory(sg);
            if (spaceGroup == null) {
                return null;
            }
        }

        return new Crystal(a, b, c, alpha, beta, gamma, sg);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof Crystal)) {
            return false;
        }
        if (this == obj) {
            return true;
        }
        return (MathUtils.equals(this.a, ((Crystal) obj).a, 0.001)
                && MathUtils.equals(this.b, ((Crystal) obj).b, 0.001)
                && MathUtils.equals(this.c, ((Crystal) obj).c, 0.001)
                && MathUtils.equals(this.alpha, ((Crystal) obj).alpha, 0.01)
                && MathUtils.equals(this.beta, ((Crystal) obj).beta, 0.01)
                && MathUtils.equals(this.gamma, ((Crystal) obj).gamma, 0.01)
                && this.spaceGroup.number == ((Crystal) obj).spaceGroup.number);
    }

    @Override
    public int hashCode() {
        int hash = HashCodeUtil.SEED;
        hash = HashCodeUtil.hash(hash, this.a);
        hash = HashCodeUtil.hash(hash, this.b);
        hash = HashCodeUtil.hash(hash, this.c);
        hash = HashCodeUtil.hash(hash, this.alpha);
        hash = HashCodeUtil.hash(hash, this.beta);
        hash = HashCodeUtil.hash(hash, this.gamma);
        hash = HashCodeUtil.hash(hash, this.spaceGroup.number);
        return hash;
    }

    /**
     * The ReplicatesCrystal over-rides this method to return the unit cell
     * rather than the ReplicateCell.
     * 
     * @return The unit cell Crystal instance.
     */
    public Crystal getUnitCell() {
        return this;
    }

    /**
     * Is this a finite system - ie. one unit cell in isolation?
     */
    public void setAperiodic(boolean aperiodic) {
        this.aperiodic = aperiodic;
    }

    public boolean aperiodic() {
        return aperiodic;
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
        if (aperiodic) {
            return x * x + y * y + z * z;
        }

        double xf = x * A00 + y * A10 + z * A20;
        double yf = x * A01 + y * A11 + z * A21;
        double zf = x * A02 + y * A12 + z * A22;

        /**
        XTEMP=INT(ABS(XTEMP)+HALF)*SIGN(ONE,-XTEMP) +XTEMP
        YTEMP=INT(ABS(YTEMP)+HALF)*SIGN(ONE,-YTEMP) +YTEMP
        ZTEMP=INT(ABS(ZTEMP)+HALF)*SIGN(ONE,-ZTEMP) +ZTEMP
         */
        xf = floor(abs(xf) + 0.5) * signum(-xf) + xf;
        yf = floor(abs(yf) + 0.5) * signum(-yf) + yf;
        zf = floor(abs(zf) + 0.5) * signum(-zf) + zf;

        x = xf * Ai00 + yf * Ai10 + zf * Ai20;
        y = xf * Ai01 + yf * Ai11 + zf * Ai21;
        z = xf * Ai02 + yf * Ai12 + zf * Ai22;

        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        return x * x + y * y + z * z;

        /*
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
        return x * x + y * y + z * z; */
    }

    public void rotateForce() {
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
        if (aperiodic) {
            return dx * dx + dy * dy + dz * dz;
        }
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

    public void averageTensor(double m[][], double r[][]) {
        int n = spaceGroup.symOps.size();
        for (int i = 0; i < n; i++) {
            SymOp symop = spaceGroup.symOps.get(i);
            double rot[][] = symop.rot;
            double rt[][] = transpose3(rot);
            double rmrt[][] = mat3mat3(mat3mat3(rot, m), rt);
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    r[j][k] += rmrt[j][k];
                }
            }
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                r[j][k] /= 6.0;
            }
        }
    }

    public void averageTensor(double v[], double r[][]) {
        int n = spaceGroup.symOps.size();
        for (int i = 0; i < n; i++) {
            SymOp symop = spaceGroup.symOps.get(i);
            double rot[][] = symop.rot;
            double rt[][] = transpose3(rot);
            double rmrt[][] = mat3mat3(mat3symvec6(rot, v), rt);
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    r[j][k] += rmrt[j][k];
                }
            }
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                r[j][k] /= 6.0;
            }
        }
    }

    /**
     * Apply a symmetry operator to an array of Cartesian coordinates. If the
     * arrays x, y or z are null or not of length n, the method returns
     * immediately. If mateX, mateY or mateZ are null or not of length n,
     * new arrays are allocated.
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
        if (x == null || y == null || z == null) {
            return;
        }
        if (x.length < n || y.length < n || z.length < n) {
            return;
        }
        if (mateX == null || mateX.length < n) {
            mateX = new double[n];
        }
        if (mateY == null || mateY.length < n) {
            mateY = new double[n];
        }
        if (mateZ == null || mateZ.length < n) {
            mateZ = new double[n];
        }

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
            double xi = xc * A00 + yc * A10 + zc * A20;
            double yi = xc * A01 + yc * A11 + zc * A21;
            double zi = xc * A02 + yc * A12 + zc * A22;
            // Apply Symmetry Operator.
            double fx = rot00 * xi + rot01 * yi + rot02 * zi + t0;
            double fy = rot10 * xi + rot11 * yi + rot12 * zi + t1;
            double fz = rot20 * xi + rot21 * yi + rot22 * zi + t2;
            // Convert back to Cartesian coordinates.
            mateX[i] = fx * Ai00 + fy * Ai10 + fz * Ai20;
            mateY[i] = fx * Ai01 + fy * Ai11 + fz * Ai21;
            mateZ[i] = fx * Ai02 + fy * Ai12 + fz * Ai22;
        }
    }

    /**
     * Apply a symmetry operator to one set of coordinates.
     *
     * @param h   Input coordinates.
     * @param k   Input coordinates.
     * @param l   Input coordinates.
     * @param mate  Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     * @param nx number of unit cell translations
     * @param ny number of unit cell translations
     * @param nz number of unit cell translations
     */
    public void applySymOp(int h, int k, int l, int mate[], SymOp symOp, int nx, int ny, int nz) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        // Apply Symmetry Operator.
        mate[0] = (int) rot[0][0] * h + (int) rot[0][1] * k + (int) rot[0][2] * l + (int) rint(nx * trans[0]);
        mate[1] = (int) rot[1][0] * h + (int) rot[1][1] * k + (int) rot[1][2] * l + (int) rint(ny * trans[1]);
        mate[2] = (int) rot[2][0] * h + (int) rot[2][1] * k + (int) rot[2][2] * l + (int) rint(nz * trans[2]);
        mate[0] = mod(mate[0], nx);
        mate[1] = mod(mate[1], ny);
        mate[2] = mod(mate[2], nz);
    }

    /**
     * Apply a fractional symmetry operator to one set of coordinates.
     *
     * @param xyz   Input coordinates.
     * @param mate  Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     */
    public void applyCartesianSymOp(double xyz[], double mate[], SymOp symOp) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        double xc = xyz[0];
        double yc = xyz[1];
        double zc = xyz[2];
        // Apply Symmetry Operator.
        mate[0] = rot[0][0] * xc + rot[0][1] * yc + rot[0][2] * zc + trans[0];
        mate[1] = rot[1][0] * xc + rot[1][1] * yc + rot[1][2] * zc + trans[1];
        mate[2] = rot[2][0] * xc + rot[2][1] * yc + rot[2][2] * zc + trans[2];
    }

    /**
     * Apply a fractional symmetry operator to one set of coordinates.
     *
     * @param xyz   Input coordinates.
     * @param mate  Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     */
    public void applySymOp(double xyz[], double mate[], SymOp symOp) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        double xc = xyz[0];
        double yc = xyz[1];
        double zc = xyz[2];
        // Convert to fractional coordinates.
        double xi = xc * A00 + yc * A10 + zc * A20;
        double yi = xc * A01 + yc * A11 + zc * A21;
        double zi = xc * A02 + yc * A12 + zc * A22;
        // Apply Symmetry Operator.
        double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi + trans[0];
        double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi + trans[1];
        double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi + trans[2];
        // Convert back to Cartesian coordinates.
        mate[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
        mate[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
        mate[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
    }

    /**
     * Apply a symmetry operator to one set of coordinates.
     *
     * @param xyz   Input coordinates.
     * @param mate  Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     */
    public void applyFracSymOp(double xyz[], double mate[], SymOp symOp) {
        double rot[][] = symOp.rot;
        double trans[] = symOp.tr;
        double xi = xyz[0];
        double yi = xyz[1];
        double zi = xyz[2];
        // Apply Symmetry Operator.
        double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi + trans[0];
        double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi + trans[1];
        double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi + trans[2];
        mate[0] = fx;
        mate[1] = fy;
        mate[2] = fz;
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
        double xi = xc * A00 + yc * A10 + zc * A20;
        double yi = xc * A01 + yc * A11 + zc * A21;
        double zi = xc * A02 + yc * A12 + zc * A22;

        // Apply Symmetry Operator.
        double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi;
        double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi;
        double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi;

        // Convert back to Cartesian coordinates.
        mate[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
        mate[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
        mate[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
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
        double hs = rot[0][0] * h + rot[0][1] * k + rot[0][2] * l;
        double ks = rot[1][0] * h + rot[1][1] * k + rot[1][2] * l;
        double ls = rot[2][0] * h + rot[2][1] * k + rot[2][2] * l;
        // Convert back to HKL
        mate.h((int) rint(hs));
        mate.k((int) rint(ks));
        mate.l((int) rint(ls));
    }

    /**
     * Apply a transpose rotation symmetry operator to one HKL.
     *
     * @param hkl
     *            Input HKL.
     * @param mate
     *            Symmetry mate HKL.
     * @param symOp
     *            The symmetry operator.
     */
    public void applyTransSymRot(HKL hkl, HKL mate, SymOp symOp) {
        double rot[][] = symOp.rot;
        double h = hkl.h();
        double k = hkl.k();
        double l = hkl.l();
        // Apply transpose Symmetry Operator.
        double hs = rot[0][0] * h + rot[1][0] * k + rot[2][0] * l;
        double ks = rot[0][1] * h + rot[1][1] * k + rot[2][1] * l;
        double ls = rot[0][2] * h + rot[1][2] * k + rot[2][2] * l;
        // Convert back to HKL
        mate.h((int) rint(hs));
        mate.k((int) rint(ks));
        mate.l((int) rint(ls));
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
            double xi = xc * A00 + yc * A10 + zc * A20;
            double yi = xc * A01 + yc * A11 + zc * A21;
            double zi = xc * A02 + yc * A12 + zc * A22;
            // Apply Symmetry Operator.
            double fx = rot00 * xi + rot01 * yi + rot02 * zi;
            double fy = rot10 * xi + rot11 * yi + rot12 * zi;
            double fz = rot20 * xi + rot21 * yi + rot22 * zi;
            // Convert back to Cartesian coordinates.
            mateX[i] = fx * Ai00 + fy * Ai10 + fz * Ai20;
            mateY[i] = fx * Ai01 + fy * Ai11 + fz * Ai21;
            mateZ[i] = fx * Ai02 + fy * Ai12 + fz * Ai22;
        }
    }

    public void toFractionalCoordinates(int n, double x[], double y[],
            double z[], double xf[], double yf[], double zf[]) {
        for (int i = 0; i < n; i++) {
            double xc = x[i];
            double yc = y[i];
            double zc = z[i];
            xf[i] = xc * A00 + yc * A10 + zc * A20;
            yf[i] = xc * A01 + yc * A11 + zc * A21;
            zf[i] = xc * A02 + yc * A12 + zc * A22;
        }
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
            frac[iX] = xc * A00 + yc * A10 + zc * A20;
            frac[iY] = xc * A01 + yc * A11 + zc * A21;
            frac[iZ] = xc * A02 + yc * A12 + zc * A22;
        }
    }

    public void toPrimaryCell(double in[], double out[]) {
        toFractionalCoordinates(in, out);
        out[0] = mod(out[0], 1.0) - 0.5;
        out[1] = mod(out[1], 1.0) - 0.5;
        out[2] = mod(out[2], 1.0) - 0.5;
        toCartesianCoordinates(out, in);
        out[0] = in[0];
        out[1] = in[1];
        out[2] = in[2];
    }

    public void toFractionalCoordinates(double x[], double xf[]) {
        double xc = x[0];
        double yc = x[1];
        double zc = x[2];
        xf[0] = xc * A00 + yc * A10 + zc * A20;
        xf[1] = xc * A01 + yc * A11 + zc * A21;
        xf[2] = xc * A02 + yc * A12 + zc * A22;
    }

    public void toCartesianCoordinates(int n, double xf[], double yf[],
            double zf[], double x[], double y[], double z[]) {
        for (int i = 0; i < n; i++) {
            double xi = xf[i];
            double yi = yf[i];
            double zi = zf[i];
            x[i] = xi * Ai00 + yi * Ai10 + zi * Ai20;
            y[i] = xi * Ai01 + yi * Ai11 + zi * Ai21;
            z[i] = xi * Ai02 + yi * Ai12 + zi * Ai22;
        }
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
            cart[iX] = xf * Ai00 + yf * Ai10 + zf * Ai20;
            cart[iY] = xf * Ai01 + yf * Ai11 + zf * Ai21;
            cart[iZ] = xf * Ai02 + yf * Ai12 + zf * Ai22;
        }
    }

    public void toCartesianCoordinates(double xf[], double x[]) {
        double fx = xf[0];
        double fy = xf[1];
        double fz = xf[2];
        x[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
        x[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
        x[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
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

    public static double sym_phase_shift(double hkl[], SymOp symOp) {
        double trans[] = symOp.tr;
        // Apply translation
        return -2.0 * PI
                * (hkl[0] * trans[0] + hkl[1] * trans[1] + hkl[2] * trans[2]);
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

    public static int mod(int a, int b) {
        int res = a % b;
        if (res < 0) {
            res += b;
        }
        return res;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n Unit Cell\n");
        sb.append(String.format(" A-axis:                               %8.3f\n", a));
        sb.append(String.format(" B-axis:                               %8.3f\n", b));
        sb.append(String.format(" C-axis:                               %8.3f\n", c));
        sb.append(String.format(" Alpha:                                %8.3f\n", alpha));
        sb.append(String.format(" Beta:                                 %8.3f\n", beta));
        sb.append(String.format(" Gamma:                                %8.3f\n", gamma));
        sb.append(String.format("\n Space group\n"));
        sb.append(String.format(" Number:                                    %3d\n", spaceGroup.number));
        sb.append(String.format(" Symbol:                               %8s\n", spaceGroup.shortName));
        sb.append(String.format(" Number of Symmetry Operators:              %3d", spaceGroup.getNumberOfSymOps()));
        return sb.toString();
    }
    private static final int XX = 0;
    private static final int YY = 1;
    private static final int ZZ = 2;
}
