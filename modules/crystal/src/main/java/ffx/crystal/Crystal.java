/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.crystal;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.rint;
import static org.apache.commons.math3.util.FastMath.signum;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;
import static org.apache.commons.math3.util.FastMath.toRadians;

import ffx.utilities.HashCodeUtil;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.mat3Mat3;
import static ffx.numerics.VectorMath.mat3SymVec6;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.transpose3;

/**
 * The Crystal class encapsulates the lattice parameters and space group that
 * describe the geometry and symmetry of a crystal. Methods are available to
 * apply the minimum image convention and space group symmetry operators.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 *
 * @see ReplicatesCrystal
 */
public class Crystal {

    private static final Logger logger = Logger.getLogger(Crystal.class.getName());
    /**
     * Length of the cell edge in the direction of the <b>a</b> basis vector.
     */
    public double a;
    /**
     * Length of the cell edge in the direction of the <b>b</b> basis vector.
     */
    public double b;
    /**
     * Length of the cell edge in the direction of the <b>c</b> basis vector.
     */
    public double c;
    /**
     * Length of the reciprocal basis vector <b>a*</b>.
     */
    public double aStar;
    /**
     * Length of the reciprocal basis vector <b>b*</b>.
     */
    public double bStar;
    /**
     * Length of the reciprocal basis vector <b>c*</b>.
     */
    public double cStar;
    /**
     * The interaxial lattice angle between <b>b</b> and <b>c</b>.
     */
    public double alpha;
    /**
     * The interaxial lattice angle between <b>a</b> and <b>c</b>.
     */
    public double beta;
    /**
     * The interaxial lattice angle between <b>a</b> and <b>b</b>.
     */
    public double gamma;
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
    private List<SymOp> symOpsCartesian;
    /**
     * The crystal unit cell volume.
     */
    public double volume;

    /**
     * Change in the volume with respect to a.
     */
    public double dVdA;
    /**
     * Change in the volume with respect to b.
     */
    public double dVdB;
    /**
     * Change in the volume with respect to c.
     */
    public double dVdC;
    /**
     * Change in the volume with respect to alpha (in Radians). This is set to
     * zero if alpha is fixed.
     */
    public double dVdAlpha;
    /**
     * Change in the volume with respect to beta (in Radians). This is set to
     * zero if beta is fixed.
     */
    public double dVdBeta;
    /**
     * Change in the volume with respect to gamma (in Radians). This is set to
     * zero if gamma is fixed.
     */
    public double dVdGamma;

    /**
     * Matrix to convert from fractional to Cartesian coordinates.
     */
    public final double Ai[][] = new double[3][3];
    /**
     * Entries in the Ai array.
     */
    public double Ai00, Ai01, Ai02, Ai10, Ai11, Ai12, Ai20, Ai21, Ai22;
    /**
     * Matrix to convert from Cartesian to fractional coordinates.
     */
    public double A[][];
    /**
     * Entries in the A array.
     */
    public double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    /**
     * The direct space metric matrix.
     */
    public final double G[][] = new double[3][3];
    /**
     * The reciprocal space metric matrix.
     */
    public double Gstar[][];
    /**
     * Interfacial radius in the direction of the A-axis.
     */
    public double interfacialRadiusA;
    /**
     * Interfacial radius in the direction of the B-axis.
     */
    public double interfacialRadiusB;
    /**
     * Interfacial radius in the direction of the C-axis.
     */
    public double interfacialRadiusC;


    private boolean aperiodic;
    public int scale_flag;
    public int scale_b[] = new int[6];
    public int scale_n;
    /**
     * An atom and one of its symmetry copies within the specialPositionCutoff
     * should be flagged to be at a special position.
     */
    public double specialPositionCutoff = 0.3;
    public double specialPositionCutoff2 = specialPositionCutoff * specialPositionCutoff;
    /**
     * Avogadro's number.
     */
    private static final double AVOGADRO = 6.02214129e23;

    /**
     * The Crystal class encapsulates the lattice parameters and space group.
     * Methods are available to apply the minimum image convention and to apply
     * space group operators.
     *
     * @param a The a-axis length.
     * @param b The b-axis length.
     * @param c The c-axis length.
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
        spaceGroup = SpaceGroup.spaceGroupFactory(sg);
        crystalSystem = spaceGroup.crystalSystem;

        if (!SpaceGroup.checkRestrictions(crystalSystem, a, b, c, alpha, beta, gamma)) {
            String message = " The lattice parameters do not satisfy the " + crystalSystem
                    + " crystal system restrictions:\n" + toString();
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

        updateCrystal();
    }

    /**
     * Update all Crystal variables that are a function of unit cell parameters.
     */
    private void updateCrystal() {

        double cos_alpha = 0.0;
        double sin_beta = 0.0;
        double cos_beta = 0.0;
        double sin_gamma = 0.0;
        double cos_gamma = 0.0;
        double beta_term = 0.0;
        double gamma_term = 0.0;

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
                dVdA = b * c;
                dVdB = a * c;
                dVdC = a * b;
                dVdAlpha = 0.0;
                dVdBeta = 0.0;
                dVdGamma = 0.0;
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
                dVdA = sin_beta * b * c;
                dVdB = sin_beta * a * c;
                dVdC = sin_beta * a * b;
                dVdAlpha = 0.0;
                dVdBeta = cos_beta * a * b * c;
                dVdGamma = 0.0;
                break;
            case HEXAGONAL:
            case TRICLINIC:
            case TRIGONAL:
            default:
                double sin_alpha = sin(toRadians(alpha));
                cos_alpha = cos(toRadians(alpha));
                sin_beta = sin(toRadians(beta));
                cos_beta = cos(toRadians(beta));
                sin_gamma = sin(toRadians(gamma));
                cos_gamma = cos(toRadians(gamma));
                beta_term = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
                gamma_term = sqrt(sin_beta * sin_beta - beta_term * beta_term);
                volume = sin_gamma * gamma_term * a * b * c;

                dVdA = sin_gamma * gamma_term * b * c;
                dVdB = sin_gamma * gamma_term * a * c;
                dVdC = sin_gamma * gamma_term * a * b;

                double dbeta = 2.0 * sin_beta * cos_beta - (2.0 * (cos_alpha - cos_beta * cos_gamma) * sin_beta * cos_gamma) / (sin_gamma * sin_gamma);
                double dgamma1 = -2.0 * (cos_alpha - cos_beta * cos_gamma) * cos_beta / sin_gamma;
                double dgamma2 = cos_alpha - cos_beta * cos_gamma;
                dgamma2 *= dgamma2 * 2.0 * cos_gamma / (sin_gamma * sin_gamma * sin_gamma);

                dVdAlpha = (cos_alpha - cos_beta * cos_gamma) * sin_alpha / (sin_gamma * gamma_term) * a * b * c;
                dVdBeta = 0.5 * sin_gamma * dbeta / gamma_term * a * b * c;
                dVdGamma = (cos_gamma * gamma_term + 0.5 * sin_gamma * (dgamma1 + dgamma2) / gamma_term) * a * b * c;

                break;
        }

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
        RealMatrix m = new Array2DRowRealMatrix(G, true);
        m = new LUDecomposition(m).getSolver().getInverse();
        Gstar = m.getData();

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
        Ai01 = Ai[0][1];
        Ai02 = Ai[0][2];
        Ai10 = Ai[1][0];
        Ai11 = Ai[1][1];
        Ai12 = Ai[1][2];
        Ai20 = Ai[2][0];
        Ai21 = Ai[2][1];
        Ai22 = Ai[2][2];

        // Invert A^-1 to get A
        m = new Array2DRowRealMatrix(Ai, true);
        m = new LUDecomposition(m).getSolver().getInverse();
        A = m.getData();

        // The columns of A are the reciprocal basis vectors
        A00 = A[0][0];
        A10 = A[1][0];
        A20 = A[2][0];
        A01 = A[0][1];
        A11 = A[1][1];
        A21 = A[2][1];
        A02 = A[0][2];
        A12 = A[1][2];
        A22 = A[2][2];

        // Reciprocal basis vector lengths
        aStar = 1.0 / sqrt(A00 * A00 + A10 * A10 + A20 * A20);
        bStar = 1.0 / sqrt(A01 * A01 + A11 * A11 + A21 * A21);
        cStar = 1.0 / sqrt(A02 * A02 + A12 * A12 + A22 * A22);
        if (logger.isLoggable(Level.FINEST)) {
            logger.finest(format(" Reciprocal Lattice Lengths: (%8.3f, %8.3f, %8.3f)",
                    aStar, bStar, cStar));
        }

        // Interfacial diameters from the dot product of the real and reciprocal vectors
        interfacialRadiusA = (Ai00 * A00 + Ai01 * A10 + Ai02 * A20) * aStar;
        interfacialRadiusB = (Ai10 * A01 + Ai11 * A11 + Ai12 * A21) * bStar;
        interfacialRadiusC = (Ai20 * A02 + Ai21 * A12 + Ai22 * A22) * cStar;

        // Divide by 2 to get radii.
        interfacialRadiusA /= 2.0;
        interfacialRadiusB /= 2.0;
        interfacialRadiusC /= 2.0;

        if (logger.isLoggable(Level.FINEST)) {
            logger.finest(format(" Interfacial radii: (%8.3f, %8.3f, %8.3f)",
                    interfacialRadiusA, interfacialRadiusB, interfacialRadiusC));
        }

        List<SymOp> symOps = spaceGroup.symOps;
        int nSymm = symOps.size();
        if (symOpsCartesian == null) {
            symOpsCartesian = new ArrayList<>(nSymm);
        } else {
            symOpsCartesian.clear();
        }

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
    }

    /**
     * This method should be called to update the unit cell parameters of a
     * crystal. The proposed parameters will only be accepted if symmetry
     * restrictions are satisfied. If so, all Crystal variables that depend on
     * the unit cell parameters will be updated.
     *
     * @param a length of the a-axis.
     * @param b length of the b-axis.
     * @param c length of the c-axis.
     * @param alpha Angle between b-axis and c-axis.
     * @param beta Angle between a-axis and c-axis.
     * @param gamma Angle between a-axis and b-axis.
     * @return The method return true if the parameters are accepted, false
     * otherwise.
     */
    public boolean changeUnitCellParameters(double a, double b, double c, double alpha, double beta,
                                            double gamma) {

        if (!SpaceGroup.checkRestrictions(crystalSystem, a, b, c, alpha, beta, gamma)) {
            if (logger.isLoggable(Level.FINE)) {
                String message = " The proposed lattice parameters do not satisfy the " + crystalSystem
                        + " crystal system restrictions and were ignored.";
                logger.fine(message);
            }
            return false;
        }

        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;

        updateCrystal();

        return true;
    }


    /**
     * Strain the unit cell vectors.
     *
     * @param dStrain
     */
    public boolean perturbCellVectors(double dStrain[][]) {

        double AA[][] = new double[3][3];
        double newAi[][] = new double[3][3];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                AA[i][j] = getUnitCell().Ai[i][j];
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    double Kronecker = 0.0;
                    if (j == k) {
                        Kronecker = 1.0;
                    }
                    // Aij' = Sum over K ( Kron_jk + dStrain_jk ) * Aij
                    newAi[i][j] += (Kronecker + dStrain[j][k]) * AA[i][k];
                }
            }
        }

        // Update a-, b-, and c-axis lengths.
        double aa = r(newAi[0]);
        double bb = r(newAi[1]);
        double cc = r(newAi[2]);

        // Update alpha, beta and gamma angles.
        double aalpha = toDegrees(acos(dot(newAi[1], newAi[2]) / (bb * cc)));
        double bbeta = toDegrees(acos(dot(newAi[0], newAi[2]) / (aa * cc)));
        double ggamma = toDegrees(acos(dot(newAi[0], newAi[1]) / (aa * bb)));

        return changeUnitCellParameters(aa, bb, cc, aalpha, bbeta, ggamma);
    }

    public double getDensity(double mass) {
        int nSymm = spaceGroup.symOps.size();
        double dens = (mass * nSymm / AVOGADRO) * (1.0e24 / volume);
        return dens;
    }

    public void setDensity(double dens, double mass) {
        double currentDensity = getDensity(mass);

        double scale = cbrt(currentDensity / dens);
        Crystal uc = getUnitCell();
        changeUnitCellParameters(uc.a * scale, uc.b * scale, uc.c * scale, alpha, beta, gamma);
        currentDensity = getDensity(mass);

        logger.info(format(" Updated density %6.3f (g/cc) with unit cell %s.",
                currentDensity, uc.toShortString()));
    }

    public void randomParameters(double dens, double mass) {
        double params[] = SpaceGroup.resetUnitCellParams(crystalSystem);
        changeUnitCellParameters(params[0], params[1], params[2], params[3], params[4], params[5]);
        setDensity(dens, mass);
    }

    /**
     * <p>
     * Setter for the field <code>specialPositionCutoff</code>.</p>
     *
     * @param cutoff a double.
     */
    public void setSpecialPositionCutoff(double cutoff) {
        specialPositionCutoff = cutoff;
        specialPositionCutoff2 = cutoff * cutoff;
    }

    /**
     * <p>
     * Getter for the field <code>specialPositionCutoff</code>.</p>
     *
     * @return a double.
     */
    public double getSpecialPositionCutoff() {
        return specialPositionCutoff;
    }

    /**
     * <p>
     * checkProperties</p>
     *
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @return a {@link ffx.crystal.Crystal} object.
     */
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

    /**
     * Two crystals are equal if all unit cell parameters are within 0.01.
     *
     * {@inheritDoc}
     */
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

        Crystal other = (Crystal) obj;

        return (abs(a - other.a) < 0.01
                && abs(b - other.b) < 0.01
                && abs(c - other.c) < 0.01
                && abs(alpha - other.alpha) < 0.01
                && abs(beta - other.beta) < 0.01
                && abs(gamma - other.gamma) < 0.01
                && spaceGroup.number == other.spaceGroup.number);
    }

    /**
     * Two crystals are equal only if all unit cell parameters are exactly the
     * same.
     *
     * @param obj the Crystal to compare to.
     * @return true if all unit cell parameters are exactly the same.
     */
    public boolean strictEquals(Object obj) {
        if (obj == null) {
            return false;
        }

        if (!(obj instanceof Crystal)) {
            return false;
        }

        if (this == obj) {
            return true;
        }

        Crystal other = (Crystal) obj;

        return (a == other.a && b == other.b && c == other.c
                && alpha == other.alpha && beta == other.beta && gamma == other.gamma
                && spaceGroup.number == other.spaceGroup.number);
    }

    /**
     * {@inheritDoc}
     */
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
     *
     * @param aperiodic a boolean.
     */
    public void setAperiodic(boolean aperiodic) {
        this.aperiodic = aperiodic;
    }

    /**
     * <p>
     * aperiodic</p>
     *
     * @return a boolean.
     */
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
        double xf = dx * A00 + dy * A10 + dz * A20;
        double yf = dx * A01 + dy * A11 + dz * A21;
        double zf = dx * A02 + dy * A12 + dz * A22;
        xf = floor(abs(xf) + 0.5) * signum(-xf) + xf;
        yf = floor(abs(yf) + 0.5) * signum(-yf) + yf;
        zf = floor(abs(zf) + 0.5) * signum(-zf) + zf;
        dx = xf * Ai00 + yf * Ai10 + zf * Ai20;
        dy = xf * Ai01 + yf * Ai11 + zf * Ai21;
        dz = xf * Ai02 + yf * Ai12 + zf * Ai22;
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * <p>
     * averageTensor</p>
     *
     * @param m an array of double.
     * @param r an array of double.
     */
    public void averageTensor(double m[][], double r[][]) {
        int n = spaceGroup.symOps.size();
        for (int i = 0; i < n; i++) {
            SymOp symop = spaceGroup.symOps.get(i);
            double rot[][] = symop.rot;
            double rt[][] = transpose3(rot);
            double rmrt[][] = mat3Mat3(mat3Mat3(rot, m), rt);
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
     * <p>
     * averageTensor</p>
     *
     * @param v an array of double.
     * @param r an array of double.
     */
    public void averageTensor(double v[], double r[][]) {
        int n = spaceGroup.symOps.size();
        for (int i = 0; i < n; i++) {
            SymOp symop = spaceGroup.symOps.get(i);
            double rot[][] = symop.rot;
            double rt[][] = transpose3(rot);
            double rmrt[][] = mat3Mat3(mat3SymVec6(rot, v), rt);
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
     * Apply a Cartesian symmetry operator to an array of Cartesian coordinates.
     * If the arrays x, y or z are null or not of length n, the method returns
     * immediately. If mateX, mateY or mateZ are null or not of length n, new
     * arrays are allocated.
     *
     * @param n Number of atoms.
     * @param x Input x coordinates.
     * @param y Input y coordinates.
     * @param z Input z coordinates.
     * @param mateX Output x coordinates.
     * @param mateY Output y coordinates.
     * @param mateZ Output z coordinates.
     * @param symOp The symmetry operator.
     */
    public void applyCartSymOp(int n, double x[], double y[], double z[],
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
            // Apply Symmetry Operator.
            mateX[i] = rot00 * xc + rot01 * yc + rot02 * zc + t0;
            mateY[i] = rot10 * xc + rot11 * yc + rot12 * zc + t1;
            mateZ[i] = rot20 * xc + rot21 * yc + rot22 * zc + t2;
        }
    }

    /**
     * Apply a fractional symmetry operator to an array of Cartesian
     * coordinates. If the arrays x, y or z are null or not of length n, the
     * method returns immediately. If mateX, mateY or mateZ are null or not of
     * length n, new arrays are allocated.
     *
     * @param n Number of atoms.
     * @param x Input x coordinates.
     * @param y Input y coordinates.
     * @param z Input z coordinates.
     * @param mateX Output x coordinates.
     * @param mateY Output y coordinates.
     * @param mateZ Output z coordinates.
     * @param symOp The symmetry operator.
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
     * @param h Input coordinates.
     * @param k Input coordinates.
     * @param l Input coordinates.
     * @param mate Symmetry mate coordinates.
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
     * @param xyz Input coordinates.
     * @param mate Symmetry mate coordinates.
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
     * @param xyz Input coordinates.
     * @param mate Symmetry mate coordinates.
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
     * @param xyz Input coordinates.
     * @param mate Symmetry mate coordinates.
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
     * @param xyz Input coordinates.
     * @param mate Symmetry mate coordinates.
     * @param symOp The symmetry operator.
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
     * Multiply coordinates by the transpose of a matrix.
     *
     * @param in input coordinates.
     * @param out output coordinates.
     * @param matrix multiply by the transpose of this matrix.
     */
    public static void applyMatrixTranspose(double in[], double out[], double matrix[][]) {
        double xc = in[0];
        double yc = in[1];
        double zc = in[2];
        out[0] = xc * matrix[0][0] + yc * matrix[1][0] + zc * matrix[2][0];
        out[1] = xc * matrix[0][1] + yc * matrix[1][1] + zc * matrix[2][1];
        out[2] = xc * matrix[0][2] + yc * matrix[1][2] + zc * matrix[2][2];
    }

    /**
     * Compute the total transformation operator R = ToCart * Rot * ToFrac.
     *
     * @param symOp Symmetry operator to apply.
     * @param rotmat Resulting transformation operator R.
     */
    public void getTransformationOperator(SymOp symOp, double rotmat[][]) {
        double rot[][] = symOp.rot;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotmat[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        rotmat[i][j] += Ai[k][i] * rot[k][l] * A[j][l];
                    }
                }
            }
        }
    }

    /**
     * Apply a symmetry operator to one set of coordinates.
     *
     * @param xyz Input coordinates.
     * @param mate Symmetry mate coordinates.
     * @param symOp The symmetry operator.
     * @param rotmat an array of double.
     */
    public void applyTransSymRot(double xyz[], double mate[], SymOp symOp, double rotmat[][]) {

        /**
         * The transformation operator R = ToCart * Rot * ToFrac
         */
        getTransformationOperator(symOp, rotmat);

        /**
         * Apply R^T (its transpose).
         */
        applyMatrixTranspose(xyz, mate, rotmat);
    }

    /**
     * Apply a symmetry operator to one HKL.
     *
     * @param hkl Input HKL.
     * @param mate Symmetry mate HKL.
     * @param symOp The symmetry operator.
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
     * @param hkl Input HKL.
     * @param mate Symmetry mate HKL.
     * @param symOp The symmetry operator.
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
     * Apply a symmetry rotation to an array of Cartesian coordinates. If the
     * arrays x, y or z are null or not of length n, the method returns
     * immediately. If mateX, mateY or mateZ are null or not of length n, new
     * arrays are allocated.
     *
     * @param n Number of atoms.
     * @param x Input x coordinates.
     * @param y Input y coordinates.
     * @param z Input z coordinates.
     * @param mateX Output x coordinates.
     * @param mateY Output y coordinates.
     * @param mateZ Output z coordinates.
     * @param symOp The symmetry operator.
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

    /**
     * Apply the transpose of a symmetry rotation to an array of Cartesian
     * coordinates. If the arrays x, y or z are null or not of length n, the
     * method returns immediately. If mateX, mateY or mateZ are null or not of
     * length n, new arrays are allocated.
     *
     * @param n Number of atoms.
     * @param x Input x coordinates.
     * @param y Input y coordinates.
     * @param z Input z coordinates.
     * @param mateX Output x coordinates.
     * @param mateY Output y coordinates.
     * @param mateZ Output z coordinates.
     * @param symOp The symmetry operator.
     * @param rotmat an array of double.
     */
    public void applyTransSymRot(int n, double x[], double y[], double z[],
                                 double mateX[], double mateY[], double mateZ[],
                                 SymOp symOp, double rotmat[][]) {

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

        /**
         * The transformation operator R = ToCart * Rot * ToFrac
         */
        getTransformationOperator(symOp, rotmat);

        for (int i = 0; i < n; i++) {
            // Apply R^T (its transpose).
            double xc = x[i];
            double yc = y[i];
            double zc = z[i];
            mateX[i] = xc * rotmat[0][0] + yc * rotmat[1][0] + zc * rotmat[2][0];
            mateY[i] = xc * rotmat[0][1] + yc * rotmat[1][1] + zc * rotmat[2][1];
            mateZ[i] = xc * rotmat[0][2] + yc * rotmat[1][2] + zc * rotmat[2][2];
        }
    }

    /**
     * <p>
     * toFractionalCoordinates</p>
     *
     * @param n a int.
     * @param x an array of double.
     * @param y an array of double.
     * @param z an array of double.
     * @param xf an array of double.
     * @param yf an array of double.
     * @param zf an array of double.
     */
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

    /**
     * <p>
     * toFractionalCoordinates</p>
     *
     * @param n a int.
     * @param cart an array of double.
     * @param frac an array of double.
     */
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

    /**
     * <p>
     * toPrimaryCell</p>
     *
     * @param in an array of double.
     * @param out an array of double.
     */
    public void toPrimaryCell(double in[], double out[]) {
        toFractionalCoordinates(in, out);
        out[0] = mod(out[0], 1.0) - 0.5;
        out[1] = mod(out[1], 1.0) - 0.5;
        out[2] = mod(out[2], 1.0) - 0.5;
        toCartesianCoordinates(out, out);
    }

    /**
     * <p>
     * toFractionalCoordinates</p>
     *
     * @param x an array of double.
     * @param xf an array of double.
     */
    public void toFractionalCoordinates(double x[], double xf[]) {
        double xc = x[0];
        double yc = x[1];
        double zc = x[2];
        xf[0] = xc * A00 + yc * A10 + zc * A20;
        xf[1] = xc * A01 + yc * A11 + zc * A21;
        xf[2] = xc * A02 + yc * A12 + zc * A22;
    }

    /**
     * <p>
     * toCartesianCoordinates</p>
     *
     * @param n a int.
     * @param xf an array of double.
     * @param yf an array of double.
     * @param zf an array of double.
     * @param x an array of double.
     * @param y an array of double.
     * @param z an array of double.
     */
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

    /**
     * <p>
     * toCartesianCoordinates</p>
     *
     * @param n a int.
     * @param frac an array of double.
     * @param cart an array of double.
     */
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

    /**
     * <p>
     * toCartesianCoordinates</p>
     *
     * @param xf an array of double.
     * @param x an array of double.
     */
    public void toCartesianCoordinates(double xf[], double x[]) {
        double fx = xf[0];
        double fy = xf[1];
        double fz = xf[2];
        x[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
        x[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
        x[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
    }

    /**
     * Minimum distance between two coordinates over all symmetry operators.
     *
     * @param xyzA Coordinate A
     * @param xyzB Coordinate B
     * @return Minimum distance in crystal
     */
    public double minDistOverSymOps(double[] xyzA, double[] xyzB) {
        double dist = 0;
        for (int i = 0; i < 3; i++) {
            double dx = xyzA[i] - xyzB[i];
            dist += (dx * dx);
        }
        double[] symB = new double[3];
        for (SymOp symOp : spaceGroup.symOps) {
            applySymOp(xyzB, symB, symOp);
            for (int i = 0; i < 3; i++) {
                symB[i] -= xyzA[i];
            }
            double d = image(symB);
            dist = d < dist ? d : dist;
        }
        return FastMath.sqrt(dist);
    }

    /**
     * <p>
     * quad_form</p>
     *
     * @param v an array of double.
     * @param mat an array of double.
     * @return a double.
     */
    public static double quad_form(double v[], double mat[][]) {
        return (v[0] * (v[0] * mat[0][0] + 2 * (v[1] * mat[0][1] + v[2] * mat[0][2]))
                + v[1] * (v[1] * mat[1][1] + 2 * (v[2] * mat[1][2]))
                + v[2] * v[2] * mat[2][2]);
    }

    /**
     * <p>
     * quad_form</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @param mat an array of double.
     * @return a double.
     */
    public static double quad_form(HKL hkl, double mat[][]) {
        return (hkl.h() * (hkl.h() * mat[0][0] + 2 * (hkl.k() * mat[0][1] + hkl.l() * mat[0][2]))
                + hkl.k() * (hkl.k() * mat[1][1] + 2 * (hkl.l() * mat[1][2]))
                + hkl.l() * hkl.l() * mat[2][2]);
    }

    /**
     * <p>
     * invressq</p>
     *
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a double.
     */
    public static double invressq(Crystal crystal, HKL hkl) {
        return quad_form(hkl, crystal.Gstar);
    }

    /**
     * <p>
     * res</p>
     *
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a double.
     */
    public static double res(Crystal crystal, HKL hkl) {
        return 1.0 / sqrt(quad_form(hkl, crystal.Gstar));
    }

    /**
     * <p>
     * sym_phase_shift</p>
     *
     * @param hkl an array of double.
     * @param symOp a {@link ffx.crystal.SymOp} object.
     * @return a double.
     */
    public static double sym_phase_shift(double hkl[], SymOp symOp) {
        double trans[] = symOp.tr;
        // Apply translation
        return -2.0 * PI
                * (hkl[0] * trans[0] + hkl[1] * trans[1] + hkl[2] * trans[2]);
    }

    /**
     * <p>
     * sym_phase_shift</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @param symOp a {@link ffx.crystal.SymOp} object.
     * @return a double.
     */
    public static double sym_phase_shift(HKL hkl, SymOp symOp) {
        double trans[] = symOp.tr;
        // Apply translation
        return -2.0 * PI
                * (hkl.h() * trans[0] + hkl.k() * trans[1] + hkl.l() * trans[2]);
    }

    /**
     * This is an atypical mod function used by crystallography methods.
     *
     * <p>
     * mod</p>
     *
     * @param a a double.
     * @param b a double.
     * @return a double.
     */
    public static double mod(double a, double b) {
        double res = a % b;
        if (res < 0.0) {
            res += b;
        }
        return res;
    }

    /**
     * This is an atypical mod function used by crystallography methods.
     *
     * <p>
     * mod</p>
     *
     * @param a a int.
     * @param b a int.
     * @return a int.
     */
    public static int mod(int a, int b) {
        int res = a % b;
        if (res < 0) {
            res += b;
        }
        return res;
    }

    /**
     * A String containing the unit cell parameters.
     *
     * @return A string with the unit cell parameters.
     */
    public String toShortString() {
        return String.format("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f", a, b, c, alpha, beta, gamma);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n Unit Cell\n");
        sb.append(String.format("  A-axis:                              %8.3f\n", a));
        sb.append(String.format("  B-axis:                              %8.3f\n", b));
        sb.append(String.format("  C-axis:                              %8.3f\n", c));
        sb.append(String.format("  Alpha:                               %8.3f\n", alpha));
        sb.append(String.format("  Beta:                                %8.3f\n", beta));
        sb.append(String.format("  Gamma:                               %8.3f\n", gamma));
        sb.append(String.format("  Space group\n"));
        sb.append(String.format("   Number:                                  %3d\n", spaceGroup.number));
        sb.append(String.format("   Symbol:                             %8s\n", spaceGroup.shortName));
        sb.append(String.format("   Number of Symmetry Operators:            %3d", spaceGroup.getNumberOfSymOps()));
        return sb.toString();
    }

    /**
     * A mask equal to 0 for X-coordinates.
     */
    private static final int XX = 0;
    /**
     * A mask equal to 1 for Y-coordinates.
     */
    private static final int YY = 1;
    /**
     * A mask equal to 2 for Z-coordinates.
     */
    private static final int ZZ = 2;
}
