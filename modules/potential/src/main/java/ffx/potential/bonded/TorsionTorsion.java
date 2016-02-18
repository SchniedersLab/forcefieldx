/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionTorsionType;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.TorsionTorsionType.units;

/**
 * The TorsionTorsion class represents two adjacent torsional angles formed by
 * five bonded atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class TorsionTorsion extends BondedTerm implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(TorsionTorsion.class.getName());
    private static final long serialVersionUID = 1L;
    public TorsionTorsionType torsionTorsionType = null;
    public final Torsion torsions[] = new Torsion[2];
    protected final boolean reversed;

    private boolean lambdaTerm = false;
    private double lambda = 1.0;
    private double dEdL = 0.0;

    /**
     * Torsion-Torsion constructor.
     *
     * @param firstBond a {@link ffx.potential.bonded.Bond} object.
     * @param angle a {@link ffx.potential.bonded.Angle} object.
     * @param lastBond a {@link ffx.potential.bonded.Bond} object.
     * @param reversed a boolean.
     */
    public TorsionTorsion(Bond firstBond, Angle angle, Bond lastBond,
            boolean reversed) {
        super();
        atoms = new Atom[5];
        bonds = new Bond[4];
        if (!reversed) {
            atoms[1] = angle.atoms[0];
            atoms[2] = angle.atoms[1];
            atoms[3] = angle.atoms[2];
            atoms[0] = firstBond.get1_2(atoms[1]);
            atoms[4] = lastBond.get1_2(atoms[3]);

            bonds[0] = firstBond;
            bonds[1] = angle.bonds[0];
            bonds[2] = angle.bonds[1];
            bonds[3] = lastBond;
        } else {
            atoms[1] = angle.atoms[2];
            atoms[2] = angle.atoms[1];
            atoms[3] = angle.atoms[0];
            atoms[0] = lastBond.get1_2(atoms[1]);
            atoms[4] = firstBond.get1_2(atoms[3]);
            bonds[0] = lastBond;
            bonds[1] = angle.bonds[1];
            bonds[2] = angle.bonds[0];
            bonds[3] = firstBond;
        }
        torsions[0] = atoms[0].getTorsion(atoms[1], atoms[2], atoms[3]);
        torsions[1] = atoms[4].getTorsion(atoms[3], atoms[2], atoms[1]);
        this.reversed = reversed;
        setID_Key(false);
    }

    /**
     * @param firstBond the first Bond.
     * @param angle the Angle.
     * @param lastBond the last Bond.
     * @param forceField the ForceField parameters to apply.
     * @return the new TorsionTorsion, or null.
     */
    public static TorsionTorsion torsionTorsionFactory(Bond firstBond,
            Angle angle, Bond lastBond, ForceField forceField) {
        int c5[] = new int[5];
        Atom atom1 = angle.atoms[0];
        Atom atom3 = angle.atoms[2];
        c5[0] = firstBond.get1_2(atom1).getAtomType().atomClass;
        c5[1] = atom1.getAtomType().atomClass;
        c5[2] = angle.atoms[1].getAtomType().atomClass;
        c5[3] = atom3.getAtomType().atomClass;
        c5[4] = lastBond.get1_2(atom3).getAtomType().atomClass;
        String key = TorsionTorsionType.sortKey(c5);
        boolean reversed = false;
        TorsionTorsionType torsionTorsionType = forceField.getTorsionTorsionType(key);
        if (torsionTorsionType == null) {
            key = TorsionTorsionType.reverseKey(c5);
            torsionTorsionType = forceField.getTorsionTorsionType(key);
            reversed = true;
        }
        if (torsionTorsionType == null) {
            return null;
        }
        TorsionTorsion torsionTorsion = new TorsionTorsion(
                firstBond, angle, lastBond, reversed);
        torsionTorsion.torsionTorsionType = torsionTorsionType;
        return torsionTorsion;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update() {
        energy(false);
    }

    protected static final double ad[] = new double[3];
    protected static final double a0[] = new double[3];
    protected static final double a1[] = new double[3];
    protected static final double a2[] = new double[3];
    protected static final double a3[] = new double[3];
    protected static final double a4[] = new double[3];
    /**
     * Constant <code>v01=new double[3]</code>
     */
    protected static final double v01[] = new double[3];
    /**
     * Constant <code>v12=new double[3]</code>
     */
    protected static final double v12[] = new double[3];
    /**
     * Constant <code>v23=new double[3]</code>
     */
    protected static final double v23[] = new double[3];
    /**
     * Constant <code>v34=new double[3]</code>
     */
    protected static final double v34[] = new double[3];
    /**
     * Constant <code>v02=new double[3]</code>
     */
    protected static final double[] v02 = new double[3];
    /**
     * Constant <code>v13=new double[3]</code>
     */
    protected static final double[] v13 = new double[3];
    /**
     * Constant <code>v24=new double[3]</code>
     */
    protected static final double[] v24 = new double[3];
    /**
     * Constant <code>t=new double[3]</code>
     */
    protected static final double t[] = new double[3];
    /**
     * Constant <code>u=new double[3]</code>
     */
    protected static final double u[] = new double[3];
    /**
     * Constant <code>v=new double[3]</code>
     */
    protected static final double v[] = new double[3];
    /**
     * Constant <code>x1=new double[3]</code>
     */
    protected static final double[] x1 = new double[3];
    /**
     * Constant <code>x2=new double[3]</code>
     */
    protected static final double[] x2 = new double[3];
    /**
     * Constant <code>tu=new double[3]</code>
     */
    protected static final double tu[] = new double[3];
    /**
     * Constant <code>uv=new double[3]</code>
     */
    protected static final double uv[] = new double[3];
    /**
     * Array of 4 spline energies surrounding the actual Torsion-Torsion
     * location.
     */
    protected static final double e[] = new double[4];
    /**
     * Array of 4 spline x-gradients surrounding the actual Torsion-Torsion
     * location.
     */
    protected static final double dx[] = new double[4];
    /**
     * Array of 4 spline y-gradients surrounding the actual Torsion-Torsion
     * location.
     */
    protected static final double dy[] = new double[4];
    /**
     * Array of 4 spline xy-gradients surrounding the actual Torsion-Torsion
     * location.
     */
    protected static final double dxy[] = new double[4];
    /**
     * Gradient on atom 0.
     */
    protected static final double g0[] = new double[3];
    /**
     * Gradient on Atom 1.
     */
    protected static final double g1[] = new double[3];
    /**
     * Gradient on Atom 2.
     */
    protected static final double g2[] = new double[3];
    /**
     * Gradient on Atom 3.
     */
    protected static final double g3[] = new double[3];
    /**
     * Constant <code>g4=new double[3]</code>
     */
    protected static final double g4[] = new double[3];

    /**
     * Evaluate the Torsion-Torsion energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        dEdL = 0.0;
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);
        atoms[4].getXYZ(a4);
        diff(a1, a0, v01);
        diff(a2, a1, v12);
        diff(a3, a2, v23);
        diff(a4, a3, v34);
        cross(v01, v12, t);
        cross(v12, v23, u);
        cross(v23, v34, v);
        cross(t, u, tu);
        cross(u, v, uv);
        double rt2 = dot(t, t);
        double ru2 = dot(u, u);
        double rv2 = dot(v, v);
        double rtru = sqrt(rt2 * ru2);
        double rurv = sqrt(ru2 * rv2);
        if (rtru != 0.0 && rurv != 0.0) {
            double r12 = r(v12);
            double cosine1 = dot(t, u) / rtru;
            cosine1 = min(1.0, max(-1.0, cosine1));
            double angle1 = toDegrees(acos(cosine1));
            double sign = dot(v01, u);
            if (sign < 0.0) {
                angle1 *= -1.0;
            }
            double r23 = r(v23);
            double cosine2 = dot(u, v) / rurv;
            cosine2 = min(1.0, max(-1.0, cosine2));
            double angle2 = toDegrees(acos(cosine2));
            sign = dot(v12, v);
            if (sign < 0.0) {
                angle2 *= -1.0;
            }
            double t1 = angle1;
            double t2 = angle2;
            sign = chktor();
            t1 *= sign;
            t2 *= sign;
            /**
             * Use bicubic interpolation to compute the spline values.
             */
            int nx = torsionTorsionType.nx;
            int ny = torsionTorsionType.ny;
            int nlow = 0;
            int nhigh = nx - 1;
            while (nhigh - nlow > 1) {
                int nt = (nhigh + nlow) / 2;
                if (torsionTorsionType.tx[nt] > t1) {
                    nhigh = nt;
                } else {
                    nlow = nt;
                }
            }
            int xlow = nlow;
            nlow = 0;
            nhigh = ny - 1;
            while (nhigh - nlow > 1) {
                int nt = (nhigh + nlow) / 2;
                if (torsionTorsionType.ty[nt] > t2) {
                    nhigh = nt;
                } else {
                    nlow = nt;
                }
            }
            int ylow = nlow;
            double x1l = torsionTorsionType.tx[xlow];
            double x1u = torsionTorsionType.tx[xlow + 1];
            double y1l = torsionTorsionType.ty[ylow];
            double y1u = torsionTorsionType.ty[ylow + 1];
            int pos2 = (ylow + 1) * nx + xlow;
            int pos1 = pos2 - nx;
            e[0] = torsionTorsionType.energy[pos1];
            e[1] = torsionTorsionType.energy[pos1 + 1];
            e[2] = torsionTorsionType.energy[pos2 + 1];
            e[3] = torsionTorsionType.energy[pos2];
            dx[0] = torsionTorsionType.dx[pos1];
            dx[1] = torsionTorsionType.dx[pos1 + 1];
            dx[2] = torsionTorsionType.dx[pos2 + 1];
            dx[3] = torsionTorsionType.dx[pos2];
            dy[0] = torsionTorsionType.dy[pos1];
            dy[1] = torsionTorsionType.dy[pos1 + 1];
            dy[2] = torsionTorsionType.dy[pos2 + 1];
            dy[3] = torsionTorsionType.dy[pos2];
            dxy[0] = torsionTorsionType.dxy[pos1];
            dxy[1] = torsionTorsionType.dxy[pos1 + 1];
            dxy[2] = torsionTorsionType.dxy[pos2 + 1];
            dxy[3] = torsionTorsionType.dxy[pos2];
            if (!gradient && !lambdaTerm) {
                energy = units * bcuint(x1l, x1u, y1l, y1u, t1, t2);
                dEdL = energy;
                energy = lambda * energy;
            } else {
                double ansy[] = new double[2];
                energy = units * bcuint1(x1l, x1u, y1l, y1u, t1, t2, ansy);
                dEdL = energy;
                energy = lambda * energy;
                double dedang1 = sign * units * toDegrees(ansy[0]);
                double dedang2 = sign * units * toDegrees(ansy[1]);
                /**
                 * Derivative components for the first angle.
                 */
                diff(a2, a0, v02);
                diff(a3, a1, v13);
                cross(t, v12, x1);
                cross(u, v12, x2);
                scalar(x1, dedang1 / (rt2 * r12), x1);
                scalar(x2, -dedang1 / (ru2 * r12), x2);
                cross(x1, v12, g0);
                cross(v02, x1, g1);
                cross(x2, v23, g2);
                sum(g1, g2, g1);
                cross(x1, v01, g2);
                cross(v13, x2, g3);
                sum(g2, g3, g2);
                cross(x2, v12, g3);
                if (lambdaTerm) {
                    atoms[0].addToLambdaXYZGradient(g0[0], g0[1], g0[2]);
                    atoms[1].addToLambdaXYZGradient(g1[0], g1[1], g1[2]);
                    atoms[2].addToLambdaXYZGradient(g2[0], g2[1], g2[2]);
                    atoms[3].addToLambdaXYZGradient(g3[0], g3[1], g3[2]);
                }
                if (gradient) {
                    atoms[0].addToXYZGradient(lambda * g0[0], lambda * g0[1], lambda * g0[2]);
                    atoms[1].addToXYZGradient(lambda * g1[0], lambda * g1[1], lambda * g1[2]);
                    atoms[2].addToXYZGradient(lambda * g2[0], lambda * g2[1], lambda * g2[2]);
                    atoms[3].addToXYZGradient(lambda * g3[0], lambda * g3[1], lambda * g3[2]);
                }
                /**
                 * Derivative components for the 2nd angle.
                 */
                diff(a4, a2, v24);
                cross(u, v23, x1);
                cross(v, v23, x2);
                scalar(x1, dedang2 / (ru2 * r23), x1);
                scalar(x2, -dedang2 / (rv2 * r23), x2);
                cross(x1, v23, g1);
                cross(v13, x1, g2);
                cross(x2, v34, g3);
                sum(g2, g3, g2);
                cross(x1, v12, g3);
                cross(v24, x2, g4);
                sum(g3, g4, g3);
                cross(x2, v23, g4);
                if (lambdaTerm) {
                    atoms[1].addToLambdaXYZGradient(g1[0], g1[1], g1[2]);
                    atoms[2].addToLambdaXYZGradient(g2[0], g2[1], g2[2]);
                    atoms[3].addToLambdaXYZGradient(g3[0], g3[1], g3[2]);
                    atoms[4].addToLambdaXYZGradient(g4[0], g4[1], g4[2]);
                }
                if (gradient) {
                    atoms[1].addToXYZGradient(g1[0] * lambda, g1[1] * lambda, g1[2] * lambda);
                    atoms[2].addToXYZGradient(g2[0] * lambda, g2[1] * lambda, g2[2] * lambda);
                    atoms[3].addToXYZGradient(g3[0] * lambda, g3[1] * lambda, g3[2] * lambda);
                    atoms[4].addToXYZGradient(g4[0] * lambda, g4[1] * lambda, g4[2] * lambda);
                }
            }
        }
        return energy;
    }

    /*
     * Log details for this Torsion-Torsion energy term.
     */
    /**
     * <p>
     * log</p>
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f",
                "Torsional-Torsion", atoms[0].getXYZIndex(), atoms[0].getAtomType().name, atoms[1].getXYZIndex(),
                atoms[1].getAtomType().name, atoms[2].getXYZIndex(), atoms[2].getAtomType().name,
                atoms[3].getXYZIndex(), atoms[3].getAtomType().name, energy));
    }

    /**
     * {@inheritDoc}
     *
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.2f,%7.2f,%7.2f)", id, torsions[0].value,
                torsions[1].value, energy);
    }
    /**
     * Vector from the central atom to site 0.
     */
    protected static double[] vc0 = new double[3];
    /**
     * Vector from the central atom to site 1.
     */
    protected static double[] vc1 = new double[3];
    /**
     * Vector from the central atom to site 2.
     */
    protected static double[] vc2 = new double[3];

    /**
     * Check for inversion of the central atom if it is chiral.
     *
     * @return The sign convention - if negative the torsion angle signs are
     * inverted.
     */
    protected double chktor() {
        ArrayList<Bond> bnds = atoms[2].getBonds();
        /**
         * To be chiral, the central atom must have 4 bonds.
         */
        if (bnds.size() == 4) {
            /**
             * Find the two atoms that are not part of the dihedral.
             */
            Atom atom1 = null;
            Atom atom2 = null;
            for (Bond b : bnds) {
                Atom a = b.get1_2(atoms[2]);
                if (a != atoms[1] && a != atoms[3]) {
                    if (atom1 == null) {
                        atom1 = a;
                    } else {
                        atom2 = a;
                    }
                }
            }
            /**
             * Choose atom1 or atom2 to use for the chiral check, depending on
             * their atom types and atomic number.
             */
            Atom atom = null;
            if (atom1.getType() > atom2.getType()) {
                atom = atom1;
            }
            if (atom2.getType() > atom1.getType()) {
                atom = atom2;
            }
            if (atom1.getAtomicNumber() > atom2.getAtomicNumber()) {
                atom = atom1;
            }
            if (atom2.getAtomicNumber() > atom1.getAtomicNumber()) {
                atom = atom2;
            }
            /**
             * Compute the signed parallelpiped volume at the central site.
             */
            if (atom != null) {
                atom.getXYZ(ad);
                atoms[1].getXYZ(a1);
                atoms[2].getXYZ(a2);
                atoms[3].getXYZ(a3);
                diff(ad, a2, vc0);
                diff(a1, a2, vc1);
                diff(a3, a2, vc2);
                double volume = vc0[0] * (vc1[1] * vc2[2] - vc1[2] * vc2[1]) + vc1[0] * (vc2[1] * vc0[2] - vc2[2] * vc0[1]) + vc2[0] * (vc0[1] * vc1[2] - vc0[2] * vc1[1]);
                if (volume < 0.0) {
                    return -1.0;
                }
            }
        }
        return 1.0;
    }
    private static double c[][] = new double[4][4];

    /**
     * <p>
     * bcuint</p>
     *
     * @param x1l a double.
     * @param x1u a double.
     * @param y1l a double.
     * @param y1u a double.
     * @param t1 a double.
     * @param t2 a double.
     * @return a double.
     */
    protected double bcuint(double x1l, double x1u, double y1l, double y1u,
            double t1, double t2) {
        double deltax = x1u - x1l;
        double deltay = y1u - y1l;
        bcucof(deltax, deltay, c);
        double tx = (t1 - x1l) / deltax;
        double ux = (t2 - y1l) / deltay;
        double ret = 0.0;
        for (int i = 3; i >= 0; i--) {
            ret = tx * ret + ((c[i][3] * ux + c[i][2]) * ux + c[i][1]) * ux + c[i][0];
        }
        return ret;
    }

    /**
     * <p>
     * bcuint1</p>
     *
     * @param x1l a double.
     * @param x1u a double.
     * @param y1l a double.
     * @param y1u a double.
     * @param t1 a double.
     * @param t2 a double.
     * @param ansy an array of double.
     * @return a double.
     */
    protected double bcuint1(double x1l, double x1u, double y1l, double y1u,
            double t1, double t2, double ansy[]) {
        double deltax = x1u - x1l;
        double deltay = y1u - y1l;
        bcucof(deltax, deltay, c);
        double tx = (t1 - x1l) / deltax;
        double ux = (t2 - y1l) / deltay;
        double ret = 0.0;
        ansy[0] = 0.0;
        ansy[1] = 0.0;
        for (int i = 3; i >= 0; i--) {
            ret = tx * ret + ((c[i][3] * ux + c[i][2]) * ux + c[i][1]) * ux + c[i][0];
            ansy[0] = ux * ansy[0] + (3.0 * c[3][i] * tx + 2.0 * c[2][i]) * tx + c[1][i];
            ansy[1] = tx * ansy[1] + (3.0 * c[i][3] * ux + 2.0 * c[i][2]) * ux + c[i][1];
        }
        ansy[0] /= deltax;
        ansy[1] /= deltay;
        return ret;
    }
    private static final double x16[] = new double[16];
    private static final double cl[] = new double[16];

    private static void bcucof(double t1, double t2, double c[][]) {
        double t1t2 = t1 * t2;
        /**
         * Pack a temporary vector of corner values.
         */
        for (int i = 0; i < 4; i++) {
            x16[i] = e[i];
            x16[i + 4] = dx[i] * t1;
            x16[i + 8] = dy[i] * t2;
            x16[i + 12] = dxy[i] * t1t2;
        }
        /**
         * Matrix multiply by the stored weight table.
         */
        for (int i = 0; i < 16; i++) {
            double xx = 0.0;
            for (int k = 0; k < 16; k++) {
                xx += wt[k][i] * x16[k];
            }
            cl[i] = xx;
        }
        /**
         * Unpack the results into the coefficient table.
         */
        int j = 0;
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                c[i][k] = cl[j++];
            }
        }
    }
    private static final double wt[][] = {
        {1.0, 0.0, -3.0, 2.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 9.0, -6.0,
            2.0, 0.0, -6.0, 4.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, -9.0, 6.0,
            -2.0, 0.0, 6.0, -4.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, -6.0, 0.0,
            0.0, -6.0, 4.0},
        {0.0, 0.0, 3.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.0, 6.0,
            0.0, 0.0, 6.0, -4.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -3.0, 2.0, -2.0, 0.0, 6.0, -4.0,
            1.0, 0.0, -3.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 3.0, -2.0,
            1.0, 0.0, -3.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 2.0, 0.0,
            0.0, 3.0, -2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -2.0, 0.0, 0.0, -6.0, 4.0,
            0.0, 0.0, 3.0, -2.0},
        {0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, -3.0,
            0.0, 2.0, -4.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -6.0, 3.0, 0.0,
            -2.0, 4.0, -2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 3.0, 0.0,
            0.0, 2.0, -2.0},
        {0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -3.0,
            0.0, 0.0, -2.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, -2.0, 4.0, -2.0,
            0.0, 1.0, -2.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, -1.0,
            0.0, 1.0, -2.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0,
            0.0, -1.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 2.0, -2.0,
            0.0, 0.0, -1.0, 1.0}};

    @Override
    public void setLambda(double lambda) {
        if (applyAllLambda()) {
            lambdaTerm = true;
            this.lambda = lambda;
        } else {
            lambdaTerm = false;
            this.lambda = 1.0;
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        if (lambdaTerm) {
            return dEdL;
        } else {
            return 0.0;
        }
    }

    @Override
    public double getd2EdL2() {
        return 0.0;
    }

    @Override
    public void getdEdXdL(double[] gradient) {
        return;
    }
}
