/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.AngleType.cubic;
import static ffx.potential.parameters.AngleType.quartic;
import static ffx.potential.parameters.AngleType.quintic;
import static ffx.potential.parameters.AngleType.sextic;
import static ffx.potential.parameters.AngleType.units;

/**
 * The Angle class represents an angle formed between three linearly bonded
 * atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class Angle extends BondedTerm implements Comparable<Angle> {

    private static final Logger logger = Logger.getLogger(Angle.class.getName());

    public enum AngleMode {

        NORMAL, IN_PLANE;
    }
    /**
     * Force field parameters to compute the angle bending energy.
     */
    public AngleType angleType;
    /**
     * AngleMode for this angle.
     */
    private AngleMode angleMode = AngleMode.NORMAL;
    /**
     * Number of hydrogens on the central atom that are not part of this Angle.
     */
    public int nh = 0;
    private double rigidScale = 1.0;
    private Atom atom4 = null;

    /**
     * Angle constructor
     *
     * @param b1 Bond that forms one leg of the angle
     * @param b2 Bond that forms the other leg of the angle
     */
    public Angle(Bond b1, Bond b2) {
        super();
        Atom a2 = b1.getCommonAtom(b2);
        Atom a1 = b1.get1_2(a2);
        Atom a3 = b2.get1_2(a2);
        b1.setAngleWith(b2);
        b2.setAngleWith(b1);
        atoms = new Atom[3];
        bonds = new Bond[2];
        atoms[1] = a2;
        atoms[0] = a1;
        atoms[2] = a3;
        bonds[0] = b1;
        bonds[1] = b2;

        a1.setAngle(this);
        a2.setAngle(this);
        a3.setAngle(this);
        setID_Key(false);
    }

    public static Angle angleFactory(Bond b1, Bond b2, ForceField forceField) {
        Angle newAngle = new Angle(b1, b2);
        Atom ac = b1.getCommonAtom(b2);
        Atom a1 = b1.get1_2(ac);
        Atom a3 = b2.get1_2(ac);
        int c[] = new int[3];
        c[0] = a1.getAtomType().atomClass;
        c[1] = ac.getAtomType().atomClass;
        c[2] = a3.getAtomType().atomClass;
        String key = AngleType.sortKey(c);
        AngleType angleType = forceField.getAngleType(key);
        if (angleType == null) {
            logger.severe("No AngleType for key: " + key);
            return null;
        }
        newAngle.setAngleType(angleType);
        return newAngle;
    }

    /**
     * <p>
     * Setter for the field <code>angleMode</code>.</p>
     *
     * @param mode a {@link ffx.potential.bonded.Angle.AngleMode} object.
     * @param a4 a {@link ffx.potential.bonded.Atom} object.
     */
    public void setAngleMode(AngleMode mode, Atom a4) {
        angleMode = mode;
        atom4 = a4;
    }

    /**
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a a {@link ffx.potential.parameters.AngleType} object.
     */
    public void setAngleType(AngleType a) {
        angleType = a;
        /**
         * Count the number of hydrogens attached to the central atom, but that
         * are not part of the angle.
         */
        ArrayList<Bond> ba = atoms[1].getBonds();
        nh = 0;
        for (Bond b1 : ba) {
            if (b1 != bonds[0] && b1 != bonds[1]) {
                Atom atom = b1.get1_2(atoms[1]);
                if (atom.getAtomType().atomicNumber == 1) {
                    nh += 1;
                }
            }
        }
        // Some angle bending parameters are generic for any number of hydrogens
        while (angleType.angle.length <= nh) {
            nh--;
        }
    }

    /**
     * <p>
     * Setter for the field <code>rigidScale</code>.</p>
     *
     * @param rigidScale a double.
     */
    public void setRigidScale(double rigidScale) {
        this.rigidScale = rigidScale;
    }

    /**
     * Get the AngleType for this angle.
     *
     * @return This angle's AngleType.
     */
    public AngleType getAngleType() {
        return angleType;
    }

    /**
     * If the specified atom is not the central atom of <b>this</b> angle, the
     * atom of the opposite leg is returned. These atoms are said to be 1-3 to
     * each other.
     *
     * @param a Atom
     * @return Atom
     */
    public Atom get1_3(Atom a) {
        if (a == atoms[0]) {
            return atoms[2];
        }
        if (a == atoms[2]) {
            return atoms[0];
        }
        return null;
    }

    public Atom getCentralAtom() {
        return atoms[1];
    }

    /**
     * Finds the common bond between <b>this</b> angle and another
     *
     * @param a An Angle that may have a common bond with <b>this</b> angle
     * @return The common Bond between this Angle and Angle a, or null if this
     * == a or no common bond exists
     */
    public Bond getCommonBond(Angle a) {
        // Comparing an angle to itself returns null
        // Comparing to a null angle return null
        if (a == this || a == null) {
            return null;
        }
        if (a.bonds[0] == bonds[0]) {
            return bonds[0];
        }
        if (a.bonds[0] == bonds[1]) {
            return bonds[1];
        }
        if (a.bonds[1] == bonds[0]) {
            return bonds[0];
        }
        if (a.bonds[1] == bonds[1]) {
            return bonds[1];
        }
        return null; // No common bond found
    }

    /**
     * Finds the other bond that makes up <b>this</b> angle
     *
     * @param b The bond to find the opposite of
     * @return The other Bond that makes up this Angle, or null if Bond b is not
     * part of this Angle
     */
    public Bond getOtherBond(Bond b) {
        if (b == bonds[0]) {
            return bonds[1];
        }
        if (b == bonds[1]) {
            return bonds[0];
        }
        return null; // b not found in angle
    }

    /**
     * If the central atom of the angle is trigonal, the 4th member of the
     * trigonal center (that is not a part of the angle) will be returned.
     *
     * @return The 4th atom of a trigonal center.
     */
    public Atom getTrigonalAtom() {
        if (atoms[1].isTrigonal()) {
            for (Bond b : atoms[1].getBonds()) {
                if (b != bonds[0] && b != bonds[1]) {
                    return b.get1_2(atoms[1]);
                }
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     *
     * Update recomputes <b>this</b> Angle's value and energy.
     */
    @Override
    public void update() {
        energy(false);
    }

    /**
     * Evaluate this Angle energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        double prefactor = units * rigidScale * angleType.forceConstant;
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        switch (angleType.angleFunction) {
            case SEXTIC:
                switch (angleMode) {
                    case NORMAL:
                        diff(a0, a1, v10);
                        diff(a2, a1, v12);
                        double rab2 = dot(v10, v10);
                        double rcb2 = dot(v12, v12);
                        if (rab2 != 0.0 && rcb2 != 0.0) {
                            cross(v12, v10, p);
                            double cosine = dot(v10, v12) / sqrt(rab2 * rcb2);
                            cosine = min(1.0, max(-1.0, cosine));
                            value = toDegrees(acos(cosine));
                            double dv = value - angleType.angle[nh];
                            double dv2 = dv * dv;
                            double dv3 = dv2 * dv;
                            double dv4 = dv2 * dv2;
                            energy = prefactor * dv2 * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4);
                            if (gradient) {
                                double deddt = prefactor * dv * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic * dv2
                                        + 5.0 * quintic * dv3 + 6.0 * sextic * dv4);
                                double rp = r(p);
                                rp = max(rp, 0.000001);
                                double terma = -deddt / (rab2 * rp);
                                double termc = deddt / (rcb2 * rp);
                                cross(v10, p, g0);
                                cross(v12, p, g2);
                                scalar(g0, terma, g0);
                                scalar(g2, termc, g2);
                                sum(g0, g2, g1);
                                scalar(g1, -1.0, g1);
                                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                            }
                            value = dv;
                        }
                        break;
                    case IN_PLANE:
                        atom4.getXYZ(a4);
                        diff(a0, a4, v10);
                        diff(a1, a4, v20);
                        diff(a2, a4, v30);
                        cross(v10, v30, p);
                        double rp2 = dot(p, p);
                        double delta = -dot(p, v20) / rp2;
                        scalar(p, delta, ip);
                        sum(ip, v20, ip);
                        diff(v10, ip, jp);
                        diff(v30, ip, kp);
                        double jp2 = dot(jp, jp);
                        double kp2 = dot(kp, kp);
                        if (jp2 != 0.0 && kp2 != 0.0) {
                            cross(kp, jp, lp);
                            double lpr = r(lp);
                            lpr = max(lpr, 0.000001);
                            double jk2 = dot(jp, kp);
                            double cosine = jk2 / Math.sqrt(jp2 * kp2);
                            cosine = min(1.0, max(-1.0, cosine));
                            value = toDegrees(acos(cosine));
                            double dv = value - angleType.angle[nh];
                            double dv2 = dv * dv;
                            double dv3 = dv2 * dv;
                            double dv4 = dv2 * dv2;
                            energy = prefactor * dv2 * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4);
                            if (gradient) {
                                double deddt = prefactor * dv * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic * dv2 + 5.0 * quintic * dv3 + 6.0 * sextic * dv4);
                                double term0 = -deddt / (jp2 * lpr);
                                double term2 = deddt / (kp2 * lpr);
                                cross(jp, lp, ded0);
                                scalar(ded0, term0, ded0);
                                cross(kp, lp, ded2);
                                scalar(ded2, term2, ded2);
                                sum(ded0, ded2, dedp);
                                scalar(dedp, -1.0, g1);
                                double delta2 = 2.0 * delta;
                                double pt2 = dot(dedp, p) / rp2;
                                cross(v20, v30, x21);
                                cross(p, v30, xp2);
                                cross(v30, g1, xd2);
                                scalar(xd2, delta, xd2);
                                scalar(xp2, delta2, xp2);
                                sum(x21, xp2, x21);
                                scalar(x21, pt2, x21);
                                sum(xd2, x21, dpd0);
                                cross(v10, v20, x01);
                                cross(p, v10, xp2);
                                cross(g1, v10, xd2);
                                scalar(xd2, delta, xd2);
                                scalar(xp2, delta2, xp2);
                                sum(x21, xp2, x21);
                                scalar(x21, pt2, x21);
                                sum(xd2, x21, dpd2);
                                sum(ded0, dpd0, g0);
                                sum(ded2, dpd2, g2);
                                sum(g0, g1, g3);
                                sum(g2, g3, g3);
                                scalar(g3, -1.0, g3);
                                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                                atom4.addToXYZGradient(g3[0], g3[1], g3[2]);
                            }
                            value = dv;
                        }
                        break;
                }
                break;
            case HARMONIC:
            default:
                switch (angleMode) {
                    case NORMAL:
                        diff(a0, a1, v10);
                        diff(a2, a1, v12);
                        double rab2 = dot(v10, v10);
                        double rcb2 = dot(v12, v12);
                        if (rab2 != 0.0 && rcb2 != 0.0) {
                            cross(v12, v10, p);
                            double cosine = dot(v10, v12) / sqrt(rab2 * rcb2);
                            cosine = min(1.0, max(-1.0, cosine));
                            value = toDegrees(acos(cosine));
                            double dv = value - angleType.angle[nh];
                            double dv2 = dv * dv;
                            energy = prefactor * dv2;
                            if (gradient) {
                                double deddt = prefactor * dv * toDegrees(2.0);
                                double rp = r(p);
                                rp = max(rp, 0.000001);
                                double terma = -deddt / (rab2 * rp);
                                double termc = deddt / (rcb2 * rp);
                                cross(v10, p, g0);
                                cross(v12, p, g2);
                                scalar(g0, terma, g0);
                                scalar(g2, termc, g2);
                                sum(g0, g2, g1);
                                scalar(g1, -1.0, g1);
                                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                            }
                            value = dv;
                        }
                        break;
                    case IN_PLANE:
                        atom4.getXYZ(a4);
                        diff(a0, a4, v10);
                        diff(a1, a4, v20);
                        diff(a2, a4, v30);
                        cross(v10, v30, p);
                        double rp2 = dot(p, p);
                        double delta = -dot(p, v20) / rp2;
                        scalar(p, delta, ip);
                        sum(ip, v20, ip);
                        diff(v10, ip, jp);
                        diff(v30, ip, kp);
                        double jp2 = dot(jp, jp);
                        double kp2 = dot(kp, kp);
                        if (jp2 != 0.0 && kp2 != 0.0) {
                            cross(kp, jp, lp);
                            double lpr = r(lp);
                            lpr = max(lpr, 0.000001);
                            double jk2 = dot(jp, kp);
                            double cosine = jk2 / Math.sqrt(jp2 * kp2);
                            cosine = min(1.0, max(-1.0, cosine));
                            value = toDegrees(acos(cosine));
                            double dv = value - angleType.angle[nh];
                            double dv2 = dv * dv;
                            energy = prefactor * dv2;
                            if (gradient) {
                                double deddt = prefactor * dv * toDegrees(2.0);
                                double term0 = -deddt / (jp2 * lpr);
                                double term2 = deddt / (kp2 * lpr);
                                cross(jp, lp, ded0);
                                scalar(ded0, term0, ded0);
                                cross(kp, lp, ded2);
                                scalar(ded2, term2, ded2);
                                sum(ded0, ded2, dedp);
                                scalar(dedp, -1.0, g1);
                                double delta2 = 2.0 * delta;
                                double pt2 = dot(dedp, p) / rp2;
                                cross(v20, v30, x21);
                                cross(p, v30, xp2);
                                cross(v30, g1, xd2);
                                scalar(xd2, delta, xd2);
                                scalar(xp2, delta2, xp2);
                                sum(x21, xp2, x21);
                                scalar(x21, pt2, x21);
                                sum(xd2, x21, dpd0);
                                cross(v10, v20, x01);
                                cross(p, v10, xp2);
                                cross(g1, v10, xd2);
                                scalar(xd2, delta, xd2);
                                scalar(xp2, delta2, xp2);
                                sum(x21, xp2, x21);
                                scalar(x21, pt2, x21);
                                sum(xd2, x21, dpd2);
                                sum(ded0, dpd0, g0);
                                sum(ded2, dpd2, g2);
                                sum(g0, g1, g3);
                                sum(g2, g3, g3);
                                scalar(g3, -1.0, g3);
                                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                                atom4.addToXYZGradient(g3[0], g3[1], g3[2]);
                            }
                            value = dv;
                        }
                        break;
                }
                break;
        }
        return energy;
    }

    /**
     * Log details for this Angle energy term.
     */
    public void log() {
        switch (angleMode) {
            case NORMAL:
                logger.info(String.format(
                        " %8s %6d-%s %6d-%s %6d-%s %7.4f  %7.4f  %10.4f", "Angle",
                        atoms[0].getXYZIndex(), atoms[0].getAtomType().name,
                        atoms[1].getXYZIndex(), atoms[1].getAtomType().name,
                        atoms[2].getXYZIndex(), atoms[2].getAtomType().name,
                        angleType.angle[nh], value, energy));
                break;
            case IN_PLANE:
                logger.info(String.format(
                        " %8s %6d-%s %6d-%s %6d-%s %7.4f  %7.4f  %10.4f", "Angle-IP",
                        atoms[0].getXYZIndex(), atoms[0].getAtomType().name,
                        atoms[1].getXYZIndex(), atoms[1].getAtomType().name,
                        atoms[2].getXYZIndex(), atoms[2].getAtomType().name,
                        angleType.angle[nh], value, energy));
                break;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compareTo(Angle a) {
        if (a == null) {
            throw new NullPointerException();
        }
        if (a == this) {
            return 0;
        }
        int this1 = atoms[1].xyzIndex;
        int a1 = a.atoms[1].xyzIndex;
        if (this1 < a1) {
            return -1;
        }
        if (this1 > a1) {
            return 1;
        }
        int this0 = atoms[0].xyzIndex;
        int a0 = a.atoms[0].xyzIndex;
        if (this0 < a0) {
            return -1;
        }
        if (this0 > a0) {
            return 1;
        }
        int this2 = atoms[2].xyzIndex;
        int a2 = a.atoms[2].xyzIndex;
        if (this2 < a2) {
            return -1;
        }
        if (this2 > a2) {
            return 1;
        }
        // There should never be duplicate, identical angle objects.
        assert (!(this0 == a0 && this1 == a1 && this2 == a2));
        return 0;
    }
    protected static final double a0[] = new double[3];
    protected static final double a1[] = new double[3];
    protected static final double a2[] = new double[3];
    protected static final double a4[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 0.
     */
    protected static final double v10[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 2.
     */
    protected static final double v12[] = new double[3];
    /**
     * Vector from Atom 3 to Atom 0.
     */
    protected static final double v30[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 0.
     */
    protected static final double v20[] = new double[3];
    /**
     * Vector v10 cross v30.
     */
    protected static final double p[] = new double[3];
    /**
     * Constant <code>ip=new double[3]</code>
     */
    protected static final double ip[] = new double[3];
    /**
     * Constant <code>jp=new double[3]</code>
     */
    protected static final double jp[] = new double[3];
    /**
     * Constant <code>kp=new double[3]</code>
     */
    protected static final double kp[] = new double[3];
    /**
     * Constant <code>lp=new double[3]</code>
     */
    protected static final double lp[] = new double[3];
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
     * Constant <code>g3=new double[3]</code>
     */
    protected static final double g3[] = new double[3];
    private static final double ded0[] = new double[3];
    private static final double ded2[] = new double[3];
    private static final double dedp[] = new double[3];
    private static final double x21[] = new double[3];
    private static final double x01[] = new double[3];
    private static final double xp2[] = new double[3];
    private static final double xd2[] = new double[3];
    private static final double dpd0[] = new double[3];
    private static final double dpd2[] = new double[3];
}
