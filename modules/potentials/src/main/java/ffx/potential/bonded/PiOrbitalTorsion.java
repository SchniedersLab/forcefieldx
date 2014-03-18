/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.potential.bonded;

import java.util.logging.Logger;

import static java.lang.Math.acos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.PiTorsionType;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.PiTorsionType.units;

/**
 * The Pi-Orbital Torsion class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class PiOrbitalTorsion extends BondedTerm {

    private static final Logger logger = Logger.getLogger(PiOrbitalTorsion.class.getName());
    private static final long serialVersionUID = 1L;
    public PiTorsionType piTorsionType = null;

    /**
     * Pi-Orbital Torsion constructor.
     *
     * @param middleBond a {@link ffx.potential.bonded.Bond} object.
     */
    public PiOrbitalTorsion(Bond middleBond) {
        super();
        atoms = new Atom[6];
        atoms[2] = middleBond.atoms[0];
        atoms[3] = middleBond.atoms[1];
        bonds = new Bond[5];
        bonds[2] = middleBond;
        int i = 0;
        for (Bond b : atoms[2].getBonds()) {
            if (b != middleBond) {
                atoms[i] = b.get1_2(atoms[2]);
                bonds[i] = b;
                i++;
            }
        }
        i = 4;
        for (Bond b : atoms[3].getBonds()) {
            if (b != middleBond) {
                atoms[i] = b.get1_2(atoms[3]);
                bonds[i - 1] = b;
                i++;
            }
        }
        setID_Key(false);
    }

    /**
     * Attempt to create a new PiOrbitalTorsion based on the supplied bond and
     * forceField.
     *
     * @param bond
     * @param forceField
     * @return a new PiOrbitalToersion, or null.
     */
    public static PiOrbitalTorsion piOrbitalTorsionFactory(Bond bond, ForceField forceField) {
        Atom atom1 = bond.getAtom(0);
        Atom atom2 = bond.getAtom(1);
        int c[] = new int[2];
        if (!atom1.isTrigonal() || !atom2.isTrigonal()) {
            return null;
        }
        c[0] = atom1.getAtomType().atomClass;
        c[1] = atom2.getAtomType().atomClass;
        String key = PiTorsionType.sortKey(c);
        PiTorsionType piTorsionType = forceField.getPiTorsionType(key);
        if (piTorsionType == null) {
            return null;
        }
        PiOrbitalTorsion piOrbitalTorsion = new PiOrbitalTorsion(
                bond);
        piOrbitalTorsion.piTorsionType = piTorsionType;
        return piOrbitalTorsion;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update() {
        energy(false);
    }
    /**
     * Vector from Atom 3 to Atom 0.
     */
    protected static final double v30[] = new double[3];
    /**
     * Vector from Atom 3 to Atom 1.
     */
    protected static final double v31[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 3.
     */
    protected static final double v23[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 4.
     */
    protected static final double v24[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 5.
     */
    protected static final double v25[] = new double[3];
    /**
     * Vector from Atom 5 to Atom 2.
     */
    protected static final double v[] = new double[3];
    /**
     * Vector v30 cross v13.
     */
    protected static final double vp[] = new double[3];
    /**
     * Vector v25 cross v52.
     */
    protected static final double vq[] = new double[3];
    /**
     * Vector xt.
     */
    protected static final double xt[] = new double[3];
    /**
     * Vector xu.
     */
    protected static final double xu[] = new double[3];
    /**
     * Vector xtu.
     */
    protected static final double xtu[] = new double[3];
    /**
     * Work array.
     */
    protected static final double vpc[] = new double[3];
    /**
     * Work array.
     */
    protected static final double vdq[] = new double[3];
    /**
     * Work array.
     */
    protected static final double vpd[] = new double[3];
    /**
     * Work array.
     */
    protected static final double vcq[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedt[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedu[] = new double[3];
    /**
     * Work array.
     */
    protected static final double temp[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedp[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedc[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedd[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedq[] = new double[3];
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
     * Gradient on Atom 4.
     */
    protected static final double g4[] = new double[3];
    /**
     * Gradient on Atom 5.
     */
    protected static final double g5[] = new double[3];

    /**
     * Evaluate the Pi-Orbital Torsion energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        diff(atoms[0].getXYZ(), atoms[3].getXYZ(), v30);
        diff(atoms[1].getXYZ(), atoms[3].getXYZ(), v31);
        diff(atoms[3].getXYZ(), atoms[2].getXYZ(), v23);
        diff(atoms[4].getXYZ(), atoms[2].getXYZ(), v24);
        diff(atoms[5].getXYZ(), atoms[2].getXYZ(), v25);
        cross(v30, v31, vp);
        cross(v24, v25, vq);
        sum(vp, atoms[2].getXYZ(), vp);
        sum(vq, atoms[3].getXYZ(), vq);
        diff(atoms[2].getXYZ(), vp, vpc);
        diff(vq, atoms[3].getXYZ(), vdq);
        cross(vpc, v23, xt);
        cross(v23, vdq, xu);
        cross(xt, xu, xtu);
        double rt2 = dot(xt, xt);
        double ru2 = dot(xu, xu);
        double rr = sqrt(rt2 * ru2);
        if (rr != 0.0) {
            double r23 = r(v23);
            double cosine = dot(xt, xu) / rr;
            double sine = dot(v23, xtu) / (r23 * rr);
            cosine = min(1.0, max(-1.0, cosine));
            value = toDegrees(acos(cosine));
            if (sine < 0.0) {
                value = -value;
            }
            if (value > 90.0) {
                value -= 180.0;
            }
            if (value < -90.0) {
                value += 180.0;
            }
            double cosine2 = cosine * cosine - sine * sine;
            double phi2 = 1.0 - cosine2;
            energy = units * piTorsionType.forceConstant * phi2;
            if (gradient) {
                double sine2 = 2.0 * cosine * sine;
                double dphi2 = 2.0 * sine2;
                double dedphi = units * piTorsionType.forceConstant * dphi2;
                diff(atoms[3].getXYZ(), vp, vpd);
                diff(vq, atoms[2].getXYZ(), vcq);
                cross(xt, v23, dedt);
                cross(xu, v23, dedu);
                scalar(dedt, dedphi / (rt2 * r23), dedt);
                scalar(dedu, -dedphi / (ru2 * r23), dedu);
                cross(dedt, v23, dedp);
                cross(vpd, dedt, dedc);
                cross(dedu, vdq, temp);
                sum(dedc, temp, dedc);
                cross(dedt, vpc, dedd);
                cross(vcq, dedu, temp);
                sum(dedd, temp, dedd);
                cross(dedu, v23, dedq);
                cross(v31, dedp, g0);
                cross(dedp, v30, g1);
                cross(v25, dedq, g4);
                cross(dedq, v24, g5);
                sum(dedc, dedp, g2);
                sum(g4, g5, temp);
                diff(g2, temp, g2);
                sum(dedd, dedq, g3);
                sum(g0, g1, temp);
                diff(g3, temp, g3);
                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
                atoms[4].addToXYZGradient(g4[0], g4[1], g4[2]);
                atoms[5].addToXYZGradient(g5[0], g5[1], g5[2]);
            }
        }
        return energy;
    }

    /**
     * Log details for this Pi-Orbital Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %10.4f %10.4f",
                "Pi-Orbital Torsion", atoms[2].getXYZIndex(), atoms[2].getAtomType().name,
                atoms[3].getXYZIndex(), atoms[3].getAtomType().name, value, energy));
    }

    /**
     * {@inheritDoc}
     *
     * Over-ridden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }
}
