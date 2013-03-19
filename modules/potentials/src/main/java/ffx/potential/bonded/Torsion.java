/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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

import static java.lang.Math.*;

import ffx.potential.parameters.TorsionType;

import static ffx.numerics.VectorMath.*;
import static ffx.potential.parameters.TorsionType.units;

/**
 * The Torsion class represents a torsional angle formed between four bonded
 * atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class Torsion extends BondedTerm {

    private static final Logger logger = Logger.getLogger(Torsion.class.getName());
    private static final long serialVersionUID = 1L;
    public TorsionType torsionType = null;

    /**
     * Torsion contructor
     *
     * @param an1 Angle that combines to form the Torsional Angle
     * @param an2 Angle that combines to form the Torsional Angle
     */
    public Torsion(Angle an1, Angle an2) {
        super();
        bonds = new Bond[3];
        bonds[1] = an1.getCommonBond(an2);
        bonds[0] = an1.getOtherBond(bonds[1]);
        bonds[2] = an2.getOtherBond(bonds[1]);
        Initialize();
    }

    /**
     * <p>compare</p>
     *
     * @param a0 a {@link ffx.potential.bonded.Atom} object.
     * @param a1 a {@link ffx.potential.bonded.Atom} object.
     * @param a2 a {@link ffx.potential.bonded.Atom} object.
     * @param a3 a {@link ffx.potential.bonded.Atom} object.
     * @return a boolean.
     */
    public boolean compare(Atom a0, Atom a1, Atom a2, Atom a3) {
        if (a0 == atoms[0] && a1 == atoms[1] && a2 == atoms[2]
                && a3 == atoms[3]) {
            return true;
        } else if (a0 == atoms[3] && a1 == atoms[2] && a2 == atoms[1]
                && a3 == atoms[0]) {
            return true;
        }
        return false;
    }

    /**
     * Torsion contructor
     *
     * @param a Angle that has one Atom in common with Bond b
     * @param b Bond that has one Atom in common with Angle A
     */
    public Torsion(Angle a, Bond b) {
        super();
        bonds = new Bond[3];
        bonds[0] = b;
        bonds[1] = a.getBond(0);
        bonds[2] = a.getBond(1);
        // See if bond 2 or bond 3 is the middle bond
        Atom atom = bonds[1].getCommonAtom(b);
        if (atom == null) {
            Bond temp = bonds[1];
            bonds[1] = bonds[2];
            bonds[2] = temp;
        }
        Initialize();
    }

    /**
     * Create a Torsion from 3 connected bonds (no error checking)
     *
     * @param b1 Bond
     * @param b2 Bond
     * @param b3 Bond
     */
    public Torsion(Bond b1, Bond b2, Bond b3) {
        super();
        bonds = new Bond[3];
        bonds[0] = b1;
        bonds[1] = b2;
        bonds[2] = b3;
        Initialize();
    }

    /**
     * Torsion Constructor
     *
     * @param n Torsion id
     */
    public Torsion(String n) {
        super(n);
    }

    // Initialization
    private void Initialize() {
        atoms = new Atom[4];
        atoms[1] = bonds[0].getCommonAtom(bonds[1]);
        atoms[0] = bonds[0].get1_2(atoms[1]);
        atoms[2] = bonds[1].get1_2(atoms[1]);
        atoms[3] = bonds[2].get1_2(atoms[2]);
        atoms[0].setTorsion(this);
        atoms[1].setTorsion(this);
        atoms[2].setTorsion(this);
        atoms[3].setTorsion(this);
        setID_Key(false);
    }

    /**
     * If the specified atom is not a central atom of <b>this</b> torsion, the
     * atom at the opposite end is returned. These atoms are said to be 1-4 to
     * each other.
     *
     * @param a Atom
     * @return Atom
     */
    public Atom get1_4(Atom a) {
        if (a == atoms[0]) {
            return atoms[3];
        }
        if (a == atoms[3]) {
            return atoms[0];
        }
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update() {
        energy(false);
    }
    /**
     * Vector from Atom 0 to Atom 1.
     */
    protected static final double v01[] = new double[3];
    /**
     * Vector from Atom 0 to Atom 2.
     */
    protected static final double v02[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 2.
     */
    protected static final double v12[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 3.
     */
    protected static final double v13[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 3.
     */
    protected static final double v23[] = new double[3];
    /**
     * Vector v01 cross v12.
     */
    protected static final double x0112[] = new double[3];
    /**
     * Vector v12 cross v23.
     */
    protected static final double x1223[] = new double[3];
    /**
     * Vector x0112 cross x12_32.
     */
    protected static final double x[] = new double[3];
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
     * Work array.
     */
    protected static final double x1[] = new double[3];
    /**
     * Work array.
     */
    protected static final double x2[] = new double[3];

    /**
     * Evaluate the Torsional Angle energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        diff(atoms[1].getXYZ(), atoms[0].getXYZ(), v01);
        diff(atoms[2].getXYZ(), atoms[1].getXYZ(), v12);
        diff(atoms[3].getXYZ(), atoms[2].getXYZ(), v23);
        cross(v01, v12, x0112);
        cross(v12, v23, x1223);
        cross(x0112, x1223, x);
        double r01_12 = dot(x0112, x0112);
        double r12_23 = dot(x1223, x1223);
        double rr = sqrt(r01_12 * r12_23);
        if (rr != 0.0) {
            double r12 = r(v12);
            double cosine = dot(x0112, x1223) / rr;
            double sine = dot(v12, x) / (r12 * rr);
            value = toDegrees(acos(cosine));
            if (sine < 0.0) {
                value = -value;
            }
            double amp[] = torsionType.amplitude;
            double tsin[] = torsionType.sine;
            double tcos[] = torsionType.cosine;
            energy = amp[0] * (1.0 + cosine * tcos[0] + sine * tsin[0]);
            double dedphi = amp[0] * (cosine * tsin[0] - sine * tcos[0]);
            double cosprev = cosine;
            double sinprev = sine;
            double n = torsionType.terms;
            for (int i = 1; i < n; i++) {
                double cosn = cosine * cosprev - sine * sinprev;
                double sinn = sine * cosprev + cosine * sinprev;
                double phi = 1.0 + cosn * tcos[i] + sinn * tsin[i];
                double dphi = (1.0 + i) * (cosn * tsin[i] - sinn * tcos[i]);
                energy = energy + amp[i] * phi;
                dedphi = dedphi + amp[i] * dphi;
                cosprev = cosn;
                sinprev = sinn;
            }
            energy = units * energy;
            if (gradient) {
                dedphi = units * dedphi;
                diff(atoms[2].getXYZ(), atoms[0].getXYZ(), v02);
                diff(atoms[3].getXYZ(), atoms[1].getXYZ(), v13);
                cross(x0112, v12, x1);
                cross(x1223, v12, x2);
                scalar(x1, dedphi / (r01_12 * r12), x1);
                scalar(x2, -dedphi / (r12_23 * r12), x2);
                cross(x1, v12, g0);
                cross(v02, x1, g1);
                cross(x2, v23, g2);
                sum(g1, g2, g1);
                cross(x1, v01, g2);
                cross(v13, x2, g3);
                sum(g2, g3, g2);
                cross(x2, v12, g3);
                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
            }
        }
        return energy;
    }

    /**
     * Log details for this Torsional Angle energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f",
                "Torsional-Angle", atoms[0].getXYZIndex(), atoms[0].getAtomType().name, atoms[1].getXYZIndex(), atoms[1].getAtomType().name, atoms[2].getXYZIndex(), atoms[2].getAtomType().name, atoms[3].getXYZIndex(), atoms[3].getAtomType().name, energy));
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
}
