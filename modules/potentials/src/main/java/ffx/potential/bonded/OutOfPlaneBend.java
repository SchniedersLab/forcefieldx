/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.potential.bonded;

import static java.lang.Math.*;

import static ffx.numerics.VectorMath.*;
import static ffx.potential.parameters.OutOfPlaneBendType.*;

import java.util.logging.Logger;

import ffx.potential.parameters.OutOfPlaneBendType;

/**
 * The OutOfPlaneBend class represents an Out-Of-Plane Bend.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class OutOfPlaneBend extends BondedTerm implements
        Comparable<OutOfPlaneBend> {

    private static final Logger logger = Logger.getLogger(OutOfPlaneBend.class.getName());

    /**
     * Force field parameters to compute the Out-of-Plane Bend energy.
     */
    public OutOfPlaneBendType outOfPlaneBendType = null;

    /**
     * OutOfPlaneBend constructor.
     *
     * @param angle
     *            Angle that contains 3 of 4 OutOfPlaneBend atoms.
     * @param atom
     *            The 4th atom of the trigonal center.
     */
    public OutOfPlaneBend(Angle angle, Atom atom) {
        super();
        atoms = new Atom[4];
        atoms[0] = angle.atoms[0];
        atoms[1] = angle.atoms[1];
        atoms[2] = angle.atoms[2];
        atoms[3] = atom;
        bonds = new Bond[3];
        bonds[0] = angle.bonds[0];
        bonds[1] = angle.bonds[1];
        bonds[2] = atoms[1].getBond(atom);
        setID_Key(false);
    }

    /**
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a
     */
    public void setAngleType(OutOfPlaneBendType a) {
        outOfPlaneBendType = a;
    }

    /**
     * Update recomputes OutOfPlaneBend value and energy.
     */
    @Override
    public void update() {
        energy(false);
    }

    /**
     * Evaluate this Out-of-Plane Bend energy.
     *
     * @param gradient
     *            Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        diff(atoms[0].getXYZ(), atoms[1].getXYZ(), v10);
        diff(atoms[2].getXYZ(), atoms[1].getXYZ(), v12);
        diff(atoms[3].getXYZ(), atoms[1].getXYZ(), v13);
        diff(atoms[0].getXYZ(), atoms[3].getXYZ(), v30);
        diff(atoms[2].getXYZ(), atoms[3].getXYZ(), v32);
        double rdb2 = dot(v13, v13);
        double rad2 = dot(v30, v30);
        double rcd2 = dot(v32, v32);
        cross(v12, v13, p);
        double ee = dot(v10, p);
        double rac2 = dot(v30, v32);
        double cc = rad2 * rcd2 - rac2 * rac2;
        if (rdb2 != 0.0 && cc != 0.0) {
            double bkk2 = rdb2 - ee * ee / cc;
            double cosine = sqrt(bkk2 / rdb2);
            cosine = min(1.0, max(-1.0, cosine));
            value = toDegrees(acos(cosine));
            double dv = value;
            double dv2 = dv * dv;
            double dv3 = dv2 * dv;
            double dv4 = dv2 * dv2;
            energy = units
                    * outOfPlaneBendType.forceConstant
                    * dv2
                    * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic
                    * dv4);
            if (gradient) {
                double deddt = units
                        * outOfPlaneBendType.forceConstant
                        * dv
                        * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic
                        * dv2 + 5.0 * quintic * dv3 + 6.0 * sextic
                        * dv4);
                double dedcos = 0.0;
                if (ee != 0.0) {
                    dedcos = -deddt * signum(ee) / sqrt(cc * bkk2);
                } else {
                    dedcos = -deddt / sqrt(cc * bkk2);
                }
                double term = ee / cc;
                scalar(v30, rcd2, sv30);
                scalar(v32, rac2, sv32);
                diff(sv30, sv32, dcda);
                scalar(dcda, term, dcda);
                scalar(v32, rad2, sv32);
                scalar(v30, rac2, sv30);
                diff(sv32, sv30, dcdc);
                scalar(dcdc, term, dcdc);
                sum(dcda, dcdc, dcdd);
                scalar(dcdd, -1.0, dcdd);
                term = ee / rdb2;
                cross(v13, v12, deda);
                cross(v10, v13, dedc);
                cross(v12, v10, dedd);
                scalar(v13, term, v13);
                sum(dedd, v13, dedd);
                sum(dcda, deda, g0);
                sum(dcdc, dedc, g2);
                sum(dcdd, dedd, g3);
                scalar(g0, dedcos, g0);
                scalar(g2, dedcos, g2);
                scalar(g3, dedcos, g3);
                sum(g0, g2, g1);
                sum(g1, g3, g1);
                scalar(g1, -1.0, g1);
                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
            }
        }
        return energy;
    }

    /**
     * Log details for this Out-of-Plane Bend energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6.4f %10.4f",
                "Out-of-Plane Bend", atoms[1].getXYZIndex(), atoms[1].getAtomType().name, atoms[3].getXYZIndex(), atoms[3].getAtomType().name, value, energy));
    }

    /**
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

    @Override
    public int compareTo(OutOfPlaneBend o) {
        if (o == null) {
            throw new NullPointerException();
        }
        if (o == this) {
            return 0;
        }
        int this1 = atoms[1].xyzIndex;
        int a1 = o.atoms[1].xyzIndex;
        if (this1 < a1) {
            return -1;
        }
        if (this1 > a1) {
            return 1;
        }
        int this3 = atoms[3].xyzIndex;
        int a3 = o.atoms[3].xyzIndex;
        if (this3 < a3) {
            return -1;
        }
        if (this3 > a3) {
            return 1;
        }
        return 0;
    }

    /**
     * Vector from Atom 1 to Atom 0.
     */
    protected static final double v10[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 2.
     */
    protected static final double v12[] = new double[3];
    /**
     * Vector from Atom 1 to Atom 3.
     */
    protected static final double v13[] = new double[3];
    /**
     * Vector from Atom 3 to Atom 0.
     */
    protected static final double v30[] = new double[3];
    /**
     * Vector from Atom 3 to Atom 2.
     */
    protected static final double v32[] = new double[3];
    /**
     * Vector v12 cross v13.
     */
    protected static final double p[] = new double[3];
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
    protected static final double sv30[] = new double[3];
    /**
     * Work array.
     */
    protected static final double sv32[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dcda[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dcdc[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dcdd[] = new double[3];
    /**
     * Work array.
     */
    protected static final double deda[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedc[] = new double[3];
    /**
     * Work array.
     */
    protected static final double dedd[] = new double[3];

}
