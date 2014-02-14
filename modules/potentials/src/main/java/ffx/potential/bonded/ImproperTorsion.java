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

import static java.lang.Math.acos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;

import ffx.potential.parameters.ImproperTorsionType;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;

/**
 * The OutOfPlaneBend class represents an Out-Of-Plane Bend.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class ImproperTorsion extends BondedTerm implements
        Comparable<ImproperTorsion> {

    private static final Logger logger = Logger.getLogger(OutOfPlaneBend.class.getName());
    /**
     * Force field parameters to compute the Out-of-Plane Bend energy.
     */
    public ImproperTorsionType improperType = null;

    /**
     * OutOfPlaneBend constructor.
     *
     * @param atom1
     * @param atom2
     * @param atom4
     * @param atom3
     */
    public ImproperTorsion(Atom atom1, Atom atom2, Atom atom3, Atom atom4) {
        super();
        atoms = new Atom[4];
        atoms[0] = atom1;
        atoms[1] = atom2;
        atoms[2] = atom3;
        atoms[3] = atom4;
        /*
         bonds = new Bond[3];
         bonds[0] = angle.bonds[0];
         bonds[1] = angle.bonds[1];
         bonds[2] = atoms[1].getBond(atom);
         */
        setID_Key(false);
    }

    /**
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     */
    public void setImproperType(ImproperTorsionType a) {
        improperType = a;
    }

    /**
     * {@inheritDoc}
     *
     * Update recomputes OutOfPlaneBend value and energy.
     */
    @Override
    public void update() {
        energy(false);
    }

    /**
     * Evaluate this Out-of-Plane Bend energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        diff(atoms[1].getXYZ(), atoms[0].getXYZ(), v10);
        diff(atoms[2].getXYZ(), atoms[1].getXYZ(), v21);
        diff(atoms[3].getXYZ(), atoms[2].getXYZ(), v32);
        cross(v10, v21, t);
        cross(v21, v32, u);
        cross(t, u, tu);
        double rt2 = dot(t, t);
        double ru2 = dot(u, u);
        double rtru = sqrt(rt2 * ru2);
        if (rtru != 0.0) {
            double rcb = r(v21);
            double cosine = dot(t, u) / rtru;
            double sine = dot(v21, tu) / (rcb * rtru);

            /* Set the improper torsional parameters for this angle */
            double v1 = 0.0;
            double c1 = 0.0;
            double s1 = 0.0;
            double v2 = 0.0;
            double c2 = 0.0;
            double s2 = 0.0;
            double v3 = 0.0;
            double c3 = 0.0;
            double s3 = 0.0;

            /* compute the multiple angle trigonometry and the phase terms */
            double cosine2 = cosine * cosine - sine * sine;
            double sine2 = 2.0 * cosine * sine;
            double cosine3 = cosine * cosine2 - sine * sine2;
            double sine3 = cosine * sine2 + sine * cosine2;
            double phi1 = 1.0 + (cosine * c1 + sine * s1);
            double phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);
            double phi3 = 1.0 + (cosine3 * c3 + sine3 * s3);
            double dphi1 = (cosine * s1 - sine * c1);
            double dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2);
            double dphi3 = 3.0 * (cosine3 * s3 - sine3 * c3);

            /* calculate improper torsion energy and master chain rule term */
            value = phi2;
            energy = ImproperTorsionType.units * (v1 * phi1 + v2 * phi2 + v3 * phi3);
            double dedphi = ImproperTorsionType.units * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);

            /* chain rule terms for first derivative components */
            if (gradient) {
                diff(atoms[2].getXYZ(), atoms[0].getXYZ(), v20);
                diff(atoms[3].getXYZ(), atoms[1].getXYZ(), v31);
                cross(t,v21,dedt);
                cross(u,v21,dedu);
                scalar(dedt,dedphi / (rt2 * rcb),dedt);
                scalar(dedu,-dedphi / (ru2 * rcb),dedu);
                /**
                 * Compute first derivative components for this angle.
                 */
                cross(v21,dedt,deda);
                //cross(v20,dedt,dedb1);
                //cross(v32,dedu,dedb2);
                //scalar(dedb2, -1.0, dedb2);
                //cross(v10,dedt,dedc1);
                //cross(v31,dedu,dedc2);
                //scalar(dedc2, -1.0, dedc2);
                cross(v21,dedu,dedd);
                /*
                double dedxia = zcb * dedyt - ycb * dedzt;
                double dedyia = xcb * dedzt - zcb * dedxt;
                double dedzia = ycb * dedxt - xcb * dedyt;

                double dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * dedzu;
                double dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * dedxu;
                double dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * dedyu;

                double dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * dedyu;
                double dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * dedzu;
                double dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * dedxu;

                double dedxid = zcb * dedyu - ycb * dedzu;
                double dedyid = xcb * dedzu - zcb * dedxu;
                double dedzid = ycb*dedxu - xcb*dedyu;
                */
                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
            }
        }
        return energy;
    }

    /**
     * Log details for this Improper Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6.4f %10.4f",
                "Improper Torsion", atoms[1].getXYZIndex(), atoms[1].getAtomType().name, atoms[3].getXYZIndex(), atoms[3].getAtomType().name, value, energy));
    }

    /**
     * {@inheritDoc}
     *
     * Overridden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compareTo(ImproperTorsion o) {
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
     * Vector from Atom 2 to Atom 1.
     */
    protected static final double v21[] = new double[3];
    /**
     * Vector from Atom 3 to Atom 2.
     */
    protected static final double v32[] = new double[3];
    protected static final double v31[] = new double[3];
    protected static final double v20[] = new double[3];
    protected static final double t[] = new double[3];
    protected static final double u[] = new double[3];
    protected static final double tu[] = new double[3];
    protected static final double dedu[] = new double[3];
    protected static final double dedt[] = new double[3];

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
