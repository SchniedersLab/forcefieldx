/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.signum;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.AtomicDoubleArray;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.OutOfPlaneBendType;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.OutOfPlaneBendType.cubic;
import static ffx.potential.parameters.OutOfPlaneBendType.quartic;
import static ffx.potential.parameters.OutOfPlaneBendType.quintic;
import static ffx.potential.parameters.OutOfPlaneBendType.sextic;
import static ffx.potential.parameters.OutOfPlaneBendType.units;

/**
 * The OutOfPlaneBend class represents an Out-Of-Plane Bend.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class OutOfPlaneBend extends BondedTerm {

    private static final Logger logger = Logger.getLogger(OutOfPlaneBend.class.getName());
    /**
     * Force field parameters to compute the Out-of-Plane Bend energy.
     */
    public OutOfPlaneBendType outOfPlaneBendType = null;

    /**
     * OutOfPlaneBend constructor.
     *
     * @param angle Angle that contains 3 of 4 OutOfPlaneBend atoms.
     * @param atom The 4th atom of the trigonal center.
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
     * Attempt to create a new OutOfPlaneBend instance for a given Angle and
     * Force Field.
     *
     * @param angle the Angle to create an OutOfPlaneBend around.
     * @param forceField the ForceField parameters to use.
     * @return a new OutOfPlaneBend if the central atom of the angle is trigonal
     * and a force field type exists.
     */
    public static OutOfPlaneBend outOfPlaneBendFactory(Angle angle, ForceField forceField) {
        Atom centralAtom = angle.atoms[1];
        if (centralAtom.isTrigonal()) {
            Atom atom4 = angle.getTrigonalAtom();
            String key = atom4.getAtomType().atomClass + " " + centralAtom.getAtomType().atomClass
                    + " 0 0";
            OutOfPlaneBendType outOfPlaneBendType = forceField.getOutOfPlaneBendType(key);
            if (outOfPlaneBendType != null) {
                angle.setAngleMode(Angle.AngleMode.IN_PLANE, atom4);
                OutOfPlaneBend newOutOfPlaneBend = new OutOfPlaneBend(
                        angle, atom4);
                newOutOfPlaneBend.outOfPlaneBendType = outOfPlaneBendType;
                return newOutOfPlaneBend;
            }
        }
        return null;
    }

    /**
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     */
    public void setAngleType(OutOfPlaneBendType a) {
        outOfPlaneBendType = a;
    }

    /**
     * Evaluate this Out-of-Plane Bend energy.
     *
     * @param gradient Evaluate the gradient.
     * @param threadID
     * @param gradX
     * @param gradY
     * @param gradZ
     * @return Returns the energy.
     */
    @Override
    public double energy(boolean gradient, int threadID,
            AtomicDoubleArray gradX,
            AtomicDoubleArray gradY,
            AtomicDoubleArray gradZ,
            AtomicDoubleArray lambdaGradX,
            AtomicDoubleArray lambdaGradY,
            AtomicDoubleArray lambdaGradZ) {

        double a0[] = new double[3];
        double a1[] = new double[3];
        double a2[] = new double[3];
        double a3[] = new double[3];
        /**
         * Vector from Atom 1 to Atom 0.
         */
        double v10[] = new double[3];
        /**
         * Vector from Atom 1 to Atom 2.
         */
        double v12[] = new double[3];
        /**
         * Vector from Atom 1 to Atom 3.
         */
        double v13[] = new double[3];
        /**
         * Vector from Atom 3 to Atom 0.
         */
        double v30[] = new double[3];
        /**
         * Vector from Atom 3 to Atom 2.
         */
        double v32[] = new double[3];
        /**
         * Vector v12 cross v13.
         */
        double p[] = new double[3];
        /**
         * Gradient on atoms 0, 1, 2 & 3.
         */
        double g0[] = new double[3];
        double g1[] = new double[3];
        double g2[] = new double[3];
        double g3[] = new double[3];
        /**
         * Work arrays.
         */
        double sv30[] = new double[3];
        double sv32[] = new double[3];
        double dcda[] = new double[3];
        double dcdc[] = new double[3];
        double dcdd[] = new double[3];
        double deda[] = new double[3];
        double dedc[] = new double[3];
        double dedd[] = new double[3];

        energy = 0.0;
        value = 0.0;
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);
        diff(a0, a1, v10);
        diff(a2, a1, v12);
        diff(a3, a1, v13);
        diff(a0, a3, v30);
        diff(a2, a3, v32);
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
                    * outOfPlaneBendType.forceConstant * dv2
                    * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4)
                    * esvLambda;
            if (gradient) {
                double deddt = units
                        * outOfPlaneBendType.forceConstant * dv
                        * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic
                            * dv2 + 5.0 * quintic * dv3 + 6.0 * sextic * dv4)
                        * esvLambda;
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
                // atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                // atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                // atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                // atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
                int i0 = atoms[0].getIndex() - 1;
                gradX.add(threadID, i0, g0[0]);
                gradY.add(threadID, i0, g0[1]);
                gradZ.add(threadID, i0, g0[2]);
                int i1 = atoms[1].getIndex() - 1;
                gradX.add(threadID, i1, g1[0]);
                gradY.add(threadID, i1, g1[1]);
                gradZ.add(threadID, i1, g1[2]);
                int i2 = atoms[2].getIndex() - 1;
                gradX.add(threadID, i2, g2[0]);
                gradY.add(threadID, i2, g2[1]);
                gradZ.add(threadID, i2, g2[2]);
                int i3 = atoms[3].getIndex() - 1;
                gradX.add(threadID, i3, g3[0]);
                gradY.add(threadID, i3, g3[1]);
                gradZ.add(threadID, i3, g3[2]);
            }
        }
        if (esvTerm) {
            final double esvLambdaInv = (esvLambda != 0.0) ? 1/esvLambda : 1.0;
            setEsvDeriv(energy * dedesvChain * esvLambdaInv);
        }
        return energy;
    }

    /**
     * Log details for this Out-of-Plane Bend energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6.4f %10.4f",
                "Out-of-Plane Bend", atoms[1].getIndex(), atoms[1].getAtomType().name, atoms[3].getIndex(), atoms[3].getAtomType().name, value, energy));
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
    public int compareTo(BondedTerm o) {
        if (o == null) {
            throw new NullPointerException();
        }
        if (o == this) {
            return 0;
        }
        if (!o.getClass().isInstance(this)) {
            return super.compareTo(o);
        }
        int this1 = atoms[1].getIndex();
        int a1 = o.atoms[1].getIndex();
        if (this1 < a1) {
            return -1;
        }
        if (this1 > a1) {
            return 1;
        }
        int this3 = atoms[3].getIndex();
        int a3 = o.atoms[3].getIndex();
        if (this3 < a3) {
            return -1;
        }
        if (this3 > a3) {
            return 1;
        }
        return 0;
    }

}
