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
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.AtomicDoubleArray;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.StretchBendType;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.StretchBendType.units;

/**
 * The StretchBend class represents a Stretch-Bend formed between three linearly
 * bonded atoms.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class StretchBend extends BondedTerm {

    private static final Logger logger = Logger.getLogger(StretchBend.class.getName());
    private static final long serialVersionUID = 1L;
    /**
     * Force field parameters to compute the Stretch-Bend energy.
     */
    private StretchBendType stretchBendType = null;
    private double rigidScale = 1.0;
    /**
     * Angle this Stretch-Bend is based on.
     */
    protected Angle angle = null;
    private double force0, force1;
    final private double angleEq;
    final private double bond0Eq;
    final private double bond1Eq;

    /**
     * <p>
     * Setter for the field <code>stretchBendType</code>.</p>
     *
     * @param stretchBendType a {@link ffx.potential.parameters.StretchBendType}
     * object.
     */
    public void setStretchBendType(StretchBendType stretchBendType) {
        this.stretchBendType = stretchBendType;
        /**
         * Match the atom class of the angle to the atom class of the
         * stretch-bend type.
         */
        if (atoms[0].getAtomType().atomClass == stretchBendType.atomClasses[0]) {
            force0 = units * stretchBendType.forceConstants[0];
            force1 = units * stretchBendType.forceConstants[1];
        } else {
            force0 = units * stretchBendType.forceConstants[1];
            force1 = units * stretchBendType.forceConstants[0];
        }
    }

    /**
     * Attempt to create a new StretchBend if a StretchBendType exists for the
     * specified Angle.
     *
     * @param angle the Angle to created the StrechBend around.
     * @param forceField the ForceField parameters to use.
     * @return a new StretchBend, or null.
     */
    public static StretchBend stretchBendFactory(Angle angle, ForceField forceField) {
        StretchBendType stretchBendType = forceField.getStretchBendType(angle.getAngleType().getKey());
        if (stretchBendType == null) {
            return null;
        }
        StretchBend stretchBend = new StretchBend(angle);
        stretchBend.setStretchBendType(stretchBendType);
        return stretchBend;
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
     * Constructor for the Stretch-Bend class.
     *
     * @param a The Angle that this stretch-bend is based on.
     */
    public StretchBend(Angle a) {
        super();
        angle = a;
        atoms = a.atoms;
        bonds = a.bonds;
        angleEq = angle.angleType.angle[angle.nh];
        bond0Eq = bonds[0].bondType.distance;
        bond1Eq = bonds[1].bondType.distance;
        setID_Key(false);
    }

    /**
     * Evaluate the Stretch-Bend energy.
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
        /**
         * Work vectors.
         */
        double v10[] = new double[3];
        double v12[] = new double[3];
        double p[] = new double[3];
        double dta[] = new double[3];
        double dtc[] = new double[3];
        /**
         * Gradient on atoms 0, 1 & 2.
         */
        double g0[] = new double[3];
        double g1[] = new double[3];
        double g2[] = new double[3];

        energy = 0.0;
        value = 0.0;
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        diff(a0, a1, v10);
        diff(a2, a1, v12);
        double rab2 = dot(v10, v10);
        double rcb2 = dot(v12, v12);
        if (rab2 != 0.0 && rcb2 != 0.0) {
            double rab = r(v10);
            double rcb = r(v12);
            cross(v12, v10, p);
            double rp = r(p);
            rp = max(rp, 0.001);
            double cosine = dot(v10, v12) / (rab * rcb);
            cosine = min(1.0, max(-1.0, cosine));
            value = toDegrees(acos(cosine));
            double e0 = rab - bond0Eq;
            double e1 = rcb - bond1Eq;
            double dt = value - angleEq;
            double dr = force0 * e0 + force1 * e1;
            double prefactor = rigidScale * esvLambda;
            energy = prefactor * dr * dt;
            if (gradient) {
                // angle chain rule terms
                double term1 = -prefactor * dr * toDegrees(1.0 / (rab2 * rp));
                double term2 = prefactor * dr * toDegrees(1.0 / (rcb2 * rp));
                cross(v10, p, dta);
                scalar(dta, term1, dta);
                cross(v12, p, dtc);
                scalar(dtc, term2, dtc);

                // bond chain rule terms
                term1 = prefactor * force0 * dt / rab;
                term2 = prefactor * force1 * dt / rcb;
                scalar(v10, term1, v10);
                scalar(v12, term2, v12);

                sum(dta, v10, g0);
                sum(dtc, v12, g2);
                sum(g0, g2, g1);
                scalar(g1, -1.0, g1);

                // atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                // atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                // atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                int i0 = atoms[0].getXYZIndex() - 1;
                gradX.add(threadID, i0, g0[0]);
                gradY.add(threadID, i0, g0[1]);
                gradZ.add(threadID, i0, g0[2]);
                int i1 = atoms[1].getXYZIndex() - 1;
                gradX.add(threadID, i1, g1[0]);
                gradY.add(threadID, i1, g1[1]);
                gradZ.add(threadID, i1, g1[2]);
                int i2 = atoms[2].getXYZIndex() - 1;
                gradX.add(threadID, i2, g2[0]);
                gradY.add(threadID, i2, g2[1]);
                gradZ.add(threadID, i2, g2[2]);
            }
        }
        if (esvTerm) {
            setEsvDeriv(energy * dedesvChain / esvLambda);
        }
        return energy;
    }

    /**
     * <p>
     * log</p>
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s"
                + "%7.4f %10.4f",
                "Stretch-Bend", atoms[0].getXYZIndex(),
                atoms[0].getAtomType().name, atoms[1].getXYZIndex(),
                atoms[1].getAtomType().name, atoms[2].getXYZIndex(),
                atoms[2].getAtomType().name, value, energy));
    }

    /**
     * {@inheritDoc}
     *
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.2f,%7.2f,%7.1f,%7.2f)", id,
                bonds[0].value, bonds[1].value, angle.value, energy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compareTo(BondedTerm sb) {
        if (!sb.getClass().isInstance(this)) {
            return super.compareTo(sb);
        }
        return angle.compareTo(((StretchBend) sb).angle);
    }
}
