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

import java.util.logging.Logger;

import ffx.numerics.AtomicDoubleArray;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.UreyBradleyType;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.potential.parameters.UreyBradleyType.cubic;
import static ffx.potential.parameters.UreyBradleyType.quartic;
import static ffx.potential.parameters.UreyBradleyType.units;

/**
 * The UreyBradley class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class UreyBradley extends BondedTerm implements Comparable<UreyBradley> {

    private static final Logger logger = Logger.getLogger(UreyBradley.class.getName());

    private static final long serialVersionUID = 1L;
    /**
     * Force field parameters to compute the Stretch-Bend energy.
     */
    public UreyBradleyType ureyBradleyType = null;
    public double rigidScale = 1.0;
    /**
     * The Angle this UreyBradley term is based on.
     */
    protected Angle angle = null;

    /**
     * <p>
     * Setter for the field <code>ureyBradleyType</code>.</p>
     *
     * @param a a {@link ffx.potential.parameters.UreyBradleyType} object.
     */
    public void setUreyBradleyType(UreyBradleyType a) {
        ureyBradleyType = a;
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
     * Constructor for the UreyBradley class.
     *
     * @param a a {@link ffx.potential.bonded.Angle} object.
     */
    public UreyBradley(Angle a) {
        super();
        angle = a;
        bonds = a.bonds;
        atoms = a.atoms;
        setID_Key(false);
    }

    /**
     * Attempt to create a new UreyBradley for the specified Angle.
     *
     * @param angle the Angle to create the UreyBradley from.
     * @param forceField the ForceField parameters to apply.
     * @return a new UreyBradley, or null.
     */
    public static UreyBradley ureyBradlyFactory(Angle angle, ForceField forceField) {
        if (angle == null) {
            return null;
        }
        UreyBradleyType ureyBradleyType = forceField.getUreyBradleyType(angle.angleType.getKey());
        if (ureyBradleyType == null) {
            return null;
        }
        UreyBradley newUreyBradley = new UreyBradley(angle);
        newUreyBradley.ureyBradleyType = ureyBradleyType;
        return newUreyBradley;
    }

    /**
     * Evaluate the Urey-Bradley energy.
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
        double a2[] = new double[3];
        /**
         * The vector from Atom 2 to Atom 0.
         */
        double v20[] = new double[3];
        /**
         * Gradient on Atoms 0 & 2.
         */
        double g0[] = new double[3];
        double g2[] = new double[3];

        atoms[0].getXYZ(a0);
        atoms[2].getXYZ(a2);
        diff(a0, a2, v20);
        value = r(v20);
        double dv = value - ureyBradleyType.distance;
        double dv2 = dv * dv;
        energy = units * rigidScale * ureyBradleyType.forceConstant * dv2 * (1.0 + cubic * dv + quartic * dv2) * esvLambda;
        if (gradient) {
            double deddt = 2.0 * units * rigidScale * ureyBradleyType.forceConstant * dv * (1.0 + 1.5 * cubic * dv + 2.0 * quartic * dv2) * esvLambda;
            double de = 0.0;
            if (value > 0.0) {
                de = deddt / value;
            }
            scalar(v20, de, g0);
            scalar(v20, -de, g2);
            // atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
            // atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
            int i0 = atoms[0].getXYZIndex() - 1;
            gradX.add(threadID, i0, g0[0]);
            gradY.add(threadID, i0, g0[1]);
            gradZ.add(threadID, i0, g0[2]);
            int i2 = atoms[2].getXYZIndex() - 1;
            gradX.add(threadID, i2, g2[0]);
            gradY.add(threadID, i2, g2[1]);
            gradZ.add(threadID, i2, g2[2]);
        }
        if (esvTerm) {
            addToEsvDeriv(energy * dedesvChain / esvLambda, UreyBradley.class);
        }
        return energy;
    }

    /*
     * Log details for this Angle energy term.
     */
    /**
     * <p>
     * log</p>
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6.4f  %6.4f  %10.4f",
                "Urey-Bradley", atoms[0].getXYZIndex(), atoms[0].getAtomType().name,
                atoms[2].getXYZIndex(), atoms[2].getAtomType().name, ureyBradleyType.distance, value,
                energy));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compareTo(UreyBradley ub) {
        return angle.compareTo(ub.angle);
    }
}
