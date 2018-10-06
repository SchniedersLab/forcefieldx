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
package ffx.potential.bonded;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.PiTorsionType;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;
import static ffx.potential.parameters.PiTorsionType.units;

/**
 * The Pi-Orbital Torsion class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PiOrbitalTorsion extends BondedTerm implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(PiOrbitalTorsion.class.getName());
    private static final long serialVersionUID = 1L;
    public PiTorsionType piTorsionType = null;
    private double lambda = 1.0;
    private double dEdL = 0.0;
    private boolean lambdaTerm = false;

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
     * @param bond       the Bond to create a PiOrbitalTorsion around.
     * @param forceField the ForceField parameters to use.
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
     * <p>
     * Evaluate the Pi-Orbital Torsion energy.
     */
    @Override
    public double energy(boolean gradient,
                         int threadID,
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
        double a4[] = new double[3];
        double a5[] = new double[3];
        /**
         * Vector from Atom 3 to Atom 0.
         */
        double v30[] = new double[3];
        /**
         * Vector from Atom 3 to Atom 1.
         */
        double v31[] = new double[3];
        /**
         * Vector from Atom 2 to Atom 3.
         */
        double v23[] = new double[3];
        /**
         * Vector from Atom 2 to Atom 4.
         */
        double v24[] = new double[3];
        /**
         * Vector from Atom 2 to Atom 5.
         */
        double v25[] = new double[3];
        /**
         * Vector from Atom 5 to Atom 2.
         */
        double v[] = new double[3];
        /**
         * Vector v30 cross v13.
         */
        double vp[] = new double[3];
        /**
         * Vector v25 cross v52.
         */
        double vq[] = new double[3];
        /**
         * Work vectors.
         */
        double xt[] = new double[3];
        double xu[] = new double[3];
        double xtu[] = new double[3];
        double vpc[] = new double[3];
        double vdq[] = new double[3];
        double vpd[] = new double[3];
        double vcq[] = new double[3];
        double dedt[] = new double[3];
        double dedu[] = new double[3];
        double temp[] = new double[3];
        double dedp[] = new double[3];
        double dedc[] = new double[3];
        double dedd[] = new double[3];
        double dedq[] = new double[3];
        /**
         * Gradient on atoms 0, 1, 2, 3, & 5.
         */
        double g0[] = new double[3];
        double g1[] = new double[3];
        double g2[] = new double[3];
        double g3[] = new double[3];
        double g4[] = new double[3];
        double g5[] = new double[3];

        energy = 0.0;
        value = 0.0;
        dEdL = 0.0;

        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);
        atoms[4].getXYZ(a4);
        atoms[5].getXYZ(a5);
        diff(a0, a3, v30);
        diff(a1, a3, v31);
        diff(a3, a2, v23);
        diff(a4, a2, v24);
        diff(a5, a2, v25);
        cross(v30, v31, vp);
        cross(v24, v25, vq);
        sum(vp, a2, vp);
        sum(vq, a3, vq);
        diff(a2, vp, vpc);
        diff(vq, a3, vdq);
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
            energy = units * piTorsionType.forceConstant * phi2 * esvLambda;
            if (esvTerm) {
                setEsvDeriv(units * piTorsionType.forceConstant * phi2 * dedesvChain);
            }
            dEdL = energy;
            energy = lambda * energy;
            if (gradient || lambdaTerm) {
                double sine2 = 2.0 * cosine * sine;
                double dphi2 = 2.0 * sine2;
                double dedphi = units * piTorsionType.forceConstant * dphi2 * esvLambda;
                diff(a3, vp, vpd);
                diff(vq, a2, vcq);
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
                if (lambdaTerm) {
                    // atoms[0].addToLambdaXYZGradient(g0[0], g0[1], g0[2]);
                    // atoms[1].addToLambdaXYZGradient(g1[0], g1[1], g1[2]);
                    // atoms[2].addToLambdaXYZGradient(g2[0], g2[1], g2[2]);
                    // atoms[3].addToLambdaXYZGradient(g3[0], g3[1], g3[2]);
                    // atoms[4].addToLambdaXYZGradient(g4[0], g4[1], g4[2]);
                    // atoms[5].addToLambdaXYZGradient(g5[0], g5[1], g5[2]);
                    int i0 = atoms[0].getIndex() - 1;
                    lambdaGradX.add(threadID, i0, g0[0]);
                    lambdaGradY.add(threadID, i0, g0[1]);
                    lambdaGradZ.add(threadID, i0, g0[2]);
                    int i1 = atoms[1].getIndex() - 1;
                    lambdaGradX.add(threadID, i1, g1[0]);
                    lambdaGradY.add(threadID, i1, g1[1]);
                    lambdaGradZ.add(threadID, i1, g1[2]);
                    int i2 = atoms[2].getIndex() - 1;
                    lambdaGradX.add(threadID, i2, g2[0]);
                    lambdaGradY.add(threadID, i2, g2[1]);
                    lambdaGradZ.add(threadID, i2, g2[2]);
                    int i3 = atoms[3].getIndex() - 1;
                    lambdaGradX.add(threadID, i3, g3[0]);
                    lambdaGradY.add(threadID, i3, g3[1]);
                    lambdaGradZ.add(threadID, i3, g3[2]);
                    int i4 = atoms[4].getIndex() - 1;
                    lambdaGradX.add(threadID, i4, g4[0]);
                    lambdaGradY.add(threadID, i4, g4[1]);
                    lambdaGradZ.add(threadID, i4, g4[2]);
                    int i5 = atoms[5].getIndex() - 1;
                    lambdaGradX.add(threadID, i5, g5[0]);
                    lambdaGradY.add(threadID, i5, g5[1]);
                    lambdaGradZ.add(threadID, i5, g5[2]);
                }
                if (gradient) {
                    // atoms[0].addToXYZGradient(lambda * g0[0], lambda * g0[1], lambda * g0[2]);
                    // atoms[1].addToXYZGradient(lambda * g1[0], lambda * g1[1], lambda * g1[2]);
                    // atoms[2].addToXYZGradient(lambda * g2[0], lambda * g2[1], lambda * g2[2]);
                    // atoms[3].addToXYZGradient(lambda * g3[0], lambda * g3[1], lambda * g3[2]);
                    // atoms[4].addToXYZGradient(lambda * g4[0], lambda * g4[1], lambda * g4[2]);
                    // atoms[5].addToXYZGradient(lambda * g5[0], lambda * g5[1], lambda * g5[2]);
                    scalar(g0, lambda, g0);
                    scalar(g1, lambda, g1);
                    scalar(g2, lambda, g2);
                    scalar(g3, lambda, g3);
                    scalar(g4, lambda, g4);
                    scalar(g5, lambda, g5);
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
                    int i4 = atoms[4].getIndex() - 1;
                    gradX.add(threadID, i4, g4[0]);
                    gradY.add(threadID, i4, g4[1]);
                    gradZ.add(threadID, i4, g4[2]);
                    int i5 = atoms[5].getIndex() - 1;
                    gradX.add(threadID, i5, g5[0]);
                    gradY.add(threadID, i5, g5[1]);
                    gradZ.add(threadID, i5, g5[2]);
                }
            }
        }
        return energy;
    }

    /**
     * Log details for this Pi-Orbital Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %10.4f %10.4f",
                "Pi-Orbital Torsion", atoms[2].getIndex(), atoms[2].getAtomType().name,
                atoms[3].getIndex(), atoms[3].getAtomType().name, value, energy));
    }

    /**
     * {@inheritDoc}
     * <p>
     * Over-ridden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (applyAllLambda()) {
            this.lambda = lambda;
            lambdaTerm = true;
        } else {
            this.lambda = 1.0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        if (lambdaTerm) {
            return dEdL;
        } else {
            return 0.0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        return;
    }
}
