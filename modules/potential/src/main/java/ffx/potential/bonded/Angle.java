//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.Constraint;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.math.Double3;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;
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
 */
public class Angle extends BondedTerm {

    private static final Logger logger = Logger.getLogger(Angle.class.getName());

    /**
     * Force field parameters to compute the angle bending energy.
     */
    public AngleType angleType;
    /**
     * Number of hydrogens on the central atom that are not part of this Angle.
     */
    public int nh = 0;
    /**
     * Scale factor to apply to angle bending.
     */
    private double rigidScale = 1.0;
    /**
     * Fourth atom for in-plane angles.
     */
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

    /**
     * {@inheritDoc}
     */
    @Override
    public int compareTo(BondedTerm a) {
        if (a == null) {
            throw new NullPointerException();
        }
        if (a == this) {
            return 0;
        }
        if (!a.getClass().isInstance(this)) {
            return super.compareTo(a);
        }
        int this1 = atoms[1].getIndex();
        int a1 = a.atoms[1].getIndex();
        if (this1 < a1) {
            return -1;
        }
        if (this1 > a1) {
            return 1;
        }

        int this0 = atoms[0].getIndex();
        int a0 = a.atoms[0].getIndex();
        if (this0 < a0) {
            return -1;
        }
        if (this0 > a0) {
            return 1;
        }
        int this2 = atoms[2].getIndex();
        int a2 = a.atoms[2].getIndex();

        return Integer.compare(this2, a2);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluate this Angle energy.
     */
    @Override
    public double energy(boolean gradient, int threadID,
                         AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
        var atomA = atoms[0];
        var atomB = atoms[1];
        var atomC = atoms[2];
        var ia = atomA.getIndex() - 1;
        var ib = atomB.getIndex() - 1;
        var ic = atomC.getIndex() - 1;
        var va = atomA.getXYZ();
        var vb = atomB.getXYZ();
        var vc = atomC.getXYZ();
        energy = 0.0;
        value = 0.0;
        var prefactor = units * rigidScale * angleType.forceConstant * esvLambda;
        switch (angleType.angleFunction) {
            case SEXTIC:
                switch (angleType.angleMode) {
                    case NORMAL:
                        var vab = va.sub(vb);
                        var vcb = vc.sub(vb);
                        var rab2 = vab.length2();
                        var rcb2 = vcb.length2();
                        if (rab2 != 0.0 && rcb2 != 0.0) {
                            var p = vcb.X(vab);
                            var cosine = min(1.0, max(-1.0, vab.dot(vcb) / sqrt(rab2 * rcb2)));
                            value = toDegrees(acos(cosine));
                            var dv = value - angleType.angle[nh];
                            var dv2 = dv * dv;
                            var dv3 = dv2 * dv;
                            var dv4 = dv2 * dv2;
                            energy = prefactor * dv2 * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4);
                            if (gradient) {
                                var deddt = prefactor * dv * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic * dv2
                                        + 5.0 * quintic * dv3 + 6.0 * sextic * dv4);
                                var rp = max(p.length(), 0.000001);
                                var terma = -deddt / (rab2 * rp);
                                var termc = deddt / (rcb2 * rp);
                                var ga = vab.X(p).scale(terma);
                                var gc = vcb.X(p).scale(termc);
                                grad.add(threadID, ia, ga);
                                grad.sub(threadID, ib, ga.add(gc));
                                grad.add(threadID, ic, gc);
                            }
                            value = dv;
                        }
                        break;
                    case IN_PLANE:
                        var vd = getAtom4XYZ();
                        int id = atom4.getIndex() - 1;
                        var vad = va.sub(vd);
                        var vbd = vb.sub(vd);
                        var vcd = vc.sub(vd);
                        var vp = vad.X(vcd);
                        var rp2 = vp.length2();
                        var delta = -vp.dot(vbd) / rp2;
                        var vip = vp.scale(delta).addI(vbd);
                        var vjp = vad.sub(vip);
                        var vkp = vcd.sub(vip);
                        var jp2 = vjp.length2();
                        var kp2 = vkp.length2();
                        if (jp2 != 0.0 && kp2 != 0.0) {
                            var cosine = min(1.0, max(-1.0, vjp.dot(vkp) / sqrt(jp2 * kp2)));
                            value = toDegrees(acos(cosine));
                            var dv = value - angleType.angle[nh];
                            var dv2 = dv * dv;
                            var dv3 = dv2 * dv;
                            var dv4 = dv2 * dv2;
                            energy = prefactor * dv2 * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4);
                            if (gradient) {
                                var deddt = prefactor * dv * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic * dv2
                                        + 5.0 * quintic * dv3 + 6.0 * sextic * dv4);
                                inPlaneGrad(threadID, grad, ia, ib, ic, id, vad, vbd, vcd,
                                        vp, rp2, vjp, jp2, vkp, kp2, delta, deddt);
                            }
                            value = dv;
                        }
                        break;
                }
                break;
            case HARMONIC:
            default:
                switch (angleType.angleMode) {
                    case NORMAL:
                        var vab = va.sub(vb);
                        var vcb = vc.sub(vb);
                        var rab2 = vab.length2();
                        var rcb2 = vcb.length2();
                        if (rab2 != 0.0 && rcb2 != 0.0) {
                            var p = vcb.X(vab);
                            var cosine = min(1.0, max(-1.0, vab.dot(vcb) / sqrt(rab2 * rcb2)));
                            value = toDegrees(acos(cosine));
                            var dv = value - angleType.angle[nh];
                            var dv2 = dv * dv;
                            energy = prefactor * dv2;
                            if (gradient) {
                                var deddt = prefactor * dv * toDegrees(2.0);
                                var rp = max(p.length(), 0.000001);
                                var terma = -deddt / (rab2 * rp);
                                var termc = deddt / (rcb2 * rp);
                                var ga = vab.X(p).scaleI(terma);
                                var gc = vcb.X(p).scaleI(termc);
                                grad.add(threadID, ia, ga);
                                grad.sub(threadID, ib, ga.add(gc));
                                grad.add(threadID, ic, gc);
                            }
                            value = dv;
                        }
                        break;
                    case IN_PLANE:
                        var vd = getAtom4XYZ();
                        var id = atom4.getIndex() - 1;
                        var vad = va.sub(vd);
                        var vbd = vb.sub(vd);
                        var vcd = vc.sub(vd);
                        var vp = vad.X(vcd);
                        var rp2 = vp.length2();
                        var delta = -vp.dot(vbd) / rp2;
                        var vip = vp.scale(delta).addI(vbd);
                        var vjp = vad.sub(vip);
                        var vkp = vcd.sub(vip);
                        var rjp2 = vjp.length2();
                        var rkp2 = vkp.length2();
                        if (rjp2 != 0.0 && rkp2 != 0.0) {
                            var cosine = min(1.0, max(-1.0, vjp.dot(vkp) / sqrt(rjp2 * rkp2)));
                            value = toDegrees(acos(cosine));
                            var dv = value - angleType.angle[nh];
                            var dv2 = dv * dv;
                            energy = prefactor * dv2;
                            if (gradient) {
                                var deddt = prefactor * dv * toDegrees(2.0);
                                inPlaneGrad(threadID, grad, ia, ib, ic, id, vad, vbd, vcd,
                                        vp, rp2, vjp, rjp2, vkp, rkp2, delta, deddt);
                            }
                            value = dv;
                        }
                        break;
                }
                break;
        }
        if (esvTerm) {
            final var esvLambdaInv = (esvLambda != 0.0) ? 1 / esvLambda : 1.0;
            setEsvDeriv(energy * dedesvChain * esvLambdaInv);
        }
        return energy;
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

    /**
     * <p>Getter for the field <code>angleMode</code>.</p>
     *
     * @return a {@link ffx.potential.parameters.AngleType.AngleMode} object.
     */
    public AngleType.AngleMode getAngleMode() {
        return angleType.angleMode;
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
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a a {@link ffx.potential.parameters.AngleType} object.
     */
    public void setAngleType(AngleType a) {
        angleType = a;

        // Count the number of hydrogens attached to the central atom, but that are not part of the angle.
        List<Bond> ba = atoms[1].getBonds();
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
     * <p>Getter for the field <code>atom4</code>.</p>
     *
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getAtom4() {
        return atom4;
    }

    /**
     * <p>getCentralAtom.</p>
     *
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getCentralAtom() {
        return atoms[1];
    }

    /**
     * Log details for this Angle energy term.
     */
    public void log() {
        switch (angleType.angleMode) {
            case NORMAL:
                logger.info(format(
                        " %-8s %6d-%s %6d-%s %6d-%s %7.4f  %7.4f  %10.4f", "Angle",
                        atoms[0].getIndex(), atoms[0].getAtomType().name,
                        atoms[1].getIndex(), atoms[1].getAtomType().name,
                        atoms[2].getIndex(), atoms[2].getAtomType().name,
                        angleType.angle[nh], value, energy));
                break;
            case IN_PLANE:
                logger.info(format(
                        " %-8s %6d-%s %6d-%s %6d-%s %7.4f  %7.4f  %10.4f", "Angle-IP",
                        atoms[0].getIndex(), atoms[0].getAtomType().name,
                        atoms[1].getIndex(), atoms[1].getAtomType().name,
                        atoms[2].getIndex(), atoms[2].getAtomType().name,
                        angleType.angle[nh], value, energy));
                break;
        }
    }

    /**
     * Log that no AngleType exists.
     *
     * @param a1  Atom 1.
     * @param ac  Atom 2.
     * @param a3  Atom 3.
     * @param key The class key.
     */
    public static void logNoAngleType(Atom a1, Atom ac, Atom a3, String key) {
        logger.severe(format("No AngleType for key: %s\n %s -> %s\n %s -> %s\n %s -> %s", key,
                a1.toString(), a1.getAtomType().toString(),
                ac.toString(), ac.getAtomType().toString(),
                a3.toString(), a3.getAtomType().toString()));
    }

    @Override
    public void setConstraint(Constraint c) {
        super.setConstraint(c);
        for (Bond b : bonds) {
            b.setConstraint(c);
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
     * {@inheritDoc}
     * <p>
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

    /**
     * <p>angleFactory.</p>
     *
     * @param b1         a {@link ffx.potential.bonded.Bond} object.
     * @param b2         a {@link ffx.potential.bonded.Bond} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @return a {@link ffx.potential.bonded.Angle} object.
     */
    static Angle angleFactory(Bond b1, Bond b2, ForceField forceField) {
        Angle newAngle = new Angle(b1, b2);
        Atom ac = b1.getCommonAtom(b2);
        Atom a1 = b1.get1_2(ac);
        Atom a3 = b2.get1_2(ac);
        int[] c = new int[3];
        c[0] = a1.getAtomType().atomClass;
        c[1] = ac.getAtomType().atomClass;
        c[2] = a3.getAtomType().atomClass;
        String key = AngleType.sortKey(c);
        AngleType angleType = forceField.getAngleType(key);
        if (angleType == null) {
            logNoAngleType(a1, ac, a3, key);
            return null;
        }
        newAngle.setAngleType(angleType);
        return newAngle;
    }

    private static void inPlaneGrad(int threadID, AtomicDoubleArray3D grad,
                                    int ia, int ib, int ic, int id, Double3 vad, Double3 vbd, Double3 vcd,
                                    Double3 vp, double rp2, Double3 vjp, double rjp2, Double3 vkp, double rkp2,
                                    double delta, double deddt) {
        // Chain rule terms for first derivative components.
        var lp = vkp.X(vjp);
        var lpr = max(lp.length(), 0.000001);
        var ded0 = vjp.X(lp).scaleI(-deddt / (rjp2 * lpr));
        var ded2 = vkp.X(lp).scaleI(deddt / (rkp2 * lpr));
        var dedp = ded0.add(ded2);
        var gb = dedp.scale(-1.0);
        var delta2 = 2.0 * delta;
        var pt2 = dedp.dot(vp) / rp2;
        var xd2 = vcd.X(gb).scaleI(delta);
        var xp2 = vp.X(vcd).scaleI(delta2);
        var x21 = vbd.X(vcd).addI(xp2).scaleI(pt2);
        var dpd0 = xd2.add(x21);
        xd2 = gb.X(vad).scaleI(delta);
        xp2 = vp.X(vad).scaleI(delta2);
        x21.addI(xp2).scaleI(pt2);
        var dpd2 = xd2.addI(x21);
        var ga = ded0.addI(dpd0);
        var gc = ded2.addI(dpd2);
        // Accumulate derivatives.
        grad.add(threadID, ia, ga);
        grad.add(threadID, ib, gb);
        grad.add(threadID, ic, gc);
        grad.sub(threadID, id, ga.addI(gb).addI(gc));
    }

    /**
     * Get the position of the 4th atom.
     *
     * @return Returns the position of the 4th atom, or throws an exception.
     */
    private Double3 getAtom4XYZ() {
        try {
            return atom4.getXYZ();
        } catch (Exception e) {
            logger.info(" Atom 4 not found for angle: " + toString());
            for (var atom : atoms) {
                logger.info(" Atom: " + atom.toString());
                logger.info(" Type: " + atom.getAtomType().toString());
            }
            throw e;
        }
    }

    /**
     * <p>
     * Setter for the field <code>InPlaneAtom</code>.</p>
     *
     * @param a4 a {@link ffx.potential.bonded.Atom} object.
     */
    void setInPlaneAtom(Atom a4) {
        if (angleType.angleMode != AngleType.AngleMode.IN_PLANE) {
            logger.severe(" Attempted to set fourth atom for a normal angle.");
        }
        atom4 = a4;
    }

    /**
     * Finds the common bond between <b>this</b> angle and another
     *
     * @param a An Angle that may have a common bond with <b>this</b> angle
     * @return The common Bond between this Angle and Angle a, or null if this
     * == a or no common bond exists
     */
    Bond getCommonBond(Angle a) {
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
    Bond getOtherBond(Bond b) {
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
    Atom getTrigonalAtom() {
        if (atoms[1].isTrigonal()) {
            for (Bond b : atoms[1].getBonds()) {
                if (b != bonds[0] && b != bonds[1]) {
                    return b.get1_2(atoms[1]);
                }
            }
        }
        return null;
    }
}
