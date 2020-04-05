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

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.signum;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.OutOfPlaneBendType;
import static ffx.potential.parameters.OutOfPlaneBendType.cubic;
import static ffx.potential.parameters.OutOfPlaneBendType.quartic;
import static ffx.potential.parameters.OutOfPlaneBendType.quintic;
import static ffx.potential.parameters.OutOfPlaneBendType.sextic;
import static ffx.potential.parameters.OutOfPlaneBendType.units;

/**
 * The OutOfPlaneBend class represents an Out-Of-Plane Bend.
 *
 * @author Michael J. Schnieders
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
     * @param atom  The 4th atom of the trigonal center.
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
        return Integer.compare(this3, a3);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluate this Out-of-Plane Bend energy.
     */
    @Override
    public double energy(boolean gradient, int threadID,
                         AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
        energy = 0.0;
        value = 0.0;
        var atomA = atoms[0];
        var atomB = atoms[1];
        var atomC = atoms[2];
        var atomD = atoms[3];
        var va = atomA.getXYZ();
        var vb = atomB.getXYZ();
        var vc = atomC.getXYZ();
        var vd = atomD.getXYZ();
        var vab = va.sub(vb);
        var vcb = vc.sub(vb);
        var vdb = vd.sub(vb);
        var vad = va.sub(vd);
        var vcd = vc.sub(vd);
        var rdb2 = vdb.length2();
        var rad2 = vad.length2();
        var rcd2 = vcd.length2();
        var vp = vcb.X(vdb);
        var ee = vab.dot(vp);
        var rac2 = vad.dot(vcd);
        var cc = rad2 * rcd2 - rac2 * rac2;
        if (rdb2 != 0.0 && cc != 0.0) {
            var bkk2 = rdb2 - ee * ee / cc;
            var cosine = min(1.0, max(-1.0, sqrt(bkk2 / rdb2)));
            value = toDegrees(acos(cosine));
            var dv = value;
            var dv2 = dv * dv;
            var dv3 = dv2 * dv;
            var dv4 = dv2 * dv2;
            energy = units * outOfPlaneBendType.forceConstant * dv2
                    * (1.0 + cubic * dv + quartic * dv2 + quintic * dv3 + sextic * dv4) * esvLambda;
            if (gradient) {
                var deddt = units * outOfPlaneBendType.forceConstant * dv
                        * toDegrees(2.0 + 3.0 * cubic * dv + 4.0 * quartic
                        * dv2 + 5.0 * quintic * dv3 + 6.0 * sextic * dv4) * esvLambda;
                var dedcos = 0.0;
                if (ee != 0.0) {
                    dedcos = -deddt * signum(ee) / sqrt(cc * bkk2);
                } else {
                    dedcos = -deddt / sqrt(cc * bkk2);
                }
                var term = ee / cc;

                // Chain rule terms for first derivative components.
                var svad = vad.scale(rcd2);
                var svcd = vcd.scale(rac2);
                var dcda = svad.sub(svcd).scaleI(term);
                svad = vad.scale(rac2);
                svcd = vcd.scale(rad2);
                var dadc = svcd.sub(svad).scaleI(term);
                var dcdd = dcda.add(dadc).scaleI(-1.0);
                var deda = vdb.X(vcb);
                var dedc = vab.X(vdb);
                var dedd = vcb.X(vab);
                dedd.addI(vdb.scaleI(ee / rdb2));

                // Atomic gradient.
                var ga = dcda.add(deda).scaleI(dedcos);
                var gc = dadc.add(dedc).scaleI(dedcos);
                var gd = dcdd.add(dedd).scaleI(dedcos);
                var gb = ga.add(gc).addI(gd).scaleI(-1.0);
                var ia = atomA.getIndex() - 1;
                var ib = atomB.getIndex() - 1;
                var ic = atomC.getIndex() - 1;
                var id = atomD.getIndex() - 1;
                grad.add(threadID, ia, ga);
                grad.add(threadID, ib, gb);
                grad.add(threadID, ic, gc);
                grad.add(threadID, id, gd);

            }
        }
        if (esvTerm) {
            var esvLambdaInv = (esvLambda != 0.0) ? 1 / esvLambda : 1.0;
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
     * Attempt to create a new OutOfPlaneBend instance for a given Angle and
     * Force Field.
     *
     * @param angle      the Angle to create an OutOfPlaneBend around.
     * @param forceField the ForceField parameters to use.
     * @return a new OutOfPlaneBend if the central atom of the angle is trigonal
     * and a force field type exists.
     */
    public static OutOfPlaneBend outOfPlaneBendFactory(Angle angle, ForceField forceField) {
        Atom centralAtom = angle.atoms[1];
        if (centralAtom.isTrigonal()) {
            Atom atom4 = angle.getTrigonalAtom();
            int class4 = atom4.getAtomType().atomClass;
            int class0 = angle.atoms[0].getAtomType().atomClass;
            int class1 = angle.atoms[1].getAtomType().atomClass;
            int class2 = angle.atoms[2].getAtomType().atomClass;

            // First check for an atom4-center-edge-edge type (also checking reversed edges).
            String key = String.format("%d %d %d %d", class4, class1, class0, class2);
            OutOfPlaneBendType oopBendType = forceField.getOutOfPlaneBendType(key);
            if (oopBendType == null) {
                key = String.format("%d %d %d %d", class4, class1, class2, class0);
                oopBendType = forceField.getOutOfPlaneBendType(key);
            }

            // Then, check for a generic OOP bend type atom4-center-any-any
            if (oopBendType == null) {
                key = String.format("%d %d 0 0", class4, class1);
                oopBendType = forceField.getOutOfPlaneBendType(key);
            }
            if (oopBendType != null) {
                if (angle.getAngleMode() == AngleType.AngleMode.IN_PLANE) {
                    angle.setInPlaneAtom(atom4);
                }
                OutOfPlaneBend newOutOfPlaneBend = new OutOfPlaneBend(angle, atom4);
                newOutOfPlaneBend.outOfPlaneBendType = oopBendType;
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
     * {@inheritDoc}
     * <p>
     * Overridden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.1f,%7.2f)", id, value, energy);
    }

}
