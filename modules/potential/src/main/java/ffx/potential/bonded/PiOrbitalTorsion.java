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
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.PiTorsionType;
import static ffx.potential.parameters.PiTorsionType.units;

/**
 * The Pi-Orbital Torsion class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PiOrbitalTorsion extends BondedTerm implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(PiOrbitalTorsion.class.getName());

    /**
     * A reference to the Pi-Torsion type in use.
     */
    public PiTorsionType piTorsionType = null;
    /**
     * Current value of lambda.
     */
    private double lambda = 1.0;
    /**
     * Current value of dE/dL.
     */
    private double dEdL = 0.0;
    /**
     * Flag to indicate use of lambda dependence.
     */
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
     * {@inheritDoc}
     * <p>
     * Evaluate the Pi-Orbital Torsion energy.
     */
    @Override
    public double energy(boolean gradient, int threadID,
                         AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
        energy = 0.0;
        value = 0.0;
        dEdL = 0.0;
        var atomA = atoms[0];
        var atomB = atoms[1];
        var atomC = atoms[2];
        var atomD = atoms[3];
        var atomE = atoms[4];
        var atomF = atoms[5];
        var va = atomA.getXYZ();
        var vb = atomB.getXYZ();
        var vc = atomC.getXYZ();
        var vd = atomD.getXYZ();
        var ve = atomE.getXYZ();
        var vf = atomF.getXYZ();
        var vad = va.sub(vd);
        var vbd = vb.sub(vd);
        var vdc = vd.sub(vc);
        var vec = ve.sub(vc);
        var vfc = vf.sub(vc);
        var vp = vad.X(vbd).addI(vc);
        var vq = vec.X(vfc).addI(vd);
        var vpc = vc.sub(vp);
        var vdq = vq.sub(vd);
        var vt = vpc.X(vdc);
        var vu = vdc.X(vdq);
        var rt2 = vt.length2();
        var ru2 = vu.length2();
        var rtru2 = rt2 * ru2;
        if (rtru2 != 0.0) {
            var rr = sqrt(rtru2);
            var rdc = vdc.length();
            var sine = vdc.dot(vt.X(vu)) / (rdc * rr);
            var cosine = min(1.0, max(-1.0, vt.dot(vu) / rr));
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
            var cosine2 = cosine * cosine - sine * sine;
            var phi2 = 1.0 - cosine2;
            energy = units * piTorsionType.forceConstant * phi2 * esvLambda;
            if (esvTerm) {
                setEsvDeriv(units * piTorsionType.forceConstant * phi2 * dedesvChain);
            }
            dEdL = energy;
            energy = lambda * energy;
            if (gradient || lambdaTerm) {
                var sine2 = 2.0 * cosine * sine;
                var dphi2 = 2.0 * sine2;
                var dedphi = units * piTorsionType.forceConstant * dphi2 * esvLambda;

                // Chain rule terms for first derivative components.
                var vdp = vd.sub(vp);
                var vqc = vq.sub(vc);
                var vdt = vt.X(vdc).scaleI(dedphi / (rt2 * rdc));
                var vdu = vu.X(vdc).scaleI(-dedphi / (ru2 * rdc));
                var dedp = vdt.X(vdc);
                var dedc = vdp.X(vdt).addI(vdu.X(vdq));
                var dedd = vdt.X(vpc).addI(vqc.X(vdu));
                var dedq = vdu.X(vdc);

                // Atomic gradient.
                var ga = vbd.X(dedp);
                var gb = dedp.X(vad);
                var ge = vfc.X(dedq);
                var gf = dedq.X(vec);
                var gc = dedc.add(dedp).subI(ge).subI(gf);
                var gd = dedd.add(dedq).subI(ga).subI(gb);
                var iA = atomA.getIndex() - 1;
                var iB = atomB.getIndex() - 1;
                var iC = atomC.getIndex() - 1;
                var iD = atomD.getIndex() - 1;
                var iE = atomE.getIndex() - 1;
                var iF = atomF.getIndex() - 1;
                if (lambdaTerm) {
                    lambdaGrad.add(threadID, iA, ga);
                    lambdaGrad.add(threadID, iB, gb);
                    lambdaGrad.add(threadID, iC, gc);
                    lambdaGrad.add(threadID, iD, gd);
                    lambdaGrad.add(threadID, iE, ge);
                    lambdaGrad.add(threadID, iF, gf);
                }
                if (gradient) {
                    grad.add(threadID, iA, ga.scaleI(lambda));
                    grad.add(threadID, iB, gb.scaleI(lambda));
                    grad.add(threadID, iC, gc.scaleI(lambda));
                    grad.add(threadID, iD, gd.scaleI(lambda));
                    grad.add(threadID, iE, ge.scaleI(lambda));
                    grad.add(threadID, iF, gf.scaleI(lambda));
                }
            }
        }
        return energy;
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
    public double getd2EdL2() {
        return 0.0;
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
    public void getdEdXdL(double[] gradient) {
        // The dEdXdL contributions are zero.
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
        int[] c = new int[2];
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
        PiOrbitalTorsion piOrbitalTorsion = new PiOrbitalTorsion(bond);
        piOrbitalTorsion.piTorsionType = piTorsionType;
        return piOrbitalTorsion;
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
}