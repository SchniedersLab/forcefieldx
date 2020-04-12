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

import java.util.List;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionTorsionType;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.potential.parameters.TorsionTorsionType.units;

/**
 * The TorsionTorsion class represents two adjacent torsional angles formed by
 * five bonded atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TorsionTorsion extends BondedTerm implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(TorsionTorsion.class.getName());
    private static final double[][] wt = {
            {1.0, 0.0, -3.0, 2.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 9.0, -6.0,
                    2.0, 0.0, -6.0, 4.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, -9.0, 6.0,
                    -2.0, 0.0, 6.0, -4.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, -6.0, 0.0,
                    0.0, -6.0, 4.0},
            {0.0, 0.0, 3.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.0, 6.0,
                    0.0, 0.0, 6.0, -4.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -3.0, 2.0, -2.0, 0.0, 6.0, -4.0,
                    1.0, 0.0, -3.0, 2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 3.0, -2.0,
                    1.0, 0.0, -3.0, 2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 2.0, 0.0,
                    0.0, 3.0, -2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -2.0, 0.0, 0.0, -6.0, 4.0,
                    0.0, 0.0, 3.0, -2.0},
            {0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, -3.0,
                    0.0, 2.0, -4.0, 2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -6.0, 3.0, 0.0,
                    -2.0, 4.0, -2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 3.0, 0.0,
                    0.0, 2.0, -2.0},
            {0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -3.0,
                    0.0, 0.0, -2.0, 2.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, -2.0, 4.0, -2.0,
                    0.0, 1.0, -2.0, 1.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, -1.0,
                    0.0, 1.0, -2.0, 1.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0,
                    0.0, -1.0, 1.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 2.0, -2.0,
                    0.0, 0.0, -1.0, 1.0}};
    /**
     * The two torsions that are coupled.
     */
    public final Torsion[] torsions = new Torsion[2];
    /**
     * The force field Torsion-Torsion type in use.
     */
    public TorsionTorsionType torsionTorsionType = null;
    /**
     * Value of lambda.
     */
    private double lambda = 1.0;
    /**
     * Value of dE/dL.
     */
    private double dEdL = 0.0;
    /**
     * Flag to indicate lambda dependence.
     */
    private boolean lambdaTerm = false;

    /**
     * Torsion-Torsion constructor.
     *
     * @param firstBond a {@link ffx.potential.bonded.Bond} object.
     * @param angle     a {@link ffx.potential.bonded.Angle} object.
     * @param lastBond  a {@link ffx.potential.bonded.Bond} object.
     * @param reversed  a boolean.
     */
    public TorsionTorsion(Bond firstBond, Angle angle, Bond lastBond,
                          boolean reversed) {
        super();
        atoms = new Atom[5];
        bonds = new Bond[4];
        if (!reversed) {
            atoms[1] = angle.atoms[0];
            atoms[2] = angle.atoms[1];
            atoms[3] = angle.atoms[2];
            atoms[0] = firstBond.get1_2(atoms[1]);
            atoms[4] = lastBond.get1_2(atoms[3]);
            bonds[0] = firstBond;
            bonds[1] = angle.bonds[0];
            bonds[2] = angle.bonds[1];
            bonds[3] = lastBond;
        } else {
            atoms[1] = angle.atoms[2];
            atoms[2] = angle.atoms[1];
            atoms[3] = angle.atoms[0];
            atoms[0] = lastBond.get1_2(atoms[1]);
            atoms[4] = firstBond.get1_2(atoms[3]);
            bonds[0] = lastBond;
            bonds[1] = angle.bonds[1];
            bonds[2] = angle.bonds[0];
            bonds[3] = firstBond;
        }
        torsions[0] = atoms[0].getTorsion(atoms[1], atoms[2], atoms[3]);
        torsions[1] = atoms[4].getTorsion(atoms[3], atoms[2], atoms[1]);
        setID_Key(false);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluate the Torsion-Torsion energy.
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
        var va = atomA.getXYZ();
        var vb = atomB.getXYZ();
        var vc = atomC.getXYZ();
        var vd = atomD.getXYZ();
        var ve = atomE.getXYZ();
        var vba = vb.sub(va);
        var vcb = vc.sub(vb);
        var vdc = vd.sub(vc);
        var ved = ve.sub(vd);
        var vt = vba.X(vcb);
        var vu = vcb.X(vdc);
        var vv = vdc.X(ved);
        var rt2 = vt.length2();
        var ru2 = vu.length2();
        var rv2 = vv.length2();
        var rtru2 = rt2 * ru2;
        var rurv2 = ru2 * rv2;
        if (rtru2 != 0.0 && rurv2 != 0.0) {
            var rtru = sqrt(rt2 * ru2);
            var rurv = sqrt(ru2 * rv2);
            var rcb = vcb.length();
            var cosine1 = min(1.0, max(-1.0, vt.dot(vu) / rtru));
            var angle1 = toDegrees(acos(cosine1));
            var sign = vba.dot(vu);
            if (sign < 0.0) {
                angle1 *= -1.0;
            }
            var rdc = vdc.length();
            var cosine2 = min(1.0, max(-1.0, vu.dot(vv) / rurv));
            var angle2 = toDegrees(acos(cosine2));
            sign = vcb.dot(vv);
            if (sign < 0.0) {
                angle2 *= -1.0;
            }
            var t1 = angle1;
            var t2 = angle2;
            sign = chktor();
            t1 *= sign;
            t2 *= sign;

            // Use bicubic interpolation to compute the spline values.
            var nx = torsionTorsionType.nx;
            var ny = torsionTorsionType.ny;
            var nlow = 0;
            var nhigh = nx - 1;
            while (nhigh - nlow > 1) {
                var nt = (nhigh + nlow) / 2;
                if (torsionTorsionType.tx[nt] > t1) {
                    nhigh = nt;
                } else {
                    nlow = nt;
                }
            }
            var xlow = nlow;
            nlow = 0;
            nhigh = ny - 1;
            while (nhigh - nlow > 1) {
                var nt = (nhigh + nlow) / 2;
                if (torsionTorsionType.ty[nt] > t2) {
                    nhigh = nt;
                } else {
                    nlow = nt;
                }
            }
            var ylow = nlow;
            var x1l = torsionTorsionType.tx[xlow];
            var x1u = torsionTorsionType.tx[xlow + 1];
            var y1l = torsionTorsionType.ty[ylow];
            var y1u = torsionTorsionType.ty[ylow + 1];
            var pos2 = (ylow + 1) * nx + xlow;
            var pos1 = pos2 - nx;

            // Array of 4 spline energies surrounding the actual Torsion-Torsion location.
            var e = new double[4];
            e[0] = torsionTorsionType.energy[pos1];
            e[1] = torsionTorsionType.energy[pos1 + 1];
            e[2] = torsionTorsionType.energy[pos2 + 1];
            e[3] = torsionTorsionType.energy[pos2];
            // Array of 4 spline x-gradients surrounding the actual Torsion-Torsion location.
            var dx = new double[4];
            dx[0] = torsionTorsionType.dx[pos1];
            dx[1] = torsionTorsionType.dx[pos1 + 1];
            dx[2] = torsionTorsionType.dx[pos2 + 1];
            dx[3] = torsionTorsionType.dx[pos2];
            // Array of 4 spline y-gradients surrounding the actual Torsion-Torsion location.
            var dy = new double[4];
            dy[0] = torsionTorsionType.dy[pos1];
            dy[1] = torsionTorsionType.dy[pos1 + 1];
            dy[2] = torsionTorsionType.dy[pos2 + 1];
            dy[3] = torsionTorsionType.dy[pos2];
            // Array of 4 spline xy-gradients surrounding the actual Torsion-Torsion location.
            var dxy = new double[4];
            dxy[0] = torsionTorsionType.dxy[pos1];
            dxy[1] = torsionTorsionType.dxy[pos1 + 1];
            dxy[2] = torsionTorsionType.dxy[pos2 + 1];
            dxy[3] = torsionTorsionType.dxy[pos2];
            if (!gradient && !lambdaTerm) {
                var bcu = bcuint(x1l, x1u, y1l, y1u, t1, t2, e, dx, dy, dxy);
                energy = units * bcu * esvLambda * lambda;
                if (esvTerm) {
                    setEsvDeriv(units * bcu * dedesvChain * lambda);
                }
                dEdL = units * bcu * esvLambda;
            } else {
                var ansy = new double[2];
                var bcu1 = bcuint1(x1l, x1u, y1l, y1u, t1, t2, e, dx, dy, dxy, ansy);
                energy = units * bcu1 * esvLambda * lambda;
                if (esvTerm) {
                    setEsvDeriv(units * bcu1 * dedesvChain * lambda);
                }
                dEdL = units * bcu1 * esvLambda;
                var dedang1 = sign * units * toDegrees(ansy[0]) * esvLambda * lambda;
                var dedang2 = sign * units * toDegrees(ansy[1]) * esvLambda * lambda;
                // Derivative components for the first angle.
                var vca = vc.sub(va);
                var vdb = vd.sub(vb);
                var vdt = vt.X(vcb).scaleI(dedang1 / (rt2 * rcb));
                var vdu = vu.X(vcb).scaleI(-dedang1 / (ru2 * rcb));
                // Gradients on each atom.
                var ga = vdt.X(vcb);
                var gb = vca.X(vdt).addI(vdu.X(vdc));
                var gc = vdb.X(vdu).addI(vdt.X(vba));
                var gd = vdu.X(vcb);
                var ia = atomA.getIndex() - 1;
                var ib = atomB.getIndex() - 1;
                var ic = atomC.getIndex() - 1;
                var id = atomD.getIndex() - 1;
                var ie = atomE.getIndex() - 1;
                if (lambdaTerm) {
                    lambdaGrad.add(threadID, ia, ga);
                    lambdaGrad.add(threadID, ib, gb);
                    lambdaGrad.add(threadID, ic, gc);
                    lambdaGrad.add(threadID, id, gd);
                }
                if (gradient) {
                    grad.add(threadID, ia, ga.scaleI(lambda));
                    grad.add(threadID, ib, gb.scaleI(lambda));
                    grad.add(threadID, ic, gc.scaleI(lambda));
                    grad.add(threadID, id, gd.scaleI(lambda));
                }
                // Derivative components for the 2nd angle.
                var vec = ve.sub(vc);
                vdt = vu.X(vdc).scaleI(dedang2 / (ru2 * rdc));
                vdu = vv.X(vdc).scaleI(-dedang2 / (rv2 * rdc));
                // Gradients on each atom.
                gb = vdt.X(vdc);
                gc = vdb.X(vdt).addI(vdu.X(ved));
                gd = vdt.X(vcb).addI(vec.X(vdu));
                var ge = vdu.X(vdc);
                if (lambdaTerm) {
                    lambdaGrad.add(threadID, ib, gb);
                    lambdaGrad.add(threadID, ic, gc);
                    lambdaGrad.add(threadID, id, gd);
                    lambdaGrad.add(threadID, ie, ge);
                }
                if (gradient) {
                    grad.add(threadID, ib, gb.scaleI(lambda));
                    grad.add(threadID, ic, gc.scaleI(lambda));
                    grad.add(threadID, id, gd.scaleI(lambda));
                    grad.add(threadID, ie, ge.scaleI(lambda));
                }
            }
        }
        return energy;
    }

    /**
     * <p>getChiralAtom.</p>
     *
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getChiralAtom() {
        Atom atom = null;
        List<Bond> bnds = atoms[2].getBonds();

        // To be chiral, the central atom must have 4 bonds.
        if (bnds.size() == 4) {
            // Find the two atoms that are not part of the dihedral.
            Atom atom1 = null;
            Atom atom2 = null;
            for (Bond b : bnds) {
                Atom a = b.get1_2(atoms[2]);
                if (a != atoms[1] && a != atoms[3]) {
                    if (atom1 == null) {
                        atom1 = a;
                    } else {
                        atom2 = a;
                    }
                }
            }
            /*
              Choose atom1 or atom2 to use for the chiral check, depending on
              their atom types and atomic number.
             */
            if (atom1.getType() > atom2.getType()) {
                atom = atom1;
            }
            if (atom2.getType() > atom1.getType()) {
                atom = atom2;
            }
            if (atom1.getAtomicNumber() > atom2.getAtomicNumber()) {
                atom = atom1;
            }
            if (atom2.getAtomicNumber() > atom1.getAtomicNumber()) {
                atom = atom2;
            }
        }
        return atom;
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
            lambdaTerm = true;
            this.lambda = lambda;
        } else {
            lambdaTerm = false;
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
        // This chain rule term is zero.
    }

    /**
     * Log details for this Torsion-Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f",
                "Torsional-Torsion", atoms[0].getIndex(), atoms[0].getAtomType().name, atoms[1].getIndex(),
                atoms[1].getAtomType().name, atoms[2].getIndex(), atoms[2].getAtomType().name,
                atoms[3].getIndex(), atoms[3].getAtomType().name, energy));
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.2f,%7.2f,%7.2f)", id, torsions[0].value,
                torsions[1].value, energy);
    }

    /**
     * <p>torsionTorsionFactory.</p>
     *
     * @param firstBond  the first Bond.
     * @param angle      the Angle.
     * @param lastBond   the last Bond.
     * @param forceField the ForceField parameters to apply.
     * @return the new TorsionTorsion, or null.
     */
    public static TorsionTorsion torsionTorsionFactory(Bond firstBond,
                                                       Angle angle, Bond lastBond, ForceField forceField) {
        int[] c5 = new int[5];
        Atom atom1 = angle.atoms[0];
        Atom atom3 = angle.atoms[2];
        c5[0] = firstBond.get1_2(atom1).getAtomType().atomClass;
        c5[1] = atom1.getAtomType().atomClass;
        c5[2] = angle.atoms[1].getAtomType().atomClass;
        c5[3] = atom3.getAtomType().atomClass;
        c5[4] = lastBond.get1_2(atom3).getAtomType().atomClass;
        String key = TorsionTorsionType.sortKey(c5);
        boolean reversed = false;
        TorsionTorsionType torsionTorsionType = forceField.getTorsionTorsionType(key);
        if (torsionTorsionType == null) {
            key = TorsionTorsionType.reverseKey(c5);
            torsionTorsionType = forceField.getTorsionTorsionType(key);
            reversed = true;
        }
        if (torsionTorsionType == null) {
            return null;
        }
        TorsionTorsion torsionTorsion = new TorsionTorsion(
                firstBond, angle, lastBond, reversed);
        torsionTorsion.torsionTorsionType = torsionTorsionType;
        return torsionTorsion;
    }

    private static void bcucof(double t1, double t2,
                               double[] e, double[] dx, double[] dy, double[] dxy,
                               double[][] c) {

        var x16 = new double[16];
        var cl = new double[16];
        var t1t2 = t1 * t2;

        // Pack a temporary vector of corner values.
        for (int i = 0; i < 4; i++) {
            x16[i] = e[i];
            x16[i + 4] = dx[i] * t1;
            x16[i + 8] = dy[i] * t2;
            x16[i + 12] = dxy[i] * t1t2;
        }

        // Matrix multiply by the stored weight table.
        for (int i = 0; i < 16; i++) {
            double xx = 0.0;
            for (int k = 0; k < 16; k++) {
                xx += wt[k][i] * x16[k];
            }
            cl[i] = xx;
        }

        // Unpack the results into the coefficient table.
        int j = 0;
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                c[i][k] = cl[j++];
            }
        }
    }

    /**
     * Check for inversion of the central atom if it is chiral.
     *
     * @return The sign convention - if negative the torsion angle signs are
     * inverted.
     */
    private double chktor() {

        // Vector from the central atom to site 0.
        var vc0 = new double[3];
        // Vector from the central atom to site 1.
        var vc1 = new double[3];
        // Vector from the central atom to site 2.
        var vc2 = new double[3];

        List<Bond> bonds = atoms[2].getBonds();
        // To be chiral, the central atom must have 4 bonds.
        if (bonds.size() == 4) {
            // Find the two atoms that are not part of the dihedral.
            Atom atom1 = null;
            Atom atom2 = null;
            for (Bond b : bonds) {
                Atom a = b.get1_2(atoms[2]);
                if (a != atoms[1] && a != atoms[3]) {
                    if (atom1 == null) {
                        atom1 = a;
                    } else {
                        atom2 = a;
                    }
                }
            }
            /*
              Choose atom1 or atom2 to use for the chiral check, depending on
              their atom types and atomic number.
             */
            Atom atom = null;
            if (atom1.getType() > atom2.getType()) {
                atom = atom1;
            }
            if (atom2.getType() > atom1.getType()) {
                atom = atom2;
            }
            if (atom1.getAtomicNumber() > atom2.getAtomicNumber()) {
                atom = atom1;
            }
            if (atom2.getAtomicNumber() > atom1.getAtomicNumber()) {
                atom = atom2;
            }

            // Compute the signed parallelpiped volume at the central site.
            if (atom != null) {
                var ad = new double[3];
                var a1 = new double[3];
                var a2 = new double[3];
                var a3 = new double[3];
                atom.getXYZ(ad);
                atoms[1].getXYZ(a1);
                atoms[2].getXYZ(a2);
                atoms[3].getXYZ(a3);
                sub(ad, a2, vc0);
                sub(a1, a2, vc1);
                sub(a3, a2, vc2);
                double volume = vc0[0] * (vc1[1] * vc2[2] - vc1[2] * vc2[1])
                        + vc1[0] * (vc2[1] * vc0[2] - vc2[2] * vc0[1])
                        + vc2[0] * (vc0[1] * vc1[2] - vc0[2] * vc1[1]);
                if (volume < 0.0) {
                    return -1.0;
                }
            }
        }
        return 1.0;
    }

    /**
     * <p>
     * bcuint</p>
     *
     * @param x1l a double.
     * @param x1u a double.
     * @param y1l a double.
     * @param y1u a double.
     * @param t1  a double.
     * @param t2  a double.
     * @param e   an array of {@link double} objects.
     * @param dx  an array of {@link double} objects.
     * @param dy  an array of {@link double} objects.
     * @param dxy an array of {@link double} objects.
     * @return a double.
     */
    private double bcuint(double x1l, double x1u, double y1l, double y1u,
                          double t1, double t2,
                          double[] e, double[] dx, double[] dy, double[] dxy) {

        var c = new double[4][4];
        var deltax = x1u - x1l;
        var deltay = y1u - y1l;
        bcucof(deltax, deltay, e, dx, dy, dxy, c);
        var tx = (t1 - x1l) / deltax;
        var ux = (t2 - y1l) / deltay;
        var ret = 0.0;
        for (int i = 3; i >= 0; i--) {
            ret = tx * ret + ((c[i][3] * ux + c[i][2]) * ux + c[i][1]) * ux + c[i][0];
        }
        return ret;
    }

    /**
     * <p>
     * bcuint1</p>
     *
     * @param x1l  a double.
     * @param x1u  a double.
     * @param y1l  a double.
     * @param y1u  a double.
     * @param t1   a double.
     * @param t2   a double.
     * @param e    an array of {@link double} objects.
     * @param dx   an array of {@link double} objects.
     * @param dy   an array of {@link double} objects.
     * @param dxy  an array of {@link double} objects.
     * @param ansy an array of {@link double} objects.
     * @return a double.
     */
    private double bcuint1(double x1l, double x1u, double y1l, double y1u,
                           double t1, double t2,
                           double[] e, double[] dx, double[] dy, double[] dxy,
                           double[] ansy) {

        var c = new double[4][4];
        var deltax = x1u - x1l;
        var deltay = y1u - y1l;
        bcucof(deltax, deltay, e, dx, dy, dxy, c);
        var tx = (t1 - x1l) / deltax;
        var ux = (t2 - y1l) / deltay;
        var ret = 0.0;
        ansy[0] = 0.0;
        ansy[1] = 0.0;
        for (int i = 3; i >= 0; i--) {
            ret = tx * ret + ((c[i][3] * ux + c[i][2]) * ux + c[i][1]) * ux + c[i][0];
            ansy[0] = ux * ansy[0] + (3.0 * c[3][i] * tx + 2.0 * c[2][i]) * tx + c[1][i];
            ansy[1] = tx * ansy[1] + (3.0 * c[i][3] * ux + 2.0 * c[i][2]) * ux + c[i][1];
        }
        ansy[0] /= deltax;
        ansy[1] /= deltay;
        return ret;
    }


}
