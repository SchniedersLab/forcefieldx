/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.StretchTorsionType;
import ffx.potential.parameters.TorsionType;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;

/**
 * The StretchTorsion class represents a coupling between a torsional angle and the three
 * bonds contained in the torsion, as defined in the 2017 AMOEBA nucleic acid force field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class StretchTorsion extends BondedTerm implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(StretchTorsion.class.getName());

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
     * Stretch Torsion force field type.
     */
    private StretchTorsionType stretchTorsionType = null;
    /**
     * Stretch Torsion force constants (may be reversed compared to storage in the StretchTorsionType instance).
     */
    private double[] constants = new double[9];
    /**
     * Torsion force field type.
     */
    private TorsionType torsionType = null;
    /**
     * First bond force field type.
     */
    public BondType bondType1 = null;
    /**
     * Second bond force field type.
     */
    public BondType bondType2 = null;
    /**
     * Third bond force field type.
     */
    public BondType bondType3 = null;

    /**
     * Create a StretchTorsion from 3 connected bonds (no error checking)
     *
     * @param b1 Bond
     * @param b2 Bond
     * @param b3 Bond
     */
    private StretchTorsion(Bond b1, Bond b2, Bond b3) {
        super();
        bonds = new Bond[3];
        bonds[0] = b1;
        bonds[1] = b2;
        bonds[2] = b3;
        initialize();
    }

    /**
     * Torsion Constructor.
     *
     * @param n Torsion id
     */
    private StretchTorsion(String n) {
        super(n);
    }

    /**
     * Attempt to create a new StretchTorsion based on the supplied torsion.
     *
     * @param torsion    the Torsion.
     * @param forceField the ForceField parameters to apply.
     * @return a new StretchTorsion, or null.
     */
    public static StretchTorsion stretchTorsionFactory(Torsion torsion, ForceField forceField) {
        TorsionType torsionType = torsion.torsionType;
        String key = torsionType.getKey();
        StretchTorsionType stretchTorsionType = forceField.getStretchTorsionType(key);
        if (stretchTorsionType != null) {
            Bond bond1 = torsion.bonds[0];
            Bond middleBond = torsion.bonds[1];
            Bond bond3 = torsion.bonds[2];
            StretchTorsion stretchTorsion = new StretchTorsion(bond1, middleBond, bond3);
            stretchTorsion.stretchTorsionType = stretchTorsionType;
            stretchTorsion.torsionType = torsion.torsionType;
            stretchTorsion.bondType1 = bond1.bondType;
            stretchTorsion.bondType2 = middleBond.bondType;
            stretchTorsion.bondType3 = bond3.bondType;
            Atom atom1 = torsion.atoms[0];
            Atom atom2 = torsion.atoms[1];
            Atom atom3 = torsion.atoms[2];
            Atom atom4 = torsion.atoms[3];
            if (atom1.getAtomType().atomClass == stretchTorsionType.atomClasses[0] &&
                    atom2.getAtomType().atomClass == stretchTorsionType.atomClasses[1] &&
                    atom3.getAtomType().atomClass == stretchTorsionType.atomClasses[2] &&
                    atom4.getAtomType().atomClass == stretchTorsionType.atomClasses[3]) {
                stretchTorsion.setFlipped(false);
            } else {
                stretchTorsion.setFlipped(true);
            }

            return stretchTorsion;
        }
        return null;
    }

    /**
     * <p>
     * compare</p>
     *
     * @param a0 a {@link ffx.potential.bonded.Atom} object.
     * @param a1 a {@link ffx.potential.bonded.Atom} object.
     * @param a2 a {@link ffx.potential.bonded.Atom} object.
     * @param a3 a {@link ffx.potential.bonded.Atom} object.
     * @return a boolean.
     */
    public boolean compare(Atom a0, Atom a1, Atom a2, Atom a3) {
        if (a0 == atoms[0] && a1 == atoms[1] && a2 == atoms[2] && a3 == atoms[3]) {
            return true;
        }
        return (a0 == atoms[3] && a1 == atoms[2] && a2 == atoms[1] && a3 == atoms[0]);
    }

    /**
     * Initialization
     */
    private void initialize() {
        atoms = new ffx.potential.bonded.Atom[4];
        atoms[1] = bonds[0].getCommonAtom(bonds[1]);
        atoms[0] = bonds[0].get1_2(atoms[1]);
        atoms[2] = bonds[1].get1_2(atoms[1]);
        atoms[3] = bonds[2].get1_2(atoms[2]);
        setID_Key(false);
        value = calculateDihedralAngle();
    }

    /**
     * <p>setFlipped.</p>
     *
     * @param flipped a boolean.
     */
    private void setFlipped(boolean flipped) {
        if (flipped) {
            constants[0] = stretchTorsionType.forceConstants[6];
            constants[1] = stretchTorsionType.forceConstants[7];
            constants[2] = stretchTorsionType.forceConstants[8];
            constants[6] = stretchTorsionType.forceConstants[0];
            constants[7] = stretchTorsionType.forceConstants[1];
            constants[8] = stretchTorsionType.forceConstants[2];
        } else {
            arraycopy(stretchTorsionType.forceConstants, 0, constants, 0, 9);
        }
    }

    /**
     * If the specified atom is not a central atom of <b>this</b> torsion, the
     * atom at the opposite end is returned. These atoms are said to be 1-4 to
     * each other.
     *
     * @param a Atom
     * @return Atom
     */
    public Atom get1_4(Atom a) {
        if (a == atoms[0]) {
            return atoms[3];
        }
        if (a == atoms[3]) {
            return atoms[0];
        }
        return null;
    }

    /**
     * Calculates the dihedral angle; useful for the constructor, when energy() has not yet been called.
     *
     * @return Value of the dihedral angle.
     */
    private double calculateDihedralAngle() {

        double theVal = 0.0;

        double[] a0 = new double[3];
        double[] a1 = new double[3];
        double[] a2 = new double[3];
        double[] a3 = new double[3];
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);

        // Vector from Atom 0 to Atom 1.
        double[] v01 = new double[3];
        // Vector from Atom 1 to Atom 2.
        double[] v12 = new double[3];
        // Vector from Atom 2 to Atom 3.
        double[] v23 = new double[3];
        diff(a1, a0, v01);
        diff(a2, a1, v12);
        diff(a3, a2, v23);

        // Vector v01 cross v12.
        double[] x0112 = new double[3];
        // Vector v12 cross v23.
        double[] x1223 = new double[3];
        // Vector x0112 cross x12_32.
        double[] x = new double[3];
        cross(v01, v12, x0112);
        cross(v12, v23, x1223);
        cross(x0112, x1223, x);
        double r01_12 = dot(x0112, x0112);
        double r12_23 = dot(x1223, x1223);
        double rr = sqrt(r01_12 * r12_23);
        if (rr != 0) {
            double r12 = r(v12);
            double cosine = dot(x0112, x1223) / rr;
            double sine = dot(v12, x) / (r12 * rr);
            theVal = toDegrees(acos(cosine));
            if (sine < 0.0) {
                theVal = -theVal;
            }
        }
        return theVal;
    }

    /**
     * Returns the array of stretch-torsion constants, in units of kcal/mol/A.
     *
     * @return Stretch-torsion constants.
     */
    public double[] getConstants() {
        return copyOf(constants, constants.length);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluate the Stretch-Torsion energy.
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
        energy = 0.0;
        value = 0.0;
        dEdL = 0.0;

        double[] a0 = new double[3];
        double[] a1 = new double[3];
        double[] a2 = new double[3];
        double[] a3 = new double[3];
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);

        // Vector from Atom 0 to Atom 1.
        double[] v01 = new double[3];
        // Vector from Atom 1 to Atom 2.
        double[] v12 = new double[3];
        // Vector from Atom 2 to Atom 3.
        double[] v23 = new double[3];
        diff(a1, a0, v01);
        diff(a2, a1, v12);
        diff(a3, a2, v23);

        double r01 = r(v01);
        double r12 = r(v12);
        double r23 = r(v23);
        if (min(min(r01, r12), r23) == 0.0) {
            return 0.0;
        }

        // Vector v01 cross v12.
        double[] x0112 = new double[3];
        // Vector v12 cross v23.
        double[] x1223 = new double[3];
        // Vector x0112 cross x12_32.
        double[] x = new double[3];
        cross(v01, v12, x0112);
        cross(v12, v23, x1223);
        cross(x0112, x1223, x);
        double r01_12 = dot(x0112, x0112);
        r01_12 = max(r01_12, 0.000001);
        double r12_23 = dot(x1223, x1223);
        r12_23 = max(r12_23, 0.000001);
        double rr = sqrt(r01_12 * r12_23);

        // Vector from Atom 0 to Atom 2.
        double[] v02 = new double[3];
        // Vector from Atom 1 to Atom 3.
        double[] v13 = new double[3];
        diff(a2, a0, v02);
        diff(a3, a1, v13);
        double cosine = dot(x0112, x1223) / rr;
        double sine = dot(v12, x) / (r12 * rr);

        value = toDegrees(acos(cosine));
        if (sine < 0.0) {
            value = -value;
        }

        // Compute multiple angle trigonometry and phase terms
        double[] tsin = torsionType.sine;
        double[] tcos = torsionType.cosine;
        double cosine2 = cosine * cosine - sine * sine;
        double sine2 = 2.0 * cosine * sine;
        double cosine3 = cosine * cosine2 - sine * sine2;
        double sine3 = cosine * sine2 + sine * cosine2;
        double phi1 = 1.0 + (cosine * tcos[0] + sine * tsin[0]);
        double phi2 = 1.0 + (cosine2 * tcos[1] + sine2 * tsin[1]);
        double phi3 = 1.0 + (cosine3 * tcos[2] + sine3 * tsin[2]);
        double dphi1 = (cosine * tsin[0] - sine * tcos[0]);
        double dphi2 = 2.0 * (cosine2 * tsin[1] - sine2 * tcos[1]);
        double dphi3 = 3.0 * (cosine3 * tsin[2] - sine3 * tcos[2]);

        // Get the stretch-torsion values for the first bond.
        double v1 = constants[0];
        double v2 = constants[1];
        double v3 = constants[2];
        double dr1 = r01 - bondType1.distance;
        double e1 = StretchTorsionType.units * dr1 * (v1 * phi1 + v2 * phi2 + v3 * phi3);

        // Get the stretch-torsion values for the second bond.
        double v4 = constants[3];
        double v5 = constants[4];
        double v6 = constants[5];
        double dr2 = r12 - bondType2.distance;
        double e2 = StretchTorsionType.units * dr2 * (v4 * phi1 + v5 * phi2 + v6 * phi3);

        // Get the stretch-torsion values for the third bond.
        double v7 = constants[6];
        double v8 = constants[7];
        double v9 = constants[8];
        double dr3 = r23 - bondType3.distance;
        double e3 = StretchTorsionType.units * dr3 * (v7 * phi1 + v8 * phi2 + v9 * phi3);

        energy = e1 + e2 + e3;

        if (esvTerm) {
            esvDerivLocal = energy * dedesvChain * lambda;
        }

        energy = energy * esvLambda * lambda;
        dEdL = energy * esvLambda;

        if (gradient || lambdaTerm) {

            // Compute derivative components for the first bond.
            double dedphi1 = StretchTorsionType.units * dr1 * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
            double ddr1 = StretchTorsionType.units * (v1 * phi1 + v2 * phi2 + v3 * phi3) / r01;

            double[] ddrdx = new double[3];
            double[] dedxt = new double[3];
            double[] dedxu = new double[3];
            scalar(v01, ddr1, ddrdx);
            cross(x0112, v12, dedxt);
            scalar(dedxt, dedphi1 / (r01_12 * r12), dedxt);
            cross(x1223, v12, dedxu);
            scalar(dedxu, -dedphi1 / (r12_23 * r12), dedxu);

            // Gradient for atoms 0, 1, 2 & 3.
            double[] g0 = new double[3];
            double[] g1 = new double[3];
            double[] g2 = new double[3];
            double[] g3 = new double[3];
            // Work vectors.
            double[] x1 = new double[3];
            double[] x2 = new double[3];

            // Determine chain rule components for the first bond.
            cross(dedxt, v12, g0);
            diff(g0, ddrdx, g0);
            cross(v02, dedxt, x1);
            cross(dedxu, v23, x2);
            sum(x1, x2, g1);
            sum(g1, ddrdx, g1);
            cross(dedxt, v01, x1);
            cross(v13, dedxu, x2);
            sum(x1, x2, g2);
            cross(dedxu, v12, g3);

            // Compute derivative components for the 2nd bond.
            double dedphi2 = StretchTorsionType.units * dr2 * (v4 * dphi1 + v5 * dphi2 + v6 * dphi3);
            double ddr2 = StretchTorsionType.units * (v4 * phi1 + v5 * phi2 + v6 * phi3) / r12;

            scalar(v12, ddr2, ddrdx);
            cross(x0112, v12, dedxt);
            scalar(dedxt, dedphi2 / (r01_12 * r12), dedxt);
            cross(x1223, v12, dedxu);
            scalar(dedxu, -dedphi2 / (r12_23 * r12), dedxu);

            //  Partial gradient on atoms 0, 1, 2 & 3.
            double[] pg0 = new double[3];
            double[] pg1 = new double[3];
            double[] pg2 = new double[3];
            double[] pg3 = new double[3];

            // Compute chain rule components for the second bond.
            cross(dedxt, v12, pg0);
            cross(v02, dedxt, x1);
            cross(dedxu, v23, x2);
            sum(x1, x2, pg1);
            diff(pg1, ddrdx, pg1);
            cross(dedxt, v01, x1);
            cross(v13, dedxu, x2);
            sum(x1, x2, pg2);
            sum(pg2, ddrdx, pg2);
            cross(dedxu, v12, pg3);

            // Accumulate derivative components.
            sum(pg0, g0, g0);
            sum(pg1, g1, g1);
            sum(pg2, g2, g2);
            sum(pg3, g3, g3);

            // Compute derivative components for the 3rd bond.
            double dedphi3 = StretchTorsionType.units * dr3 * (v7 * dphi1 + v8 * dphi2 + v9 * dphi3);
            double ddr3 = StretchTorsionType.units * (v7 * phi1 + v8 * phi2 + v9 * phi3) / r23;
            scalar(v23, ddr3, ddrdx);
            cross(x0112, v12, dedxt);
            scalar(dedxt, dedphi3 / (r01_12 * r12), dedxt);
            cross(x1223, v12, dedxu);
            scalar(dedxu, -dedphi3 / (r12_23 * r12), dedxu);

            // Compute chain rule components for the third bond.
            cross(dedxt, v12, pg0);
            cross(v02, dedxt, x1);
            cross(dedxu, v23, x2);
            sum(x1, x2, pg1);
            cross(dedxt, v01, x1);
            cross(v13, dedxu, x2);
            sum(x1, x2, pg2);
            diff(pg2, ddrdx, pg2);
            cross(dedxu, v12, pg3);
            sum(pg3, ddrdx, pg3);

            // Accumulate derivative components.
            sum(pg0, g0, g0);
            sum(pg1, g1, g1);
            sum(pg2, g2, g2);
            sum(pg3, g3, g3);

            scalar(g0, esvLambda, g0);
            scalar(g1, esvLambda, g1);
            scalar(g2, esvLambda, g2);
            scalar(g3, esvLambda, g3);

            int i0 = atoms[0].getIndex() - 1;
            int i1 = atoms[1].getIndex() - 1;
            int i2 = atoms[2].getIndex() - 1;
            int i3 = atoms[3].getIndex() - 1;

            if (lambdaTerm) {
                lambdaGradX.add(threadID, i0, g0[0]);
                lambdaGradY.add(threadID, i0, g0[1]);
                lambdaGradZ.add(threadID, i0, g0[2]);

                lambdaGradX.add(threadID, i1, g1[0]);
                lambdaGradY.add(threadID, i1, g1[1]);
                lambdaGradZ.add(threadID, i1, g1[2]);

                lambdaGradX.add(threadID, i2, g2[0]);
                lambdaGradY.add(threadID, i2, g2[1]);
                lambdaGradZ.add(threadID, i2, g2[2]);

                lambdaGradX.add(threadID, i3, g3[0]);
                lambdaGradY.add(threadID, i3, g3[1]);
                lambdaGradZ.add(threadID, i3, g3[2]);
            }
            if (gradient) {
                scalar(g0, lambda, g0);
                scalar(g1, lambda, g1);
                scalar(g2, lambda, g2);
                scalar(g3, lambda, g3);

                gradX.add(threadID, i0, g0[0]);
                gradY.add(threadID, i0, g0[1]);
                gradZ.add(threadID, i0, g0[2]);

                gradX.add(threadID, i1, g1[0]);
                gradY.add(threadID, i1, g1[1]);
                gradZ.add(threadID, i1, g1[2]);

                gradX.add(threadID, i2, g2[0]);
                gradY.add(threadID, i2, g2[1]);
                gradZ.add(threadID, i2, g2[2]);

                gradX.add(threadID, i3, g3[0]);
                gradY.add(threadID, i3, g3[1]);
                gradZ.add(threadID, i3, g3[2]);
            }
        }

        return energy;
    }

    /**
     * Log details for this Torsional Angle energy term.
     */
    public void log() {
        logger.info(String.format(" %-8s %6d-%s %6d-%s %6d-%s %6d-%s %10.4f",
                "Stretch-Torsion",
                atoms[0].getIndex(), atoms[0].getAtomType().name,
                atoms[1].getIndex(), atoms[1].getAtomType().name,
                atoms[2].getIndex(), atoms[2].getAtomType().name,
                atoms[3].getIndex(), atoms[3].getAtomType().name, energy));
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden toString Method returns the Term's id.
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

    /**
     * Returns the mathematical form of a stretch-torsion as an OpenMM-parsable String.
     *
     * @return Mathematical form of the stretch-torsion coupling.
     */
    public static String stretchTorsionForm() {
        return mathForm;
    }

    /**
     * Functional form for OpenMM.
     */
    private static final String mathForm;

    static {
        /*
          Defined constants:
          p1-p4 are particles 1-4.
          m is a bond number, from 1-3, representing bonds p1-p2, p2-p3, p3-p4.
          n is a periodicity, from 1-3.

          k[m][n] is a set of 9 energy constants defined in the parameter file for this stretch-torsion.

          bVal[m] is the current bond distance for bond m.
          b[m] is the equilibrium distance for bond m.

          tVal is the current value of the 1-2-3-4 dihedral angle.

          phi[m] is a phase offset constant; phi1 = phi3 = 0, phi2 = pi.
         */

        StringBuilder mathFormBuilder = new StringBuilder();

        for (int m = 1; m < 4; m++) {
            for (int n = 1; n < 4; n++) {
                // kmn * (bm - bm(equil)) * (1 + cos(n*tors + phi(n)))
                mathFormBuilder.append(String.format("k%d%d*(bVal%d-b%d)*(1+cos(%d*tVal+phi%d))+", m, n, m, m, n, n));
            }
        }
        int lenStr = mathFormBuilder.length();
        mathFormBuilder.replace(lenStr - 1, lenStr, ";");

        for (int m = 1; m < 4; m++) {
            mathFormBuilder.append(String.format("bVal%d=distance(p%d,p%d);", m, m, (m + 1)));
        }

        mathFormBuilder.append("tVal=dihedral(p1,p2,p3,p4)");
        mathForm = mathFormBuilder.toString();
    }
}
