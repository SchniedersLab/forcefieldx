/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ImproperTorsionType;

import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;

/**
 * The ImproperTorsion class represents an Improper Torsion.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class ImproperTorsion extends BondedTerm implements
        Comparable<ImproperTorsion> {

    private static final Logger logger = Logger.getLogger(ImproperTorsion.class.getName());
    /**
     * Force field parameters to compute the ImproperTorsion energy.
     */
    public ImproperTorsionType improperType = null;

    public double scaleFactor = 1.0;

    /**
     * Convert Torsional Angle energy to kcal/mole.
     *
     * @since 1.0
     */
    public double units = 0.5;

    /**
     * ImproperTorsion constructor.
     *
     * @param atom1 Atom number 1.
     * @param atom2 Atom number 2.
     * @param atom3 Atom number 3.
     * @param atom4 Atom number 4.
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
//        print();
    }

    /**
     * ImproperTorsion factory method.
     *
     * @param atom Create the improper torsion around this atom.
     * @param forceField retrieve parameters from this ForceField.
     * @return the ImproperTorsion if created, or null otherwise.
     */
    public static ArrayList<ImproperTorsion> improperTorsionFactory(Atom atom, ForceField forceField) {

        if (atom == null) {
            return null;
        }

        Atom atoms[] = new Atom[4];
        atoms[2] = atom;
        ArrayList<Bond> bonds = atom.getBonds();
        if (bonds == null || bonds.size() != 3) {
            return null;
        }

        for (int i = 0; i < 3; i++) {
            Bond bond = bonds.get(i);
            Atom atom2 = bond.get1_2(atom);
            if (i == 2) {
                atoms[3] = atom2;
            } else {
                atoms[i] = atom2;
            }
        }

        ArrayList<ImproperTorsion> improperTorsions = new ArrayList<>();
        Collection<ImproperTorsionType> types = forceField.getImproperTypes();

        double units = forceField.getDouble(ForceField.ForceFieldDouble.IMPTORUNIT, 0.5);

        for (ImproperTorsionType type : types) {
            int classes[] = new int[4];
            classes[0] = atoms[0].getAtomType().atomClass;
            classes[1] = atoms[1].getAtomType().atomClass;
            classes[2] = atoms[2].getAtomType().atomClass;
            classes[3] = atoms[3].getAtomType().atomClass;
            boolean assigned = type.assigned(classes);
            if (assigned) {
                // Finalize atom ordering.
                if (classes[3] == atoms[3].getAtomType().atomClass || type.atomClasses[3] == 0) {
                    // do nothing.
                } else if (classes[3] == atoms[1].getAtomType().atomClass) {
                    Atom temp = atoms[3];
                    atoms[3] = atoms[1];
                    atoms[1] = temp;
                } else {
                    Atom temp = atoms[0];
                    atoms[0] = atoms[3];
                    atoms[3] = temp;
                }
                if (classes[1] == atoms[1].getAtomType().atomClass || type.atomClasses[1] == 0) {
                    // do nothing.
                } else if (classes[1] == atoms[0].getAtomType().atomClass) {
                    Atom temp = atoms[1];
                    atoms[1] = atoms[0];
                    atoms[0] = temp;
                }
                ImproperTorsion improperTorsion = new ImproperTorsion(atoms[0], atoms[1], atoms[2], atoms[3]);
                improperTorsion.setImproperType(type);
                improperTorsion.units = units;
                improperTorsions.add(improperTorsion);
                int c0 = type.atomClasses[0];
                int c1 = type.atomClasses[1];
                int c3 = type.atomClasses[3];
                if (c0 == c1 && c1 == c3) {
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsion = new ImproperTorsion(atoms[0], atoms[3], atoms[2], atoms[1]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsions.add(improperTorsion);
                    improperTorsion = new ImproperTorsion(atoms[1], atoms[0], atoms[2], atoms[3]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsions.add(improperTorsion);
                    improperTorsion = new ImproperTorsion(atoms[1], atoms[3], atoms[2], atoms[0]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsions.add(improperTorsion);
                    improperTorsion = new ImproperTorsion(atoms[3], atoms[0], atoms[2], atoms[1]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsions.add(improperTorsion);
                    improperTorsion = new ImproperTorsion(atoms[3], atoms[1], atoms[2], atoms[0]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 1.0 / 6.0;
                    improperTorsions.add(improperTorsion);
                } else if (c0 == c1) {
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsion = new ImproperTorsion(atoms[1], atoms[0], atoms[2], atoms[3]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsions.add(improperTorsion);
                } else if (c0 == c3) {
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsion = new ImproperTorsion(atoms[3], atoms[1], atoms[2], atoms[0]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsions.add(improperTorsion);
                } else if (c1 == c3) {
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsion = new ImproperTorsion(atoms[0], atoms[3], atoms[2], atoms[1]);
                    improperTorsion.setImproperType(type);
                    improperTorsion.units = units;
                    improperTorsion.scaleFactor = 0.5;
                    improperTorsions.add(improperTorsion);
                }
            }
        }

        if (improperTorsions.isEmpty()) {
            return null;
        }

        return improperTorsions;
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
     * Update recomputes ImproperTorsion value and energy.
     */
    @Override
    public void update() {
        energy(false);
    }

    /**
     * Evaluate this Improper Torsion energy.
     *
     * @param gradient Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        energy = 0.0;
        value = 0.0;
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);
        diff(a1, a0, v10);
        diff(a2, a1, v21);
        diff(a3, a2, v32);
        cross(v10, v21, t);
        cross(v21, v32, u);
        cross(t, u, tu);
        double rt2 = dot(t, t);
        double ru2 = dot(u, u);
        double rtru = sqrt(rt2 * ru2);
        if (rtru != 0.0) {

            /* Set the improper torsional parameters for this angle */
            //double v1 = 0.0;
            //double c1 = 0.0;
            //double s1 = 0.0;
            double v2 = improperType.k;
            double c2 = improperType.cos;
            double s2 = improperType.sin;
            //double v3 = 0.0;
            // double c3 = 0.0;
            // double s3 = 0.0;

            /* compute the multiple angle trigonometry and the phase terms */
            double rcb = r(v21);
            double cosine = dot(t, u) / rtru;
            double sine = dot(v21, tu) / (rcb * rtru);
            double cosine2 = cosine * cosine - sine * sine;
            double sine2 = 2.0 * cosine * sine;
            //double cosine3 = cosine * cosine2 - sine * sine2;
            //double sine3 = cosine * sine2 + sine * cosine2;
            //double phi1 = 1.0 + (cosine * c1 + sine * s1);
            double phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);
            //double phi3 = 1.0 + (cosine3 * c3 + sine3 * s3);
            //double dphi1 = (cosine * s1 - sine * c1);
            double dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2);
            //double dphi3 = 3.0 * (cosine3 * s3 - sine3 * c3);

            /* calculate improper torsion energy and master chain rule term */
            value = Math.toDegrees(Math.acos(cosine));
            //energy = ImproperTorsionType.units * (v1 * phi1 + v2 * phi2 + v3 * phi3);
            //double dedphi = ImproperTorsionType.units * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
            energy = units * (v2 * phi2) * scaleFactor;
            double dedphi = units * (v2 * dphi2) * scaleFactor;

            if (gradient) {
                /**
                 * Chain rule terms for first derivative components.
                 */
                diff(a2, a0, v20);
                diff(a3, a1, v31);
                cross(t, v21, dedt);
                cross(u, v21, dedu);
                scalar(dedt, dedphi / (rt2 * rcb), dedt);
                scalar(dedu, -dedphi / (ru2 * rcb), dedu);
                /**
                 * Compute first derivative components for this angle.
                 */
                cross(dedt, v21, g0);
                cross(dedt, v20, g1a);
                cross(dedu, v32, g1);
                scalar(g1a, -1.0, g1a);
                sum(g1a, g1, g1);
                cross(dedt, v10, g2a);
                cross(dedu, v31, g2);
                scalar(g2, -1.0, g2);
                sum(g2a, g2, g2);
                cross(dedu, v21, g3);
                /**
                 * Accumulate derivatives.
                 */
                atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
                atoms[1].addToXYZGradient(g1[0], g1[1], g1[2]);
                atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
                atoms[3].addToXYZGradient(g3[0], g3[1], g3[2]);
            }
        }

        // log();
        return energy;
    }

    /**
     * Log details for this Improper Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s %6d-%s %6.4f %10.4f",
                "Improper Torsion", atoms[0].getXYZIndex(), atoms[0].getAtomType().name, atoms[1].getXYZIndex(), atoms[1].getAtomType().name, atoms[2].getXYZIndex(), atoms[2].getAtomType().name, atoms[3].getXYZIndex(), atoms[3].getAtomType().name,
                value, energy));
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
    protected static final double a0[] = new double[3];
    protected static final double a1[] = new double[3];
    protected static final double a2[] = new double[3];
    protected static final double a3[] = new double[3];
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
    /**
     * Vector from Atom 3 to Atom 1.
     */
    protected static final double v31[] = new double[3];
    /**
     * Vector from Atom 2 to Atom 0.
     */
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
    protected static final double g1a[] = new double[3];
    /**
     * Work array.
     */
    protected static final double g2a[] = new double[3];
}
