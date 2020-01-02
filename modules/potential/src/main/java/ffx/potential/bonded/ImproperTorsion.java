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
import java.util.Collection;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ImproperTorsionType;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;

/**
 * The ImproperTorsion class represents an Improper Torsion.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ImproperTorsion extends BondedTerm {

    private static final Logger logger = Logger.getLogger(ImproperTorsion.class.getName());

    /**
     * Force field parameters to compute the ImproperTorsion energy.
     */
    public ImproperTorsionType improperType = null;
    /**
     * Scale factor.
     */
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
    private ImproperTorsion(Atom atom1, Atom atom2, Atom atom3, Atom atom4) {
        super();
        atoms = new Atom[4];
        atoms[0] = atom1;
        atoms[1] = atom2;
        atoms[2] = atom3;
        atoms[3] = atom4;
        setID_Key(false);
    }

    /**
     * ImproperTorsion factory method.
     *
     * @param atom       Create the improper torsion around this atom.
     * @param forceField retrieve parameters from this ForceField.
     * @return the ImproperTorsion if created, or null otherwise.
     */
    static ArrayList<ImproperTorsion> improperTorsionFactory(Atom atom, ForceField forceField) {

        if (atom == null) {
            return null;
        }

        Atom[] atoms = new Atom[4];
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
        double units = forceField.getDouble("IMPTORUNIT", 1.0);
        boolean done = false;

        // No wild card matches.
        for (ImproperTorsionType type : types) {
            int[] classes = new int[4];
            classes[0] = atoms[0].getAtomType().atomClass;
            classes[1] = atoms[1].getAtomType().atomClass;
            classes[2] = atoms[2].getAtomType().atomClass;
            classes[3] = atoms[3].getAtomType().atomClass;
            boolean assigned = type.assigned(classes, false, false);
            if (assigned) {
                done = true;

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

            if (done) {
                break;
            }
        }

        // Wild card matches for the first or second class.
        if (!done) {
            for (ImproperTorsionType type : types) {
                int[] classes = new int[4];
                classes[0] = atoms[0].getAtomType().atomClass;
                classes[1] = atoms[1].getAtomType().atomClass;
                classes[2] = atoms[2].getAtomType().atomClass;
                classes[3] = atoms[3].getAtomType().atomClass;
                boolean assigned = type.assigned(classes, true, false);
                if (assigned) {
                    done = true;
                    createWildCardImproperTorsion(atoms, classes, type, units, improperTorsions);
                    break;
                }
            }
        }

        // Wild card matches for the first, second and third classes.
        if (!done) {
            for (ImproperTorsionType type : types) {
                int[] classes = new int[4];
                classes[0] = atoms[0].getAtomType().atomClass;
                classes[1] = atoms[1].getAtomType().atomClass;
                classes[2] = atoms[2].getAtomType().atomClass;
                classes[3] = atoms[3].getAtomType().atomClass;
                boolean assigned = type.assigned(classes, true, true);
                if (assigned) {
                    createWildCardImproperTorsion(atoms, classes, type, units, improperTorsions);
                    break;
                }
            }
        }

        if (improperTorsions.isEmpty()) {
            return null;
        }

        return improperTorsions;
    }

    private static void createWildCardImproperTorsion(Atom[] atoms, int[] classes, ImproperTorsionType type,
                                                      double units, ArrayList<ImproperTorsion> improperTorsions) {
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

        // One or more zero classes in the improper torsion type
        ImproperTorsion improperTorsion = new ImproperTorsion(atoms[0], atoms[1], atoms[2], atoms[3]);
        improperTorsion.setImproperType(type);
        improperTorsion.units = units;
        improperTorsions.add(improperTorsion);
        improperTorsion.scaleFactor = 1.0 / 3.0;
        improperTorsions.add(improperTorsion);

        improperTorsion = new ImproperTorsion(atoms[1], atoms[3], atoms[2], atoms[0]);
        improperTorsion.setImproperType(type);
        improperTorsion.units = units;
        improperTorsions.add(improperTorsion);
        improperTorsion.scaleFactor = 1.0 / 3.0;
        improperTorsions.add(improperTorsion);

        improperTorsion = new ImproperTorsion(atoms[3], atoms[0], atoms[2], atoms[1]);
        improperTorsion.setImproperType(type);
        improperTorsion.units = units;
        improperTorsions.add(improperTorsion);
        improperTorsion.scaleFactor = 1.0 / 3.0;
        improperTorsions.add(improperTorsion);
    }

    /**
     * Set a reference to the force field parameters for <b>this</b> Angle.
     *
     * @param a a {@link ffx.potential.parameters.ImproperTorsionType} object.
     */
    private void setImproperType(ImproperTorsionType a) {
        improperType = a;
    }

    /**
     * Log details for this Improper Torsion energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6d-%s %6d-%s %6.4f %10.4f",
                "Improper Torsion", atoms[0].getIndex(), atoms[0].getAtomType().name, atoms[1].getIndex(), atoms[1].getAtomType().name, atoms[2].getIndex(), atoms[2].getAtomType().name, atoms[3].getIndex(), atoms[3].getAtomType().name,
                value, energy));
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluate this Improper Torsion energy.
     */
    @Override
    public double energy(boolean gradient, int threadID,
                         AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
        energy = 0.0;
        value = 0.0;
        double[] a0 = new double[3];
        double[] a1 = new double[3];
        double[] a2 = new double[3];
        double[] a3 = new double[3];
        atoms[0].getXYZ(a0);
        atoms[1].getXYZ(a1);
        atoms[2].getXYZ(a2);
        atoms[3].getXYZ(a3);

        // Vector from Atom 1 to Atom 0.
        double[] v10 = new double[3];
        // Vector from Atom 2 to Atom 1.
        double[] v21 = new double[3];
        // Vector from Atom 3 to Atom 2.
        double[] v32 = new double[3];
        diff(a1, a0, v10);
        diff(a2, a1, v21);
        diff(a3, a2, v32);

        double[] t = new double[3];
        double[] u = new double[3];
        double[] tu = new double[3];
        double[] dedu = new double[3];
        double[] dedt = new double[3];
        cross(v10, v21, t);
        cross(v21, v32, u);
        cross(t, u, tu);
        double rt2 = dot(t, t);
        double ru2 = dot(u, u);
        double rtru = sqrt(rt2 * ru2);
        if (rtru != 0.0) {

            // Set the improper torsional parameters for this angle
            double v2 = improperType.k;
            double c2 = improperType.cos;
            double s2 = improperType.sin;

            // Compute the multiple angle trigonometry and the phase terms
            double rcb = r(v21);
            double cosine = dot(t, u) / rtru;
            double sine = dot(v21, tu) / (rcb * rtru);
            double cosine2 = cosine * cosine - sine * sine;
            double sine2 = 2.0 * cosine * sine;
            double phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);
            double dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2);

            // Calculate improper torsion energy and master chain rule term
            value = toDegrees(acos(cosine));
            final double desvPrefactor = units * scaleFactor;
            final double prefactor = units * scaleFactor * esvLambda;
            energy = prefactor * (v2 * phi2);
            double dedphi = prefactor * (v2 * dphi2);
            if (esvTerm) {
                setEsvDeriv(desvPrefactor * (v2 * phi2) * dedesvChain);
            }

            if (gradient) {
                // Vector from Atom 3 to Atom 1.
                double[] v31 = new double[3];
                // Vector from Atom 2 to Atom 0.
                double[] v20 = new double[3];
                // Chain rule terms for first derivative components.
                diff(a2, a0, v20);
                diff(a3, a1, v31);
                cross(t, v21, dedt);
                cross(u, v21, dedu);
                scalar(dedt, dedphi / (rt2 * rcb), dedt);
                scalar(dedu, -dedphi / (ru2 * rcb), dedu);
                // Compute first derivative components for this angle.
                // Gradient on atoms 0, 1, 2 & 3.
                double[] g0 = new double[3];
                double[] g1 = new double[3];
                double[] g2 = new double[3];
                double[] g3 = new double[3];
                // Work arrays.
                double[] g1a = new double[3];
                double[] g2a = new double[3];
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
                // Atom indices
                int i0 = atoms[0].getIndex() - 1;
                int i1 = atoms[1].getIndex() - 1;
                int i2 = atoms[2].getIndex() - 1;
                int i3 = atoms[3].getIndex() - 1;
                // Accumulate derivatives.
                grad.add(threadID, i0, g0[0], g0[1], g0[2]);
                grad.add(threadID, i1, g1[0], g1[1], g1[2]);
                grad.add(threadID, i2, g2[0], g2[1], g2[2]);
                grad.add(threadID, i3, g3[0], g3[1], g3[2]);
            }
        }

        return energy;
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

}
