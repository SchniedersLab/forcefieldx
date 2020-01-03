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
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import org.apache.commons.math3.util.FastMath;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BioType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.norm;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;

/**
 * Utilities for placing atoms.
 *
 * @author Michael Schnieders
 * @since 1.0
 */
public class BondedUtils {

    private static final Logger logger = Logger.getLogger(BondedUtils.class.getName());
    private static final double eps = 0.0000001d;

    /**
     * This routine was derived from a similar routine in TINKER.
     *
     * @param atom   a {@link ffx.potential.bonded.Atom} object.
     * @param ia     a {@link ffx.potential.bonded.Atom} object.
     * @param bond   a double.
     * @param ib     a {@link ffx.potential.bonded.Atom} object.
     * @param angle1 a double.
     * @param ic     a {@link ffx.potential.bonded.Atom} object.
     * @param angle2 a double.
     * @param chiral a int.
     */
    public static void intxyz(Atom atom, Atom ia, double bond, Atom ib, double angle1, Atom ic, double angle2, int chiral) {
        double[] xa = new double[3];
        xa = (ia == null) ? null : ia.getXYZ(xa);
        double[] xb = new double[3];
        xb = (ib == null) ? null : ib.getXYZ(xb);
        double[] xc = new double[3];
        xc = (ic == null) ? null : ic.getXYZ(xc);
        atom.moveTo(determineIntxyz(xa, bond, xb, angle1, xc, angle2, chiral));
    }

    /**
     * This routine was derived from a similar routine in TINKER. It determines
     * at what coordinates an atom would be placed without moving or calling any
     * atoms, relying solely upon coordinates. Passed arrays are copied into
     * local arrays to avoid any over-writing of the passed arrays.
     * <p>
     * The chiral argument is 0 if angle2 is a dihedral.
     * Else, if angle2 is the atom-ia-ic angle:
     * -1 indicates left-hand-rule placement.
     * +1 indicates right-hand-rule placement.
     * +3 indicates trigonal planar placement.
     * <p>
     * Chiral +3 replaces the angle1 and angle2 constraints with a planarity
     * constraint, and minimized, equipartitioned deviation from angle1 and angle2.
     *
     * @param ia     a double[] of atomic coordinates.
     * @param bond   a double.
     * @param ib     a double[] of atomic coordinates.
     * @param angle1 a double.
     * @param ic     a double[] of atomic coordinates.
     * @param angle2 a double.
     * @param chiral 0, 1, -1, or 3.
     * @return A double[] with XYZ coordinates at which an atom would be placed.
     */
    static double[] determineIntxyz(double[] ia, double bond, double[] ib, double angle1, double[] ic, double angle2, int chiral) {
        if (ia != null && !Double.isFinite(bond)) {
            throw new IllegalArgumentException(String.format(" Passed bond length is non-finite %f", bond));
        } else if (ib != null && !Double.isFinite(angle1)) {
            throw new IllegalArgumentException(String.format(" Passed angle is non-finite %f", angle1));
        } else if (ic != null && !Double.isFinite(angle2)) {
            throw new IllegalArgumentException(String.format(" Passed dihedral/improper is non-finite %f", angle2));
        }

        if (chiral == 3) {
            double[] negChiral = determineIntxyz(ia, bond, ib, angle1, ic, angle2, -1);
            double[] posChiral = determineIntxyz(ia, bond, ib, angle1, ic, angle2, 1);
            double[] displacement = new double[3];
            double dispMag = 0;

            for (int i = 0; i < 3; i++) {
                // First, draw the midpoint between the positive- and negative- chiral solutions.
                displacement[i] = 0.5 * (posChiral[i] + negChiral[i]);
                // Second, take the displacement from a2 to this midpoint.
                displacement[i] -= ia[i];
                // Third, accumulate into the vector magnitude.
                dispMag += (displacement[i] * displacement[i]);
            }

            dispMag = FastMath.sqrt(dispMag);
            double extend = bond / dispMag;
            assert extend > 0.999; // Should be >= 1.0, with slight machine-precision tolerance.

            double[] outXYZ = new double[3];
            for (int i = 0; i < 3; i++) {
                displacement[i] *= extend;
                outXYZ[i] = displacement[i] + ia[i];
            }
            return outXYZ;
        }

        angle1 = toRadians(angle1);
        angle2 = toRadians(angle2);
        double zcos0 = cos(angle1);
        double zcos1 = cos(angle2);
        double zsin0 = sin(angle1);
        double zsin1 = sin(angle2);
        double[] ret = new double[3];
        double x[] = new double[3];

        // No partners
        if (ia == null) {
            x[0] = x[1] = x[2] = 0.0;
        } else if (ib == null) {
            double xa[] = new double[3];
            arraycopy(ia, 0, xa, 0, ia.length);
            // One partner - place on the z-axis
            x[0] = xa[0];
            x[1] = xa[1];
            x[2] = xa[2] + bond;
        } else if (ic == null) {
            double[] xa = new double[3];
            double[] xb = new double[3];
            double[] xab = new double[3];
            for (int i = 0; i < ia.length; i++) {
                xa[i] = ia[i];
                xb[i] = ib[i];
            }
            // Two partners - place in the xz-plane
            diff(xa, xb, xab);
            double rab = r(xab);
            norm(xab, xab);
            double cosb = xab[2];
            double sinb = sqrt(xab[0] * xab[0] + xab[1] * xab[1]);
            double cosg, sing;
            if (sinb == 0.0d) {
                cosg = 1.0d;
                sing = 0.0d;
            } else {
                cosg = xab[1] / sinb;
                sing = xab[0] / sinb;
            }
            double xtmp = bond * zsin0;
            double ztmp = rab - bond * zcos0;
            x[0] = xb[0] + xtmp * cosg + ztmp * sing * sinb;
            x[1] = xb[1] - xtmp * sing + ztmp * cosg * sinb;
            x[2] = xb[2] + ztmp * cosb;
        } else if (chiral == 0) {
            double[] xa = new double[3];
            double[] xb = new double[3];
            double[] xc = new double[3];
            double[] xab = new double[3];
            double[] xbc = new double[3];
            double[] xt = new double[3];
            double[] xu = new double[3];
            for (int i = 0; i < ia.length; i++) {
                xa[i] = ia[i];
                xb[i] = ib[i];
                xc[i] = ic[i];
            }
            // General case - with a dihedral
            diff(xa, xb, xab);
            norm(xab, xab);
            diff(xb, xc, xbc);
            norm(xbc, xbc);
            xt[0] = xab[2] * xbc[1] - xab[1] * xbc[2];
            xt[1] = xab[0] * xbc[2] - xab[2] * xbc[0];
            xt[2] = xab[1] * xbc[0] - xab[0] * xbc[1];
            double cosine = xab[0] * xbc[0] + xab[1] * xbc[1] + xab[2] * xbc[2];
            double sine = sqrt(max(1.0d - cosine * cosine, eps));
            if (abs(cosine) >= 1.0d) {
                logger.warning("Undefined Dihedral");
            }
            scalar(xt, 1.0d / sine, xt);
            xu[0] = xt[1] * xab[2] - xt[2] * xab[1];
            xu[1] = xt[2] * xab[0] - xt[0] * xab[2];
            xu[2] = xt[0] * xab[1] - xt[1] * xab[0];
            x[0] = xa[0] + bond * (xu[0] * zsin0 * zcos1 + xt[0] * zsin0 * zsin1 - xab[0] * zcos0);
            x[1] = xa[1] + bond * (xu[1] * zsin0 * zcos1 + xt[1] * zsin0 * zsin1 - xab[1] * zcos0);
            x[2] = xa[2] + bond * (xu[2] * zsin0 * zcos1 + xt[2] * zsin0 * zsin1 - xab[2] * zcos0);
        } else if (abs(chiral) == 1) {
            double[] xa = new double[3];
            double[] xb = new double[3];
            double[] xc = new double[3];
            double[] xba = new double[3];
            double[] xac = new double[3];
            double[] xt = new double[3];
            for (int i = 0; i < ia.length; i++) {
                xa[i] = ia[i];
                xb[i] = ib[i];
                xc[i] = ic[i];
            }
            diff(xb, xa, xba);
            norm(xba, xba);
            diff(xa, xc, xac);
            norm(xac, xac);
            xt[0] = xba[2] * xac[1] - xba[1] * xac[2];
            xt[1] = xba[0] * xac[2] - xba[2] * xac[0];
            xt[2] = xba[1] * xac[0] - xba[0] * xac[1];
            double cosine = xba[0] * xac[0] + xba[1] * xac[1] + xba[2] * xac[2];
            double sine2 = max(1.0d - cosine * cosine, eps);
            if (abs(cosine) >= 1.0d) {
                logger.warning("Defining Atom Colinear");
            }
            double a = (-zcos1 - cosine * zcos0) / sine2;
            double b = (zcos0 + cosine * zcos1) / sine2;
            double c = (1.0d + a * zcos1 - b * zcos0) / sine2;
            if (c > eps) {
                c = chiral * sqrt(c);
            } else if (c < -eps) {
                c = sqrt((a * xac[0] + b * xba[0]) * (a * xac[0] + b * xba[0])
                        + (a * xac[1] + b * xba[1]) * (a * xac[1] + b * xba[1])
                        + (a * xac[2] + b * xba[2]) * (a * xac[2] + b * xba[2]));
                a /= c;
                b /= c;
                c = 0.0d;
            } else {
                c = 0.0d;
            }
            x[0] = xa[0] + bond * (a * xac[0] + b * xba[0] + c * xt[0]);
            x[1] = xa[1] + bond * (a * xac[1] + b * xba[1] + c * xt[1]);
            x[2] = xa[2] + bond * (a * xac[2] + b * xba[2] + c * xt[2]);
        }
        arraycopy(x, 0, ret, 0, ret.length);
        return ret;
    }

    /**
     * <p>findAtomType.</p>
     *
     * @param key        a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @return a {@link ffx.potential.parameters.AtomType} object.
     */
    public static AtomType findAtomType(int key, ForceField forceField) {
        BioType bioType = forceField.getBioType(Integer.toString(key));
        if (bioType != null) {
            AtomType atomType = forceField.getAtomType(Integer.toString(bioType.atomType));
            if (atomType != null) {
                return atomType;
            } else {
                logger.severe(format("The atom type %s was not found for biotype %s.", bioType.atomType,
                        bioType.toString()));
            }
        }
        return null;
    }

    /**
     * <p>buildBond.</p>
     *
     * @param a1         a {@link ffx.potential.bonded.Atom} object.
     * @param a2         a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Bond} object.
     */
    public static Bond buildBond(Atom a1, Atom a2, ForceField forceField, ArrayList<Bond> bondList) {
        Bond bond = new Bond(a1, a2);
        int[] c = new int[2];
        c[0] = a1.getAtomType().atomClass;
        c[1] = a2.getAtomType().atomClass;
        String key = BondType.sortKey(c);
        BondType bondType = forceField.getBondType(key);
        if (bondType == null) {
            logger.severe(format("No BondType for key: %s\n %s\n %s\n %s\n %s", key,
                    a1.toString(), a1.getAtomType().toString(),
                    a2.toString(), a2.getAtomType().toString()));
        } else {
            bond.setBondType(bondType);
        }
        if (bondList != null) {
            bondList.add(bond);
        }
        return bond;
    }

    /**
     * <p>buildHeavy.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.MSGroup} object.
     * @param atomName   a {@link java.lang.String} object.
     * @param bondedTo   a {@link ffx.potential.bonded.Atom} object.
     * @param key        a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    public static Atom buildHeavy(MSGroup residue, String atomName, Atom bondedTo, int key, ForceField forceField, ArrayList<Bond> bondList)
            throws MissingHeavyAtomException {
        Atom atom = (Atom) residue.getAtomNode(atomName);
        AtomType atomType = findAtomType(key, forceField);
        if (atom == null) {
            throw new MissingHeavyAtomException(atomName, atomType, bondedTo);
        }
        atom.setAtomType(atomType);
        if (bondedTo != null) {
            buildBond(atom, bondedTo, forceField, bondList);
        }
        return atom;
    }

    /**
     * <p>buildHeavyAtom.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.MSGroup} object.
     * @param atomName   a {@link java.lang.String} object.
     * @param ia         a {@link ffx.potential.bonded.Atom} object.
     * @param bond       a double.
     * @param ib         a {@link ffx.potential.bonded.Atom} object.
     * @param angle1     a double.
     * @param ic         a {@link ffx.potential.bonded.Atom} object.
     * @param angle2     a double.
     * @param chiral     a int.
     * @param atomType   a {@link ffx.potential.parameters.AtomType} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    private static Atom buildHeavyAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
                                       Atom ic, double angle2, int chiral, AtomType atomType, ForceField forceField, ArrayList<Bond> bondList) {
        Atom atom = (Atom) residue.getAtomNode(atomName);
        if (atomType == null) {
            return null;
        }
        if (atom == null) {
            String resName = ia.getResidueName();
            int resSeq = ia.getResidueNumber();
            Character chainID = ia.getChainID();
            Character altLoc = ia.getAltLoc();
            String segID = ia.getSegID();
            double occupancy = ia.getOccupancy();
            double tempFactor = ia.getTempFactor();
            atom = new Atom(0, atomName, altLoc, new double[3], resName, resSeq, chainID,
                    occupancy, tempFactor, segID, true);
            residue.addMSNode(atom);
            intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
        }
        atom.setAtomType(atomType);
        buildBond(ia, atom, forceField, bondList);
        return atom;
    }

    /**
     * <p>buildHeavy.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.MSGroup} object.
     * @param atomName   a {@link java.lang.String} object.
     * @param ia         a {@link ffx.potential.bonded.Atom} object.
     * @param bond       a double.
     * @param ib         a {@link ffx.potential.bonded.Atom} object.
     * @param angle1     a double.
     * @param ic         a {@link ffx.potential.bonded.Atom} object.
     * @param angle2     a double.
     * @param chiral     a int.
     * @param lookUp     a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public static Atom buildHeavy(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
                                  Atom ic, double angle2, int chiral, int lookUp, ForceField forceField, ArrayList<Bond> bondList) {
        AtomType atomType = findAtomType(lookUp, forceField);
        return buildHeavyAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType, forceField, bondList);
    }

    /**
     * <p>buildHydrogen.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.MSGroup} object.
     * @param atomName   a {@link java.lang.String} object.
     * @param ia         a {@link ffx.potential.bonded.Atom} object.
     * @param bond       a double.
     * @param ib         a {@link ffx.potential.bonded.Atom} object.
     * @param angle1     a double.
     * @param ic         a {@link ffx.potential.bonded.Atom} object.
     * @param angle2     a double.
     * @param chiral     a int.
     * @param lookUp     a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public static Atom buildHydrogen(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
                                     Atom ic, double angle2, int chiral, int lookUp, ForceField forceField, ArrayList<Bond> bondList) {
        AtomType atomType = findAtomType(lookUp, forceField);
        return buildHydrogenAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType, forceField, bondList);
    }

    /**
     * <p>buildHydrogenAtom.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.MSGroup} object.
     * @param atomName   a {@link java.lang.String} object.
     * @param ia         a {@link ffx.potential.bonded.Atom} object.
     * @param bond       a double.
     * @param ib         a {@link ffx.potential.bonded.Atom} object.
     * @param angle1     a double.
     * @param ic         a {@link ffx.potential.bonded.Atom} object.
     * @param angle2     a double.
     * @param chiral     a int.
     * @param atomType   a {@link ffx.potential.parameters.AtomType} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public static Atom buildHydrogenAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
                                         Atom ic, double angle2, int chiral, AtomType atomType, ForceField forceField, ArrayList<Bond> bondList) {
        if (atomType == null) {
            return null;
        }
        Atom atom = (Atom) residue.getAtomNode(atomName);
        // It may be a Deuterium
        if (atom == null) {
            String dAtomName = atomName.replaceFirst("H", "D");
            atom = (Atom) residue.getAtomNode(dAtomName);
        }

        // Basic checking for unspecified H atoms attached to waters.
        if (residue instanceof Molecule && atom == null) {
            Molecule molec = (Molecule) residue;
            String molName = molec.getName().toUpperCase();
            if (molName.startsWith("HOH") || molName.startsWith("WAT") || molName.startsWith("TIP3")) {
                atom = (Atom) molec.getAtomNode("H");
                if (atom == null) {
                    atom = (Atom) molec.getAtomNode("D");
                }
                if (atom != null) {
                    atom.setName(atomName);
                }
            }
        }
        if (atom == null) {
            String resName = ia.getResidueName();
            int resSeq = ia.getResidueNumber();
            Character chainID = ia.getChainID();
            Character altLoc = ia.getAltLoc();
            String segID = ia.getSegID();
            double occupancy = ia.getOccupancy();
            double tempFactor = ia.getTempFactor();
            atom = new Atom(0, atomName, altLoc, new double[3], resName, resSeq, chainID,
                    occupancy, tempFactor, segID, true);
            residue.addMSNode(atom);
            intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
        }
        atom.setAtomType(atomType);
        buildBond(ia, atom, forceField, bondList);
        return atom;
    }

    /**
     * This exception is thrown when an atom type could not be assigned.
     */
    public static class MissingAtomTypeException extends Exception {

        public final Residue residue;
        public final Atom atom;

        public MissingAtomTypeException(Residue residue, Atom atom) {
            this.residue = residue;
            this.atom = atom;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(format(" Atom %s", atom.toString()));
            if (residue != null) {
                sb.append(format("\n of residue %c-%s", residue.getChainID(), residue.toString()));
            }
            sb.append("\n could not be assigned an atom type.\n");
            return sb.toString();
        }
    }

    /**
     * This exception is thrown when a heavy atom is not found.
     */
    public static class MissingHeavyAtomException extends Exception {

        public final String atomName;
        public final AtomType atomType;
        final Atom bondedTo;

        public MissingHeavyAtomException(String atomName, AtomType atomType, Atom bondedTo) {
            this.atomName = atomName;
            this.atomType = atomType;
            this.bondedTo = bondedTo;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            if (atomType != null) {
                sb.append(format("\n An atom of type\n %s\n", atomType.toString()));
            } else {
                sb.append(format("\n Atom %s", atomName));
            }
            sb.append(" was not found");
            if (bondedTo != null) {
                sb.append(format(" bonded to atom %s ", bondedTo));
            }
            sb.append(".\n");
            return sb.toString();
        }
    }
}
