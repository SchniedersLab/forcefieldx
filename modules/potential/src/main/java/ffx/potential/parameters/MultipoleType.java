/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential.parameters;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.norm;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;

/**
 * The MultipoleType class defines a multipole in its local frame.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class MultipoleType extends BaseType implements Comparator<String> {

    private static final Logger logger = Logger.getLogger(MultipoleType.class.getName());

    /**
     * The local multipole frame is defined by the Z-then-X or Bisector
     * convention.
     */
    public enum MultipoleFrameDefinition {

        ZONLY, ZTHENX, BISECTOR, ZTHENBISECTOR, TRISECTOR
    }
    /**
     * Conversion from electron-Angstroms to Debyes
     */
    public static final double DEBYE = 4.80321;
    /**
     * Conversion from electron-Angstroms^2 to Buckinghams
     */
    public static final double BUCKINGHAM = DEBYE * DEBYE;
    /**
     * Conversion from Bohr to Angstroms
     */
    public static final double BOHR = 0.52917720859;
    /**
     * Conversion from Bohr^2 to Angstroms^2
     */
    public static final double BOHR2 = BOHR * BOHR;
    /**
     * Conversion from electron**2/Ang to kcal/mole.
     */
    public static final double ELECTRIC = 332.063709;
    /**
     * Partial atomic charge (e).
     */
    public final double charge;
    /**
     * Atomic dipole. 1 x 3 (e Angstroms).
     */
    public final double dipole[];
    /**
     * Atomic quadrupole. 3 x 3 (e Angstroms^2).
     */
    public final double quadrupole[][];
    /**
     * Local frame definition method.
     */
    public final MultipoleFrameDefinition frameDefinition;
    /**
     * Atom types that define the local frame of this multipole.
     */
    public final int[] frameAtomTypes;
    /**
     * Charge, dipole, and quadrupole packed into tensor notation: c, dx, dy,
     * dz, qxx, qyy, qzz, qxy, qxz, qyz
     */
    private final double[] multipole;

    public static final double[] zeroM = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public static final double[] zeroD = new double[]{0.0, 0.0, 0.0};

    /**
     * Multipole Constructor. Conversion to electron Angstroms should be
     * requested only when reading multipole values from the force field file.
     */
    public MultipoleType(double[] multipole, int[] frameAtomTypes,
            MultipoleFrameDefinition frameDefinition, boolean convertFromBohr) {
        super(ForceField.ForceFieldType.MULTIPOLE, frameAtomTypes);
        this.multipole = (convertFromBohr) ? bohrToElectronAngstroms(multipole) : multipole;
        this.frameAtomTypes = frameAtomTypes;
        this.frameDefinition = frameDefinition;
        charge = multipole[t000];
        dipole = unpackDipole(multipole);
        quadrupole = unpackQuad(multipole);
        checkMultipole();
    }

    public MultipoleType(double charge, double[] dipole, double[][] quadrupole,
            int[] frameAtomTypes, MultipoleFrameDefinition frameDefinition, boolean convertFromBohr) {
        super(ForceField.ForceFieldType.MULTIPOLE, frameAtomTypes);
        this.charge = charge;
        if (convertFromBohr) {
            this.multipole = bohrToElectronAngstroms(pack(charge, dipole, quadrupole));
            this.dipole = unpackDipole(multipole);
            this.quadrupole = unpackQuad(multipole);
        } else {
            this.multipole = pack(charge, dipole, quadrupole);
            this.dipole = dipole;
            this.quadrupole = quadrupole;
        }
        this.frameAtomTypes = frameAtomTypes;
        this.frameDefinition = frameDefinition;
        checkMultipole();
    }

    /**
     * @return An uneditable copy of this type's multipole. To make changes, use
     * getMultipoleReference().
     */
    public double[] getMultipole() {
        return new double[]{
            multipole[t000],
            multipole[t100], multipole[t010], multipole[t001],
            multipole[t200], multipole[t020], multipole[t002], multipole[t110], multipole[t101], multipole[t011]};
    }

    /**
     * Exposes multipole array for editing; makes this process explicit.
     */
    public double[] getMultipoleReference() {
        return multipole;
    }

    /**
     * @return An uneditable copy of this type's charge. To make changes, use
     * getMultipoleReference().
     */
    public double getCharge() {
        return multipole[t000];
    }

    /**
     * @return An uneditable copy of this type's dipole. To make changes, use
     * getMultipoleReference().
     */
    public double[] getDipole() {
        return new double[]{multipole[t100], multipole[t010], multipole[t001]};
    }

    /**
     * @return An uneditable copy of this type's quadrupole. To make changes,
     * use getMultipoleReference().
     */
    public double[][] getQuadrupole() {
        return new double[][]{
            {multipole[t200], multipole[t110], multipole[t101]},
            {multipole[t110], multipole[t020], multipole[t011]},
            {multipole[t101], multipole[t011], multipole[t002]}};
    }

    /**
     * <p>
     * incrementType</p>
     *
     * @param increment a int.
     */
    public void incrementType(int increment) {
        for (int i = 0; i < frameAtomTypes.length; i++) {
            // Frame atom types of 0 are unchanged.
            if (frameAtomTypes[i] > 0) {
                frameAtomTypes[i] += increment;
            } else if (frameAtomTypes[i] < 0) {
                frameAtomTypes[i] -= increment;
            }
        }
        setKey(frameAtomTypes);
    }

    /**
     * Remap new atom types to known internal ones.
     *
     * @param typeMap a lookup between new atom types and known atom types.
     *
     * @return
     */
    public MultipoleType patchTypes(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        int len = frameAtomTypes.length;
        /**
         * Look for a MultipoleType that contain a mapped atom class.
         */
        for (AtomType newType : typeMap.keySet()) {
            for (int i = 0; i < len; i++) {
                if (frameAtomTypes[i] == newType.type || frameAtomTypes[i] == 0) {
                    count++;
                }
            }
        }
        /**
         * If found, create a new MultipoleType that bridges to known classes.
         */
        if (count > 0 && count < len) {
            int newFrame[] = Arrays.copyOf(frameAtomTypes, len);
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < len; i++) {
                    if (frameAtomTypes[i] == newType.type) {
                        AtomType knownType = typeMap.get(newType);
                        newFrame[i] = knownType.type;
                    }
                }
            }
            return new MultipoleType(multipole, newFrame, frameDefinition, false);
        }
        return null;
    }

    private void checkMultipole() {
        double[][] quadrupole = unpackQuad(multipole);
        // Check symmetry.
        if (Math.abs(quadrupole[0][1] - quadrupole[1][0]) > 1.0e-6) {
            logger.warning("Multipole component Qxy != Qyx");
            logger.info(this.toString());
        }
        if (Math.abs(quadrupole[0][2] - quadrupole[2][0]) > 1.0e-6) {
            logger.warning("Multipole component Qxz != Qzx");
            logger.info(this.toString());
        }
        if (Math.abs(quadrupole[1][2] - quadrupole[2][1]) > 1.0e-6) {
            logger.warning("Multipole component Qyz != Qzy");
            logger.info(this.toString());
        }
        // Warn if the multipole is not traceless.
        if (Math.abs(quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2]) > 1.0e-5) {
            logger.log(Level.WARNING, format("Multipole is not traceless: %12.8f",
                    Math.abs(quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2])));
            logger.info(this.toString());
        }
    }

    public static boolean assignMultipole(Atom atom, ForceField forceField,
            double multipole[], int i, int axisAtom[][], MultipoleFrameDefinition frame[]) {
        MultipoleType type = multipoleTypeFactory(atom, forceField);
        if (type == null) {
            return false;
        }
        System.arraycopy(type.getMultipole(), 0, multipole, 0, 10);
        axisAtom[i] = atom.getAxisAtomIndices();
        frame[i] = atom.getMultipoleType().frameDefinition;
        return true;
    }

    public static MultipoleType multipoleTypeFactory(Atom atom, ForceField forceField) {
        AtomType atomType = atom.getAtomType();
        if (atomType == null) {
            String message = " Multipoles can only be assigned to atoms that have been typed.";
            logger.severe(message);
            return null;
        }

        PolarizeType polarizeType = forceField.getPolarizeType(atomType.getKey());
        if (polarizeType != null) {
            atom.setPolarizeType(polarizeType);
        } else {
            String message = " No polarization type was found for " + atom.toString();
            logger.fine(message);
            double polarizability = 0.0;
            double thole = 0.0;
            int polarizationGroup[] = null;
            polarizeType = new PolarizeType(atomType.type,
                    polarizability, thole, polarizationGroup);
            forceField.addForceFieldType(polarizeType);
            atom.setPolarizeType(polarizeType);
        }

        String key;
        // No reference atoms.
        key = atomType.getKey() + " 0 0";
        MultipoleType multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
            atom.setMultipoleType(multipoleType);
            Atom axisAtom = null;
            atom.setAxisAtoms(axisAtom);
            return multipoleType;
        }

        // No bonds.
        List<Bond> bonds = atom.getBonds();
        if (bonds == null || bonds.size() < 1) {
            String message = "Multipoles can only be assigned after bonded relationships are defined.\n";
            logger.severe(message);
        }

        // 1 reference atom.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            key = atomType.getKey() + " " + atom2.getAtomType().getKey() + " 0";
            multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
                atom.setMultipoleType(multipoleType);
                atom.setAxisAtoms(atom2);
                return multipoleType;
            }
        }

        // 2 reference atoms.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                key = atomType.getKey() + " " + key2 + " " + key3;
                multipoleType = forceField.getMultipoleType(key);
                if (multipoleType != null) {
                    atom.setMultipoleType(multipoleType);
                    atom.setAxisAtoms(atom2, atom3);
                    return multipoleType;
                }
            }
        }

        /**
         * 3 reference atoms.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                for (Bond b3 : bonds) {
                    if (b == b3 || b2 == b3) {
                        continue;
                    }
                    Atom atom4 = b3.get1_2(atom);
                    String key4 = atom4.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        atom.setMultipoleType(multipoleType);
                        atom.setAxisAtoms(atom2, atom3, atom4);
                        return multipoleType;
                    }
                }
                List<Angle> angles = atom.getAngles();
                for (Angle angle : angles) {
                    Atom atom4 = angle.get1_3(atom);
                    if (atom4 != null) {
                        String key4 = atom4.getAtomType().getKey();
                        key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                        multipoleType = forceField.getMultipoleType(key);
                        if (multipoleType != null) {
                            atom.setMultipoleType(multipoleType);
                            atom.setAxisAtoms(atom2, atom3, atom4);
                            return multipoleType;
                        }
                    }
                }
            }
        }

        /**
         * Revert to a 2 reference atom definition that may include a 1-3 site.
         * For example a hydrogen on water.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            List<Angle> angles = atom.getAngles();
            for (Angle angle : angles) {
                Atom atom3 = angle.get1_3(atom);
                if (atom3 != null) {
                    String key3 = atom3.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        atom.setMultipoleType(multipoleType);
                        atom.setAxisAtoms(atom2, atom3);
                        return multipoleType;
                    }
                    for (Angle angle2 : angles) {
                        Atom atom4 = angle2.get1_3(atom);
                        if (atom4 != null && atom4 != atom3) {
                            String key4 = atom4.getAtomType().getKey();
                            key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                            multipoleType = forceField.getMultipoleType(key);
                            if (multipoleType != null) {
                                atom.setMultipoleType(multipoleType);
                                atom.setAxisAtoms(atom2, atom3, atom4);
                                return multipoleType;
                            }
                        }
                    }
                }
            }
        }
        return null;
    }

    private static double[] bohrToElectronAngstroms(double[] multipole) {
        return new double[]{
            multipole[t000],
            multipole[t100] *= BOHR,
            multipole[t010] *= BOHR,
            multipole[t001] *= BOHR,
            multipole[t200] *= BOHR2,
            multipole[t020] *= BOHR2,
            multipole[t002] *= BOHR2,
            multipole[t110] *= BOHR2,
            multipole[t101] *= BOHR2,
            multipole[t011] *= BOHR2};
    }

    /**
     * @return Whether this multipole underwent chiral inversion.
     */
    public static boolean checkMultipoleChirality(double[] multipole, MultipoleFrameDefinition frame,
            double localOrigin[], double frameCoords[][]) {
        if (frame != MultipoleFrameDefinition.ZTHENX) {
            return false;
        }
        double zAxis[] = new double[3];
        double xAxis[] = new double[3];
        double yAxis[] = new double[3];
        double yMinOrigin[] = new double[3];
        zAxis[0] = frameCoords[0][0];
        zAxis[1] = frameCoords[0][1];
        zAxis[2] = frameCoords[0][2];
        xAxis[0] = frameCoords[1][0];
        xAxis[1] = frameCoords[1][1];
        xAxis[2] = frameCoords[1][2];
        yAxis[0] = frameCoords[2][0];
        yAxis[1] = frameCoords[2][1];
        yAxis[2] = frameCoords[2][2];
        diff(localOrigin, yAxis, yMinOrigin);
        diff(zAxis, yAxis, zAxis);
        diff(xAxis, yAxis, xAxis);
        double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
        double c2 = xAxis[1] * yMinOrigin[2] - xAxis[2] * yMinOrigin[1];
        double c3 = yMinOrigin[1] * zAxis[2] - yMinOrigin[2] * zAxis[1];
        double vol = yMinOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
        return (vol < 0.0);
    }

    public static void invertMultipoleChirality(double[] mpole) {
        mpole[t010] = -mpole[t010];
        mpole[t110] = -mpole[t110];
        mpole[t011] = -mpole[t011];
    }

    /**
     * Return the rotation matrix for the local to lab frame.
     * @param frame the multipole frame definition
     * @param localOrigin the local origin of the frame
     * @param frameCoords the coordinates of the frame atoms
     * @return the rotation matrix
     */
    public static double[][] getRotationMatrix(MultipoleFrameDefinition frame,
            double localOrigin[], double frameCoords[][]) {
        double[][] rotmat = new double[3][3];
        double zAxis[] = new double[3];
        double xAxis[] = new double[3];
        switch (frame) {
            case BISECTOR:
                zAxis[0] = frameCoords[0][0];
                zAxis[1] = frameCoords[0][1];
                zAxis[2] = frameCoords[0][2];
                xAxis[0] = frameCoords[1][0];
                xAxis[1] = frameCoords[1][1];
                xAxis[2] = frameCoords[1][2];
                diff(zAxis, localOrigin, zAxis);
                norm(zAxis, zAxis);
                diff(xAxis, localOrigin, xAxis);
                norm(xAxis, xAxis);
                sum(xAxis, zAxis, zAxis);
                norm(zAxis, zAxis);
                rotmat[0][2] = zAxis[0];
                rotmat[1][2] = zAxis[1];
                rotmat[2][2] = zAxis[2];
                double dot = dot(xAxis, zAxis);
                scalar(zAxis, dot, zAxis);
                diff(xAxis, zAxis, xAxis);
                norm(xAxis, xAxis);
                break;
            case ZTHENBISECTOR:
                double yAxis[] = new double[3];
                zAxis[0] = frameCoords[0][0];
                zAxis[1] = frameCoords[0][1];
                zAxis[2] = frameCoords[0][2];
                xAxis[0] = frameCoords[1][0];
                xAxis[1] = frameCoords[1][1];
                xAxis[2] = frameCoords[1][2];
                yAxis[0] = frameCoords[2][0];
                yAxis[1] = frameCoords[2][1];
                yAxis[2] = frameCoords[2][2];
                diff(zAxis, localOrigin, zAxis);
                norm(zAxis, zAxis);
                rotmat[0][2] = zAxis[0];
                rotmat[1][2] = zAxis[1];
                rotmat[2][2] = zAxis[2];
                diff(xAxis, localOrigin, xAxis);
                norm(xAxis, xAxis);
                diff(yAxis, localOrigin, yAxis);
                norm(yAxis, yAxis);
                sum(xAxis, yAxis, xAxis);
                norm(xAxis, xAxis);
                dot = dot(xAxis, zAxis);
                scalar(zAxis, dot, zAxis);
                diff(xAxis, zAxis, xAxis);
                norm(xAxis, xAxis);
                break;
            case ZONLY:
                zAxis[0] = frameCoords[0][0];
                zAxis[1] = frameCoords[0][1];
                zAxis[2] = frameCoords[0][2];
                diff(zAxis, localOrigin, zAxis);
                norm(zAxis, zAxis);
                rotmat[0][2] = zAxis[0];
                rotmat[1][2] = zAxis[1];
                rotmat[2][2] = zAxis[2];
                xAxis[0] = random();
                xAxis[1] = random();
                xAxis[2] = random();
                dot = dot(xAxis, zAxis);
                scalar(zAxis, dot, zAxis);
                diff(xAxis, zAxis, xAxis);
                norm(xAxis, xAxis);
                break;
            case ZTHENX:
            default:
                zAxis[0] = frameCoords[0][0];
                zAxis[1] = frameCoords[0][1];
                zAxis[2] = frameCoords[0][2];
                xAxis[0] = frameCoords[1][0];
                xAxis[1] = frameCoords[1][1];
                xAxis[2] = frameCoords[1][2];
                diff(zAxis, localOrigin, zAxis);
                norm(zAxis, zAxis);
                rotmat[0][2] = zAxis[0];
                rotmat[1][2] = zAxis[1];
                rotmat[2][2] = zAxis[2];
                diff(xAxis, localOrigin, xAxis);
                dot = dot(xAxis, zAxis);
                scalar(zAxis, dot, zAxis);
                diff(xAxis, zAxis, xAxis);
                norm(xAxis, xAxis);
        }
        // Set the X elements.
        rotmat[0][0] = xAxis[0];
        rotmat[1][0] = xAxis[1];
        rotmat[2][0] = xAxis[2];
        // Set the Y elements.
        rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0] * rotmat[2][2];
        rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0] * rotmat[0][2];
        rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0] * rotmat[1][2];
        return rotmat;
    }

    public static void rotateDipole(double rotmat[][], double dipole[], double rotatedDipole[]) {
        for (int i = 0; i < 3; i++) {
            double[] rotmati = rotmat[i];
            for (int j = 0; j < 3; j++) {
                rotatedDipole[i] += rotmati[j] * dipole[j];
            }
        }
    }

    public static void rotateMultipole(double rotmat[][], double dipole[],
            double quadrupole[][], double rotatedDipole[], double rotatedQuadrupole[][]) {
        for (int i = 0; i < 3; i++) {
            double[] rotmati = rotmat[i];
            double[] quadrupolei = rotatedQuadrupole[i];
            for (int j = 0; j < 3; j++) {
                double[] rotmatj = rotmat[j];
                rotatedDipole[i] += rotmati[j] * dipole[j];
                if (j < i) {
                    quadrupolei[j] = rotatedQuadrupole[j][i];
                } else {
                    for (int k = 0; k < 3; k++) {
                        double[] localQuadrupolek = quadrupole[k];
                        quadrupolei[j] += rotmati[k]
                                * (rotmatj[0] * localQuadrupolek[0]
                                + rotmatj[1] * localQuadrupolek[1]
                                + rotmatj[2] * localQuadrupolek[2]);
                    }
                }
            }
        }
    }

    /**
     * Pack charge, dipole, quad into 1d tensor array form.
     */
    public static double[] pack(double charge, double[] dipl, double[][] quad) {
        return new double[]{
            charge,
            dipl[0], dipl[1], dipl[2],
            quad[0][0], quad[1][1], quad[2][2], quad[0][1], quad[0][2], quad[1][2]};
    }

    /**
     * Unpack dipole from 1d tensor-form multipole.
     */
    private static double[] unpackDipole(double[] mpole) {
        return new double[]{mpole[t100], mpole[t010], mpole[t001]};
    }

    /**
     * Unpack quadrupole from 1d tensor-form multipole.
     */
    private static double[][] unpackQuad(double[] mpole) {
        return new double[][]{
            {mpole[t200], mpole[t110], mpole[t101]},
            {mpole[t110], mpole[t020], mpole[t011]},
            {mpole[t101], mpole[t011], mpole[t002]}};
    }

    public static double[] scale(MultipoleType type, double[] cdtScales) {
        return scale(type.getMultipole(), cdtScales);
    }

    public static double[] scale(double[] multipole, double[] cdtScales) {
        double chargeScale = cdtScales[0];
        double dipoleScale = cdtScales[1];
        double quadScale = cdtScales[2];
        return new double[]{
            multipole[t000] * chargeScale,
            multipole[t100] * dipoleScale,
            multipole[t010] * dipoleScale,
            multipole[t001] * dipoleScale,
            multipole[t200] * quadScale,
            multipole[t020] * quadScale,
            multipole[t002] * quadScale,
            multipole[t110] * quadScale,
            multipole[t101] * quadScale,
            multipole[t011] * quadScale};
    }

    /**
     * Nicely formatted multipole string. Dipole and qaudrupole are in
     * electron-Bohrs and electron-Bohrs^2, respectively.
     *
     * @return String
     */
    public String toBohrString() {
        StringBuilder multipoleBuffer = new StringBuilder("multipole");
        if (frameDefinition == MultipoleFrameDefinition.BISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
        } else if (frameDefinition == MultipoleFrameDefinition.ZTHENBISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[3]));
        } else if (frameDefinition == MultipoleFrameDefinition.TRISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[3]));
        } else {
            for (int i : frameAtomTypes) {
                multipoleBuffer.append(format("  %5d", i));
            }
        }
        if (frameAtomTypes.length == 3) {
            multipoleBuffer.append("       ");
        }
        multipoleBuffer.append(format("  % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f % 7.5f \\\n"
                + "%11$s % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f % 7.5f",
                multipole[t000],
                multipole[t100] / BOHR, multipole[t010] / BOHR, multipole[t001] / BOHR,
                multipole[t200] / BOHR2,
                multipole[t110] / BOHR2, multipole[t020] / BOHR2,
                multipole[t101] / BOHR2, multipole[t011] / BOHR2, multipole[t002] / BOHR2,
                "                                      "));
        return multipoleBuffer.toString();
    }

    /**
     * Nicely formatted multipole string. Dipole and qaudrupole are in
     * electron-Bohrs and electron-Bohrs^2, respectively.
     *
     * @return String
     */
    public String toCompactBohrString() {
        StringBuilder multipoleBuffer = new StringBuilder("mpol ");
        if (frameDefinition == MultipoleFrameDefinition.BISECTOR) {
            multipoleBuffer.append(format("(%3d,%3d,%3d): ",
                    frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2]));
        } else if (frameDefinition == MultipoleFrameDefinition.ZTHENBISECTOR) {
            multipoleBuffer.append(format("(%3d,%3d,%3d,%3d): ",
                    frameAtomTypes[0], frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
        } else if (frameDefinition == MultipoleFrameDefinition.TRISECTOR) {
            multipoleBuffer.append(format("(%3d,%3d,%3d,%3d): ",
                    frameAtomTypes[0], -frameAtomTypes[1], -frameAtomTypes[2], -frameAtomTypes[3]));
        } else {
            multipoleBuffer.append(format("("));
            for (int i : frameAtomTypes) {
                multipoleBuffer.append(format("%3d,", i));
            }
            int comma = multipoleBuffer.lastIndexOf(",");
            multipoleBuffer.replace(comma, comma + 1, "");
            multipoleBuffer.append(format("): "));
        }
        multipoleBuffer.append(format("[%6.3f / %6.3f %6.3f %6.3f / %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]",
                charge,
                dipole[0] / BOHR, dipole[1] / BOHR, dipole[2] / BOHR,
                quadrupole[0][0] / BOHR2,
                quadrupole[1][0] / BOHR2, quadrupole[1][1] / BOHR2,
                quadrupole[2][0] / BOHR2, quadrupole[2][1] / BOHR2, quadrupole[2][2] / BOHR2));
        return multipoleBuffer.toString();
    }

    /**
     * Nicely formatted multipole string. Dipole and qaudrupole are in units of
     * Debye and Buckinghams, respectively.
     *
     * @return String
     */
    public String toDebyeString() {
        StringBuilder multipoleBuffer = new StringBuilder("multipole");
        if (frameDefinition == MultipoleFrameDefinition.BISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
        } else if (frameDefinition == MultipoleFrameDefinition.ZTHENBISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[3]));
        } else if (frameDefinition == MultipoleFrameDefinition.TRISECTOR) {
            multipoleBuffer.append(format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(format("  %5d", -frameAtomTypes[3]));
        } else {
            for (int i : frameAtomTypes) {
                multipoleBuffer.append(format("  %5d", i));
            }
        }
        if (frameAtomTypes.length == 3) {
            multipoleBuffer.append("       ");
        }
        multipoleBuffer.append(format(
                "  %7.5f \\\n"
                + "%11s %7.5f %7.5f %7.5f \\\n"
                + "%11s %7.5f \\\n"
                + "%11s %7.5f %7.5f \\\n"
                + "%11s %7.5f %7.5f %7.5f",
                multipole[t000],
                "", multipole[t100] * DEBYE, multipole[t010] * DEBYE, multipole[t001] * DEBYE,
                "", multipole[t200] * BUCKINGHAM,
                "", multipole[t110] * BUCKINGHAM, multipole[t020] * BUCKINGHAM,
                "", multipole[t101] * BUCKINGHAM, multipole[t011] * BUCKINGHAM, multipole[t002] * BUCKINGHAM));
        return multipoleBuffer.toString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return toBohrString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {
        String keys1[] = s1.split(" ");
        String keys2[] = s2.split(" ");

        int len = keys1.length;
        if (keys1.length > keys2.length) {
            len = keys2.length;
        }
        int c1[] = new int[len];
        int c2[] = new int[len];
        for (int i = 0; i < len; i++) {
            c1[i] = abs(Integer.parseInt(keys1[i]));
            c2[i] = abs(Integer.parseInt(keys2[i]));
            if (c1[i] < c2[i]) {
                return -1;
            } else if (c1[i] > c2[i]) {
                return 1;
            }
        }

        if (keys1.length < keys2.length) {
            return -1;
        } else if (keys1.length > keys2.length) {
            return 1;
        }

        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null || !(other instanceof MultipoleType)) {
            return false;
        }
        MultipoleType multipoleType = (MultipoleType) other;
        int c[] = multipoleType.frameAtomTypes;
        if (c.length != frameAtomTypes.length) {
            return false;
        }
        for (int i = 0; i < c.length; i++) {
            if (c[i] != this.frameAtomTypes[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Arrays.hashCode(frameAtomTypes);
        return hash;
    }

    /**
     * Average two MultipoleType instances. The atom types that define the frame
     * of the new type must be supplied.
     */
    public static MultipoleType averageTypes(MultipoleType multipoleType1, MultipoleType multipoleType2, int[] multipoleFrameTypes) {
        if (multipoleType1 == null || multipoleType2 == null || multipoleFrameTypes != null) {
            return null;
        }
        if (multipoleType1.frameDefinition != multipoleType2.frameDefinition) {
            return null;
        }
        MultipoleType[] types = {multipoleType1, multipoleType2};
        double[] weights = {0.5, 0.5};
        double[] averagedMultipole = weightMultipole(types, weights);
        if (averagedMultipole == null) {
            return null;
        }
        return new MultipoleType(averagedMultipole, multipoleFrameTypes, multipoleType1.frameDefinition, false);
    }

    public static MultipoleType weightMultipoleTypes(MultipoleType[] types, double[] weights, int[] frameAtomTypes) {
        double[] weightedMultipole = weightMultipole(types, weights);
        if (weightedMultipole == null) {
            return null;
        }
        return new MultipoleType(weightedMultipole, frameAtomTypes, types[0].frameDefinition, false);
    }

    /**
     * Create a new multipole representing a weighted average.
     */
    public static double[] weightMultipole(MultipoleType[] types, double[] weights) {
        if (types == null || weights == null || types.length != weights.length) {
            throw new IllegalArgumentException();
        }
        if (Arrays.asList(types).contains(null)) {
            // Multipoles have not yet been assigned.
            return null;
        }
        for (MultipoleType type : types) {
            if (type.frameDefinition != types[0].frameDefinition) {
                logger.warning(String.format("Multipole frame definition mismatch during weighting:\n\t%s->%s,\n\t%s->%s",
                        types[0].toString(), types[0].frameDefinition.toString(),
                        type.toString(), type.frameDefinition.toString()));
                throw new IllegalArgumentException();
            }
        }
        double[] weightedMultipole = new double[10];
        fill(weightedMultipole, 0.0);
        for (int idx = 0; idx < types.length; idx++) {
            double[] multipole = types[idx].getMultipole();
            for (int comp = 0; comp < 10; comp++) {
                weightedMultipole[comp] += weights[idx] * multipole[comp];
            }
        }
        return weightedMultipole;
    }

    /**
     * Indices into a 1D tensor array based on compressed tensor notation. This
     * makes multipole code much easier to read.
     */
    public static final int t000 = 0, chrg = t000;
    /**
     * Constant <code>t100=1</code>
     */
    public static final int t100 = 1, diplx = t100;
    /**
     * Constant <code>t010=2</code>
     */
    public static final int t010 = 2, diply = t010;
    /**
     * Constant <code>t001=3</code>
     */
    public static final int t001 = 3, diplz = t100;
    /**
     * Constant <code>t200=4</code>
     */
    public static final int t200 = 4, quadxx = t200;
    /**
     * Constant <code>t020=5</code>
     */
    public static final int t020 = 5, quadyy = t020;
    /**
     * Constant <code>t002=6</code>
     */
    public static final int t002 = 6, quadzz = t002;
    /**
     * Constant <code>t110=7</code>
     */
    public static final int t110 = 7, quadxy = t110;
    /**
     * Constant <code>t101=8</code>
     */
    public static final int t101 = 8, quadxz = t101;
    /**
     * Constant <code>t011=9</code>
     */
    public static final int t011 = 9, quadyz = t011;
    /**
     * Constant <code>t300=10</code>
     */
    public static final int t300 = 10;
    /**
     * Constant <code>t030=11</code>
     */
    public static final int t030 = 11;
    /**
     * Constant <code>t003=12</code>
     */
    public static final int t003 = 12;
    /**
     * Constant <code>t210=13</code>
     */
    public static final int t210 = 13;
    /**
     * Constant <code>t201=14</code>
     */
    public static final int t201 = 14;
    /**
     * Constant <code>t120=15</code>
     */
    public static final int t120 = 15;
    /**
     * Constant <code>t021=16</code>
     */
    public static final int t021 = 16;
    /**
     * Constant <code>t102=17</code>
     */
    public static final int t102 = 17;
    /**
     * Constant <code>t012=18</code>
     */
    public static final int t012 = 18;
    /**
     * Constant <code>t111=19</code>
     */
    public static final int t111 = 19;
}
