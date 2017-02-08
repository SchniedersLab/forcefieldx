/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
import java.util.logging.Logger;

import static java.lang.String.format;

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
     * Charge, dipole, and quadrupole packed into 1d.
     */
    public final double[] packedMultipole;

    /**
     * Multipole Constructor.
     *
     * @param charge double
     * @param dipole double[]
     * @param quadrupole double[]
     * @param multipoleFrameTypes int[]
     * @param frameDefinition a
     * {@link ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition}
     * object.
     */
    public MultipoleType(double charge, double dipole[], double quadrupole[][],
            int[] multipoleFrameTypes, MultipoleFrameDefinition frameDefinition,
            boolean convertFromBohr) {
        super(ForceField.ForceFieldType.MULTIPOLE, multipoleFrameTypes);
        this.charge = charge;
        this.dipole = dipole;
        this.quadrupole = quadrupole;
        this.frameAtomTypes = multipoleFrameTypes;
        this.frameDefinition = frameDefinition;
        if (convertFromBohr) {
            convertBohrToElectronAngstroms();
        }
        checkMultipole();
        packedMultipole = new double[]{charge, dipole[0], dipole[1], dipole[2],
                quadrupole[0][0], quadrupole[1][1], quadrupole[2][2],
                quadrupole[0][1], quadrupole[0][2], quadrupole[1][2]};
    }
    
    /**
     * This assumes the dipole and quadrupole are in
     * units of Bohr, and are converted to electron-Angstroms and
     * electron-Angstroms^2, respectively, before the constructor returns.
     */
    public MultipoleType(double charge, double dipole[], double quadrupole[][],
            int[] multipoleFrameTypes, MultipoleFrameDefinition frameDefinition) {
        this(charge, dipole, quadrupole, multipoleFrameTypes, frameDefinition, true);
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
            return new MultipoleType(charge, dipole, quadrupole, newFrame, frameDefinition);
        }
        return null;
    }

    private void checkMultipole() {
        // Check symmetry.
        double check = Math.abs(quadrupole[0][1] - quadrupole[1][0]);
        if (check > 1.0e-6) {
            logger.warning("Multipole component Qxy != Qyx");
            print();
        }
        check = Math.abs(quadrupole[0][2] - quadrupole[2][0]);
        if (check > 1.0e-6) {
            logger.warning("Multipole component Qxz != Qzx");
            print();
        }
        check = Math.abs(quadrupole[1][2] - quadrupole[2][1]);
        if (check > 1.0e-6) {
            logger.warning("Multipole component Qyz != Qzy");
            print();
        }
        // Warn if the multipole is not traceless.
        double sum = quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2];
        if (Math.abs(sum) > 1.0e-5) {
            String message = format("Multipole is not traceless: %7.5f \n%s",
                    sum, toBohrString());
            logger.warning(message);
        }
    }
    
    public static boolean assignMultipole(Atom atom, ForceField forceField,
            double multipole[], int i, int axisAtom[][], MultipoleFrameDefinition frame[]) {
        AtomType atomType = atom.getAtomType();
        if (atomType == null) {
            String message = " Multipoles can only be assigned to atoms that have been typed.";
            logger.severe(message);
            return false;
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
            atom.setMultipoleType(multipoleType, null);
            multipole[t000] = multipoleType.charge;
            multipole[t100] = multipoleType.dipole[0];
            multipole[t010] = multipoleType.dipole[1];
            multipole[t001] = multipoleType.dipole[2];
            multipole[t200] = multipoleType.quadrupole[0][0];
            multipole[t020] = multipoleType.quadrupole[1][1];
            multipole[t002] = multipoleType.quadrupole[2][2];
            multipole[t110] = multipoleType.quadrupole[0][1];
            multipole[t101] = multipoleType.quadrupole[0][2];
            multipole[t011] = multipoleType.quadrupole[1][2];
            axisAtom[i] = null;
            frame[i] = multipoleType.frameDefinition;
            return true;
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
            multipoleType = multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
                int multipoleReferenceAtoms[] = new int[1];
                multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                atom.setMultipoleType(multipoleType, null);
                multipole[t000] = multipoleType.charge;
                multipole[t100] = multipoleType.dipole[0];
                multipole[t010] = multipoleType.dipole[1];
                multipole[t001] = multipoleType.dipole[2];
                multipole[t200] = multipoleType.quadrupole[0][0];
                multipole[t020] = multipoleType.quadrupole[1][1];
                multipole[t002] = multipoleType.quadrupole[2][2];
                multipole[t110] = multipoleType.quadrupole[0][1];
                multipole[t101] = multipoleType.quadrupole[0][2];
                multipole[t011] = multipoleType.quadrupole[1][2];
                axisAtom[i] = multipoleReferenceAtoms;
                frame[i] = multipoleType.frameDefinition;
                return true;
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
                    int multipoleReferenceAtoms[] = new int[2];
                    multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                    multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                    atom.setMultipoleType(multipoleType, null);
                    multipole[t000] = multipoleType.charge;
                    multipole[t100] = multipoleType.dipole[0];
                    multipole[t010] = multipoleType.dipole[1];
                    multipole[t001] = multipoleType.dipole[2];
                    multipole[t200] = multipoleType.quadrupole[0][0];
                    multipole[t020] = multipoleType.quadrupole[1][1];
                    multipole[t002] = multipoleType.quadrupole[2][2];
                    multipole[t110] = multipoleType.quadrupole[0][1];
                    multipole[t101] = multipoleType.quadrupole[0][2];
                    multipole[t011] = multipoleType.quadrupole[1][2];
                    axisAtom[i] = multipoleReferenceAtoms;
                    frame[i] = multipoleType.frameDefinition;
                    return true;
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
                        int multipoleReferenceAtoms[] = new int[3];
                        multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                        multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                        multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                        atom.setMultipoleType(multipoleType, null);
                        multipole[t000] = multipoleType.charge;
                        multipole[t100] = multipoleType.dipole[0];
                        multipole[t010] = multipoleType.dipole[1];
                        multipole[t001] = multipoleType.dipole[2];
                        multipole[t200] = multipoleType.quadrupole[0][0];
                        multipole[t020] = multipoleType.quadrupole[1][1];
                        multipole[t002] = multipoleType.quadrupole[2][2];
                        multipole[t110] = multipoleType.quadrupole[0][1];
                        multipole[t101] = multipoleType.quadrupole[0][2];
                        multipole[t011] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
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
                            int multipoleReferenceAtoms[] = new int[3];
                            multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                            multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                            multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                            atom.setMultipoleType(multipoleType, null);
                            multipole[t000] = multipoleType.charge;
                            multipole[t100] = multipoleType.dipole[0];
                            multipole[t010] = multipoleType.dipole[1];
                            multipole[t001] = multipoleType.dipole[2];
                            multipole[t200] = multipoleType.quadrupole[0][0];
                            multipole[t020] = multipoleType.quadrupole[1][1];
                            multipole[t002] = multipoleType.quadrupole[2][2];
                            multipole[t110] = multipoleType.quadrupole[0][1];
                            multipole[t101] = multipoleType.quadrupole[0][2];
                            multipole[t011] = multipoleType.quadrupole[1][2];
                            axisAtom[i] = multipoleReferenceAtoms;
                            frame[i] = multipoleType.frameDefinition;
                            return true;
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
                        int multipoleReferenceAtoms[] = new int[2];
                        multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                        multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                        atom.setMultipoleType(multipoleType, null);
                        multipole[t000] = multipoleType.charge;
                        multipole[t100] = multipoleType.dipole[0];
                        multipole[t010] = multipoleType.dipole[1];
                        multipole[t001] = multipoleType.dipole[2];
                        multipole[t200] = multipoleType.quadrupole[0][0];
                        multipole[t020] = multipoleType.quadrupole[1][1];
                        multipole[t002] = multipoleType.quadrupole[2][2];
                        multipole[t110] = multipoleType.quadrupole[0][1];
                        multipole[t101] = multipoleType.quadrupole[0][2];
                        multipole[t011] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
                    }
                    for (Angle angle2 : angles) {
                        Atom atom4 = angle2.get1_3(atom);
                        if (atom4 != null && atom4 != atom3) {
                            String key4 = atom4.getAtomType().getKey();
                            key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                            multipoleType = forceField.getMultipoleType(key);
                            if (multipoleType != null) {
                                int multipoleReferenceAtoms[] = new int[3];
                                multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                                multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                                multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                                atom.setMultipoleType(multipoleType, null);
                                multipole[t000] = multipoleType.charge;
                                multipole[t100] = multipoleType.dipole[0];
                                multipole[t010] = multipoleType.dipole[1];
                                multipole[t001] = multipoleType.dipole[2];
                                multipole[t200] = multipoleType.quadrupole[0][0];
                                multipole[t020] = multipoleType.quadrupole[1][1];
                                multipole[t002] = multipoleType.quadrupole[2][2];
                                multipole[t110] = multipoleType.quadrupole[0][1];
                                multipole[t101] = multipoleType.quadrupole[0][2];
                                multipole[t011] = multipoleType.quadrupole[1][2];
                                axisAtom[i] = multipoleReferenceAtoms;
                                frame[i] = multipoleType.frameDefinition;
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }
    
    private void convertBohrToElectronAngstroms() {
        for (int i = 0; i < 3; i++) {
            dipole[i] *= BOHR;
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                quadrupole[i][j] *= BOHR2;
            }
        }
    }

    public static void checkMultipoleChirality(MultipoleFrameDefinition frame,
            double localOrigin[], double frameCoords[][],
            double dipole[], double quadrupole[][]) {
        if (frame != MultipoleFrameDefinition.ZTHENX) {
            return;
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
        if (vol < 0.0) {
            dipole[1] = -dipole[1];
            quadrupole[0][1] = -quadrupole[0][1];
            quadrupole[1][0] = -quadrupole[1][0];
            quadrupole[1][2] = -quadrupole[1][2];
            quadrupole[2][1] = -quadrupole[2][1];
        }
    }

    public static void getRotationMatrix(MultipoleFrameDefinition frame,
            double localOrigin[], double frameCoords[][], double rotmat[][]) {
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
    }

    public static void rotateDipole(double rotmat[][], double dipole[],
            double rotatedDipole[]) {
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
                + "%11$s % 7.5f % 7.5f % 7.5f \\\n" + "%11$s % 7.5f \\\n"
                + "%11$s % 7.5f % 7.5f \\\n" + "%11$s % 7.5f % 7.5f % 7.5f",
                charge, dipole[0] / BOHR, dipole[1] / BOHR, dipole[2] / BOHR,
                quadrupole[0][0] / BOHR2, quadrupole[1][0] / BOHR2,
                quadrupole[1][1] / BOHR2, quadrupole[2][0] / BOHR2,
                quadrupole[2][1] / BOHR2, quadrupole[2][2] / BOHR2,
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
            multipoleBuffer.replace(comma, comma+1, "");
            multipoleBuffer.append(format("): "));
        }
        multipoleBuffer.append(format("[%6.3f / %6.3f %6.3f %6.3f / %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]",
                charge, dipole[0] / BOHR, dipole[1] / BOHR, dipole[2] / BOHR,
                quadrupole[0][0] / BOHR2, quadrupole[1][0] / BOHR2,
                quadrupole[1][1] / BOHR2, quadrupole[2][0] / BOHR2,
                quadrupole[2][1] / BOHR2, quadrupole[2][2] / BOHR2));
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
        multipoleBuffer.append(format("  % 7.5f\\\n"
                + "%11$s % 7.5f % 7.5f % 7.5f\\\n" + "%11$s % 7.5f\\\n"
                + "%11$s % 7.5f % 7.5f\\\n" + "%11$s % 7.5f % 7.5f % 7.5f",
                charge, dipole[0] * DEBYE, dipole[1] * DEBYE,
                dipole[2] * DEBYE, quadrupole[0][0] * BUCKINGHAM,
                quadrupole[1][0] * BUCKINGHAM, quadrupole[1][1] * BUCKINGHAM,
                quadrupole[2][0] * BUCKINGHAM, quadrupole[2][1] * BUCKINGHAM,
                quadrupole[2][2] * BUCKINGHAM,
                "                                      "));
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
     *
     * @param multipoleType1
     * @param multipoleType2
     * @param multipoleFrameTypes
     * @return
     */
    public static MultipoleType average(MultipoleType multipoleType1, MultipoleType multipoleType2, int[] multipoleFrameTypes) {
        if (multipoleType1 == null || multipoleType2 == null || multipoleFrameTypes != null) {
            return null;
        }
        if (multipoleType1.frameDefinition != multipoleType2.frameDefinition) {
            return null;
        }
        MultipoleFrameDefinition frame = multipoleType1.frameDefinition;
        double charge = (multipoleType1.charge + multipoleType2.charge) / 2.0;
        double[] dipole = new double[3];
        double[][] quadrupole = new double[3][3];
        for (int i = 0; i < 3; i++) {
            dipole[i] = (multipoleType1.dipole[i] + multipoleType2.dipole[i]) / 2.0;
            for (int j = 0; j < 3; j++) {
                quadrupole[i][j] = (multipoleType1.quadrupole[i][j] + multipoleType2.quadrupole[i][j]) / 2.0;
            }
        }
        return new MultipoleType(charge, dipole, quadrupole, multipoleFrameTypes, frame);
    }
    
    /**
     * Create a new MultipoleType representing the weighted average of 
     * @param types
     * @param weights
     * @param frameTypes
     * @return 
     */
    public static MultipoleType scale(MultipoleType[] types, double[] weights, int[] frameTypes) {
        if (Arrays.asList(types).contains(null)) {
            // Multipoles have not yet been assigned.
            return null;
        }
        if (types == null || weights == null || types.length != weights.length) {
            throw new IllegalArgumentException();
        }
        MultipoleFrameDefinition frame = types[0].frameDefinition;
        for (MultipoleType type : types) {
            if (type.frameDefinition != frame) {
                logger.severe("All frame definitions must match.");
            }
        }
//        double denominator = Arrays.stream(weights).sum();
//        if (denominator != 1.0) {
//            logger.warning("Input multipole weights did not sum to unity; normalizing.");
//            for (int w = 0; w < weights.length; w++) {
//                weights[w] = weights[w] / denominator;
//            }
//        }
        double weightedCharge = 0.0;
        double weightedDipole[] = new double[3];
        double weightedQuadrupole[][] = new double[3][3];
        for (int d = 0; d < 3; d++) {
            weightedDipole[d] = 0.0;
            for (int q = 0; q < 3; q++) {
                weightedQuadrupole[d][q] = 0.0;
            }
        }
        for (int i = 0; i < types.length; i++) {
            MultipoleType type = types[i];
            double weight = weights[i];
            weightedCharge += type.charge * weight;
            for (int d = 0; d < 3; d++) {
                weightedDipole[d] += type.dipole[d] * weight;
                for (int q = 0; q < 3; q++) {
                    weightedQuadrupole[d][q] += type.quadrupole[d][q] * weight;
                }
            }
        }
        return new MultipoleType(weightedCharge, weightedDipole, weightedQuadrupole,
                frameTypes, frame, false);
    }

    /**
     * Indices into a 1D tensor array based on compressed tensor notation. This
     * makes multipole code much easier to read.
     */
    public static final int t000 = 0;
    /**
     * Constant <code>t100=1</code>
     */
    public static final int t100 = 1;
    /**
     * Constant <code>t010=2</code>
     */
    public static final int t010 = 2;
    /**
     * Constant <code>t001=3</code>
     */
    public static final int t001 = 3;
    /**
     * Constant <code>t200=4</code>
     */
    public static final int t200 = 4;
    /**
     * Constant <code>t020=5</code>
     */
    public static final int t020 = 5;
    /**
     * Constant <code>t002=6</code>
     */
    public static final int t002 = 6;
    /**
     * Constant <code>t110=7</code>
     */
    public static final int t110 = 7;
    /**
     * Constant <code>t101=8</code>
     */
    public static final int t101 = 8;
    /**
     * Constant <code>t011=9</code>
     */
    public static final int t011 = 9;
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
