/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.parameters;

import static java.lang.Math.abs;

import java.util.Arrays;
import java.util.Comparator;
import java.util.logging.Logger;

/**
 * The MultipoleType class defines a multipole in its local frame.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class MultipoleType extends BaseType implements Comparator<String> {

    private static final Logger logger = Logger.getLogger(MultipoleType.class.getName());

    /**
     * The local multipole frame is defined by the Z-then-X or Bisector
     * convention.
     */
    public enum MultipoleFrameDefinition {

        ZTHENX, BISECTOR, ZTHENBISECTOR, TRISECTOR
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
     * Multipole Constructor. This assumes the dipole and quadrupole are in
     * units of Bohr, and are converted to electron-Angstroms and
     * electron-Angstroms^2, respectively, before the constructor returns.
     *
     * @param charge
     *            double
     * @param dipole
     *            double[]
     * @param quadrupole
     *            double[]
     * @param multipoleFrameTypes
     *            int[]
     */
    public MultipoleType(double charge, double dipole[], double quadrupole[][],
            int[] multipoleFrameTypes, MultipoleFrameDefinition frameDefinition) {
        super(ForceField.ForceFieldType.MULTIPOLE, multipoleFrameTypes);
        this.charge = charge;
        this.dipole = dipole;
        this.quadrupole = quadrupole;
        this.frameAtomTypes = multipoleFrameTypes;
        this.frameDefinition = frameDefinition;
        initMultipole();
    }

    private void initMultipole() {
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
        // Convert to electron-Angstroms
        for (int i = 0; i < 3; i++) {
            dipole[i] *= BOHR;
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                quadrupole[i][j] *= BOHR * BOHR;
            }
        }
        // Warn if the multipole is not traceless.
        double sum = quadrupole[0][0] + quadrupole[1][1] + quadrupole[2][2];
        if (Math.abs(sum) > 1.0e-5) {
            String message = String.format("Multipole is not traceless: %7.5f",
                    sum);
            logger.warning(message + "\n" + toBohrString());
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
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
        } else if (frameDefinition == MultipoleFrameDefinition.ZTHENBISECTOR) {
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[3]));
        } else if (frameDefinition == MultipoleFrameDefinition.TRISECTOR){
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[3]));
        } else {
            for (int i : frameAtomTypes) {
                multipoleBuffer.append(String.format("  %5d", i));
            }
        }
        if (frameAtomTypes.length == 3) {
            multipoleBuffer.append("       ");
        }
        multipoleBuffer.append(String.format("  % 7.5f \\\n"
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
     * Nicely formatted multipole string. Dipole and qaudrupole are in units of
     * Debye and Buckinghams, respectively.
     *
     * @return String
     */
    public String toDebyeString() {
        StringBuilder multipoleBuffer = new StringBuilder("multipole");
        if (frameDefinition == MultipoleFrameDefinition.BISECTOR) {
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
        } else if (frameDefinition == MultipoleFrameDefinition.ZTHENBISECTOR) {
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[3]));
        } else if (frameDefinition == MultipoleFrameDefinition.TRISECTOR){
            multipoleBuffer.append(String.format("  %5d", frameAtomTypes[0]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[1]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[2]));
            multipoleBuffer.append(String.format("  %5d", -frameAtomTypes[3]));
        } else {
            for (int i : frameAtomTypes) {
                multipoleBuffer.append(String.format("  %5d", i));
            }
        }
        if (frameAtomTypes.length == 3) {
            multipoleBuffer.append("       ");
        }
        multipoleBuffer.append(String.format("  % 7.5f\\\n"
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

    @Override
    public String toString() {
        return toBohrString();
    }

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

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Arrays.hashCode(frameAtomTypes);
        return hash;
    }

    /**
     * Indices into a 1D tensor array based on compressed tensor notation. This
     * makes multipole code much easier to read.
     */
    public static final int t000 = 0;
    public static final int t100 = 1;
    public static final int t010 = 2;
    public static final int t001 = 3;
    public static final int t200 = 4;
    public static final int t020 = 5;
    public static final int t002 = 6;
    public static final int t110 = 7;
    public static final int t101 = 8;
    public static final int t011 = 9;
    public static final int t300 = 10;
    public static final int t030 = 11;
    public static final int t003 = 12;
    public static final int t210 = 13;
    public static final int t201 = 14;
    public static final int t120 = 15;
    public static final int t021 = 16;
    public static final int t102 = 17;
    public static final int t012 = 18;
    public static final int t111 = 19;
}
