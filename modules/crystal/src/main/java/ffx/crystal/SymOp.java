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
package ffx.crystal;

/**
 * The SymOp class defines the rotation and translation of a single symmetry
 * operator.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class SymOp {

    /**
     * The rotation matrix in fractional coordinates.
     */
    public final double[][] rot;
    /**
     * The translation vector in fractional coordinates.
     */
    public final double[] tr;

    /**
     * The SymOp constructor.
     *
     * @param rot The rotation matrix.
     * @param tr The translation vector.
     */
    public SymOp(double[][] rot, double[] tr) {
        this.rot = rot;
        this.tr = tr;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(" Rotation operator:\n");
        sb.append(String.format(" [[%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]]\n",
                rot[0][0], rot[0][1], rot[0][2],
                rot[1][0], rot[1][1], rot[1][2],
                rot[2][0], rot[2][1], rot[2][2]));
        sb.append(" Translation:\n");
        sb.append(String.format(" [%4.2f,%4.2f,%4.2f]", tr[0], tr[1], tr[2]));
        return sb.toString();
    }

    public static String trtoString(double tr) {
        if (tr == f12) {
            return "+1/2";
        }
        if (tr == f13) {
            return "+1/3";
        }
        if (tr == f23) {
            return "+2/3";
        }
        if (tr == f14) {
            return "+1/4";
        }
        if (tr == f34) {
            return "+3/4";
        }
        if (tr == f16) {
            return "+1/6";
        }
        if (tr == f56) {
            return "+5/6";
        }

        return "";
    }

    public String toXYZString() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < 3; i++) {
            boolean s = false;
            if (rot[i][0] < 0.0) {
                sb.append("-X");
                if (tr[0] > 0.0) {
                    sb.append(trtoString(tr[0]));
                }
                s = true;
            } else if (rot[i][0] > 0.0) {
                sb.append("X");
                if (tr[0] > 0.0) {
                    sb.append(trtoString(tr[0]));
                }
                s = true;
            }
            if (rot[i][1] < 0.0) {
                sb.append("-Y");
                if (tr[1] > 0.0) {
                    sb.append(trtoString(tr[1]));
                }
                s = true;
            } else if (rot[i][1] > 0.0) {
                sb.append(s ? "+Y" : "Y");
                if (tr[1] > 0.0) {
                    sb.append(trtoString(tr[1]));
                }
                s = true;
            }
            if (rot[i][2] < 0.0) {
                sb.append("-Z");
                if (tr[2] > 0.0) {
                    sb.append(trtoString(tr[2]));
                }
                s = true;
            } else if (rot[i][2] > 0.0) {
                sb.append(s ? "+Z" : "Z");
                if (tr[2] > 0.0) {
                    sb.append(trtoString(tr[2]));
                }
                s = true;
            }

            if (i < 2) {
                sb.append(", ");
            }
        }
        return sb.toString();
    }
    private static final double zero = 0.0;
    private static final double f12 = 1.0 / 2.0;
    private static final double f13 = 1.0 / 3.0;
    private static final double f23 = 2.0 / 3.0;
    private static final double f14 = 1.0 / 4.0;
    private static final double f34 = 3.0 / 4.0;
    private static final double f16 = 1.0 / 6.0;
    private static final double f56 = 5.0 / 6.0;
    private static final double[] X = {1.0, zero, zero};
    private static final double[] Y = {zero, 1.0, zero};
    private static final double[] Z = {zero, zero, 1.0};
    private static final double[] mX = {-1.0, zero, zero};
    private static final double[] mY = {zero, -1.0, zero};
    private static final double[] mZ = {zero, zero, -1.0};
    private static final double[] XmY = {1.0, -1.0, zero};
    private static final double[] mXY = {-1.0, 1.0, zero};
    public static final double[][] Rot_Z_mY_X = {Z, mY, X};
    public static final double[][] Rot_Y_mX_mZ = {Y, mX, mZ};
    public static final double[][] Rot_XmY_X_mZ = {XmY, X, mZ};
    public static final double[][] Rot_mX_Y_mZ = {mX, Y, mZ};
    public static final double[][] Rot_X_mZ_Y = {X, mZ, Y};
    public static final double[][] Rot_Y_mXY_Z = {Y, mXY, Z};
    public static final double[][] Rot_Y_mX_Z = {Y, mX, Z};
    public static final double[][] Rot_XmY_X_Z = {XmY, X, Z};
    public static final double[][] Rot_mX_mXY_mZ = {mX, mXY, mZ};
    public static final double[][] Rot_Y_Z_X = {Y, Z, X};
    public static final double[][] Rot_mY_mZ_X = {mY, mZ, X};
    public static final double[][] Rot_X_Z_mY = {X, Z, mY};
    public static final double[][] Rot_XmY_mY_Z = {XmY, mY, Z};
    public static final double[][] Rot_Y_X_mZ = {Y, X, mZ};
    public static final double[][] Rot_Y_mZ_X = {Y, mZ, X};
    public static final double[][] Rot_mXY_Y_Z = {mXY, Y, Z};
    public static final double[][] Rot_mX_mY_mZ = {mX, mY, mZ};
    public static final double[][] Rot_X_Y_mZ = {X, Y, mZ};
    public static final double[][] Rot_mXY_mX_Z = {mXY, mX, Z};
    public static final double[][] Rot_mZ_mY_mX = {mZ, mY, mX};
    public static final double[][] Rot_X_mZ_mY = {X, mZ, mY};
    public static final double[][] Rot_X_Y_Z = {X, Y, Z};
    public static final double[][] Rot_mY_mX_mZ = {mY, mX, mZ};
    public static final double[][] Rot_mY_X_Z = {mY, X, Z};
    public static final double[][] Rot_Z_X_Y = {Z, X, Y};
    public static final double[][] Rot_X_XmY_Z = {X, XmY, Z};
    public static final double[][] Rot_mY_X_mZ = {mY, X, mZ};
    public static final double[][] Rot_mY_Z_mX = {mY, Z, mX};
    public static final double[][] Rot_mY_Z_X = {mY, Z, X};
    public static final double[][] Rot_mX_mZ_mY = {mX, mZ, mY};
    public static final double[][] Rot_mX_Z_Y = {mX, Z, Y};
    public static final double[][] Rot_mZ_mX_mY = {mZ, mX, mY};
    public static final double[][] Rot_X_XmY_mZ = {X, XmY, mZ};
    public static final double[][] Rot_mY_XmY_mZ = {mY, XmY, mZ};
    public static final double[][] Rot_Z_X_mY = {Z, X, mY};
    public static final double[][] Rot_mZ_mY_X = {mZ, mY, X};
    public static final double[][] Rot_X_Z_Y = {X, Z, Y};
    public static final double[][] Rot_Z_mX_mY = {Z, mX, mY};
    public static final double[][] Rot_mX_Z_mY = {mX, Z, mY};
    public static final double[][] Rot_X_mY_Z = {X, mY, Z};
    public static final double[][] Rot_mY_mX_Z = {mY, mX, Z};
    public static final double[][] Rot_Z_mY_mX = {Z, mY, mX};
    public static final double[][] Rot_mX_mY_Z = {mX, mY, Z};
    public static final double[][] Rot_Z_Y_X = {Z, Y, X};
    public static final double[][] Rot_mZ_Y_mX = {mZ, Y, mX};
    public static final double[][] Rot_Y_Z_mX = {Y, Z, mX};
    public static final double[][] Rot_mY_XmY_Z = {mY, XmY, Z};
    public static final double[][] Rot_mXY_Y_mZ = {mXY, Y, mZ};
    public static final double[][] Rot_mZ_mX_Y = {mZ, mX, Y};
    public static final double[][] Rot_mX_mZ_Y = {mX, mZ, Y};
    public static final double[][] Rot_mX_Y_Z = {mX, Y, Z};
    public static final double[][] Rot_X_mY_mZ = {X, mY, mZ};
    public static final double[][] Rot_mZ_X_Y = {mZ, X, Y};
    public static final double[][] Rot_Y_mZ_mX = {Y, mZ, mX};
    public static final double[][] Rot_mY_mZ_mX = {mY, mZ, mX};
    public static final double[][] Rot_mZ_Y_X = {mZ, Y, X};
    public static final double[][] Rot_Z_Y_mX = {Z, Y, mX};
    public static final double[][] Rot_mXY_mX_mZ = {mXY, mX, mZ};
    public static final double[][] Rot_XmY_mY_mZ = {XmY, mY, mZ};
    public static final double[][] Rot_Z_mX_Y = {Z, mX, Y};
    public static final double[][] Rot_mX_mXY_Z = {mX, mXY, Z};
    public static final double[][] Rot_Y_mXY_mZ = {Y, mXY, mZ};
    public static final double[][] Rot_mZ_X_mY = {mZ, X, mY};
    public static final double[][] Rot_Y_X_Z = {Y, X, Z};
    public static final double[] Tr_0_0_34 = {zero, zero, f34};
    public static final double[] Tr_12_0_34 = {f12, zero, f34};
    public static final double[] Tr_0_0_56 = {zero, zero, f56};
    public static final double[] Tr_12_0_12 = {f12, zero, f12};
    public static final double[] Tr_0_12_12 = {zero, f12, f12};
    public static final double[] Tr_12_0_14 = {f12, zero, f14};
    public static final double[] Tr_0_12_14 = {zero, f12, f14};
    public static final double[] Tr_14_14_14 = {f14, f14, f14};
    public static final double[] Tr_0_12_34 = {zero, f12, f34};
    public static final double[] Tr_34_14_14 = {f34, f14, f14};
    public static final double[] Tr_0_0_0 = {zero, zero, zero};
    public static final double[] Tr_23_13_56 = {f23, f13, f56};
    public static final double[] Tr_14_14_34 = {f14, f14, f34};
    public static final double[] Tr_12_12_0 = {f12, f12, zero};
    public static final double[] Tr_23_13_13 = {f23, f13, f13};
    public static final double[] Tr_13_23_23 = {f13, f23, f23};
    public static final double[] Tr_12_12_12 = {f12, f12, f12};
    public static final double[] Tr_12_12_14 = {f12, f12, f14};
    public static final double[] Tr_14_34_14 = {f14, f34, f14};
    public static final double[] Tr_12_12_34 = {f12, f12, f34};
    public static final double[] Tr_0_0_23 = {zero, zero, f23};
    public static final double[] Tr_0_12_0 = {zero, f12, zero};
    public static final double[] Tr_14_34_34 = {f14, f34, f34};
    public static final double[] Tr_34_34_14 = {f34, f34, f14};
    public static final double[] Tr_12_0_0 = {f12, zero, zero};
    public static final double[] Tr_34_34_34 = {f34, f34, f34};
    public static final double[] Tr_0_0_13 = {zero, zero, f13};
    public static final double[] Tr_0_0_12 = {zero, zero, f12};
    public static final double[] Tr_13_23_16 = {f13, f23, f16};
    public static final double[] Tr_0_0_14 = {zero, zero, f14};
    public static final double[] Tr_0_0_16 = {zero, zero, f16};
    public static final double[] Tr_34_14_34 = {f34, f14, f34};
}
