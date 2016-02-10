/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.crystal;

/**
 * The SymOp class defines the rotation and translation of a single symmetry
 * operator.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 *
 * @see SpaceGroup
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

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(" Rotation operator:\n");
        sb.append(String.format(" [[%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]]\n",
                rot[0][0], rot[0][1], rot[0][2],
                rot[1][0], rot[1][1], rot[1][2],
                rot[2][0], rot[2][1], rot[2][2]));
        sb.append(" Translation:\n");
        sb.append(String.format(" [%4.2f,%4.2f,%4.2f]", tr[0], tr[1], tr[2]));
        return sb.toString();
    }

    /**
     * <p>
     * trtoString</p>
     *
     * @param tr a double.
     * @return a {@link java.lang.String} object.
     */
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

    /**
     * <p>
     * toXYZString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toXYZString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 3; i++) {
            boolean s = false;
            if (rot[i][0] < 0.0) {
                sb.append("-X");
                s = true;
            } else if (rot[i][0] > 0.0) {
                sb.append("X");
                s = true;
            }
            if (rot[i][1] < 0.0) {
                sb.append("-Y");
                s = true;
            } else if (rot[i][1] > 0.0) {
                sb.append(s ? "+Y" : "Y");
                s = true;
            }
            if (rot[i][2] < 0.0) {
                sb.append("-Z");
            } else if (rot[i][2] > 0.0) {
                sb.append(s ? "+Z" : "Z");
            }

            if (tr[i] > 0.0) {
                sb.append(trtoString(tr[i]));
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
    /**
     * Constant <code>Rot_Z_mY_X={Z, mY, X}</code>
     */
    public static final double[][] Rot_Z_mY_X = {Z, mY, X};
    /**
     * Constant <code>Rot_Y_mX_mZ={Y, mX, mZ}</code>
     */
    public static final double[][] Rot_Y_mX_mZ = {Y, mX, mZ};
    /**
     * Constant <code>Rot_XmY_X_mZ={XmY, X, mZ}</code>
     */
    public static final double[][] Rot_XmY_X_mZ = {XmY, X, mZ};
    /**
     * Constant <code>Rot_mX_Y_mZ={mX, Y, mZ}</code>
     */
    public static final double[][] Rot_mX_Y_mZ = {mX, Y, mZ};
    /**
     * Constant <code>Rot_X_mZ_Y={X, mZ, Y}</code>
     */
    public static final double[][] Rot_X_mZ_Y = {X, mZ, Y};
    /**
     * Constant <code>Rot_Y_mXY_Z={Y, mXY, Z}</code>
     */
    public static final double[][] Rot_Y_mXY_Z = {Y, mXY, Z};
    /**
     * Constant <code>Rot_Y_mX_Z={Y, mX, Z}</code>
     */
    public static final double[][] Rot_Y_mX_Z = {Y, mX, Z};
    /**
     * Constant <code>Rot_XmY_X_Z={XmY, X, Z}</code>
     */
    public static final double[][] Rot_XmY_X_Z = {XmY, X, Z};
    /**
     * Constant <code>Rot_mX_mXY_mZ={mX, mXY, mZ}</code>
     */
    public static final double[][] Rot_mX_mXY_mZ = {mX, mXY, mZ};
    /**
     * Constant <code>Rot_Y_Z_X={Y, Z, X}</code>
     */
    public static final double[][] Rot_Y_Z_X = {Y, Z, X};
    /**
     * Constant <code>Rot_mY_mZ_X={mY, mZ, X}</code>
     */
    public static final double[][] Rot_mY_mZ_X = {mY, mZ, X};
    /**
     * Constant <code>Rot_X_Z_mY={X, Z, mY}</code>
     */
    public static final double[][] Rot_X_Z_mY = {X, Z, mY};
    /**
     * Constant <code>Rot_XmY_mY_Z={XmY, mY, Z}</code>
     */
    public static final double[][] Rot_XmY_mY_Z = {XmY, mY, Z};
    /**
     * Constant <code>Rot_Y_X_mZ={Y, X, mZ}</code>
     */
    public static final double[][] Rot_Y_X_mZ = {Y, X, mZ};
    /**
     * Constant <code>Rot_Y_mZ_X={Y, mZ, X}</code>
     */
    public static final double[][] Rot_Y_mZ_X = {Y, mZ, X};
    /**
     * Constant <code>Rot_mXY_Y_Z={mXY, Y, Z}</code>
     */
    public static final double[][] Rot_mXY_Y_Z = {mXY, Y, Z};
    /**
     * Constant <code>Rot_mX_mY_mZ={mX, mY, mZ}</code>
     */
    public static final double[][] Rot_mX_mY_mZ = {mX, mY, mZ};
    /**
     * Constant <code>Rot_X_Y_mZ={X, Y, mZ}</code>
     */
    public static final double[][] Rot_X_Y_mZ = {X, Y, mZ};
    /**
     * Constant <code>Rot_mXY_mX_Z={mXY, mX, Z}</code>
     */
    public static final double[][] Rot_mXY_mX_Z = {mXY, mX, Z};
    /**
     * Constant <code>Rot_mZ_mY_mX={mZ, mY, mX}</code>
     */
    public static final double[][] Rot_mZ_mY_mX = {mZ, mY, mX};
    /**
     * Constant <code>Rot_X_mZ_mY={X, mZ, mY}</code>
     */
    public static final double[][] Rot_X_mZ_mY = {X, mZ, mY};
    /**
     * Constant <code>Rot_X_Y_Z={X, Y, Z}</code>
     */
    public static final double[][] Rot_X_Y_Z = {X, Y, Z};
    /**
     * Constant <code>Rot_mY_mX_mZ={mY, mX, mZ}</code>
     */
    public static final double[][] Rot_mY_mX_mZ = {mY, mX, mZ};
    /**
     * Constant <code>Rot_mY_X_Z={mY, X, Z}</code>
     */
    public static final double[][] Rot_mY_X_Z = {mY, X, Z};
    /**
     * Constant <code>Rot_Z_X_Y={Z, X, Y}</code>
     */
    public static final double[][] Rot_Z_X_Y = {Z, X, Y};
    /**
     * Constant <code>Rot_X_XmY_Z={X, XmY, Z}</code>
     */
    public static final double[][] Rot_X_XmY_Z = {X, XmY, Z};
    /**
     * Constant <code>Rot_mY_X_mZ={mY, X, mZ}</code>
     */
    public static final double[][] Rot_mY_X_mZ = {mY, X, mZ};
    /**
     * Constant <code>Rot_mY_Z_mX={mY, Z, mX}</code>
     */
    public static final double[][] Rot_mY_Z_mX = {mY, Z, mX};
    /**
     * Constant <code>Rot_mY_Z_X={mY, Z, X}</code>
     */
    public static final double[][] Rot_mY_Z_X = {mY, Z, X};
    /**
     * Constant <code>Rot_mX_mZ_mY={mX, mZ, mY}</code>
     */
    public static final double[][] Rot_mX_mZ_mY = {mX, mZ, mY};
    /**
     * Constant <code>Rot_mX_Z_Y={mX, Z, Y}</code>
     */
    public static final double[][] Rot_mX_Z_Y = {mX, Z, Y};
    /**
     * Constant <code>Rot_mZ_mX_mY={mZ, mX, mY}</code>
     */
    public static final double[][] Rot_mZ_mX_mY = {mZ, mX, mY};
    /**
     * Constant <code>Rot_X_XmY_mZ={X, XmY, mZ}</code>
     */
    public static final double[][] Rot_X_XmY_mZ = {X, XmY, mZ};
    /**
     * Constant <code>Rot_mY_XmY_mZ={mY, XmY, mZ}</code>
     */
    public static final double[][] Rot_mY_XmY_mZ = {mY, XmY, mZ};
    /**
     * Constant <code>Rot_Z_X_mY={Z, X, mY}</code>
     */
    public static final double[][] Rot_Z_X_mY = {Z, X, mY};
    /**
     * Constant <code>Rot_mZ_mY_X={mZ, mY, X}</code>
     */
    public static final double[][] Rot_mZ_mY_X = {mZ, mY, X};
    /**
     * Constant <code>Rot_X_Z_Y={X, Z, Y}</code>
     */
    public static final double[][] Rot_X_Z_Y = {X, Z, Y};
    /**
     * Constant <code>Rot_Z_mX_mY={Z, mX, mY}</code>
     */
    public static final double[][] Rot_Z_mX_mY = {Z, mX, mY};
    /**
     * Constant <code>Rot_mX_Z_mY={mX, Z, mY}</code>
     */
    public static final double[][] Rot_mX_Z_mY = {mX, Z, mY};
    /**
     * Constant <code>Rot_X_mY_Z={X, mY, Z}</code>
     */
    public static final double[][] Rot_X_mY_Z = {X, mY, Z};
    /**
     * Constant <code>Rot_mY_mX_Z={mY, mX, Z}</code>
     */
    public static final double[][] Rot_mY_mX_Z = {mY, mX, Z};
    /**
     * Constant <code>Rot_Z_mY_mX={Z, mY, mX}</code>
     */
    public static final double[][] Rot_Z_mY_mX = {Z, mY, mX};
    /**
     * Constant <code>Rot_mX_mY_Z={mX, mY, Z}</code>
     */
    public static final double[][] Rot_mX_mY_Z = {mX, mY, Z};
    /**
     * Constant <code>Rot_Z_Y_X={Z, Y, X}</code>
     */
    public static final double[][] Rot_Z_Y_X = {Z, Y, X};
    /**
     * Constant <code>Rot_mZ_Y_mX={mZ, Y, mX}</code>
     */
    public static final double[][] Rot_mZ_Y_mX = {mZ, Y, mX};
    /**
     * Constant <code>Rot_Y_Z_mX={Y, Z, mX}</code>
     */
    public static final double[][] Rot_Y_Z_mX = {Y, Z, mX};
    /**
     * Constant <code>Rot_mY_XmY_Z={mY, XmY, Z}</code>
     */
    public static final double[][] Rot_mY_XmY_Z = {mY, XmY, Z};
    /**
     * Constant <code>Rot_mXY_Y_mZ={mXY, Y, mZ}</code>
     */
    public static final double[][] Rot_mXY_Y_mZ = {mXY, Y, mZ};
    /**
     * Constant <code>Rot_mZ_mX_Y={mZ, mX, Y}</code>
     */
    public static final double[][] Rot_mZ_mX_Y = {mZ, mX, Y};
    /**
     * Constant <code>Rot_mX_mZ_Y={mX, mZ, Y}</code>
     */
    public static final double[][] Rot_mX_mZ_Y = {mX, mZ, Y};
    /**
     * Constant <code>Rot_mX_Y_Z={mX, Y, Z}</code>
     */
    public static final double[][] Rot_mX_Y_Z = {mX, Y, Z};
    /**
     * Constant <code>Rot_X_mY_mZ={X, mY, mZ}</code>
     */
    public static final double[][] Rot_X_mY_mZ = {X, mY, mZ};
    /**
     * Constant <code>Rot_mZ_X_Y={mZ, X, Y}</code>
     */
    public static final double[][] Rot_mZ_X_Y = {mZ, X, Y};
    /**
     * Constant <code>Rot_Y_mZ_mX={Y, mZ, mX}</code>
     */
    public static final double[][] Rot_Y_mZ_mX = {Y, mZ, mX};
    /**
     * Constant <code>Rot_mY_mZ_mX={mY, mZ, mX}</code>
     */
    public static final double[][] Rot_mY_mZ_mX = {mY, mZ, mX};
    /**
     * Constant <code>Rot_mZ_Y_X={mZ, Y, X}</code>
     */
    public static final double[][] Rot_mZ_Y_X = {mZ, Y, X};
    /**
     * Constant <code>Rot_Z_Y_mX={Z, Y, mX}</code>
     */
    public static final double[][] Rot_Z_Y_mX = {Z, Y, mX};
    /**
     * Constant <code>Rot_mXY_mX_mZ={mXY, mX, mZ}</code>
     */
    public static final double[][] Rot_mXY_mX_mZ = {mXY, mX, mZ};
    /**
     * Constant <code>Rot_XmY_mY_mZ={XmY, mY, mZ}</code>
     */
    public static final double[][] Rot_XmY_mY_mZ = {XmY, mY, mZ};
    /**
     * Constant <code>Rot_Z_mX_Y={Z, mX, Y}</code>
     */
    public static final double[][] Rot_Z_mX_Y = {Z, mX, Y};
    /**
     * Constant <code>Rot_mX_mXY_Z={mX, mXY, Z}</code>
     */
    public static final double[][] Rot_mX_mXY_Z = {mX, mXY, Z};
    /**
     * Constant <code>Rot_Y_mXY_mZ={Y, mXY, mZ}</code>
     */
    public static final double[][] Rot_Y_mXY_mZ = {Y, mXY, mZ};
    /**
     * Constant <code>Rot_mZ_X_mY={mZ, X, mY}</code>
     */
    public static final double[][] Rot_mZ_X_mY = {mZ, X, mY};
    /**
     * Constant <code>Rot_Y_X_Z={Y, X, Z}</code>
     */
    public static final double[][] Rot_Y_X_Z = {Y, X, Z};
    /**
     * Constant <code>Tr_0_0_34={zero, zero, f34}</code>
     */
    public static final double[] Tr_0_0_34 = {zero, zero, f34};
    /**
     * Constant <code>Tr_12_0_34={f12, zero, f34}</code>
     */
    public static final double[] Tr_12_0_34 = {f12, zero, f34};
    /**
     * Constant <code>Tr_0_0_56={zero, zero, f56}</code>
     */
    public static final double[] Tr_0_0_56 = {zero, zero, f56};
    /**
     * Constant <code>Tr_12_0_12={f12, zero, f12}</code>
     */
    public static final double[] Tr_12_0_12 = {f12, zero, f12};
    /**
     * Constant <code>Tr_0_12_12={zero, f12, f12}</code>
     */
    public static final double[] Tr_0_12_12 = {zero, f12, f12};
    /**
     * Constant <code>Tr_12_0_14={f12, zero, f14}</code>
     */
    public static final double[] Tr_12_0_14 = {f12, zero, f14};
    /**
     * Constant <code>Tr_0_12_14={zero, f12, f14}</code>
     */
    public static final double[] Tr_0_12_14 = {zero, f12, f14};
    /**
     * Constant <code>Tr_14_14_14={f14, f14, f14}</code>
     */
    public static final double[] Tr_14_14_14 = {f14, f14, f14};
    /**
     * Constant <code>Tr_0_12_34={zero, f12, f34}</code>
     */
    public static final double[] Tr_0_12_34 = {zero, f12, f34};
    /**
     * Constant <code>Tr_34_14_14={f34, f14, f14}</code>
     */
    public static final double[] Tr_34_14_14 = {f34, f14, f14};
    /**
     * Constant <code>Tr_0_0_0={zero, zero, zero}</code>
     */
    public static final double[] Tr_0_0_0 = {zero, zero, zero};
    /**
     * Constant <code>Tr_23_13_56={f23, f13, f56}</code>
     */
    public static final double[] Tr_23_13_56 = {f23, f13, f56};
    /**
     * Constant <code>Tr_14_14_34={f14, f14, f34}</code>
     */
    public static final double[] Tr_14_14_34 = {f14, f14, f34};
    /**
     * Constant <code>Tr_12_12_0={f12, f12, zero}</code>
     */
    public static final double[] Tr_12_12_0 = {f12, f12, zero};
    /**
     * Constant <code>Tr_23_13_13={f23, f13, f13}</code>
     */
    public static final double[] Tr_23_13_13 = {f23, f13, f13};
    /**
     * Constant <code>Tr_13_23_23={f13, f23, f23}</code>
     */
    public static final double[] Tr_13_23_23 = {f13, f23, f23};
    /**
     * Constant <code>Tr_12_12_12={f12, f12, f12}</code>
     */
    public static final double[] Tr_12_12_12 = {f12, f12, f12};
    /**
     * Constant <code>Tr_12_12_14={f12, f12, f14}</code>
     */
    public static final double[] Tr_12_12_14 = {f12, f12, f14};
    /**
     * Constant <code>Tr_14_34_14={f14, f34, f14}</code>
     */
    public static final double[] Tr_14_34_14 = {f14, f34, f14};
    /**
     * Constant <code>Tr_12_12_34={f12, f12, f34}</code>
     */
    public static final double[] Tr_12_12_34 = {f12, f12, f34};
    /**
     * Constant <code>Tr_0_0_23={zero, zero, f23}</code>
     */
    public static final double[] Tr_0_0_23 = {zero, zero, f23};
    /**
     * Constant <code>Tr_0_12_0={zero, f12, zero}</code>
     */
    public static final double[] Tr_0_12_0 = {zero, f12, zero};
    /**
     * Constant <code>Tr_14_34_34={f14, f34, f34}</code>
     */
    public static final double[] Tr_14_34_34 = {f14, f34, f34};
    /**
     * Constant <code>Tr_34_34_14={f34, f34, f14}</code>
     */
    public static final double[] Tr_34_34_14 = {f34, f34, f14};
    /**
     * Constant <code>Tr_12_0_0={f12, zero, zero}</code>
     */
    public static final double[] Tr_12_0_0 = {f12, zero, zero};
    /**
     * Constant <code>Tr_34_34_34={f34, f34, f34}</code>
     */
    public static final double[] Tr_34_34_34 = {f34, f34, f34};
    /**
     * Constant <code>Tr_0_0_13={zero, zero, f13}</code>
     */
    public static final double[] Tr_0_0_13 = {zero, zero, f13};
    /**
     * Constant <code>Tr_0_0_12={zero, zero, f12}</code>
     */
    public static final double[] Tr_0_0_12 = {zero, zero, f12};
    /**
     * Constant <code>Tr_13_23_16={f13, f23, f16}</code>
     */
    public static final double[] Tr_13_23_16 = {f13, f23, f16};
    /**
     * Constant <code>Tr_0_0_14={zero, zero, f14}</code>
     */
    public static final double[] Tr_0_0_14 = {zero, zero, f14};
    /**
     * Constant <code>Tr_0_0_16={zero, zero, f16}</code>
     */
    public static final double[] Tr_0_0_16 = {zero, zero, f16};
    /**
     * Constant <code>Tr_34_14_34={f34, f14, f34}</code>
     */
    public static final double[] Tr_34_14_34 = {f34, f14, f34};
}
