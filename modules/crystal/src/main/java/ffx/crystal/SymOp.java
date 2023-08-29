// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.crystal;

import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.MatrixMath.mat3Inverse;
import static ffx.numerics.math.MatrixMath.mat4Mat4;
import static ffx.numerics.math.ScalarMath.mod;
import static java.lang.Double.parseDouble;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.rint;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The SymOp class defines the rotation and translation of a single symmetry operator.
 *
 * @author Michael J. Schnieders
 * @see SpaceGroup
 * @since 1.0
 */
public class SymOp {

  /** Constant <code>ZERO = 0.0</code> */
  private static final double ZERO = 0.0;
  /** Constant <code>ZERO_ROTATION = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}</code> */
  public static final double[][] ZERO_ROTATION = {
      {1.0, ZERO, ZERO},
      {ZERO, 1.0, ZERO},
      {ZERO, ZERO, 1.0}};
  /** Constant <code>Tr_0_0_0={ZERO, ZERO, ZERO}</code> */
  public static final double[] Tr_0_0_0 = {ZERO, ZERO, ZERO};
  /** Constant <code>f12 = 1.0 / 2.0</code> */
  private static final double f12 = 1.0 / 2.0;
  /** Constant <code>Tr_12_0_12={f12, ZERO, f12}</code> */
  static final double[] Tr_12_0_12 = {f12, ZERO, f12};
  /** Constant <code>Tr_0_12_12={ZERO, f12, f12}</code> */
  static final double[] Tr_0_12_12 = {ZERO, f12, f12};
  /** Constant <code>Tr_12_12_0={f12, f12, ZERO}</code> */
  static final double[] Tr_12_12_0 = {f12, f12, ZERO};
  /** Constant <code>Tr_12_12_12={f12, f12, f12}</code> */
  static final double[] Tr_12_12_12 = {f12, f12, f12};
  /** Constant <code>Tr_0_12_0={ZERO, f12, ZERO}</code> */
  static final double[] Tr_0_12_0 = {ZERO, f12, ZERO};
  /** Constant <code>Tr_12_0_0={f12, ZERO, ZERO}</code> */
  static final double[] Tr_12_0_0 = {f12, ZERO, ZERO};
  /** Constant <code>Tr_0_0_12={ZERO, ZERO, f12}</code> */
  static final double[] Tr_0_0_12 = {ZERO, ZERO, f12};
  /** Constant <code>f13 = 1.0 / 3.0</code> */
  private static final double f13 = 1.0 / 3.0;
  /** Constant <code>Tr_0_0_13={ZERO, ZERO, f13}</code> */
  static final double[] Tr_0_0_13 = {ZERO, ZERO, f13};
  /** Constant <code>f23 = 2.0 / 3.0</code> */
  private static final double f23 = 2.0 / 3.0;
  /** Constant <code>Tr_23_13_13={f23, f13, f13}</code> */
  static final double[] Tr_23_13_13 = {f23, f13, f13};
  /** Constant <code>Tr_13_23_23={f13, f23, f23}</code> */
  static final double[] Tr_13_23_23 = {f13, f23, f23};
  /** Constant <code>Tr_0_0_23={ZERO, ZERO, f23}</code> */
  static final double[] Tr_0_0_23 = {ZERO, ZERO, f23};
  /** Constant <code>f14 = 1.0 / 4.0</code> */
  private static final double f14 = 1.0 / 4.0;
  /** Constant <code>Tr_12_0_14={f12, ZERO, f14}</code> */
  static final double[] Tr_12_0_14 = {f12, ZERO, f14};
  /** Constant <code>Tr_0_12_14={ZERO, f12, f14}</code> */
  static final double[] Tr_0_12_14 = {ZERO, f12, f14};
  /** Constant <code>Tr_14_14_14={f14, f14, f14}</code> */
  static final double[] Tr_14_14_14 = {f14, f14, f14};
  /** Constant <code>Tr_12_12_14={f12, f12, f14}</code> */
  static final double[] Tr_12_12_14 = {f12, f12, f14};
  /** Constant <code>Tr_0_0_14={ZERO, ZERO, f14}</code> */
  static final double[] Tr_0_0_14 = {ZERO, ZERO, f14};
  /** Constant <code>f34 = 3.0 / 4.0</code> */
  private static final double f34 = 3.0 / 4.0;
  /** Constant <code>Tr_0_0_34={ZERO, ZERO, f34}</code> */
  static final double[] Tr_0_0_34 = {ZERO, ZERO, f34};
  /** Constant <code>Tr_12_0_34={f12, ZERO, f34}</code> */
  static final double[] Tr_12_0_34 = {f12, ZERO, f34};
  /** Constant <code>Tr_0_12_34={ZERO, f12, f34}</code> */
  static final double[] Tr_0_12_34 = {ZERO, f12, f34};
  /** Constant <code>Tr_34_14_14={f34, f14, f14}</code> */
  static final double[] Tr_34_14_14 = {f34, f14, f14};
  /** Constant <code>Tr_14_14_34={f14, f14, f34}</code> */
  static final double[] Tr_14_14_34 = {f14, f14, f34};
  /** Constant <code>Tr_14_34_14={f14, f34, f14}</code> */
  static final double[] Tr_14_34_14 = {f14, f34, f14};
  /** Constant <code>Tr_12_12_34={f12, f12, f34}</code> */
  static final double[] Tr_12_12_34 = {f12, f12, f34};
  /** Constant <code>Tr_14_34_34={f14, f34, f34}</code> */
  static final double[] Tr_14_34_34 = {f14, f34, f34};
  /** Constant <code>Tr_34_34_14={f34, f34, f14}</code> */
  static final double[] Tr_34_34_14 = {f34, f34, f14};
  /** Constant <code>Tr_34_34_34={f34, f34, f34}</code> */
  static final double[] Tr_34_34_34 = {f34, f34, f34};
  /** Constant <code>Tr_34_14_34={f34, f14, f34}</code> */
  static final double[] Tr_34_14_34 = {f34, f14, f34};
  /** Constant <code>f16 = 1.0 / 6.0</code> */
  private static final double f16 = 1.0 / 6.0;
  /** Constant <code>Tr_13_23_16={f13, f23, f16}</code> */
  static final double[] Tr_13_23_16 = {f13, f23, f16};
  /** Constant <code>Tr_0_0_16={ZERO, ZERO, f16}</code> */
  static final double[] Tr_0_0_16 = {ZERO, ZERO, f16};
  /** Constant <code>f56 = 5.0 / 6.0</code> */
  private static final double f56 = 5.0 / 6.0;
  /** Constant <code>Tr_0_0_56={ZERO, ZERO, f56}</code> */
  static final double[] Tr_0_0_56 = {ZERO, ZERO, f56};
  /** Constant <code>Tr_23_13_56={f23, f13, f56}</code> */
  static final double[] Tr_23_13_56 = {f23, f13, f56};
  /** Constant <code>X = {1.0, ZERO, ZERO}</code> */
  private static final double[] X = {1.0, ZERO, ZERO};
  /** Constant <code>Y = {ZERO, 1.0, ZERO}</code> */
  private static final double[] Y = {ZERO, 1.0, ZERO};
  /** Constant <code>Z = {ZERO, ZERO, 1.0}</code> */
  private static final double[] Z = {ZERO, ZERO, 1.0};
  /** Constant <code>Rot_Y_Z_X={Y, Z, X}</code> */
  static final double[][] Rot_Y_Z_X = {Y, Z, X};
  /** Constant <code>Rot_X_Y_Z={X, Y, Z}</code> */
  static final double[][] Rot_X_Y_Z = {X, Y, Z};
  /** Constant <code>Rot_Z_X_Y={Z, X, Y}</code> */
  static final double[][] Rot_Z_X_Y = {Z, X, Y};
  /** Constant <code>Rot_X_Z_Y={X, Z, Y}</code> */
  static final double[][] Rot_X_Z_Y = {X, Z, Y};
  /** Constant <code>Rot_Z_Y_X={Z, Y, X}</code> */
  static final double[][] Rot_Z_Y_X = {Z, Y, X};
  /** Constant <code>Rot_Y_X_Z={Y, X, Z}</code> */
  static final double[][] Rot_Y_X_Z = {Y, X, Z};
  /** Constant <code>mX = {-1.0, ZERO, ZERO}</code> */
  private static final double[] mX = {-1.0, ZERO, ZERO};
  /** Constant <code>Rot_Y_mX_Z={Y, mX, Z}</code> */
  static final double[][] Rot_Y_mX_Z = {Y, mX, Z};
  /** Constant <code>Rot_mX_Z_Y={mX, Z, Y}</code> */
  static final double[][] Rot_mX_Z_Y = {mX, Z, Y};
  /** Constant <code>Rot_Y_Z_mX={Y, Z, mX}</code> */
  static final double[][] Rot_Y_Z_mX = {Y, Z, mX};
  /** Constant <code>Rot_mX_Y_Z={mX, Y, Z}</code> */
  static final double[][] Rot_mX_Y_Z = {mX, Y, Z};
  /** Constant <code>Rot_Z_Y_mX={Z, Y, mX}</code> */
  static final double[][] Rot_Z_Y_mX = {Z, Y, mX};
  /** Constant <code>Rot_Z_mX_Y={Z, mX, Y}</code> */
  static final double[][] Rot_Z_mX_Y = {Z, mX, Y};
  /** Constant <code>mY = {ZERO, -1.0, ZERO}</code> */
  private static final double[] mY = {ZERO, -1.0, ZERO};
  /** Constant <code>Rot_Z_mY_X={Z, mY, X}</code> */
  static final double[][] Rot_Z_mY_X = {Z, mY, X};
  /** Constant <code>Rot_X_Z_mY={X, Z, mY}</code> */
  static final double[][] Rot_X_Z_mY = {X, Z, mY};
  /** Constant <code>Rot_mY_X_Z={mY, X, Z}</code> */
  static final double[][] Rot_mY_X_Z = {mY, X, Z};
  /** Constant <code>Rot_mY_Z_mX={mY, Z, mX}</code> */
  static final double[][] Rot_mY_Z_mX = {mY, Z, mX};
  /** Constant <code>Rot_mY_Z_X={mY, Z, X}</code> */
  static final double[][] Rot_mY_Z_X = {mY, Z, X};
  /** Constant <code>Rot_Z_X_mY={Z, X, mY}</code> */
  static final double[][] Rot_Z_X_mY = {Z, X, mY};
  /** Constant <code>Rot_Z_mX_mY={Z, mX, mY}</code> */
  static final double[][] Rot_Z_mX_mY = {Z, mX, mY};
  /** Constant <code>Rot_mX_Z_mY={mX, Z, mY}</code> */
  static final double[][] Rot_mX_Z_mY = {mX, Z, mY};
  /** Constant <code>Rot_X_mY_Z={X, mY, Z}</code> */
  static final double[][] Rot_X_mY_Z = {X, mY, Z};
  /** Constant <code>Rot_mY_mX_Z={mY, mX, Z}</code> */
  static final double[][] Rot_mY_mX_Z = {mY, mX, Z};
  /** Constant <code>Rot_Z_mY_mX={Z, mY, mX}</code> */
  static final double[][] Rot_Z_mY_mX = {Z, mY, mX};
  /** Constant <code>Rot_mX_mY_Z={mX, mY, Z}</code> */
  static final double[][] Rot_mX_mY_Z = {mX, mY, Z};
  /** Constant <code>mZ = {ZERO, ZERO, -1.0}</code> */
  private static final double[] mZ = {ZERO, ZERO, -1.0};
  /** Constant <code>Rot_Y_mX_mZ={Y, mX, mZ}</code> */
  static final double[][] Rot_Y_mX_mZ = {Y, mX, mZ};
  /** Constant <code>Rot_mX_Y_mZ={mX, Y, mZ}</code> */
  static final double[][] Rot_mX_Y_mZ = {mX, Y, mZ};
  /** Constant <code>Rot_X_mZ_Y={X, mZ, Y}</code> */
  static final double[][] Rot_X_mZ_Y = {X, mZ, Y};
  /** Constant <code>Rot_mY_mZ_X={mY, mZ, X}</code> */
  static final double[][] Rot_mY_mZ_X = {mY, mZ, X};
  /** Constant <code>Rot_Y_X_mZ={Y, X, mZ}</code> */
  static final double[][] Rot_Y_X_mZ = {Y, X, mZ};
  /** Constant <code>Rot_Y_mZ_X={Y, mZ, X}</code> */
  static final double[][] Rot_Y_mZ_X = {Y, mZ, X};
  /** Constant <code>Rot_mX_mY_mZ={mX, mY, mZ}</code> */
  static final double[][] Rot_mX_mY_mZ = {mX, mY, mZ};
  /** Constant <code>Rot_X_Y_mZ={X, Y, mZ}</code> */
  static final double[][] Rot_X_Y_mZ = {X, Y, mZ};
  /** Constant <code>Rot_mZ_mY_mX={mZ, mY, mX}</code> */
  static final double[][] Rot_mZ_mY_mX = {mZ, mY, mX};
  /** Constant <code>Rot_X_mZ_mY={X, mZ, mY}</code> */
  static final double[][] Rot_X_mZ_mY = {X, mZ, mY};
  /** Constant <code>Rot_mY_mX_mZ={mY, mX, mZ}</code> */
  static final double[][] Rot_mY_mX_mZ = {mY, mX, mZ};
  /** Constant <code>Rot_mY_X_mZ={mY, X, mZ}</code> */
  static final double[][] Rot_mY_X_mZ = {mY, X, mZ};
  /** Constant <code>Rot_mX_mZ_mY={mX, mZ, mY}</code> */
  static final double[][] Rot_mX_mZ_mY = {mX, mZ, mY};
  /** Constant <code>Rot_mZ_mX_mY={mZ, mX, mY}</code> */
  static final double[][] Rot_mZ_mX_mY = {mZ, mX, mY};
  /** Constant <code>Rot_mZ_mY_X={mZ, mY, X}</code> */
  static final double[][] Rot_mZ_mY_X = {mZ, mY, X};
  /** Constant <code>Rot_mZ_Y_mX={mZ, Y, mX}</code> */
  static final double[][] Rot_mZ_Y_mX = {mZ, Y, mX};
  /** Constant <code>Rot_mZ_mX_Y={mZ, mX, Y}</code> */
  static final double[][] Rot_mZ_mX_Y = {mZ, mX, Y};
  /** Constant <code>Rot_mX_mZ_Y={mX, mZ, Y}</code> */
  static final double[][] Rot_mX_mZ_Y = {mX, mZ, Y};
  /** Constant <code>Rot_X_mY_mZ={X, mY, mZ}</code> */
  static final double[][] Rot_X_mY_mZ = {X, mY, mZ};
  /** Constant <code>Rot_mZ_X_Y={mZ, X, Y}</code> */
  static final double[][] Rot_mZ_X_Y = {mZ, X, Y};
  /** Constant <code>Rot_Y_mZ_mX={Y, mZ, mX}</code> */
  static final double[][] Rot_Y_mZ_mX = {Y, mZ, mX};
  /** Constant <code>Rot_mY_mZ_mX={mY, mZ, mX}</code> */
  static final double[][] Rot_mY_mZ_mX = {mY, mZ, mX};
  /** Constant <code>Rot_mZ_Y_X={mZ, Y, X}</code> */
  static final double[][] Rot_mZ_Y_X = {mZ, Y, X};
  /** Constant <code>Rot_mZ_X_mY={mZ, X, mY}</code> */
  static final double[][] Rot_mZ_X_mY = {mZ, X, mY};
  /** Constant <code>XmY = {1.0, -1.0, ZERO}</code> */
  private static final double[] XmY = {1.0, -1.0, ZERO};
  /** Constant <code>Rot_XmY_X_mZ={XmY, X, mZ}</code> */
  static final double[][] Rot_XmY_X_mZ = {XmY, X, mZ};
  /** Constant <code>Rot_XmY_X_Z={XmY, X, Z}</code> */
  static final double[][] Rot_XmY_X_Z = {XmY, X, Z};
  /** Constant <code>Rot_XmY_mY_Z={XmY, mY, Z}</code> */
  static final double[][] Rot_XmY_mY_Z = {XmY, mY, Z};
  /** Constant <code>Rot_X_XmY_Z={X, XmY, Z}</code> */
  static final double[][] Rot_X_XmY_Z = {X, XmY, Z};
  /** Constant <code>Rot_X_XmY_mZ={X, XmY, mZ}</code> */
  static final double[][] Rot_X_XmY_mZ = {X, XmY, mZ};
  /** Constant <code>Rot_mY_XmY_mZ={mY, XmY, mZ}</code> */
  static final double[][] Rot_mY_XmY_mZ = {mY, XmY, mZ};
  /** Constant <code>Rot_mY_XmY_Z={mY, XmY, Z}</code> */
  static final double[][] Rot_mY_XmY_Z = {mY, XmY, Z};
  /** Constant <code>Rot_XmY_mY_mZ={XmY, mY, mZ}</code> */
  static final double[][] Rot_XmY_mY_mZ = {XmY, mY, mZ};
  /** Constant <code>mXY = {-1.0, 1.0, ZERO}</code> */
  private static final double[] mXY = {-1.0, 1.0, ZERO};
  /** Constant <code>Rot_Y_mXY_Z={Y, mXY, Z}</code> */
  static final double[][] Rot_Y_mXY_Z = {Y, mXY, Z};
  /** Constant <code>Rot_mX_mXY_mZ={mX, mXY, mZ}</code> */
  static final double[][] Rot_mX_mXY_mZ = {mX, mXY, mZ};
  /** Constant <code>Rot_mXY_Y_Z={mXY, Y, Z}</code> */
  static final double[][] Rot_mXY_Y_Z = {mXY, Y, Z};
  /** Constant <code>Rot_mXY_mX_Z={mXY, mX, Z}</code> */
  static final double[][] Rot_mXY_mX_Z = {mXY, mX, Z};
  /** Constant <code>Rot_mXY_Y_mZ={mXY, Y, mZ}</code> */
  static final double[][] Rot_mXY_Y_mZ = {mXY, Y, mZ};
  /** Constant <code>Rot_mXY_mX_mZ={mXY, mX, mZ}</code> */
  static final double[][] Rot_mXY_mX_mZ = {mXY, mX, mZ};
  /** Constant <code>Rot_mX_mXY_Z={mX, mXY, Z}</code> */
  static final double[][] Rot_mX_mXY_Z = {mX, mXY, Z};
  /** Constant <code>Rot_Y_mXY_mZ={Y, mXY, mZ}</code> */
  static final double[][] Rot_Y_mXY_mZ = {Y, mXY, mZ};
  /** The rotation matrix in fractional coordinates. */
  public final double[][] rot;
  /** The translation vector in fractional coordinates. */
  public final double[] tr;
  /** A mask equal to 0 for X-coordinates. */
  private static final int XX = 0;
  /** A mask equal to 1 for Y-coordinates. */
  private static final int YY = 1;
  /** A mask equal to 2 for Z-coordinates. */
  private static final int ZZ = 2;
  /** Replicates position for the given symmetry operator (LxMxN). */
  public final int[] replicatesVector;

  /**
   * The SymOp constructor using a rotation matrix and translation vector.
   *
   * @param rot The rotation matrix.
   * @param tr The translation vector.
   */
  public SymOp(double[][] rot, double[] tr) {
    this.rot = rot;
    this.tr = tr;
    replicatesVector = new int[]{0, 0, 0};
  }

  /**
   * The SymOp constructor using a rotation matrix and translation vector.
   *
   * @param rot The rotation matrix.
   * @param tr The translation vector.
   * @param replicatesVector Describes symmetry operators location within replicates crystal.
   */
  public SymOp(double[][] rot, double[] tr, int[] replicatesVector) {
    this.rot = rot;
    this.tr = tr;
    this.replicatesVector = replicatesVector;
  }


  /**
   * The SymOp constructor using a 4x4 matrix.
   *
   * @param m The rotation matrix and translation vector as a 4x4 matrix.
   */
  public SymOp(double[][] m) {
    this.rot = new double[3][3];
    rot[0][0] = m[0][0];
    rot[0][1] = m[0][1];
    rot[0][2] = m[0][2];
    rot[1][0] = m[1][0];
    rot[1][1] = m[1][1];
    rot[1][2] = m[1][2];
    rot[2][0] = m[2][0];
    rot[2][1] = m[2][1];
    rot[2][2] = m[2][2];

    this.tr = new double[3];
    tr[0] = m[0][3] / m[3][3];
    tr[1] = m[1][3] / m[3][3];
    tr[2] = m[2][3] / m[3][3];

    replicatesVector = new int[]{0, 0, 0};
  }

  /**
   * symPhaseShift
   *
   * @param hkl an array of double.
   * @return a double.
   */
  public double symPhaseShift(double[] hkl) {
    // Apply translation
    return -2.0 * PI * (hkl[0] * tr[0] + hkl[1] * tr[1] + hkl[2] * tr[2]);
  }

  /**
   * symPhaseShift
   *
   * @param hkl a {@link HKL} object.
   * @return a double.
   */
  public double symPhaseShift(HKL hkl) {
    // Apply translation
    return -2.0 * PI * (hkl.getH() * tr[0] + hkl.getK() * tr[1] + hkl.getL() * tr[2]);
  }

  /**
   * Return the SymOp as a 4x4 matrix.
   *
   * @return A 4x4 matrix representation of the SymOp.
   */
  public double[][] asMatrix() {
    return new double[][] {
        {rot[0][0], rot[0][1], rot[0][2], tr[0]},
        {rot[1][0], rot[1][1], rot[1][2], tr[1]},
        {rot[2][0], rot[2][1], rot[2][2], tr[2]},
        {0.0, 0.0, 0.0, 1.0}};
  }

  /**
   * Return the combined SymOp that is equivalent to first applying <code>this</code> SymOp and then
   * the argument. Note: Applied as rotation then translation.
   * <code>X' = S_arg(S_this(X))</code>
   * <code>X' = S_combined(X)</code>
   *
   * @param symOp The SymOp to append to <code>this</code> SymOp.
   * @return The combined SymOp.
   */
  public SymOp append(SymOp symOp) {
    return new SymOp(mat4Mat4(symOp.asMatrix(), asMatrix()));
  }

  /**
   * Return the combined SymOp that is equivalent to first applying the argument and then
   * <code>this</code> SymOp. Note: Applied as rotation then translation.
   * <code>X' = S_this(S_arg(X))</code>
   * <code>X' = S_combined(X)</code>
   *
   * @param symOp The SymOp to prepend to <code>this</code> SymOp.
   * @return The combined SymOp.
   */
  public SymOp prepend(SymOp symOp) {
    return new SymOp(mat4Mat4(asMatrix(), symOp.asMatrix()));
  }

  /**
   * Return the combined SymOp that is equivalent to first applying symOp1 and then SymOp2. Note:
   * Applied as rotation then translation.
   * <p>
   * <code>X' = S_2(S_1(X))</code>
   * <code>X' = S_combined(X)</code>
   *
   * @param symOp1 The fist SymOp.
   * @param symOp2 The second SymOp.
   * @return The combined SymOp.
   */
  public static SymOp combineSymOps(SymOp symOp1, SymOp symOp2) {
    return new SymOp(mat4Mat4(symOp2.asMatrix(), symOp1.asMatrix()));
  }

  /**
   * Apply a Cartesian symmetry operator to an array of Cartesian coordinates. If the arrays x, y or
   * z are null or not of length n, the method returns immediately. If mateX, mateY or mateZ are null
   * or not of length n, the method returns immediately.
   *
   * @param n Number of atoms.
   * @param x Input cartesian x-coordinates.
   * @param y Input cartesian y-coordinates.
   * @param z Input cartesian z-coordinates.
   * @param mateX Output cartesian x-coordinates.
   * @param mateY Output cartesian y-coordinates.
   * @param mateZ Output cartesian z-coordinates.
   * @param symOp The cartesian symmetry operator.
   */
  public static void applyCartSymOp(int n, double[] x, double[] y, double[] z,
      double[] mateX, double[] mateY, double[] mateZ, SymOp symOp) {
    if (x == null || y == null || z == null) {
      throw new IllegalArgumentException("The input arrays x, y and z must not be null.");
    }
    if (x.length < n || y.length < n || z.length < n) {
      throw new IllegalArgumentException("The input arrays x, y and z must be of length n: " + n);
    }
    if (mateX == null || mateY == null || mateZ == null) {
      throw new IllegalArgumentException("The output arrays mateX, mateY and mateZ must not be null.");
    }
    if (mateX.length < n || mateY.length < n || mateZ.length < n) {
      throw new IllegalArgumentException("The output arrays mateX, mateY and mateZ must be of length n: " + n);
    }

    final double[][] rot = symOp.rot;
    final double[] trans = symOp.tr;

    final double rot00 = rot[0][0];
    final double rot10 = rot[1][0];
    final double rot20 = rot[2][0];
    final double rot01 = rot[0][1];
    final double rot11 = rot[1][1];
    final double rot21 = rot[2][1];
    final double rot02 = rot[0][2];
    final double rot12 = rot[1][2];
    final double rot22 = rot[2][2];
    final double t0 = trans[0];
    final double t1 = trans[1];
    final double t2 = trans[2];
    for (int i = 0; i < n; i++) {
      double xc = x[i];
      double yc = y[i];
      double zc = z[i];
      // Apply Symmetry Operator.
      mateX[i] = rot00 * xc + rot01 * yc + rot02 * zc + t0;
      mateY[i] = rot10 * xc + rot11 * yc + rot12 * zc + t1;
      mateZ[i] = rot20 * xc + rot21 * yc + rot22 * zc + t2;
    }
  }

  /**
   * Apply a cartesian symmetry operator to an array of coordinates.
   *
   * @param xyz Input  cartesian coordinates.
   * @param mate Symmetry mate  cartesian coordinates.
   * @param symOp The cartesian symmetry operator.
   */
  public static void applyCartesianSymOp(double[] xyz, double[] mate, SymOp symOp) {
    applyCartesianSymOp(xyz, mate, symOp, null);
  }

  /**
   * Apply a cartesian symmetry operator to an array of coordinates.
   *
   * @param xyz Input  cartesian coordinates.
   * @param mate Symmetry mate  cartesian coordinates.
   * @param symOp The cartesian symmetry operator.
   * @param mask Only apply the SymOp if the per atom mask is true.
   */
  public static void applyCartesianSymOp(double[] xyz, double[] mate, SymOp symOp, boolean[] mask) {
    var rot = symOp.rot;
    var trans = symOp.tr;

    assert (xyz.length % 3 == 0);
    assert (xyz.length == mate.length);
    int len = xyz.length / 3;

    // Load the SymOp into local variables.
    var rot00 = rot[0][0];
    var rot10 = rot[1][0];
    var rot20 = rot[2][0];
    var rot01 = rot[0][1];
    var rot11 = rot[1][1];
    var rot21 = rot[2][1];
    var rot02 = rot[0][2];
    var rot12 = rot[1][2];
    var rot22 = rot[2][2];
    var tx = trans[0];
    var ty = trans[1];
    var tz = trans[2];

    if (mask == null) {
      for (int i = 0; i < len; i++) {
        int index = i * 3;
        var xc = xyz[index + XX];
        var yc = xyz[index + YY];
        var zc = xyz[index + ZZ];
        // Apply Symmetry Operator.
        mate[index + XX] = rot00 * xc + rot01 * yc + rot02 * zc + tx;
        mate[index + YY] = rot10 * xc + rot11 * yc + rot12 * zc + ty;
        mate[index + ZZ] = rot20 * xc + rot21 * yc + rot22 * zc + tz;
      }
    } else {
      for (int i = 0; i < len; i++) {
        int index = i * 3;
        var xc = xyz[index + XX];
        var yc = xyz[index + YY];
        var zc = xyz[index + ZZ];
        if (mask[i]) {
          // Apply Symmetry Operator.
          mate[index + XX] = rot00 * xc + rot01 * yc + rot02 * zc + tx;
          mate[index + YY] = rot10 * xc + rot11 * yc + rot12 * zc + ty;
          mate[index + ZZ] = rot20 * xc + rot21 * yc + rot22 * zc + tz;
        } else {
          mate[index + XX] = xc;
          mate[index + YY] = yc;
          mate[index + ZZ] = zc;
        }
      }
    }
  }

  /**
   * Apply a fractional symmetry operator to one set of coordinates.
   *
   * @param xyz Input fractional coordinates.
   * @param mate Symmetry mate fractional coordinates.
   * @param symOp The fractional symmetry operator.
   */
  public static void applyFracSymOp(double[] xyz, double[] mate, SymOp symOp) {
    var rot = symOp.rot;
    var trans = symOp.tr;
    var xf = xyz[0];
    var yf = xyz[1];
    var zf = xyz[2];
    // Apply Symmetry Operator.
    mate[0] = rot[0][0] * xf + rot[0][1] * yf + rot[0][2] * zf + trans[0];
    mate[1] = rot[1][0] * xf + rot[1][1] * yf + rot[1][2] * zf + trans[1];
    mate[2] = rot[2][0] * xf + rot[2][1] * yf + rot[2][2] * zf + trans[2];
  }

  /**
   * Apply a symmetry operator to one set of coordinates.
   *
   * @param h Input coordinates.
   * @param k Input coordinates.
   * @param l Input coordinates.
   * @param mate Symmetry mate coordinates.
   * @param symOp The symmetry operator.
   * @param nx number of unit cell translations
   * @param ny number of unit cell translations
   * @param nz number of unit cell translations
   */
  public static void applySymOp(int h, int k, int l, int[] mate, SymOp symOp, int nx, int ny,
      int nz) {
    var rot = symOp.rot;
    var trans = symOp.tr;
    // Apply Symmetry Operator.
    mate[0] =
        (int) rot[0][0] * h + (int) rot[0][1] * k + (int) rot[0][2] * l + (int) rint(nx * trans[0]);
    mate[1] =
        (int) rot[1][0] * h + (int) rot[1][1] * k + (int) rot[1][2] * l + (int) rint(ny * trans[1]);
    mate[2] =
        (int) rot[2][0] * h + (int) rot[2][1] * k + (int) rot[2][2] * l + (int) rint(nz * trans[2]);
    mate[0] = mod(mate[0], nx);
    mate[1] = mod(mate[1], ny);
    mate[2] = mod(mate[2], nz);
  }

  /**
   * Apply a symmetry operator to one HKL.
   *
   * @param hkl Input HKL.
   * @param mate Symmetry mate HKL.
   * @param symOp The symmetry operator.
   */
  public static void applySymRot(HKL hkl, HKL mate, SymOp symOp) {
    var rot = symOp.rot;
    double h = hkl.getH();
    double k = hkl.getK();
    double l = hkl.getL();
    double hs = rot[0][0] * h + rot[0][1] * k + rot[0][2] * l;
    double ks = rot[1][0] * h + rot[1][1] * k + rot[1][2] * l;
    double ls = rot[2][0] * h + rot[2][1] * k + rot[2][2] * l;
    // Convert back to HKL
    mate.setH((int) rint(hs));
    mate.setK((int) rint(ks));
    mate.setL((int) rint(ls));
  }

  /**
   * Apply a Cartesian symmetry rotation to an array of Cartesian coordinates. The length of xyz must
   * be divisible by 3 and mate must have the same length.
   *
   * @param xyz Input cartesian x, y, z-coordinates.
   * @param mate Output cartesian x, y, z-coordinates.
   * @param symOp The fractional symmetry operator.
   */
  public static void applyCartesianSymRot(double[] xyz, double[] mate, SymOp symOp) {
    applyCartesianSymRot(xyz, mate, symOp, null);
  }

  /**
   * Apply a Cartesian symmetry rotation to an array of Cartesian coordinates. The length of xyz must
   * be divisible by 3 and mate must have the same length.
   *
   * @param xyz Input cartesian x, y, z-coordinates.
   * @param mate Output cartesian x, y, z-coordinates.
   * @param symOp The fractional symmetry operator.
   * @param mask Only apply the SymOp if the per atom mask is true.
   */
  public static void applyCartesianSymRot(double[] xyz, double[] mate, SymOp symOp, boolean[] mask) {
    int l = xyz.length;
    int n = l / 3;
    assert (l % 3 == 0);
    assert (mate.length == l);

    // Load the rotation matrix
    var rot = symOp.rot;
    var rot00 = rot[0][0];
    var rot10 = rot[1][0];
    var rot20 = rot[2][0];
    var rot01 = rot[0][1];
    var rot11 = rot[1][1];
    var rot21 = rot[2][1];
    var rot02 = rot[0][2];
    var rot12 = rot[1][2];
    var rot22 = rot[2][2];

    if (mask == null) {
      for (int i = 0; i < n; i++) {
        int index = i * 3;
        var xi = xyz[index + XX];
        var yi = xyz[index + YY];
        var zi = xyz[index + ZZ];
        // Apply Symmetry Operator.
        mate[index + XX] = rot00 * xi + rot01 * yi + rot02 * zi;
        mate[index + YY] = rot10 * xi + rot11 * yi + rot12 * zi;
        mate[index + ZZ] = rot20 * xi + rot21 * yi + rot22 * zi;
      }
    } else {
      for (int i = 0; i < n; i++) {
        int index = i * 3;
        var xi = xyz[index + XX];
        var yi = xyz[index + YY];
        var zi = xyz[index + ZZ];
        if (mask[i]) {
          // Apply Symmetry Operator.
          mate[index + XX] = rot00 * xi + rot01 * yi + rot02 * zi;
          mate[index + YY] = rot10 * xi + rot11 * yi + rot12 * zi;
          mate[index + ZZ] = rot20 * xi + rot21 * yi + rot22 * zi;
        } else {
          mate[index + XX] = xi;
          mate[index + YY] = yi;
          mate[index + ZZ] = zi;
        }
      }
    }
  }

  /**
   * Apply a transpose rotation symmetry operator to one HKL.
   *
   * @param hkl Input HKL.
   * @param mate Symmetry mate HKL.
   * @param symOp The symmetry operator.
   */
  public static void applyTransSymRot(HKL hkl, HKL mate, SymOp symOp) {
    double[][] rot = symOp.rot;
    double h = hkl.getH();
    double k = hkl.getK();
    double l = hkl.getL();
    // Apply transpose Symmetry Operator.
    double hs = rot[0][0] * h + rot[1][0] * k + rot[2][0] * l;
    double ks = rot[0][1] * h + rot[1][1] * k + rot[2][1] * l;
    double ls = rot[0][2] * h + rot[1][2] * k + rot[2][2] * l;
    // Convert back to HKL
    mate.setH((int) rint(hs));
    mate.setK((int) rint(ks));
    mate.setL((int) rint(ls));
  }

  /**
   * Print a Sym Op matrix as a continued line string.
   *
   * @param symOp Symmetry operation to print.
   * @return Continued line string.
   */
  public static String asMatrixString(SymOp symOp) {
    double[][] values = symOp.asMatrix();
    return format("""
                     %14.8f %14.8f %14.8f \\
                               %14.8f %14.8f %14.8f \\
                               %14.8f %14.8f %14.8f \\
                               %14.8f %14.8f %14.8f \
                    """, values[0][0],
            values[0][1], values[0][2], values[1][0], values[1][1], values[1][2], values[2][0],
            values[2][1], values[2][2], values[0][3], values[1][3], values[2][3]);
  }

  /**
   * Invert a symmetry operator.
   *
   * @param symOp Original symmetry operator of which the inverse is desired.
   * @return SymOp The inverse symmetry operator of the one supplied.
   */
  public static SymOp invertSymOp(SymOp symOp) {
    var tr = symOp.tr;
    var rot = symOp.rot;
    var inv = mat3Inverse(rot);
    return new SymOp(inv,
        new double[] {-dot(inv[0], tr), -dot(inv[1], tr), -dot(inv[2], tr)});
  }

  /**
   * Create a SymOp from an input String.
   *
   * @param s Input <code>String</code> containing 12 white-space delimited double values.
   * @return The SymOp.
   */
  public static SymOp parse(String s) {
    String[] tokens = s.split(" +");
    if (tokens.length < 12) {
      return null;
    }
    return new SymOp(new double[][] {
        {parseDouble(tokens[0]), parseDouble(tokens[1]), parseDouble(tokens[2])},
        {parseDouble(tokens[3]), parseDouble(tokens[4]), parseDouble(tokens[5])},
        {parseDouble(tokens[6]), parseDouble(tokens[7]), parseDouble(tokens[8])}},
        new double[] {parseDouble(tokens[9]), parseDouble(tokens[10]), parseDouble(tokens[11])});
  }

  /**
   * Generate a random Cartesian Symmetry Operator.
   *
   * @param scalar The range of translations will be from -scalar/2 .. scalar/2.
   * @return A Cartesian SymOp with a random rotation and translation.
   */
  public static SymOp randomSymOpFactory(double scalar) {
    double[] tr = {scalar * (random() - 0.5), scalar * (random() - 0.5), scalar * (random() - 0.5)};
    return randomSymOpFactory(tr);
  }

  /**
   * Generate a random Cartesian Symmetry Operator.
   *
   * <p>The random rotation matrix is derived from: Arvo, James (1992), "Fast random rotation
   * matrices", in David Kirk, Graphics Gems III, San Diego: Academic Press Professional, pp.
   * 117–120, ISBN 978-0-12-409671-4
   *
   * @param tr The translations to apply.
   * @return A Cartesian SymOp with a random rotation and translation.
   */
  public static SymOp randomSymOpFactory(double[] tr) {
    double[][] rot = new double[3][3];
    double PI2 = 2.0 * PI;
    double[] x = new double[3];
    x[0] = random();
    x[1] = random();
    x[2] = random();
    /* Rotation about the pole (Z).      */
    double theta = x[0] * PI2;
    /* For direction of pole deflection. */
    double phi = x[1] * PI2;
    /* For magnitude of pole deflection. */
    double z = x[2] * 2.0;
    /*
     Compute a vector V used for distributing points over the sphere via
     the reflection I - V Transpose(V). This formulation of V will
     guarantee that if x[1] and x[2] are uniformly distributed, the
     reflected points will be uniform on the sphere. Note that V has
     length sqrt(2) to eliminate the 2 in the Householder matrix.
    */
    double r = sqrt(z);
    double Vx = sin(phi) * r;
    double Vy = cos(phi) * r;
    double Vz = sqrt(2.0 - z);
    /*
     Compute the row vector S = Transpose(V) * R, where R is a simple
     rotation by theta about the z-axis. No need to compute Sz since it's
     just Vz.
    */
    double st = sin(theta);
    double ct = cos(theta);
    double Sx = Vx * ct - Vy * st;
    double Sy = Vx * st + Vy * ct;

    // Construct the rotation matrix ( V Transpose(V) - I ) R, which is equivalent to V S - R.
    rot[0][0] = Vx * Sx - ct;
    rot[0][1] = Vx * Sy - st;
    rot[0][2] = Vx * Vz;
    rot[1][0] = Vy * Sx + st;
    rot[1][1] = Vy * Sy - ct;
    rot[1][2] = Vy * Vz;
    rot[2][0] = Vz * Sx;
    rot[2][1] = Vz * Sy;
    // This equals Vz * Vz - 1.0
    rot[2][2] = 1.0 - z;

    return new SymOp(rot, tr);
  }

  /**
   * trtoString
   *
   * @param tr a double.
   * @return a {@link java.lang.String} object.
   */
  private static String trtoString(double tr) {
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

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder(" Rotation operator:\n");
    sb.append(format(" [[%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]\n  [%4.1f,%4.1f,%4.1f]]\n",
        rot[0][0], rot[0][1], rot[0][2],
        rot[1][0], rot[1][1], rot[1][2],
        rot[2][0], rot[2][1], rot[2][2]));
    sb.append(" Translation:\n");
    sb.append(format(" [%4.2f,%4.2f,%4.2f]", tr[0], tr[1], tr[2]));
    return sb.toString();
  }

  /**
   * Print the symmetry operator with double precision.
   *
   * @return String of rotation/translation with double precision.
   */
  public String toStringPrecise() {
    StringBuilder sb = new StringBuilder(" Rotation operator:\n");
    sb.append(format(
        " [[%18.16e,%18.16e,%18.16e]\n  [%18.16e,%18.16e,%18.16e]\n  [%18.16e,%18.16e,%18.16e]]\n",
        rot[0][0], rot[0][1], rot[0][2],
        rot[1][0], rot[1][1], rot[1][2],
        rot[2][0], rot[2][1], rot[2][2]));
    sb.append(" Translation:\n");
    sb.append(format(" [%18.16e,%18.16e,%18.16e]", tr[0], tr[1], tr[2]));
    return sb.toString();
  }

  /**
   * toXYZString
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
}
