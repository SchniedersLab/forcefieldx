// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.numerics.multipole;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;

import static java.util.Arrays.copyOf;

/**
 * The QIFrame class defines a quasi-internal frame between two atoms.
 * <p>
 * The Z-axis of the QI frame is defined as the vector between the two atoms. The X- and Y-axes are
 * then defined to create a right-handed coordinate system.
 * <p>
 * A rotation matrix from the global frame to the QI frame is constructed, and vice versa. Using the
 * rotation matrices, methods are provided to rotate vectors and multipoles between the two frames.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class QIFrameSIMD {

  // Rotation Matrix from Global to QI.
  private DoubleVector r00Vec, r01Vec, r02Vec;
  private DoubleVector r10Vec, r11Vec, r12Vec;
  private DoubleVector r20Vec, r21Vec, r22Vec;

  // Rotation Matrix from QI to Global.
  private DoubleVector ir00Vec, ir01Vec, ir02Vec;
  private DoubleVector ir10Vec, ir11Vec, ir12Vec;
  private DoubleVector ir20Vec, ir21Vec, ir22Vec;

  /**
   * QIFrame constructor
   * <p>
   * (dx = 0, dy = 0, dz = 1).
   */
  public QIFrameSIMD() {
    DoubleVector zero = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 0.0);
    setQIVector(zero, zero, zero);
  }

  /**
   * QIFrame constructor.
   */
  public QIFrameSIMD(DoubleVector dx, DoubleVector dy, DoubleVector dz) {
    setQIVector(dx, dy, dz);
  }

  /**
   * QIFrame constructor.
   *
   * @param r Separation along each axis.
   */
  public QIFrameSIMD(DoubleVector[] r) {
    setQIVector(r[0], r[1], r[2]);
  }

  /**
   * Update the QIFrame rotation matrix.
   *
   * @param r Separation along each axis.
   */
  public void setQIVector(DoubleVector[] r) {
    setQIVector(r[0], r[1], r[2]);
  }

  /**
   * Update the QIFrame rotation matrix.
   *
   * @param dx Separation along the x-axis.
   * @param dy Separation along the y-axis.
   * @param dz Separation along the z-axis.
   */
  public final void setQIVector(DoubleVector dx, DoubleVector dy, DoubleVector dz) {
    // The QI Z-axis is along the separation vector.
    DoubleVector[] zAxis = {dx, dy, dz};

    // The "guess" for the QI X-axis cannot be along the QI separation vector.
    DoubleVector[] xAxis = copyOf(zAxis, 3);

    // Check if the QI separation vector is the global X-axis.
    VectorMask<Double> qiVectorIsXAxis = dz.eq(0.0).and(dy.eq(0.0));
    // If the QI separation vector is the global X-axis, add 1 to the Y-component of the QI X-axis.
    xAxis[1] = xAxis[1].add(1.0, qiVectorIsXAxis);
    // Otherwise, add 1 to the X-component of the QI X-axis.
    xAxis[0] = xAxis[0].add(1.0, qiVectorIsXAxis.not());

    // Normalize the QI Z-axis.
    normalizeVec(zAxis, zAxis);
    ir02Vec = zAxis[0];
    ir12Vec = zAxis[1];
    ir22Vec = zAxis[2];

    // Finalize the QI X-axis.
    DoubleVector dot = dotVec(xAxis, zAxis);
    scaleVec(zAxis, dot, zAxis);
    subVec(xAxis, zAxis, xAxis);
    normalizeVec(xAxis, xAxis);
    ir00Vec = xAxis[0];
    ir10Vec = xAxis[1];
    ir20Vec = xAxis[2];

    // Cross the QI X-axis and Z-axis to get the QI Y-axis.
    ir01Vec = ir20Vec.mul(ir12Vec).sub(ir10Vec.mul(ir22Vec));
    ir11Vec = ir00Vec.mul(ir22Vec).sub(ir20Vec.mul(ir02Vec));
    ir21Vec = ir10Vec.mul(ir02Vec).sub(ir00Vec.mul(ir12Vec));

    // Set the forward elements as the transpose of the inverse matrix.
    r00Vec = ir00Vec;
    r11Vec = ir11Vec;
    r22Vec = ir22Vec;
    r01Vec = ir10Vec;
    r02Vec = ir20Vec;
    r10Vec = ir01Vec;
    r12Vec = ir21Vec;
    r20Vec = ir02Vec;
    r21Vec = ir12Vec;
  }

  /**
   * Compute the dot product of two vectors.
   *
   * @param a The first vector.
   * @param b The second vector.
   * @return The dot product of the two vectors.
   */
  public static DoubleVector dotVec(DoubleVector[] a, DoubleVector[] b) {
    return a[0].fma(b[0], a[1].fma(b[1], a[2].mul(b[2])));
  }

  /**
   * Scale a vector by a scalar.
   *
   * @param a     The vector to scale.
   * @param scale The scalar to multiply the vector by.
   * @param ret   The vector to store the result in.
   * @return The scaled vector.
   */
  public static DoubleVector[] scaleVec(DoubleVector[] a, DoubleVector scale, DoubleVector[] ret) {
    ret[0] = a[0].mul(scale);
    ret[1] = a[1].mul(scale);
    ret[2] = a[2].mul(scale);
    return ret;
  }

  /**
   * Subtract two vectors.
   *
   * @param a   The first vector.
   * @param b   The second vector.
   * @param ret The vector to store the result in.
   * @return The difference of the two vectors.
   */
  public static DoubleVector[] subVec(DoubleVector[] a, DoubleVector[] b, DoubleVector[] ret) {
    ret[0] = a[0].sub(b[0]);
    ret[1] = a[1].sub(b[1]);
    ret[2] = a[2].sub(b[2]);
    return ret;
  }

  /**
   * Normalize a vector.
   *
   * @param n   The vector to normalize.
   * @param ret The vector to store the result in.
   * @return The normalized vector.
   */
  public static DoubleVector[] normalizeVec(DoubleVector[] n, DoubleVector[] ret) {
    return scaleVec(n, DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 1.0).div(dotVec(n, n).sqrt()), ret);
  }

  /**
   * Update the QIFrame rotation matrix and rotate the multipoles.
   *
   * @param r  Separation along each axis.
   * @param mI PolarizableMultipole for site I.
   * @param mK PolarizableMultipole for site K.
   */
  public void setAndRotate(DoubleVector[] r, PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    setAndRotate(r[0], r[1], r[2], mI, mK);
  }

  /**
   * Update the QIFrame rotation matrix and rotate the multipoles.
   *
   * @param dx Separation along the x-axis.
   * @param dy Separation along the y-axis.
   * @param dz Separation along the z-axis.
   * @param mI PolarizableMultipole for site I.
   * @param mK PolarizableMultipole for site K.
   */
  public void setAndRotate(DoubleVector dx, DoubleVector dy, DoubleVector dz,
                           PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    setQIVector(dx, dy, dz);
    rotatePolarizableMultipole(mI);
    rotatePolarizableMultipole(mK);
  }

  /**
   * Rotate the permanent multipole and induced dipole.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotatePolarizableMultipole(PolarizableMultipoleSIMD m) {
    rotatePermanentMultipole(m);
    rotateInducedDipoles(m);
  }

  /**
   * Rotate the permanent multipole.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotatePermanentMultipole(PolarizableMultipoleSIMD m) {
    // Rotate the permanent dipole.
    DoubleVector dx = m.dx;
    DoubleVector dy = m.dy;
    DoubleVector dz = m.dz;
    m.dx = r00Vec.mul(dx).add(r01Vec.mul(dy)).add(r02Vec.mul(dz));
    m.dy = r10Vec.mul(dx).add(r11Vec.mul(dy)).add(r12Vec.mul(dz));
    m.dz = r20Vec.mul(dx).add(r21Vec.mul(dy)).add(r22Vec.mul(dz));

    // Rotate the quadrupole.
    DoubleVector qxx = m.qxx;
    DoubleVector qyy = m.qyy;
    DoubleVector qzz = m.qzz;
    // The Multipole class stores 2.0 times the off-diagonal components.
    DoubleVector qxy = m.qxy.mul(0.5);
    DoubleVector qxz = m.qxz.mul(0.5);
    DoubleVector qyz = m.qyz.mul(0.5);

    m.qxx = r00Vec.mul(r00Vec.mul(qxx).add(r01Vec.mul(qxy)).add(r02Vec.mul(qxz)))
        .add(r01Vec.mul((r00Vec.mul(qxy).add(r01Vec.mul(qyy)).add(r02Vec.mul(qyz)))))
        .add(r02Vec.mul((r00Vec.mul(qxz).add(r01Vec.mul(qyz)).add(r02Vec.mul(qzz)))));

    m.qxy = r00Vec.mul(r10Vec.mul(qxx).add(r11Vec.mul(qxy)).add(r12Vec.mul(qxz)))
        .add(r01Vec.mul((r10Vec.mul(qxy).add(r11Vec.mul(qyy)).add(r12Vec.mul(qyz)))))
        .add(r02Vec.mul((r10Vec.mul(qxz).add(r11Vec.mul(qyz)).add(r12Vec.mul(qzz)))));
    m.qxy = m.qxy.mul(2.0);

    m.qxz = r00Vec.mul(r20Vec.mul(qxx).add(r21Vec.mul(qxy)).add(r22Vec.mul(qxz)))
        .add(r01Vec.mul((r20Vec.mul(qxy).add(r21Vec.mul(qyy)).add(r22Vec.mul(qyz)))))
        .add(r02Vec.mul((r20Vec.mul(qxz).add(r21Vec.mul(qyz)).add(r22Vec.mul(qzz)))));
    m.qxz = m.qxz.mul(2.0);

    m.qyy = r10Vec.mul(r10Vec.mul(qxx).add(r11Vec.mul(qxy)).add(r12Vec.mul(qxz)))
        .add(r11Vec.mul((r10Vec.mul(qxy).add(r11Vec.mul(qyy)).add(r12Vec.mul(qyz)))))
        .add(r12Vec.mul((r10Vec.mul(qxz).add(r11Vec.mul(qyz)).add(r12Vec.mul(qzz)))));

    m.qyz = r10Vec.mul(r20Vec.mul(qxx).add(r21Vec.mul(qxy)).add(r22Vec.mul(qxz)))
        .add(r11Vec.mul((r20Vec.mul(qxy).add(r21Vec.mul(qyy)).add(r22Vec.mul(qyz)))))
        .add(r12Vec.mul((r20Vec.mul(qxz).add(r21Vec.mul(qyz)).add(r22Vec.mul(qzz)))));
    m.qyz = m.qyz.mul(2.0);

    m.qzz = r20Vec.mul(r20Vec.mul(qxx).add(r21Vec.mul(qxy)).add(r22Vec.mul(qxz)))
        .add(r21Vec.mul((r20Vec.mul(qxy).add(r21Vec.mul(qyy)).add(r22Vec.mul(qyz)))))
        .add(r22Vec.mul((r20Vec.mul(qxz).add(r21Vec.mul(qyz)).add(r22Vec.mul(qzz)))));

  }

  /**
   * Rotate the induced dipoles components.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotateInducedDipoles(PolarizableMultipoleSIMD m) {
    DoubleVector dx = m.ux;
    DoubleVector dy = m.uy;
    DoubleVector dz = m.uz;
    m.ux = r00Vec.mul(dx).add(r01Vec.mul(dy)).add(r02Vec.mul(dz));
    m.uy = r10Vec.mul(dx).add(r11Vec.mul(dy)).add(r12Vec.mul(dz));
    m.uz = r20Vec.mul(dx).add(r21Vec.mul(dy)).add(r22Vec.mul(dz));
    dx = m.px;
    dy = m.py;
    dz = m.pz;
    m.px = r00Vec.mul(dx).add(r01Vec.mul(dy)).add(r02Vec.mul(dz));
    m.py = r10Vec.mul(dx).add(r11Vec.mul(dy)).add(r12Vec.mul(dz));
    m.pz = r20Vec.mul(dx).add(r21Vec.mul(dy)).add(r22Vec.mul(dz));
    // Set update the averaged induced and chain rule terms.
    m.sx = m.ux.add(m.px).mul(0.5);
    m.sy = m.uy.add(m.py).mul(0.5);
    m.sz = m.uz.add(m.pz).mul(0.5);
  }

  /**
   * Rotate a vector in the QI frame into the global frame.
   *
   * @param v The vector to rotate (in-place).
   */
  public void toGlobal(DoubleVector[] v) {
    DoubleVector vx = v[0];
    DoubleVector vy = v[1];
    DoubleVector vz = v[2];
    // X-component.
    v[0] = ir00Vec.mul(vx).add(ir01Vec.mul(vy)).add(ir02Vec.mul(vz));
    // Y-component.
    v[1] = ir10Vec.mul(vx).add(ir11Vec.mul(vy)).add(ir12Vec.mul(vz));
    // Z-component.
    v[2] = ir20Vec.mul(vx).add(ir21Vec.mul(vy)).add(ir22Vec.mul(vz));
  }
}
