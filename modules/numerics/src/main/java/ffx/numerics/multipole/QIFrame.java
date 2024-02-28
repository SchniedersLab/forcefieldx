// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.DoubleMath.normalize;
import static ffx.numerics.math.DoubleMath.scale;
import static ffx.numerics.math.DoubleMath.sub;
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
public class QIFrame {

  // Rotation Matrix from Global to QI.
  private double r00, r01, r02;
  private double r10, r11, r12;
  private double r20, r21, r22;

  // Rotation Matrix from QI to Global.
  private double ir00, ir01, ir02;
  private double ir10, ir11, ir12;
  private double ir20, ir21, ir22;

  /**
   * QIFrame constructor
   * <p>
   * (dx = 0, dy = 0, dz = 1).
   */
  public QIFrame() {
    setQIVector(0.0, 0.0, 1.0);
  }

  /**
   * QIFrame constructor.
   *
   * @param dx Separation along the x-axis.
   * @param dy Separation along the y-axis.
   * @param dz Separation along the z-axis.
   */
  public QIFrame(double dx, double dy, double dz) {
    setQIVector(dx, dy, dz);
  }

  /**
   * QIFrame constructor.
   *
   * @param r Separation along each axis.
   */
  public QIFrame(double[] r) {
    setQIVector(r[0], r[1], r[2]);
  }

  /**
   * Update the QIFrame rotation matrix.
   *
   * @param r Separation along each axis.
   */
  public void setQIVector(double[] r) {
    setQIVector(r[0], r[1], r[2]);
  }

  /**
   * Update the QIFrame rotation matrix.
   *
   * @param dx Separation along the x-axis.
   * @param dy Separation along the y-axis.
   * @param dz Separation along the z-axis.
   */
  public final void setQIVector(double dx, double dy, double dz) {
    // The QI Z-axis is along the separation vector.
    double[] zAxis = {dx, dy, dz};

    // The "guess" for the QI X-axis cannot be along the QI separation vector.
    double[] xAxis = copyOf(zAxis, 3);
    if (dy != 0.0 || dz != 0.0) {
      // The QI separation vector is not along the global X-axis.
      // Set the QI X-axis "guess" to the QI Z-axis plus add 1 to X-component.
      xAxis[0] += 1.0;
    } else {
      // The QI separation vector is the global X-axis.
      // Set the QI X-axis "guess" to the QI Z-axis plus add 1 to Y-component.
      xAxis[1] += 1.0;
    }

    // Normalize the QI Z-axis.
    normalize(zAxis, zAxis);
    ir02 = zAxis[0];
    ir12 = zAxis[1];
    ir22 = zAxis[2];

    // Finalize the QI X-axis.
    double dot = dot(xAxis, zAxis);
    scale(zAxis, dot, zAxis);
    sub(xAxis, zAxis, xAxis);
    normalize(xAxis, xAxis);
    ir00 = xAxis[0];
    ir10 = xAxis[1];
    ir20 = xAxis[2];

    // Cross the QI X-axis and Z-axis to get the QI Y-axis.
    ir01 = ir20 * ir12 - ir10 * ir22;
    ir11 = ir00 * ir22 - ir20 * ir02;
    ir21 = ir10 * ir02 - ir00 * ir12;

    // Set the forward elements as the transpose of the inverse matrix.
    r00 = ir00;
    r11 = ir11;
    r22 = ir22;
    r01 = ir10;
    r02 = ir20;
    r10 = ir01;
    r12 = ir21;
    r20 = ir02;
    r21 = ir12;
  }

  /**
   * Update the QIFrame rotation matrix and rotate the multipoles.
   *
   * @param r Separation along each axis.
   * @param mI PolarizableMultipole for site I.
   * @param mK PolarizableMultipole for site K.
   */
  public void setAndRotate(double[] r, PolarizableMultipole mI, PolarizableMultipole mK) {
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
  public void setAndRotate(double dx, double dy, double dz,
      PolarizableMultipole mI, PolarizableMultipole mK) {
    setQIVector(dx, dy, dz);
    rotatePolarizableMultipole(mI);
    rotatePolarizableMultipole(mK);
  }

  /**
   * Rotate the permanent multipole and induced dipole.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotatePolarizableMultipole(PolarizableMultipole m) {
    rotatePermanentMultipole(m);
    rotateInducedDipoles(m);
  }

  /**
   * Rotate the permanent multipole.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotatePermanentMultipole(PolarizableMultipole m) {
    // Rotate the permanent dipole.
    double dx = m.dx;
    double dy = m.dy;
    double dz = m.dz;
    m.dx = r00 * dx + r01 * dy + r02 * dz;
    m.dy = r10 * dx + r11 * dy + r12 * dz;
    m.dz = r20 * dx + r21 * dy + r22 * dz;

    // Rotate the quadrupole.
    double qxx = m.qxx;
    double qyy = m.qyy;
    double qzz = m.qzz;
    // The Multipole class stores 2.0 times the off-diagonal components.
    double qxy = m.qxy * 0.5;
    double qxz = m.qxz * 0.5;
    double qyz = m.qyz * 0.5;

    m.qxx = r00 * (r00 * qxx + r01 * qxy + r02 * qxz)
        + r01 * (r00 * qxy + r01 * qyy + r02 * qyz)
        + r02 * (r00 * qxz + r01 * qyz + r02 * qzz);

    m.qxy = r00 * (r10 * qxx + r11 * qxy + r12 * qxz)
        + r01 * (r10 * qxy + r11 * qyy + r12 * qyz)
        + r02 * (r10 * qxz + r11 * qyz + r12 * qzz);
    m.qxy *= 2.0;

    m.qxz = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
        + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
        + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);
    m.qxz *= 2.0;

    m.qyy = r10 * (r10 * qxx + r11 * qxy + r12 * qxz)
        + r11 * (r10 * qxy + r11 * qyy + r12 * qyz)
        + r12 * (r10 * qxz + r11 * qyz + r12 * qzz);

    m.qyz = r10 * (r20 * qxx + r21 * qxy + r22 * qxz)
        + r11 * (r20 * qxy + r21 * qyy + r22 * qyz)
        + r12 * (r20 * qxz + r21 * qyz + r22 * qzz);
    m.qyz *= 2.0;

    m.qzz = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
        + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
        + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);
  }

  /**
   * Rotate the induced dipoles components.
   *
   * @param m PolarizableMultipole to rotate.
   */
  public void rotateInducedDipoles(PolarizableMultipole m) {
    double dx = m.ux;
    double dy = m.uy;
    double dz = m.uz;
    m.ux = r00 * dx + r01 * dy + r02 * dz;
    m.uy = r10 * dx + r11 * dy + r12 * dz;
    m.uz = r20 * dx + r21 * dy + r22 * dz;
    dx = m.px;
    dy = m.py;
    dz = m.pz;
    m.px = r00 * dx + r01 * dy + r02 * dz;
    m.py = r10 * dx + r11 * dy + r12 * dz;
    m.pz = r20 * dx + r21 * dy + r22 * dz;
    // Set update the averaged induced and chain rule terms.
    m.sx = (m.ux + m.px) * 0.5;
    m.sy = (m.uy + m.py) * 0.5;
    m.sz = (m.uz + m.pz) * 0.5;
  }

  /**
   * Rotate a vector in the QI frame into the global frame.
   *
   * @param v The vector to rotate (in-place).
   */
  public void toGlobal(double[] v) {
    double vx = v[0];
    double vy = v[1];
    double vz = v[2];
    // X-component.
    v[0] = ir00 * vx + ir01 * vy + ir02 * vz;
    // Y-component.
    v[1] = ir10 * vx + ir11 * vy + ir12 * vz;
    // Z-component.
    v[2] = ir20 * vx + ir21 * vy + ir22 * vz;
  }
}
