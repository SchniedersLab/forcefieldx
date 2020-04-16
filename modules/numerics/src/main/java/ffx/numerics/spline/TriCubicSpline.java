// ******************************************************************************
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
// ******************************************************************************
package ffx.numerics.spline;

/**
 * TriCubicSpline class.
 *
 * @author Timothy D. Fenn
 * @see <a href="http://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf"
 *     target="_blank">Catmull-Rom splines</a>
 * @since 1.0
 */
public class TriCubicSpline {

  /** Tau for the smoothing matrix. */
  private static final double tau = 0.25;
  /** Smoothing matrix: Catmull-Rom spline with tau=0.25 */
  private static final double[][] catmullRomMat =
      new double[][] {
        {0.0, 1.0, 0.0, 0.0},
        {-tau, 0.0, tau, 0.0},
        {2.0 * tau, tau - 3.0, 3.0 - 2.0 * tau, -tau},
        {-tau, 2.0 - tau, tau - 2.0, tau}
      };

  private final double[] p;
  private final double[] q;
  private final double[] r;
  private final double[] u;
  private final double[] v;
  private final double[] w;
  private final double[] dp;
  private final double[] dq;
  private final double[] dr;
  private final double[] du;
  private final double[] dv;
  private final double[] dw;

  /** Initialize Spline function. */
  public TriCubicSpline() {
    dw = new double[4];
    dv = new double[4];
    du = new double[4];
    dr = new double[4];
    dq = new double[4];
    dp = new double[4];
    w = new double[4];
    v = new double[4];
    u = new double[4];
    r = new double[4];
    q = new double[4];
    p = new double[4];
  }

  /**
   * Determine the spline value at a given point.
   *
   * @param dx delta between point and previous grid point in X
   * @param dy delta between point and previous grid point in Y
   * @param dz delta between point and previous grid point in Z
   * @param scalar 3d array in x,y,z order of 3D scalar data
   * @param g gradient array (can be null)
   * @return the interpolated scalar value at the requested point
   */
  public double spline(double dx, double dy, double dz, double[][][] scalar, double[] g) {

    // p(s) = u . catmull-rom matrix . p^T applied in 3 dimensions (u, v, w)
    u[0] = 1.0;
    v[0] = 1.0;
    w[0] = 1.0;
    for (int i = 1; i < 4; i++) {
      u[i] = u[i - 1] * dx;
      v[i] = v[i - 1] * dy;
      w[i] = w[i - 1] * dz;
    }

    // Derivatives
    du[0] = dv[0] = dw[0] = 0.0;
    du[1] = dv[1] = dw[1] = 1.0;
    du[2] = 2.0 * dx;
    du[3] = 3.0 * dx * dx;
    dv[2] = 2.0 * dy;
    dv[3] = 3.0 * dy * dy;
    dw[2] = 2.0 * dz;
    dw[3] = 3.0 * dz * dz;

    // vec4mat4 - could put this in VectorMath class
    for (int i = 0; i < 4; i++) {
      p[i] = q[i] = r[i] = 0.0;
      dp[i] = dq[i] = dr[i] = 0.0;
      for (int j = 0; j < 4; j++) {
        p[i] += u[j] * catmullRomMat[j][i];
        q[i] += v[j] * catmullRomMat[j][i];
        r[i] += w[j] * catmullRomMat[j][i];
        dp[i] += du[j] * catmullRomMat[j][i];
        dq[i] += dv[j] * catmullRomMat[j][i];
        dr[i] += dw[j] * catmullRomMat[j][i];
      }
    }

    // Tensor products
    double sum = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          sum += p[i] * q[j] * r[k] * scalar[i][j][k];
          gx += dp[i] * q[j] * r[k] * scalar[i][j][k];
          gy += p[i] * dq[j] * r[k] * scalar[i][j][k];
          gz += p[i] * q[j] * dr[k] * scalar[i][j][k];
        }
      }
    }

    if (g != null) {
      g[0] = gx;
      g[1] = gy;
      g[2] = gz;
    }

    return sum;
  }
}
