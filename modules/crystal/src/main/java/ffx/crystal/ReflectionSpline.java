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

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;

/**
 * The ReflectionSpline class represents a reflection spline basis.
 *
 * @author Timothy D. Fenn<br>
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank"> K. Cowtan, J. Appl.
 *     Cryst. (2002). 35, 655-663 </a>
 * @since 1.0
 */
public class ReflectionSpline {

  private final ReflectionList reflectionList;
  private final int nParams;
  private double f;
  private int i0, i1, i2;
  private double dfi0, dfi1, dfi2;

  /**
   * Constructor for ReflectionSpline.
   *
   * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
   * @param nParams an int.
   */
  public ReflectionSpline(ReflectionList reflectionList, int nParams) {
    this.reflectionList = reflectionList;
    this.nParams = nParams;
  }

  /**
   * dfi0
   *
   * @return a double.
   */
  public double dfi0() {
    return dfi0;
  }

  /**
   * dfi1
   *
   * @return a double.
   */
  public double dfi1() {
    return dfi1;
  }

  /**
   * dfi2
   *
   * @return a double.
   */
  public double dfi2() {
    return dfi2;
  }

  /**
   * f
   *
   * @return a double.
   */
  public double f() {
    return f;
  }

  /**
   * Evaluate basis function and derivative at a given resolution
   * <p>
   * Equations 24 and 25 in Cowtan et al.
   *
   * @param invResSq resolution of desired spline interpolation
   * @param params current spline parameters
   * @return value at invResSq
   */
  public double f(double invResSq, double[] params) {
    double s = nParams * reflectionList.ordinal(invResSq);
    int i = (int) floor(s);
    double ds = s - i - 0.5;
    i0 = min(max(0, i - 1), nParams - 1);
    i1 = min(max(0, i), nParams - 1);
    i2 = min(max(0, i + 1), nParams - 1);

    f =
        params[i0] * 0.5 * (ds - 0.5) * (ds - 0.5)
            + params[i1] * (0.75 - ds * ds)
            + params[i2] * 0.5 * (ds + 0.5) * (ds + 0.5);

    dfi0 = 0.5 * (ds - 0.5) * (ds - 0.5);
    dfi1 = 0.75 - ds * ds;
    dfi2 = 0.5 * (ds + 0.5) * (ds + 0.5);

    return f;
  }

  /**
   * i0
   *
   * @return an int.
   */
  public int i0() {
    return i0;
  }

  /**
   * i1
   *
   * @return an int.
   */
  public int i1() {
    return i1;
  }

  /**
   * i2
   *
   * @return an int.
   */
  public int i2() {
    return i2;
  }
}
