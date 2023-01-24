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
package ffx.numerics.spline;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collection;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Parameterized test of the UniformBSpline class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class UniformBSplineTest {

  private static double tolerance = 1.0e-6;
  private final String info;
  private final int order;
  private final double x;
  private final double[] coefficients;
  private final double[] expected;
  private final double[] expectedDiff;
  private final int deriveOrder;
  private final double[][] theta;
  private final double[][] bSplineWork;

  /**
   * Parameterized test of the UniformBSpline class.
   *
   * @param info
   * @param order
   * @param x
   * @param expected
   * @param expectedDiff
   * @since 1.0
   */
  public UniformBSplineTest(
      String info, int order, double x, double[] expected, double[] expectedDiff) {
    this.info = info;
    this.order = order;
    this.x = x;
    this.expected = expected;
    this.expectedDiff = expectedDiff;

    deriveOrder = order - 1;
    coefficients = new double[order];
    theta = new double[order][deriveOrder + 1];
    bSplineWork = new double[order][order];
  }

  @Parameters
  public static Collection<Object[]> data() {
    // Order 3, x = 0.0 (
    double[] c30 = {0.500000e0, 0.5000000e0, 0.000000e0};
    double[] d30 = {-1.0000000e0, 1.0000000e0, 0.0000000e0};
    // Order 3, x = 0.5 (At a grid point).
    double[] c31 = {0.1250000e0, 0.7500000e0, 0.1250000e0};
    double[] d31 = {-0.5000000e0, 0.0000000e0, 0.5000000e0};

    // Order 5, x = 0.00
    double[] c50 = {0.0416666e0, 0.4583333e0, 0.4583333e0, 0.0416666e0, 0.0000000e0};
    double[] d50 = {-0.1666667e0, -0.5000000e0, 0.5000000e0, 0.1666667e0, 0.0000000e0};
    // Order 5, x = 0.25
    double[] c51 = {0.0131836e0, 0.3248700e0, 0.5608720e0, 0.1009110e0, 0.0001627e0};
    double[] d51 = {-0.0703125e0, -0.5416670e0, 0.2968750e0, 0.3125000e0, 0.0026042e0};
    // Order 5, x = 0.50
    double[] c52 = {0.0026041e0, 0.1979167e0, 0.5989583e0, 0.1979167e0, 0.0026041e0};
    double[] d52 = {-0.0208333e0, -0.4583333e0, 0.0000000e0, 0.4583333e0, 0.0208333e0};
    // Order 5, x = 0.75
    double[] c53 = {0.0001627e0, 0.1009110e0, 0.5608720e0, 0.3248700e0, 0.0131836e0};
    double[] d53 = {-0.0026042e0, -0.3125000e0, -0.2968750e0, 0.5416670e0, 0.0703125e0};
    // Order 5, x = 1.0
    double[] c54 = {0.0000000e0, 0.0416666e0, 0.4583333e0, 0.4583333e0, 0.0416666e0};
    double[] d54 = {0.0000000e0, -0.1666667e0, -0.5000000e0, 0.5000000e0, 0.1666667e0};

    return Arrays.asList(
        new Object[][] {
            {"Order 3, x = 0.00", 3, 0.00e0, c30, d30},
            {"Order 3, x = 0.50", 3, 0.50e0, c31, d31},
            {"Order 5, x = 0.00", 5, 0.00e0, c50, d50},
            {"Order 5, x = 0.25", 5, 0.25e0, c51, d51},
            {"Order 5, x = 0.50", 5, 0.50e0, c52, d52},
            {"Order 5, x = 0.75", 5, 0.75e0, c53, d53},
            {"Order 5, x = 1.00", 5, 1.00e0, c54, d54}
        });
  }

  /**
   * Test of bSpline method, of class UniformBSpline.
   *
   * @since 1.0
   */
  @Test
  public void testBSpline() {
    UniformBSpline.bSpline(x, order, coefficients);
    for (int i = 0; i < order; i++) {
      assertEquals(info, expected[i], coefficients[i], tolerance);
    }
  }

  /**
   * Test of bSplineDerivatives method, of class UniformBSpline.
   *
   * @since 1.0
   */
  @Test
  public void testBSplineDerivatives() {
    UniformBSpline.bSplineDerivatives(x, order, deriveOrder, theta, bSplineWork);
    for (int i = 0; i < order; i++) {
      double actual = theta[i][0];
      assertEquals(info, expected[i], actual, tolerance);
      actual = theta[i][1];
      assertEquals(info, expectedDiff[i], actual, tolerance);
    }
  }
}
