/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
package ffx.numerics;

import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * Parameterized test of the UniformBSpline class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class UniformBSplineTest {

    @Parameters
    public static Collection<Object[]> data() {
        // Order 5, x = 0.00
        double c1[] = {0.0416666e0, 0.4583333e0, 0.4583333e0, 0.0416666e0, 0.0000000e0};
        double d1[] = {-0.1666667e0, -0.5000000e0, 0.5000000e0, 0.1666667e0, 0.0000000e0};
        // Order 5, x = 0.25
        double c2[] = {0.0131836e0, 0.3248700e0, 0.5608720e0, 0.1009110e0, 0.0001627e0};
        double d2[] = {-0.0703125e0, -0.5416670e0, 0.2968750e0, 0.3125000e0, 0.0026042e0};
        // Order 5, x = 0.50
        double c3[] = {0.0026041e0, 0.1979167e0, 0.5989583e0, 0.1979167e0, 0.0026041e0};
        double d3[] = {-0.0208333e0, -0.4583333e0, 0.0000000e0, 0.4583333e0, 0.0208333e0};

        return Arrays.asList(new Object[][] {
                    {"Order 5, x = 0.00", 5, 0.00e0, c1, d1},
                    {"Order 5, x = 0.25", 5, 0.25e0, c2, d2},
                    {"Order 5, x = 0.50", 5, 0.50e0, c3, d3}
                });
    }
    private static double tolerance = 1.0e-6;
    private final String info;
    private final int order;
    private final double x;
    private final double coefficients[];
    private final double expected[];
    private final double expectedDiff[];
    private final int deriveOrder;
    private final double theta[][];
    private final double bSplineWork[][];

    /**
     * Parameterized test of the UniformBSpline class.
     *
     * @param info
     * @param order
     * @param x
     * @param expected
     * @param expectedDiff
     *
     * @since 1.0
     */
    public UniformBSplineTest(String info, int order, double x, double expected[], double expectedDiff[]) {
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
