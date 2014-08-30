/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.numerics;

import java.util.Arrays;
import java.util.Collection;

import static java.lang.Math.sqrt;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class InverseSqrtTest {

    /**
     * Java double precision follows the IEEE 754 Binary Floating-Point
     * Arithmetic standard. Each double consumes 8 bytes of storage and offers
     * 52 binary digits of precision (14-15 decimal digits). This implementation
     * of Inverse SQRT passes for a tolerance of 1.0e-12.
     */
    private static final double tolerance = 1.0e-12;

    /**
     * The expected values are found using java.lang.Math.sqrt.
     */
    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {"Test 1.0e-6", 1.0e-6},
            {"Test 1.0e-5", 1.0e-5},
            {"Test 1.0e-4", 1.0e-4},
            {"Test 1.0e-3", 1.0e-3},
            {"Test 1.0e-2", 1.0e-2},
            {"Test 1.0e-1", 1.0e-1},
            {"Test 1.0", 1.0},
            {"Test 1.0e1", 1.0e1},
            {"Test 1.0e2", 1.0e2},
            {"Test 1.0e3", 1.0e3},
            {"Test 1.0e4", 1.0e4},
            {"Test 1.0e5", 1.0e5},
            {"Test 1.0e6", 1.0e6},
            {"Test 1.0e7", 1.0e7},
            {"Test 1.0e8", 1.0e8},
            {"Test 2.0", 2.0},
            {"Test 2.12345", 2.12345e0},
            {"Test 4.0", 4.0},
            {"Test 6.1", 6.1},
            {"Test 12345.12345", 12345.12345},});
    }
    private final String info;
    private final double x;

    public InverseSqrtTest(String info, double x) {
        this.info = info;
        this.x = x;
    }

    /**
     * Test inverse sqrt method, of class InverseSqrt.
     */
    @Test
    public void testInverseSqrt() {
        double actual = InverseSqrt.inverseSQRT(x);
        double expected = 1.0 / sqrt(x);
        assertEquals(info, expected, actual, tolerance);
    }

}
