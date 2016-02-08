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
public class SquareRootTest {

    private static final double tolerance = 1.0e-13;

    /**
     * The expected values are found using java.lang.Math.sqrt.
     */
    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {"Test 1.0e-2", 1.0e-2},
            {"Test 1.0e-1", 1.0e-1},
            {"Test 1.0e0", 1.0e0},
            {"Test 1.0e1", 1.0e1},
            {"Test 1.0e2", 1.0e2},});
    }
    private final String info;
    private final double x;

    public SquareRootTest(String info, double x) {
        this.info = info;
        this.x = x;
    }

    /**
     * Test sqrt method, of class InverseSqrt.
     */
    @Test
    public void testSqrt() {
        for (int i = 0; i < 1000; i++) {
            double increment = i * x / 100.0;
            double actual = SquareRoot.sqrt(x + increment);
            double expected = sqrt(x + increment);
            assertEquals(info + " + " + increment, expected, actual, tolerance);
        }
    }

    /**
     * Test inverse sqrt method, of class InverseSqrt.
     */
    @Test
    public void testInverseSqrt() {
        for (int i = 0; i < 100; i++) {
            double increment = i * x / 10.0;
            double actual = SquareRoot.isqrt(x + increment);
            double expected = 1.0 / sqrt(x + increment);
            assertEquals(info + " + " + increment, expected, actual, tolerance);
        }
    }

}
