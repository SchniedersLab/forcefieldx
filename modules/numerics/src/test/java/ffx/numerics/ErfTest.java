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

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class ErfTest {

    /**
     * Java double precision follows the IEEE 754 Binary Floating-Point
     * Arithmetic standard. Each double consumes 8 bytes of storage and offers
     * 52 binary digits of precision (14-15 decimal digits). This implementation
     * of Erf passes for a tolerance of 1.0e-15 and (as one might expect) fails
     * using 1.0e-16.
     */
    private static final double tolerance = 1.0e-15;

    /**
     * The expected values were found to 20 decimal points of precision using
     * Mathematica: Erf[SetPrecision[x, 20]] Erfc[SetPrecision[x, 20]]
     */
    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Test 0.0", 0.0e0, 0.0e0},
        {"Test 0.1; below the first branch point.", 0.1e0, 0.1124629160182848984e0},
        {"Test 0.46875; at the first branch point.", 0.46875e0, 0.4926134732179379916e0},
        {"Test 1.0; between the branch points.", 1.0e0, 0.842700792949714869e0},
        {"Test 4.0; at the second branch point.", 4.0e0, 1.0e0 - 1.5417257900280018852e-8},
        {"Test 5.0; above the second branch point.", 5.0e0, 1.0e0 - 1.5374597944280348502e-12}
        });
    }
    private final String info;
    private final double x;
    private final double expected;

    public ErfTest(String info, double x, double expected) {
        this.info = info;
        this.x = x;
        this.expected = expected;
    }

    /**
     * Test of erf method, of class Erf.
     */
    @Test
    public void testErf() {
        double actual = Erf.erf(x);
        assertEquals(info, expected, actual, tolerance);
    }

    /**
     * Test of erfc method, of class Erf.
     */
    @Test
    public void testErfc() {
        double actual = Erf.erfc(x);
        assertEquals(info, 1.0 - expected, actual, tolerance);
    }
}
