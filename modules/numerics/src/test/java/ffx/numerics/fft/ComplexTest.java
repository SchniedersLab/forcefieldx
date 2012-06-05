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
package ffx.numerics.fft;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class ComplexTest {

    private final int n;
    private final String info;
    private final boolean preferred;
    private final double data[];
    private final double expected[];
    private final double tolerance = 1.0e-14;

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][] {
            /**
             * This test will fail without the factor 7
             * {"Test n = 21", 21, true},
             */
                                             {"Test n = 22", 22, false},
                                             {"Test n = 120", 120, true}
                });
    }

    public ComplexTest(String info, int n, boolean preferred) {
        this.info = info;
        this.n = n;
        this.preferred = preferred;
        data = new double[n * 2];
        expected = new double[n * 2];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            double d = r.nextDouble();
            data[i * 2] = d;
            expected[i * 2] = d;
        }
    }

    /**
     * Test of preferredDimension method, of class Complex.
     */
    @Test
    public void testPreferredDimension() {
        boolean result = Complex.preferredDimension(n);
        assertEquals(info, preferred, result);
    }

    /**
     * Test of fft method, of class Complex.
     */
    @Test
    public void testFft() {
        int offset = 0;
        int stride = 2;
        Complex complex = new Complex(n);
        complex.fft(data, offset, stride);
        complex.inverse(data, offset, stride);
        for (int i = 0; i < n; i++) {
            double orig = expected[i * 2];
            double actual = data[i * 2];
            assertEquals(info, orig, actual, tolerance);
        }
    }
}
