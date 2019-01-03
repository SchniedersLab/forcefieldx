/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
public class RealTest {

    private final int n;
    private final int paddedN;
    private final String info;
    private final double data[];
    private final double complexData[];
    private final double expected[];
    private final double tolerance = 1.0e-13;

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Test n = 20", 20},
                {"Test n = 22", 22},
                {"Test n = 120", 120}
        });
    }

    public RealTest(String info, int n) {
        this.info = info;
        this.n = n;
        assert (n % 2 == 0);
        paddedN = n + 2;
        data = new double[paddedN];
        complexData = new double[n * 2];
        expected = new double[n];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            double d = r.nextDouble();
            data[i] = d;
            complexData[i * 2] = d;
            expected[i] = d;
        }
    }

    /**
     * Test of fft method, of class Complex.
     */
    @Test
    public void testFft() {
        int offset = 0;
        int stride = 2;
        Real real = new Real(n);
        Complex complex = new Complex(n);
        real.fft(data, offset);
        complex.fft(complexData, offset, stride);
        for (int i = 0; i < n; i++) {
            double expect = complexData[i];
            double actual = data[i];
            assertEquals(info + " @ " + i, expect, actual, tolerance);
        }
        real.inverse(data, offset);
        for (int i = 0; i < n; i++) {
            double orig = expected[i];
            double actual = data[i];
            assertEquals(info, orig, actual, tolerance);
        }
    }
}
