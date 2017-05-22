/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.numerics.fft;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class Complex3DTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Test nx=32, ny=32, nz=32}", 32, 32, 32},
        {"Test nx=32, ny=45, nz=21}", 32, 45, 21}
        });
    }
    private final String info;
    private final int nx;
    private final int ny;
    private final int nz;
    private final int tot;
    private final double data[];
    private final double expected[];
    private final double recip[];
    private final double tolerance = 1.0e-14;

    public Complex3DTest(String info, int nx, int ny, int nz) {
        this.info = info;
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        tot = nx * ny * nz;
        data = new double[tot * 2];
        expected = new double[tot];
        recip = new double[tot];
    }

    @Before
    public void setUp() {
        Random random = new Random();
        for (int i = 0; i < tot; i++) {
            int index = i * 2;
            double r = random.nextDouble();
            data[index] = r;
            expected[i] = r;
            recip[i] = 1.0e0;
        }
    }

    /**
     * Test of the fft and ifft methods, of class Complex3D.
     */
    @Test
    public void testFft() {
        Complex3D complex3D = new Complex3D(nx, ny, nz);
        complex3D.fft(data);
        complex3D.ifft(data);
        for (int i = 0; i < tot; i++) {
            int index = i * 2;
            double actual = data[index] / tot;
            double orig = expected[i];
            assertEquals(info, orig, actual, tolerance);
        }
    }

    /**
     * Test of convolution method, of class Complex3D.
     */
    @Test
    public void testConvolution() {
        Complex3D complex3D = new Complex3D(nx, ny, nz);
        complex3D.setRecip(recip);
        complex3D.convolution(data);
        for (int i = 0; i < tot; i++) {
            int index = i * 2;
            double actual = data[index] / tot;
            double orig = expected[i];
            assertEquals(info, orig, actual, tolerance);
        }
    }

    /**
     * Disable quasi-internal frame for all tests from this class.
     */
    @org.junit.BeforeClass
    public static void setUpClass() {
        qiPrevious = Boolean.valueOf(System.getProperty("pme.qi", "false"));
        System.setProperty("pme.qi", "false");
    }
    private static boolean qiPrevious;
    @org.junit.AfterClass
    public static void tearDownClass() {
        System.setProperty("pme.qi", String.valueOf(qiPrevious));
    }
}
