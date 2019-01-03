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
public class Real3DTest {

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
    private final int paddedTot;
    private final int tot;
    private final double data[];
    private final double expected[];
    private final double recip[];
    private final double tolerance = 1.0e-14;

    public Real3DTest(String info, int nx, int ny, int nz) {
        this.info = info;
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        paddedTot = (nx + 2) * ny * nz;
        tot = nx * ny * nz;
        data = new double[paddedTot];
        expected = new double[tot];
        recip = new double[tot];
    }

    @Before
    public void setUp() {
        Random random = new Random();
        int paddedIndex = 0;
        int index = 0;
        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    double r = random.nextDouble();
                    data[paddedIndex++] = r;
                    expected[index] = r;
                    recip[index++] = 1.0e0;
                }
                paddedIndex += 2;
            }
        }
    }

    /**
     * Test of the fft and ifft methods, of class Real3D.
     */
    @Test
    public void testFft() {
        Real3D real3D = new Real3D(nx, ny, nz);
        real3D.fft(data);
        real3D.ifft(data);
        int paddedIndex = 0;
        int index = 0;
        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    double actual = data[paddedIndex++] / tot;
                    double orig = expected[index++];
                    assertEquals(info, orig, actual, tolerance);
                }
                paddedIndex += 2;
            }
        }
    }

    /**
     * Test of convolution method, of class Real3D.
     */
    @Test
    public void testConvolution() {
        Real3D real3D = new Real3D(nx, ny, nz);
        real3D.setRecip(recip);
        real3D.convolution(data);
        int paddedIndex = 0;
        int index = 0;
        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    double actual = data[paddedIndex++] / tot;
                    double orig = expected[index++];
                    assertEquals(info, orig, actual, tolerance);
                }
                paddedIndex += 2;
            }
        }
    }
}
