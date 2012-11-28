/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * Parameterized Test of the TensorRecursion class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class TensorRecursionTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Test {1.1,1.2,1.3} for order 5", 1.1e0, 1.2e0, 1.3e0, 5, 56}
                });
    }
    private final double tolerance = 1.0e-15;
    private final double r[] = new double[3];
    private final double tensors[];
    private final double noStorageTensors[];
    private final int order;
    private final int tensorCount;
    private final String info;

    public TensorRecursionTest(String info, double x, double y, double z, int order, int tensorCount) {
        this.info = info;
        this.order = order;
        r[0] = x;
        r[1] = y;
        r[2] = z;
        this.tensorCount = tensorCount;
        tensors = new double[tensorCount];
        noStorageTensors = new double[tensorCount];
    }

    /**
     * Test of tensorIndex method, of class TensorRecursion.
     */
    @Test
    public void testTensorIndex() {
        int dx = 1;
        int dy = 0;
        int dz = 0;
        int expResult = 1;
        int result = TensorRecursion.tensorIndex(dx, dy, dz, order);
        assertEquals(info, expResult, result);
    }

    /**
     * Test of tensorCount method, of class TensorRecursion.
     *
     * @since 1.0
     */
    @Test
    public void testTensorCount() {
        int result = TensorRecursion.tensorCount(order);
        assertEquals(info, tensorCount, result);
    }

    /**
     * Test of tensorRecursion and noStorageTensorRecursion methods, of class
     * TensorRecursion.
     *
     * @since 1.0
     */
    @Test
    public void testTensorRecursion() {
        TensorRecursion tensorRecursion = new TensorRecursion(order);
        tensorRecursion.tensorRecursion(r, tensors);
        tensorRecursion.noStorageTensorRecursion(r, noStorageTensors);
        for (int i = 0; i < tensorCount; i++) {
            double expect = noStorageTensors[i];
            double actual = tensors[i];
            //System.out.println(String.format("%d: %10.5f %10.5f", i, expect, actual));
            assertEquals(info + " @ " + i, expect, actual, tolerance);
        }
    }
}
