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
package ffx.numerics;

import java.util.Arrays;
import java.util.Collection;

import static java.util.Arrays.fill;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

import ffx.numerics.MultipoleTensor.COORDINATES;
import ffx.numerics.MultipoleTensor.OPERATOR;

/**
 * Parameterized Test of the MultipoleTensor class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class MultipoleTensorTest {

    private static final double xyz[] = {1.1, 1.2, 1.3, 0.5};

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {"Order 5 Coulomb", xyz[0], xyz[1], xyz[2], xyz[3], 5, 56, OPERATOR.COULOMB, 0.0, 0.0, 0.0},
            {"Order 5 Screened Coulomb", xyz[0], xyz[1], xyz[2], xyz[3], 5, 56, OPERATOR.SCREENED_COULOMB, 0.45, 0.0, 0.0},
            {"Order 4 Thole Field", xyz[0], xyz[1], xyz[2], xyz[3], 4, 35, OPERATOR.THOLE_FIELD, 0.0, 1.0, 1.0}
        });
    }
    private final double tolerance = 1.0e-14;
    private final double fdTolerance = 1.0e-6;
    private final double r[] = new double[3];
    private final double tensor[];
    private final double noStorageTensor[];
    private final double fastTensor[];
    private final int order;
    private final int tensorCount;
    private final String info;
    private final OPERATOR operator;
    private final double beta;
    private final double damp;
    private final double aiak;
    private double lambdaFunction;
    
    private static final boolean testAllQI;
    static {
        String testQIstring = System.getProperty("testAllQI");
        if (testQIstring != null) {
            testAllQI = Boolean.parseBoolean(testQIstring);
        } else {
            testAllQI = false;
        }
    }

    private final double Qi[] = {0.11,
        0.21, 0.31, 0.41,
        -0.51, -0.61, 1.12, 0.71, 0.81, 0.91};
    private final double Qk[] = {0.11,
        0.21, 0.31, 0.41,
        -0.51, -0.61, 1.12, 0.71, 0.81, 0.91};

    public MultipoleTensorTest(String info, double x, double y, double z, double lambdaFunction,
            int order, int tensorCount, OPERATOR operator,
            double beta, double damp, double aiak) {
        this.info = info;
        this.order = order;
        r[0] = x;
        r[1] = y;
        r[2] = z;
        this.lambdaFunction = lambdaFunction;
        this.tensorCount = tensorCount;
        tensor = new double[tensorCount];
        noStorageTensor = new double[tensorCount];
        fastTensor = new double[tensorCount];
        this.operator = operator;
        this.beta = beta;
        this.damp = damp;
        this.aiak = aiak;
    }

    /**
     * Test of ti method, of class MultipoleTensor.
     */
    @Test
    public void tensorIndexTest() {
        int dx = 1;
        int dy = 0;
        int dz = 0;
        int expResult = 1;
        int result = MultipoleTensor.ti(dx, dy, dz, order);
        assertEquals(info, expResult, result);
    }

    /**
     * Test of tensorCount method, of class MultipoleTensor.
     *
     * @since 1.0
     */
    @Test
    public void tensorCountTest() {
        int result = MultipoleTensor.tensorCount(order);
        assertEquals(info, tensorCount, result);
    }

    /**
     * Test of recursion and noStorageRecursion methods, of class
     * MultipoleTensor.
     *
     * @since 1.0
     */
    @Test
    public void multipoleTensorTest() {
        MultipoleTensor multipoleTensor = new MultipoleTensor(
                operator, COORDINATES.GLOBAL, order, beta);
		multipoleTensor.setTholeDamping(damp, aiak);
		
		/**
         * Check Cartesian Tensors in the Global frame.
         */
        multipoleTensor.noStorageRecursion(r, noStorageTensor);
        multipoleTensor.recursion(r, tensor);
        multipoleTensor.generateTensor();
        multipoleTensor.getTensor(fastTensor);
        
        for (int i = 0; i < tensorCount; i++) {
            double expect = noStorageTensor[i];
            double actual = tensor[i];
            assertEquals(info + " @ " + i, expect, actual, tolerance);
            if (order == 4 || order == 5) {
                expect = noStorageTensor[i];
                actual = fastTensor[i];
                assertEquals(info + " @ " + i, expect, actual, tolerance);
            }
        }

        /**
         * Check QI Tensors in a quasi-internal frame.
         */
        // Set x and y = 0.0
        r[0] = 0.0;
        r[1] = 0.0;
        fill(noStorageTensor, 0.0);
        fill(tensor, 0.0);
        fill(fastTensor, 0.0);
        multipoleTensor.noStorageRecursion(r, noStorageTensor);
        multipoleTensor.recursion(r, tensor);
        multipoleTensor.setTensor(fastTensor);
        multipoleTensor.generateTensor();
        multipoleTensor.getTensor(fastTensor);

        for (int i = 0; i < tensorCount; i++) {
            double expect = noStorageTensor[i];
            double actual = tensor[i];
            assertEquals(info + " @ " + i, expect, actual, tolerance);
            if (order == 4 || order == 5) {
                expect = noStorageTensor[i];
                actual = fastTensor[i];
                assertEquals(info + " @ " + i, expect, actual, tolerance);
            }
        }
    }

    @Test
    public void finiteDifferenceTest() {
        MultipoleTensor multipoleTensor = new MultipoleTensor(
                operator, COORDINATES.GLOBAL, order, beta);
        multipoleTensor.setTholeDamping(damp, aiak);

        multipoleTensor.recursion(r, tensor);

        double tensorsPx[] = new double[tensorCount];
        double tensorsNx[] = new double[tensorCount];
        double tensorsPy[] = new double[tensorCount];
        double tensorsNy[] = new double[tensorCount];
        double tensorsPz[] = new double[tensorCount];
        double tensorsNz[] = new double[tensorCount];
        double delta = 1.0e-5;
        double delta2 = delta * 2;
        r[0] += delta;
        multipoleTensor.noStorageRecursion(r, tensorsPx);
        r[0] -= delta2;
        multipoleTensor.noStorageRecursion(r, tensorsNx);
        r[0] += delta;

        r[1] += delta;
        multipoleTensor.noStorageRecursion(r, tensorsPy);
        r[1] -= delta2;
        multipoleTensor.noStorageRecursion(r, tensorsNy);
        r[1] += delta;

        r[2] += delta;
        multipoleTensor.noStorageRecursion(r, tensorsPz);
        r[2] -= delta2;
        multipoleTensor.noStorageRecursion(r, tensorsNz);
        r[2] += delta;

        tensorFiniteDifference(multipoleTensor, delta2,
                tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);

        multipoleTensor.recursion(r, tensor);
        tensorFiniteDifference(multipoleTensor, delta2,
                tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);

        multipoleTensor.generateTensor();
        multipoleTensor.getTensor(tensor);
        tensorFiniteDifference(multipoleTensor, delta2,
                tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);
    }

    @Test
    public void finiteDifferenceQITest() {
        MultipoleTensor multipoleTensor = new MultipoleTensor(
                operator, COORDINATES.QI, order, beta);
        multipoleTensor.setTholeDamping(damp, aiak);

        // Set x and y = 0.0
        r[0] = 0.0;
        r[1] = 0.0;

        multipoleTensor.noStorageRecursion(r, tensor);
        double tensorsPz[] = new double[tensorCount];
        double tensorsNz[] = new double[tensorCount];
        double delta = 1.0e-5;
        double delta2 = delta * 2;

        r[2] += delta;
        multipoleTensor.noStorageRecursion(r, tensorsPz);
        r[2] -= delta2;
        multipoleTensor.noStorageRecursion(r, tensorsNz);
        r[2] += delta;

        tensorFiniteDifferenceQI(multipoleTensor, delta2, tensorsPz, tensorsNz);

        multipoleTensor.recursion(r, tensor);
        tensorFiniteDifferenceQI(multipoleTensor, delta2, tensorsPz, tensorsNz);

        multipoleTensor.generateTensor();
        multipoleTensor.getTensor(tensor);
        tensorFiniteDifferenceQI(multipoleTensor, delta2, tensorsPz, tensorsNz);
    }

    @Test
    public void energyAndForceTest() {

        if (operator == OPERATOR.THOLE_FIELD) {
            return;
        }

        MultipoleTensor multipoleTensor = new MultipoleTensor(
                operator, COORDINATES.GLOBAL, order, beta);
        r[0] = 1.1;
        r[1] = 1.2;
        r[2] = 1.3;
        double delta = 1.0e-5;
        double delta2 = 2.0 * 1.0e-5;
        double Fi[] = new double[3];
        double Ti[] = new double[3];
        double Tk[] = new double[3];
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);

        double energy = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);

        double aX = -Fi[0];
        double aY = -Fi[1];
        double aZ = -Fi[2];

        double analyticdEdF = multipoleTensor.getdEdZbuff();

        r[0] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posX = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[0] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negX = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[0] += delta;
        r[1] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posY = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[1] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negY = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[1] += delta;
        r[2] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posZ = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[2] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negZ = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[2] += delta;

        lambdaFunction += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posF = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        lambdaFunction -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negF = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        lambdaFunction += delta;

        double expect = analyticdEdF;
        double actual = (posF - negF) / delta2;
        assertEquals(info + " Global dE/dF", expect, actual, fdTolerance);
        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " Force Z", expect, actual, fdTolerance);
    }

    @Test
    public void energyAndForceQITest() {
        if (!testAllQI) {
            return;
        }

        if (operator == OPERATOR.THOLE_FIELD) {
            return;
        }

        MultipoleTensor multipoleTensor = new MultipoleTensor(
                operator, COORDINATES.QI, order, beta);
        double delta = 1.0e-5;
        double delta2 = 2.0 * 1.0e-5;
        double Fi[] = new double[3];
        double Ti[] = new double[3];
        double Tk[] = new double[3];

        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double energy = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);

        double aX = -Fi[0];
        double aY = -Fi[1];
        double aZ = -Fi[2];

        double analyticdEdF = multipoleTensor.getdEdZbuff();

        r[0] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posX = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[0] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negX = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[0] += delta;

        r[1] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posY = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[1] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negY = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[1] += delta;

        r[2] += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posZ = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[2] -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negZ = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        r[2] += delta;

        lambdaFunction += delta;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double posF = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        lambdaFunction -= delta2;
        multipoleTensor.generateTensor(r, lambdaFunction, Qi, Qk);
        double negF = multipoleTensor.multipoleEnergy(Fi, Ti, Tk);
        lambdaFunction += delta;

        double expect = analyticdEdF;
        double actual = (posF - negF) / delta2;
        assertEquals(info + " QI dE/dF", expect, actual, fdTolerance);

        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " QI Force X", expect, actual, fdTolerance);

        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " QI Force Y", expect, actual, fdTolerance);

        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " QI Force Z", expect, actual, fdTolerance);
    }

    private void tensorFiniteDifference(MultipoleTensor multipoleTensor, double delta2,
            double tensorsPx[], double tensorsNx[],
            double tensorsPy[], double tensorsNy[],
            double tensorsPz[], double tensorsNz[]) {

        int start = 0;
        /**
         * We do not calculate the zeroth term for Thole damping.
         */
        if (operator == OPERATOR.THOLE_FIELD) {
            start = 1;
        }

        /**
         * Test the partial derivatives for all tensor components.
         */
        for (int l = start; l < order; l++) {
            // Test X derivative
            double expect = tensor[multipoleTensor.ti(l + 1, 0, 0)];
            double actual = (tensorsPx[multipoleTensor.ti(l, 0, 0)] - tensorsNx[multipoleTensor.ti(l, 0, 0)]) / delta2;
            assertEquals(info + " @ " + l, expect, actual, fdTolerance);
            // Test Y derivative
            expect = tensor[multipoleTensor.ti(l, 1, 0)];
            actual = (tensorsPy[multipoleTensor.ti(l, 0, 0)] - tensorsNy[multipoleTensor.ti(l, 0, 0)]) / delta2;
            assertEquals(info + " @ " + l, expect, actual, fdTolerance);
            // Test Z derivative
            expect = tensor[multipoleTensor.ti(l, 0, 1)];
            actual = (tensorsPz[multipoleTensor.ti(l, 0, 0)] - tensorsNz[multipoleTensor.ti(l, 0, 0)]) / delta2;
            assertEquals(info + " @ " + l, expect, actual, fdTolerance);
        }
        for (int l = 0; l < order; l++) {
            for (int m = 1; m < order - l; m++) {
                // Test X derivative
                double expect = tensor[multipoleTensor.ti(l + 1, m, 0)];
                double actual = (tensorsPx[multipoleTensor.ti(l, m, 0)] - tensorsNx[multipoleTensor.ti(l, m, 0)]) / delta2;
                assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                // Test Y derivative
                expect = tensor[multipoleTensor.ti(l, m + 1, 0)];
                actual = (tensorsPy[multipoleTensor.ti(l, m, 0)] - tensorsNy[multipoleTensor.ti(l, m, 0)]) / delta2;
                assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                // Test Z derivative
                expect = tensor[multipoleTensor.ti(l, m, 1)];
                actual = (tensorsPz[multipoleTensor.ti(l, m, 0)] - tensorsNz[multipoleTensor.ti(l, m, 0)]) / delta2;
                assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m < order - l; m++) {
                for (int n = 1; n < order - l - m; n++) {
                    // Test X derivative
                    double expect = tensor[multipoleTensor.ti(l + 1, m, n)];
                    double actual = (tensorsPx[multipoleTensor.ti(l, m, n)] - tensorsNx[multipoleTensor.ti(l, m, n)]) / delta2;
                    assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                    // Test Y derivative
                    expect = tensor[multipoleTensor.ti(l, m + 1, n)];
                    actual = (tensorsPy[multipoleTensor.ti(l, m, n)] - tensorsNy[multipoleTensor.ti(l, m, n)]) / delta2;
                    assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                    // Test Z derivative
                    expect = tensor[multipoleTensor.ti(l, m, n + 1)];
                    actual = (tensorsPz[multipoleTensor.ti(l, m, n)] - tensorsNz[multipoleTensor.ti(l, m, n)]) / delta2;
                    assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                }
            }
        }
    }

    private void tensorFiniteDifferenceQI(MultipoleTensor multipoleTensor, double delta2,
            double tensorsPz[], double tensorsNz[]) {

        int start = 0;
        /**
         * We do not calculate the zeroth term for Thole damping.
         */
        if (operator == OPERATOR.THOLE_FIELD) {
            start = 1;
        }

        /**
         * Test the partial derivatives for all tensor components.
         */
        for (int l = start; l < order; l++) {
            // Test Z derivative
            double expect = tensor[multipoleTensor.ti(l, 0, 1)];
            double actual = (tensorsPz[multipoleTensor.ti(l, 0, 0)] - tensorsNz[multipoleTensor.ti(l, 0, 0)]) / delta2;
            assertEquals(info + " @ " + l, expect, actual, fdTolerance);
        }
        for (int l = 0; l < order; l++) {
            for (int m = 1; m < order - l; m++) {
                // Test Z derivative
                double expect = tensor[multipoleTensor.ti(l, m, 1)];
                double actual = (tensorsPz[multipoleTensor.ti(l, m, 0)] - tensorsNz[multipoleTensor.ti(l, m, 0)]) / delta2;
                assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
            }
        }
        // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l < order; l++) {
            for (int m = 0; m < order - l; m++) {
                for (int n = 1; n < order - l - m; n++) {
                    // Test Z derivative
                    double expect = tensor[multipoleTensor.ti(l, m, n + 1)];
                    double actual = (tensorsPz[multipoleTensor.ti(l, m, n)] - tensorsNz[multipoleTensor.ti(l, m, n)]) / delta2;
                    assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
                }
            }
        }
    }

}
