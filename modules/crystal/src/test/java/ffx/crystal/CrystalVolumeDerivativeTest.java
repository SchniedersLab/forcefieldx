/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.crystal;

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
public class CrystalVolumeDerivativeTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {"Triclinic (3TRW)",
                28.38, 31.73, 36.75, 90.12, 99.61, 96.52, "P-1"},
            {"Cubic (3RN4)",
                92.69, 92.69, 92.69, 90.0, 90.0, 90.0, "P23"},
            {"Orthorhombic (2WLD)",
                79.09, 94.81, 100.85, 90.0, 90.0, 90.0, "P222"},
            {"Monoclinic (3V0E)",
                50.85, 38.60, 89.83, 90.0, 103.99, 90.0, "P2"},
            {"Hexagonal (4DAC)",
                63.67, 63.67, 40.40, 90.0, 90.0, 120.0, "P6"},
            {"Tetragonal (3TI8)",
                112.56, 112.56, 66.81, 90.0, 90.0, 90.0, "P4"},
            {"Trigonal (3TNZ)",
                108.99, 108.99, 49.40, 90.0, 90.0, 120.0, "P3"}
        });
    }
    private final String info;
    private final Crystal crystal;

    public CrystalVolumeDerivativeTest(
            String info, double a, double b, double c,
            double alpha, double beta, double gamma, String sg) {
        this.info = info;
        this.crystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
    }
    /**
     * Finite-Difference parameters.
     */
    private double eps = 0.00001;
    private double epsD = Math.toDegrees(eps);
    private double eps2 = 2.0 * eps;
    private double tolerance = eps * 10.0;

    @Test
    public void finiteDifferenceTest() {

        /**
         * Current unit cell parameters.
         */
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;

        /**
         * Analytic volume derivatives with respect to unit parameters.
         */
        double dVdA = crystal.dVdA;
        double dVdB = crystal.dVdB;
        double dVdC = crystal.dVdC;
        double dVdAlpha = crystal.dVdAlpha;
        double dVdBeta = crystal.dVdBeta;
        double dVdGamma = crystal.dVdGamma;
        double dV;

        switch (crystal.spaceGroup.crystalSystem) {
            case TRICLINIC:
                testdVdA();
                testdVdB();
                testdVdC();
                testdVdAlpha();
                testdVdBeta();
                testdVdGamma();
                break;
            case MONOCLINIC:
                // alpha == gamma == 90.0
                testdVdA();
                testdVdB();
                testdVdC();
                testdVdBeta();
                break;
            case ORTHORHOMBIC:
                // alpha == beta == gamma == 90.0
                testdVdA();
                testdVdB();
                testdVdC();
                break;
            case TETRAGONAL:
                // a == b && alpha == beta == gamma == 90.0
                crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A & B should agree.", dVdA + dVdB, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                testdVdC();
                break;
            case TRIGONAL:
                if (a == b && b == c && alpha == beta && beta == gamma) {
                    // Rhombohedral axes, primitive cell.
                    crystal.changeUnitCellParameters(a + eps, b + eps, c + eps, alpha, beta, gamma);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a - eps, b - eps, c - eps, alpha, beta, gamma);
                    dV -= crystal.volume;
                    double dTot = dVdA + dVdB + dVdC;
                    assertEquals(info + " analytic and FD volume derivatives with respect to the sides should agree.", dTot, dV / eps2, tolerance);
                    assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                    crystal.changeUnitCellParameters(a, b, c, alpha + epsD, beta + epsD, gamma + epsD);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a, b, c, alpha - epsD, beta - epsD, gamma - epsD);
                    dV -= crystal.volume;
                    dTot = dVdAlpha + dVdBeta + dVdGamma;
                    assertEquals(info + " analytic and FD volume derivatives with respect to the angles should agree.", dTot, dV / eps2, tolerance);
                    assertEquals(info + " analytic volume derivatives with respect to Alpha and Beta should be equal.", dVdAlpha, dVdBeta, tolerance);
                } else if (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0) {
                    // Hexagonal axes, triple obverse cell.
                    crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                    dV -= crystal.volume;
                    assertEquals(info + " analytic and FD volume derivatives with respect to A & B should agree.", dVdA + dVdB, dV / eps2, tolerance);
                    assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                    testdVdC();
                }
                break;
            case HEXAGONAL:
                // a == b && alpha == beta == 90.0 && gamma == 120.0
                crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A & B should agree.", dVdA + dVdB, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                testdVdC();
                break;
            case CUBIC:
                // a == b == c && alpha == beta == gamma == 90.0
                crystal.changeUnitCellParameters(a + eps, b + eps, c + eps, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c - eps, alpha, beta, gamma);
                dV -= crystal.volume;
                double dTot = dVdA + dVdB + dVdC;
                assertEquals(info + " analytic and FD volume derivatives with respect to A, B & C should agree.", dTot, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                break;
        }
    }

    private void testdVdA() {
        double dVdA = crystal.dVdA;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a + eps, b, c, alpha, beta, gamma);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a - eps, b, c, alpha, beta, gamma);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
    }

    private void testdVdB() {
        double dVdB = crystal.dVdB;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a, b + eps, c, alpha, beta, gamma);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a, b - eps, c, alpha, beta, gamma);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to B should agree.", dVdB, dV / eps2, tolerance);
    }

    private void testdVdC() {
        double dVdC = crystal.dVdC;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a, b, c + eps, alpha, beta, gamma);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a, b, c - eps, alpha, beta, gamma);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to C should agree.", dVdC, dV / eps2, tolerance);
    }

    private void testdVdAlpha() {
        double dVdAlpha = crystal.dVdAlpha;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a, b, c, alpha + epsD, beta, gamma);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a, b, c, alpha - epsD, beta, gamma);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to Alpha should agree.", dVdAlpha, dV / eps2, tolerance);
    }

    private void testdVdBeta() {
        double dVdBeta = crystal.dVdBeta;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a, b, c, alpha, beta + epsD, gamma);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a, b, c, alpha, beta - epsD, gamma);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to Beta should agree.", dVdBeta, dV / eps2, tolerance);
    }

    private void testdVdGamma() {
        double dVdGamma = crystal.dVdGamma;
        double a = crystal.a;
        double b = crystal.b;
        double c = crystal.c;
        double alpha = crystal.alpha;
        double beta = crystal.beta;
        double gamma = crystal.gamma;
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma + epsD);
        double dV = crystal.volume;
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma - epsD);
        dV -= crystal.volume;
        assertEquals(info + " analytic and FD volume derivatives with respect to Gamma should agree.", dVdGamma, dV / eps2, tolerance);
    }
}
