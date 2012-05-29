/**
 * Title: Force Field X Description: Force Field X - Software for Molecular
 * Biophysics Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
                    {"Cubic",
                        12.0, 12.0, 12.0, 90.0, 90.0, 90.0, "P1"},
                    {"Monoclinic",
                        50.85, 38.60, 89.83, 90.0, 103.99, 90.0, "P2"},
                    {"Hexagonal",
                        63.67, 63.67, 40.40, 90.0, 90.0, 120.0, "P6"}
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

    @Test
    public void finiteDifferenceTest() {

        /**
         * Finite-Difference parameters.
         */
        double eps = 0.00001;
        double epsD = Math.toDegrees(eps);
        double eps2 = 2.0 * eps;
        double tolerance = eps * 10.0;
        
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
        
        switch (crystal.spaceGroup.crystalSystem) {
            case TRICLINIC:
                crystal.changeUnitCellParameters(a + eps, b, c, alpha, beta, gamma);
                double dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                break;
            case MONOCLINIC:
                // (alpha == beta || alpha == gamma);
                crystal.changeUnitCellParameters(a + eps, b, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                if (alpha == beta) {
                    crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma + epsD);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma - epsD);
                    dV -= crystal.volume;
                    assertEquals(info + " analytic and FD volume derivatives with respect to Gamma should agree.", dVdGamma, dV / eps2, tolerance);
                } else if (alpha == gamma) {
                    crystal.changeUnitCellParameters(a, b, c, alpha, beta + epsD, gamma);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a, b, c, alpha, beta - epsD, gamma);
                    dV -= crystal.volume;
                    assertEquals(info + " analytic and FD volume derivatives with respect to Beta should agree.", dVdBeta, dV / eps2, tolerance);
                }
                break;
            case ORTHORHOMBIC:
                // (alpha == 90.0 && beta == 90.0 && gamma == 90.0);
                crystal.changeUnitCellParameters(a + eps, b, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                break;
            case TETRAGONAL:
                // (a == b && alpha == 90.0 && beta == 90.0 && gamma == 90.0);
                crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                break;
            case TRIGONAL:
                // Rombohedral axes, primitive cell.
                if (a == b && b == c && alpha == beta && beta == gamma) {
                    crystal.changeUnitCellParameters(a + eps, b + eps, c + eps, alpha, beta, gamma);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a - eps, b - eps, c - eps, alpha, beta, gamma);
                    dV -= crystal.volume;
                    assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                    assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                    // Hexagonal axes, triple obverse cell.
                } else if (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0) {
                    crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                    dV = crystal.volume;
                    crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                    dV -= crystal.volume;
                    assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                    assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                } else {
                    // Programming error.
                }
                break;
            case HEXAGONAL:
                // (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0);
                crystal.changeUnitCellParameters(a + eps, b + eps, c, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                break;
            case CUBIC:
                // (a == b && b == c && alpha == 90.0 && beta == 90.0 && gamma == 90.0);
                crystal.changeUnitCellParameters(a + eps, b + eps, c + eps, alpha, beta, gamma);
                dV = crystal.volume;
                crystal.changeUnitCellParameters(a - eps, b - eps, c - eps, alpha, beta, gamma);
                dV -= crystal.volume;
                assertEquals(info + " analytic and FD volume derivatives with respect to A should agree.", dVdA, dV / eps2, tolerance);
                assertEquals(info + " analytic volume derivatives with respect to A and B should be equal.", dVdA, dVdB, tolerance);
                break;
        }
    }
}
