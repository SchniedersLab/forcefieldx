// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.numerics.integrate;

import static ffx.numerics.integrate.Integration.generateTestData_v1;
import static ffx.numerics.integrate.Integration.halfBinComposite;
import static ffx.numerics.integrate.Integration.leftBoole;
import static ffx.numerics.integrate.Integration.leftRectangularMethod;
import static ffx.numerics.integrate.Integration.leftSimpsons;
import static ffx.numerics.integrate.Integration.leftTrapInput;
import static ffx.numerics.integrate.Integration.rightBoole;
import static ffx.numerics.integrate.Integration.rightRectangularMethod;
import static ffx.numerics.integrate.Integration.rightSimpsons;
import static ffx.numerics.integrate.Integration.rightTrapInput;
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

/**
 * The IntegrationTest is a JUnit test for the Integration program that ensures that the integrals
 * are calculated correctly. This test is run using known integrals calculated with the equation:
 *
 * <p>y = 10 sin(6x) - 7 cos(5x) + 11 sin(8x).
 *
 * @author Claire O'Connell
 */
public class IntegrationTest {

  /** Create array with pointers to doubles that will contain known integrals. */
  private double[] knownIntegral;

  /** Compares the calculated integrals with the known values. */
  @Test
  public void integrationTest() {

    /*
     *  Calculate the integrals using the left hand trapezoidal, Simpson's,
     * and Boole's methods using data generated with the bounds 1 and 201
     * with an interval of .1. The second four are the right handed integrals
     * in the same order.
     */
    double[] calculatedIntegral = new double[8];

    calculatedIntegral[0] = leftTrapInput(generateTestData_v1());
    calculatedIntegral[1] =
        leftSimpsons(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 1, "left");
    calculatedIntegral[2] =
        leftBoole(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 2, "left");
    calculatedIntegral[3] = leftRectangularMethod(generateTestData_v1());

    calculatedIntegral[4] = rightTrapInput(generateTestData_v1());
    calculatedIntegral[5] =
        rightSimpsons(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 1, "right");
    calculatedIntegral[6] =
        rightBoole(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 2, "right");
    calculatedIntegral[7] = rightRectangularMethod(generateTestData_v1());

    // Set the delta value for the assertEquals comparison.
    double DELTA = 1e-8;

    // Assert that the known integrals and calculated integrals are the same.
    for (int i = 0; i < 8; i++) {
      assertEquals(knownIntegral[i], calculatedIntegral[i], DELTA);
    }
  }

  /** Initializes the array before testing. */
  @Before
  public void setUp() {
    // Instantiate the knownIntegral array.
    knownIntegral = new double[8];

    /*The answers are in the order of the trapezoidal integral first, the
    Simpson's second, Boole's third, and rectangular method last. The
    integrals are calculated with the bounds 1 and 201 with an interval of
    .1. The first four are using left hand integrals and the second four
    use right hand integrals.
    */
    knownIntegral[0] = 2.9684353512887753;
    knownIntegral[1] = 2.9687126459508564;
    knownIntegral[2] = 2.968712622691869;
    knownIntegral[3] = 2.936215172510247;

    knownIntegral[4] = 3.0006927996084642;
    knownIntegral[5] = 3.000977174918476;
    knownIntegral[6] = 3.000977149598709;
    knownIntegral[7] = 2.968898509509485;
  }
}
