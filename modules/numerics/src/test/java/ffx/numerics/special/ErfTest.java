// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.numerics.special;

import ffx.utilities.FFXTest;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class ErfTest extends FFXTest {

  /**
   * Java double precision follows the IEEE 754 Binary Floating-Point Arithmetic standard. Each
   * double consumes 8 bytes of storage and offers 52 binary digits of precision (14-15 decimal
   * digits). This implementation of Erf passes for a tolerance of 1.0e-15 and (as one might expect)
   * fails using 1.0e-16.
   */
  private static final double tolerance = 1.0e-15;

  private final String info;
  private final double x;
  private final double expected;

  public ErfTest(String info, double x, double expected) {
    this.info = info;
    this.x = x;
    this.expected = expected;
  }

  /**
   * The expected values were found to 20 decimal points of precision using Mathematica:
   * Erf[SetPrecision[x, 20]] Erfc[SetPrecision[x, 20]]
   */
  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            {"Test 0.0", 0.0e0, 0.0e0},

            // Very small value tests
            {"Test 1.0e-16; near xSmall threshold.", 1.0e-16, 1.128379167095513e-16},
            {"Test 1.0e-10; small value.", 1.0e-10, 1.128379167095513e-10},

            {"Test 0.1; below the first branch point.", 0.1e0, 0.1124629160182849},
            {"Test 0.46; just below first branch point.", 0.46, 0.4846553900016797},
            {"Test 0.46875; at the first branch point.", 0.46875, 0.492613473217938},
            {"Test 0.47; just above first branch point.", 0.47, 0.4937450508860821},
            {"Test 1.0; between the branch points.", 1.0e0, 0.842700792949715},
            {"Test 3.9; just below second branch point.", 3.9, 1.0 - 3.479224859723177e-8},
            {"Test 4.0; at the second branch point.", 4.0, 1.0 - 1.5417257900280018852e-8},
            {"Test 4.1; just above second branch point.", 4.1, 1.0 - 6.700027654084919e-9},
            {"Test 5.0; above the second branch point.", 5.0e0, 1.0e0 - 1.5374597944280348502e-12},

            // Very large value tests
            {"Test 10.0; large value.", 10.0e0, 1.0 - 2.088487583762545e-45},
            {"Test 26.0; near xBig threshold.", 26.0, 1.0e0},

            // Negative value tests (erf is an odd function: erf(-x) = -erf(x))
            {"Test -0.1; negative value below first branch point.", -0.1e0, -0.1124629160182849},
            {"Test -1.0; negative value between branch points.", -1.0e0, -0.842700792949715},
            {"Test -5.0; negative value above second branch point.", -5.0e0, -1.0e0 + 1.5374597944280348502e-12}
        });
  }

  /**
   * Test of erf method, of class Erf.
   */
  @Test
  public void testErf() {
    double actual = Erf.erf(x);
    Assert.assertEquals(info, expected, actual, tolerance);
  }

  /**
   * Test of erfc method, of class Erf.
   */
  @Test
  public void testErfc() {
    double actual = Erf.erfc(x);
    Assert.assertEquals(info, 1.0 - expected, actual, tolerance);
  }

  /**
   * Test of erf method with special values.
   */
  @Test
  public void testErfSpecialCases() {
    // Test NaN - not using parameterized test to avoid issues with NaN comparison
    double nanResult = Erf.erf(Double.NaN);
    Assert.assertTrue("erf(NaN) should be NaN", Double.isNaN(nanResult));

    // Test positive infinity
    double posInfResult = Erf.erf(Double.POSITIVE_INFINITY);
    Assert.assertEquals("erf(+Infinity) should be 1.0", 1.0, posInfResult, 0.0);

    // Test negative infinity
    double negInfResult = Erf.erf(Double.NEGATIVE_INFINITY);
    Assert.assertEquals("erf(-Infinity) should be -1.0", -1.0, negInfResult, 0.0);
  }

  /**
   * Test of erfc method with special values.
   */
  @Test
  public void testErfcSpecialCases() {
    // Test NaN - not using parameterized test to avoid issues with NaN comparison
    double nanResult = Erf.erfc(Double.NaN);
    Assert.assertTrue("erfc(NaN) should be NaN", Double.isNaN(nanResult));

    // Test positive infinity
    double posInfResult = Erf.erfc(Double.POSITIVE_INFINITY);
    Assert.assertEquals("erfc(+Infinity) should be 0.0", 0.0, posInfResult, 0.0);

    // Test negative infinity
    double negInfResult = Erf.erfc(Double.NEGATIVE_INFINITY);
    Assert.assertEquals("erfc(-Infinity) should be 2.0", 2.0, negInfResult, 0.0);
  }

}
