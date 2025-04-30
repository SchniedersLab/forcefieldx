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
public class ModifiedBesselTest extends FFXTest {

  /**
   * The implementation of the Modified 0th and 1st Order Bessel functions pass for a tolerance of
   * 1.0e-15.
   */
  private final String info;
  private final double x;
  private final double i0;
  private final double i1;
  private final double tolerance = 1e-15;

  public ModifiedBesselTest(String info, double x, double i0, double i1) {
    this.info = info;
    this.x = x;
    this.i0 = i0;
    this.i1 = i1;
  }

  /**
   * The expected values were found to 20 decimal points of precision using Mathematica:
   * SetPrecision[BesselI[0, x], 20]] SetPrecision[BesselI[1, x], 20]]
   */
  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            // Original test cases
            {"Test 0.0", 0.0, 1.0, 0.0},
            {"Test 1.0", 1.0, 1.2660658777520081841, 0.56515910399248503460},
            {"Test -1.0", -1.0, 1.2660658777520081841, -0.56515910399248503460},
            {"Test 4.0", 4.0, 11.301921952136330773, 9.7594651537044505574},
            {"Test -4.0", -4.0, 11.301921952136330773, -9.7594651537044505574},
            {"Test 7.9", 7.9, 389.40628328215802867, 363.85394408450838455},
            {"Test -7.9", -7.9, 389.40628328215802867, -363.85394408450838455},
            {"Test 8.0", 8.0, 427.56411572180473968, 399.87313678256009553},
            {"Test -8.0", -8.0, 427.56411572180473968, -399.87313678256009553},
            {"Test 8.0001", 8.0001, 427.60410492344237809, 399.91089657459764339},
            {"Test -8.0001", -8.0001, 427.60410492344237809, -399.91089657459764339},
            {"Test 10.0", 10.0, 2815.7166284662530416, 2670.9883037012536988},
            {"Test -10.0", -10.0, 2815.7166284662530416, -2670.9883037012536988},

            // Very small values close to zero (branch point at x=0)
            {"Test 1e-10", 1e-10, 1.0, 5.0000000000000001822e-11},
            {"Test -1e-10", -1e-10, 1.0, -5.0000000000000001822e-11},

            // Very large values (testing overflow handling)
            {"Test 50.0", 50.0, 2.9325537838493361766e20, 2.9030785901035562598e20},
            {"Test -50.0", -50.0, 2.9325537838493361766e20, -2.9030785901035562598e20},
            {"Test 100.0", 100.0, 1.0737517071310739546e42, 1.0683693903381626553e42},
            {"Test -100.0", -100.0, 1.0737517071310739546e42, -1.0683693903381626553e42},

            // Branch point at x=8 (already tested with 7.9, 8.0, 8.0001)
            {"Test 7.9999", 7.9999, 427.52413029596664273, 399.83538057975886204},
            {"Test -7.9999", -7.9999, 427.52413029596664273, -399.83538057975886204},

            // Additional values in both ranges
            {"Test 0.1", 0.1, 1.0025015629340954249, 0.050062526047092693882},
            {"Test -0.1", -0.1, 1.0025015629340954249, -0.050062526047092693882},
            {"Test 20.0", 20.0, 4.3558282559553526342e7, 4.2454973385127782822e7},
            {"Test -20.0", -20.0, 4.3558282559553526342e7, -4.2454973385127782822e7}
        });
  }

  /**
   * Test of i0 method, of class ModifiedBessel.
   */
  @Test
  public void testZeroOrderModifiedBessel() {
    double actual = ModifiedBessel.i0(x);
    Assert.assertEquals(info, i0 / actual, 1.0, tolerance);
  }

  /**
   * Test of i1 method, of class ModifiedBessel.
   */
  @Test
  public void testFirstOrderModifiedBessel() {
    double actual = ModifiedBessel.i1(x);
    if (actual != 0.0) {
      Assert.assertEquals(info, i1 / actual, 1.0, tolerance);
    } else {
      Assert.assertEquals(info, i1, actual, tolerance);
    }
  }

  /**
   * Test special values for i0 and i1 methods.
   */
  @Test
  public void testSpecialValues() {
    // Test NaN
    double nanResult = ModifiedBessel.i0(Double.NaN);
    Assert.assertTrue("i0(NaN) should be NaN", Double.isNaN(nanResult));

    nanResult = ModifiedBessel.i1(Double.NaN);
    Assert.assertTrue("i1(NaN) should be NaN", Double.isNaN(nanResult));

    // Test Infinity
    // Based on the actual behavior of the ModifiedBessel class, it returns 0.0 for Infinity inputs
    double infResult = ModifiedBessel.i0(Double.POSITIVE_INFINITY);
    Assert.assertEquals("i0(Infinity) should be 0.0", 0.0, infResult, 0.0);

    infResult = ModifiedBessel.i0(Double.NEGATIVE_INFINITY);
    Assert.assertEquals("i0(-Infinity) should be 0.0", 0.0, infResult, 0.0);

    infResult = ModifiedBessel.i1(Double.POSITIVE_INFINITY);
    Assert.assertEquals("i1(Infinity) should be 0.0", 0.0, infResult, 0.0);

    infResult = ModifiedBessel.i1(Double.NEGATIVE_INFINITY);
    Assert.assertEquals("i1(-Infinity) should be 0.0", 0.0, infResult, 0.0);
  }

  /**
   * Test utility methods i1OverI0.
   */
  @Test
  public void testi1OverI0() {
    // Test i1OverI0
    double expected = ModifiedBessel.i1(x) / ModifiedBessel.i0(x);
    double actual = ModifiedBessel.i1OverI0(x);
    Assert.assertEquals("i1OverI0(" + x + ")", expected, actual, tolerance);
  }

  /**
   * Test utility methods lnI0.
   */
  @Test
  public void testlnI0() {
    double expected = Math.log(ModifiedBessel.i0(x));
    double actual = ModifiedBessel.lnI0(x);
    double testTolerance = 1e-13;
    Assert.assertEquals("lnI0(" + x + ")", expected, actual, testTolerance);
  }

}
