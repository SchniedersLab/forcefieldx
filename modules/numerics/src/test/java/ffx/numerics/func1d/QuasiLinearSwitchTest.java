/**
 * Title: Force Field X.
 *
 * <p>Description: Force Field X - Software for Molecular Biophysics.
 *
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2001-2021.
 *
 * <p>This file is part of Force Field X.
 *
 * <p>Force Field X is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * <p>Force Field X is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * <p>You should have received a copy of the GNU General Public License along with Force Field X; if
 * not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 * <p>Linking this library statically or dynamically with other modules is making a combined work
 * based on this library. Thus, the terms and conditions of the GNU General Public License cover the
 * whole combination.
 *
 * <p>As a special exception, the copyright holders of this library give you permission to link this
 * library with independent modules to produce an executable, regardless of the license terms of
 * these independent modules, and to copy and distribute the resulting executable under terms of
 * your choice, provided that you also meet, for each linked independent module, the terms and
 * conditions of the license of that module. An independent module is a module which is not derived
 * from or based on this library. If you modify this library, you may extend this exception to your
 * version of the library, but you are not obligated to do so. If you do not wish to do so, delete
 * this exception statement from your version.
 */
package ffx.numerics.func1d;

import static ffx.numerics.func1d.QuasiLinearThetaMap.Branch;
import static ffx.numerics.func1d.QuasiLinearThetaMap.Branch.A;
import static ffx.numerics.func1d.QuasiLinearThetaMap.Branch.B;
import static ffx.numerics.func1d.QuasiLinearThetaMap.Branch.C;
import static ffx.numerics.func1d.QuasiLinearThetaMap.Branch.D;
import static java.lang.String.format;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.OptionalDouble;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

@RunWith(Parameterized.class)
public class QuasiLinearSwitchTest {
  private static final double TOL = 1E-11;
  private final String info;
  private final QuasiLinearThetaMap mapF;
  private final double[] expectedConsts;
  private final double[][] expecteds;
  private final Branch[] primaryBranch;
  private final Branch[] secondaryBranch;
  private final int nVals;

  public QuasiLinearSwitchTest(
      String info,
      OptionalDouble theta0,
      double[] eConsts,
      double[][] eVals,
      Branch[] primary,
      Branch[] secondary) {
    this.info = info;
    mapF =
        theta0.isPresent()
            ? new QuasiLinearThetaMap(theta0.getAsDouble())
            : new QuasiLinearThetaMap();
    expectedConsts = eConsts;
    expecteds = eVals;
    primaryBranch = primary;
    secondaryBranch = secondary;
    nVals = eVals.length;
    assertTrue(
        "Inequal number of arguments for expected vals & branches!",
        nVals == primary.length && nVals == secondary.length);
  }

  @Parameterized.Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
          {
            // Test each endpoint & each midpoint.
            "Default Linear Theta Map",
            OptionalDouble.empty(),
            // r, a, b, c
            new double[] {
              6.586293806435631, 0.32876610686815605, -0.016424593043157933, -5.586293806435631
            },
            new double[][] {
              {-Math.PI, 1.0, 0.0, -3.2931469032178153},
              {(0.05 - Math.PI), 0.995884423889855, -0.16458874650913316, -3.28903132710767},
              {(0.1 - Math.PI), 0.9835479823563423, -0.32876610686815605, -3.2766948855741576},
              {0.5 * (0.1 - Math.PI), 0.4835616946565922, -0.32876610686815605, 0.0},
              {-0.1, 0.0164520176436577, -0.32876610686815605, 3.276694885574158},
              {-0.05, 0.004115576110145532, -0.16458874650913274, 3.28903132710767},
              {0, 0, 0, 3.2931469032178153},
              {0.05, 0.004115576110145532, 0.16458874650913274, 3.28903132710767},
              {0.1, 0.0164520176436577, 0.32876610686815605, 3.276694885574158},
              {0.5 * (Math.PI - 0.1), 0.48356169465659227, 0.32876610686815605, 0.0},
              {Math.PI - 0.1, 0.9835479823563427, 0.3287661068681553, -3.2766948855741576},
              {Math.PI - 0.05, 0.995884423889855, 0.16458874650913316, -3.28903132710767},
              {Math.PI, 1.0, 0.0, -3.2931469032178153}
            },
            new QuasiLinearThetaMap.Branch[] {D, D, D, C, A, A, A, A, A, B, D, D, D},
            new QuasiLinearThetaMap.Branch[] {D, D, C, C, C, A, A, A, B, B, B, D, D}
          },
          {
            "Linear Theta Map: t0 = 0.1 (manual)",
            OptionalDouble.of(0.1),
            new double[] {
              6.586293806435631, 0.32876610686815605, -0.016424593043157933, -5.586293806435631
            },
            new double[][] {
              {-Math.PI, 1.0, 0.0, -3.2931469032178153},
              {(0.05 - Math.PI), 0.995884423889855, -0.16458874650913316, -3.28903132710767},
              {(0.1 - Math.PI), 0.9835479823563423, -0.32876610686815605, -3.2766948855741576},
              {0.5 * (0.1 - Math.PI), 0.4835616946565922, -0.32876610686815605, 0.0},
              {-0.1, 0.0164520176436577, -0.32876610686815605, 3.276694885574158},
              {-0.05, 0.004115576110145532, -0.16458874650913274, 3.28903132710767},
              {0, 0, 0, 3.2931469032178153},
              {0.05, 0.004115576110145532, 0.16458874650913274, 3.28903132710767},
              {0.1, 0.0164520176436577, 0.32876610686815605, 3.276694885574158},
              {0.5 * (Math.PI - 0.1), 0.48356169465659227, 0.32876610686815605, 0.0},
              {Math.PI - 0.1, 0.9835479823563427, 0.3287661068681553, -3.2766948855741576},
              {Math.PI - 0.05, 0.995884423889855, 0.16458874650913316, -3.28903132710767},
              {Math.PI, 1.0, 0.0, -3.2931469032178153}
            },
            new QuasiLinearThetaMap.Branch[] {D, D, D, C, A, A, A, A, A, B, D, D, D},
            new QuasiLinearThetaMap.Branch[] {D, D, C, C, C, A, A, A, B, B, B, D, D}
          },
        });
  }

  @Test
  public void testPoints() {
    for (int i = 0; i < 4; i++) {
      assertEquals(
          "Incorrect constant " + i + ": ", mapF.getConstants()[i], expectedConsts[i], TOL);
    }
    for (int i = 0; i < nVals; i++) {
      double[] exp = expecteds[i];
      Branch prim = primaryBranch[i];
      Branch secon = secondaryBranch[i];
      if (prim == secon) {
        testMidPoint(exp, prim);
      } else {
        testJoint(exp, prim, secon);
      }
    }
  }

  @Override
  public String toString() {
    return info;
  }

  private void testMidPoint(double[] expect, Branch br) {
    double t = expect[0];
    assertEquals("Incorrect branching!", br, mapF.getBranch(t));
    double eVal = expect[1];
    double val = mapF.valueAt(t);
    assertEquals("A midpoint value was incorrect!", val, eVal, TOL);
    val = mapF.val(t, br);
    assertEquals("A midpoint value was incorrect when the branch was specified!", val, eVal, TOL);

    eVal = expect[2];
    val = mapF.firstDerivative(t);
    assertEquals("A midpoint first derivative was incorrect!", val, eVal, TOL);
    val = mapF.fd(t, br);
    assertEquals(
        "A midpoint first derivative was incorrect when the branch was specified!", val, eVal, TOL);

    eVal = expect[3];
    val = mapF.secondDerivative(t);
    assertEquals("A midpoint second derivative was incorrect!", val, eVal, TOL);
    val = mapF.sd(t, br);
    assertEquals(
        "A midpoint second derivative was incorrect when the branch was specified!",
        val,
        eVal,
        TOL);
  }

  private void testJoint(double[] expect, Branch prim, Branch second) {
    double t = expect[0];
    Branch received = mapF.getBranch(t);
    boolean rightBranch = (prim == received) || (second == received);
    assertTrue(
        format(
            " Did not find any expected branch! Found: %s. "
                + "Primary expected: %s. Secondary expected: %s",
            received, prim, second),
        rightBranch);

    double eVal = expect[1];
    double val = mapF.valueAt(t);
    assertEquals("A joint value was incorrect!", val, eVal, TOL);
    val = mapF.val(t, prim);
    assertEquals("A joint value was incorrect when the branch was specified!", val, eVal, TOL);
    val = mapF.val(t, second);
    assertEquals("A joint value was incorrect for the alternate branch!", val, eVal, TOL);

    eVal = expect[2];
    val = mapF.firstDerivative(t);
    assertEquals("A joint first derivative was incorrect!", val, eVal, TOL);
    val = mapF.fd(t, prim);
    assertEquals(
        "A joint first derivative was incorrect when the branch was specified!", val, eVal, TOL);
    val = mapF.fd(t, second);
    assertEquals(
        "A joint first derivative was incorrect for the alternate branch!", val, eVal, TOL);

    eVal = expect[3];
    val = mapF.sd(t, prim);
    assertEquals("A joint second derivative was incorrect!", val, eVal, TOL);
  }
}
