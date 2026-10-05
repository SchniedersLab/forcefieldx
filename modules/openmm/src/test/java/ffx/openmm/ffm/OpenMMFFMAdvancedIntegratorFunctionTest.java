// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.openmm.ffm;

import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native tests for advanced integrators and tabulated functions. */
public class OpenMMFFMAdvancedIntegratorFunctionTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testTabulatedFunctionRoundTrips() {
    double[] values1 = {0.0, 0.25, 0.5, 0.75};
    double[] values3 = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    try (Continuous1DFunction continuous = new Continuous1DFunction(values1, -1.0, 2.0, false);
         Continuous2DFunction continuous2 =
             new Continuous2DFunction(values1, 2, 2, -1.0, 1.0, -2.0, 2.0, false);
         Continuous3DFunction continuous3 =
             new Continuous3DFunction(values3, 2, 2, 2, -1.0, 1.0, -2.0, 2.0, 0.0, 1.0, false);
         Discrete1DFunction discrete = new Discrete1DFunction(values1);
         Discrete2DFunction discrete2 = new Discrete2DFunction(2, 2, values1);
         Discrete3DFunction discrete3 = new Discrete3DFunction(2, 2, 2, values3)) {
      Continuous1DFunction.Parameters c1 = continuous.getFunctionParameters();
      assertArrayEquals(values1, c1.values(), 0.0);
      assertEquals(-1.0, c1.min(), 0.0);
      assertEquals(2.0, c1.max(), 0.0);
      assertArrayEquals(values1, continuous2.getFunctionParameters().values(), 0.0);
      assertEquals(2, continuous2.getFunctionParameters().xsize());
      assertArrayEquals(values3, continuous3.getFunctionParameters().values(), 0.0);
      assertEquals(2, continuous3.getFunctionParameters().zsize());
      assertArrayEquals(values1, discrete.getFunctionParameters(), 0.0);
      assertArrayEquals(values1, discrete2.getFunctionParameters().values(), 0.0);
      assertArrayEquals(values3, discrete3.getFunctionParameters().values(), 0.0);

      double[] replacement = {0.1, 0.2, 0.3, 0.4};
      double[] replacement3 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      continuous.setFunctionParameters(replacement, 0.0, 3.0);
      continuous2.setFunctionParameters(replacement, 2, 2, 0.0, 2.0, -1.0, 1.0);
      continuous3.setFunctionParameters(
          replacement3, 2, 2, 2, 0.0, 2.0, -1.0, 1.0, 0.0, 1.0);
      discrete.setFunctionParameters(replacement);
      discrete2.setFunctionParameters(2, 2, replacement);
      discrete3.setFunctionParameters(2, 2, 2, replacement3);
      assertArrayEquals(replacement, continuous.getFunctionParameters().values(), 0.0);
      assertArrayEquals(replacement, continuous2.getFunctionParameters().values(), 0.0);
      assertArrayEquals(replacement3, continuous3.getFunctionParameters().values(), 0.0);
      assertArrayEquals(replacement, discrete.getFunctionParameters(), 0.0);
      assertArrayEquals(replacement, discrete2.getFunctionParameters().values(), 0.0);
      assertArrayEquals(replacement3, discrete3.getFunctionParameters().values(), 0.0);
    }
  }

  @Test
  public void testCustomAndCompoundIntegrators() {
    try (CustomIntegrator custom = new CustomIntegrator(0.001)) {
      custom.addGlobalVariable("counter", 2.5);
      assertEquals(2.5, custom.getGlobalVariable(0), 0.0);
      assertEquals(2.5, custom.getGlobalVariableByName("counter"), 0.0);
      custom.setGlobalVariableByName("counter", 4.5);
      assertEquals(4.5, custom.getGlobalVariable(0), 0.0);
      assertEquals("counter", custom.getGlobalVariableName(0));
      custom.addComputeGlobal("counter", "counter+1");
      assertEquals(1, custom.getNumComputations());
      custom.setKineticEnergyExpression("0");
      assertEquals("0", custom.getKineticEnergyExpression());
      custom.setRandomNumberSeed(31);
      assertEquals(31, custom.getRandomNumberSeed());
    }

    CompoundIntegrator compound = new CompoundIntegrator();
    VerletIntegrator child = new VerletIntegrator(0.001);
    int index = compound.addIntegrator(child);
    assertEquals(0, index);
    assertTrue(child.isDestroyed());
    assertEquals(1, compound.getNumIntegrators());
    assertEquals(0.001, compound.getStepSize(), 0.0);
    compound.destroy();
  }

  @Test
  public void testNoseHooverAndReporterLifecycle() {
    try (NoseHooverIntegrator nose = new NoseHooverIntegrator(0.001, 300.0, 1.0, 3, 3, 2);
         MinimizationReporter reporter = new MinimizationReporter()) {
      assertEquals(300.0, nose.getTemperature(0), 0.0);
      assertEquals(1.0, nose.getCollisionFrequency(0), 0.0);
      assertEquals(1, nose.getNumThermostats());
      nose.setTemperature(305.0, 0);
      assertEquals(305.0, nose.getTemperature(0), 0.0);
      assertTrue(reporter.getPointer().address() != 0);
    }
  }
}
