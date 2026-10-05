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
package ffx.openmm.ffm.drude;

import ffx.openmm.ffm.OpenMMRuntime;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native construction and parameter tests for Drude FFM façades. */
public class OpenMMFFMDrudeTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testDrudeForceParameters() {
    try (DrudeForce force = new DrudeForce()) {
      force.addParticle(1, 0, -1, -1, -1, -0.5, 0.001, 1.1, 1.2);
      force.addParticle(3, 2, 1, 4, 5, -0.4, 0.002, 0.8, 0.9);
      force.addScreenedPair(0, 1, 1.3);
      assertEquals(2, force.getNumParticles());
      assertEquals(1, force.getNumScreenedPairs());
      assertEquals(new DrudeForce.ParticleParameters(
          3, 2, 1, 4, 5, -0.4, 0.002, 0.8, 0.9), force.getParticleParameters(1));
      assertEquals(new DrudeForce.ScreenedPairParameters(0, 1, 1.3),
          force.getScreenedPairParameters(0));

      force.setParticleParameters(1, 4, 3, 2, 1, 0, -0.3, 0.003, 0.7, 0.6);
      force.setScreenedPairParameters(0, 1, 0, 1.4);
      force.setUsesPeriodicBoundaryConditions(true);
      assertEquals(new DrudeForce.ParticleParameters(
          4, 3, 2, 1, 0, -0.3, 0.003, 0.7, 0.6), force.getParticleParameters(1));
      assertEquals(new DrudeForce.ScreenedPairParameters(1, 0, 1.4),
          force.getScreenedPairParameters(0));
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testBaseDrudeIntegrator() {
    try (DrudeIntegrator integrator = new DrudeIntegrator(0.001)) {
      assertEquals(0.0, integrator.getDrudeTemperature(), 0.0);
      assertEquals(0.0, integrator.getMaxDrudeDistance(), 0.0);
      integrator.setDrudeTemperature(1.0);
      integrator.setMaxDrudeDistance(0.03);
      integrator.setRandomNumberSeed(14);
      assertEquals(1.0, integrator.getDrudeTemperature(), 0.0);
      assertEquals(0.03, integrator.getMaxDrudeDistance(), 0.0);
      assertEquals(14, integrator.getRandomNumberSeed());
    }
  }

  @Test
  public void testDrudeLangevinIntegrator() {
    try (DrudeLangevinIntegrator integrator =
             new DrudeLangevinIntegrator(0.001, 300.0, 1.0, 1.0, 20.0)) {
      assertEquals(0.001, integrator.getStepSize(), 0.0);
      assertEquals(300.0, integrator.getTemperature(), 0.0);
      assertEquals(1.0, integrator.getFriction(), 0.0);
      assertEquals(20.0, integrator.getDrudeFriction(), 0.0);
      integrator.setTemperature(310.0);
      integrator.setFriction(2.0);
      integrator.setDrudeFriction(30.0);
      integrator.setDrudeTemperature(2.0);
      assertEquals(310.0, integrator.getTemperature(), 0.0);
      assertEquals(2.0, integrator.getFriction(), 0.0);
      assertEquals(30.0, integrator.getDrudeFriction(), 0.0);
      assertEquals(2.0, integrator.getDrudeTemperature(), 0.0);
    }
  }

  @Test
  public void testDrudeScfIntegrator() {
    try (DrudeSCFIntegrator integrator = new DrudeSCFIntegrator(0.001)) {
      assertEquals(0.001, integrator.getStepSize(), 0.0);
      double initialTolerance = integrator.getMinimizationErrorTolerance();
      integrator.setMinimizationErrorTolerance(2.5);
      assertEquals(2.5, integrator.getMinimizationErrorTolerance(), 0.0);
      assertTrue(initialTolerance > 0.0);
    }
  }

  @Test
  public void testDrudeNoseHooverIntegrator() {
    try (DrudeNoseHooverIntegrator integrator =
             new DrudeNoseHooverIntegrator(0.001, 300.0, 1.0, 1.0, 20.0, 3, 3, 1)) {
      assertEquals(0.001, integrator.getStepSize(), 0.0);
      assertEquals(0.02, integrator.getMaxDrudeDistance(), 1.0e-12);
      integrator.setMaxDrudeDistance(0.025);
      assertEquals(0.025, integrator.getMaxDrudeDistance(), 0.0);
      assertEquals(0.001, integrator.getStepSize(), 0.0);
    }
  }
}
