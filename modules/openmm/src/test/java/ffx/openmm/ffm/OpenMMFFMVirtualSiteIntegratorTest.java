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
import static org.junit.Assume.assumeTrue;

/** Native construction and parameter tests for virtual sites and standard integrators. */
public class OpenMMFFMVirtualSiteIntegratorTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testAverageAndOutOfPlaneSites() {
    try (TwoParticleAverageSite two = new TwoParticleAverageSite(0, 1, 0.25, 0.75);
         ThreeParticleAverageSite three = new ThreeParticleAverageSite(0, 1, 2, 0.2, 0.3, 0.5);
         OutOfPlaneSite out = new OutOfPlaneSite(0, 1, 2, 0.1, 0.2, 0.3)) {
      assertEquals(2, two.getNumParticles());
      assertEquals(1, two.getParticle(1));
      assertEquals(0.25, two.getWeight(0), 0.0);
      assertEquals(3, three.getNumParticles());
      assertEquals(2, three.getParticle(2));
      assertEquals(0.5, three.getWeight(2), 0.0);
      assertEquals(0.1, out.getWeight12(), 0.0);
      assertEquals(0.2, out.getWeight13(), 0.0);
      assertEquals(0.3, out.getWeightCross(), 0.0);
    }
  }

  @Test
  public void testLocalCoordinatesSite() {
    try (IntArray particles = new IntArray(0);
         DoubleArray origin = new DoubleArray(3);
         DoubleArray x = new DoubleArray(3);
         DoubleArray y = new DoubleArray(3)) {
      particles.append(0);
      particles.append(1);
      particles.append(2);
      origin.set(0, 1.0);
      origin.set(1, 0.0);
      origin.set(2, 0.0);
      x.set(0, -1.0);
      x.set(1, 0.5);
      x.set(2, 0.5);
      y.set(0, 0.0);
      y.set(1, -1.0);
      y.set(2, 1.0);

      try (LocalCoordinatesSite site =
               new LocalCoordinatesSite(particles, origin, x, y, new Vec3(0.1, 0.2, 0.3))) {
        assertEquals(3, site.getNumParticles());
        assertArrayEquals(new double[]{1.0, 0.0, 0.0}, site.getOriginWeights(), 0.0);
        assertArrayEquals(new double[]{-1.0, 0.5, 0.5}, site.getXWeights(), 0.0);
        assertArrayEquals(new double[]{0.0, -1.0, 1.0}, site.getYWeights(), 0.0);
        assertEquals(new Vec3(0.1, 0.2, 0.3), site.getLocalPosition());
      }
    }
    try (LocalCoordinatesSite site = new LocalCoordinatesSite(
        0, 1, 2, new Vec3(1.0, 0.0, 0.0), new Vec3(-1.0, 0.5, 0.5),
        new Vec3(0.0, -1.0, 1.0), new Vec3(0.1, 0.2, 0.3))) {
      assertEquals(3, site.getNumParticles());
      assertArrayEquals(new double[]{1.0, 0.0, 0.0}, site.getOriginWeights(), 0.0);
    }
  }

  @Test
  public void testLangevinAndBrownianIntegratorParameters() {
    try (LangevinMiddleIntegrator langevin = new LangevinMiddleIntegrator(0.002, 300.0, 1.5);
         LangevinIntegrator alias = new LangevinIntegrator(0.003, 310.0, 2.0);
         BrownianIntegrator brownian = new BrownianIntegrator(290.0, 3.0, 0.004)) {
      assertEquals(300.0, langevin.getTemperature(), 0.0);
      assertEquals(1.5, langevin.getFriction(), 0.0);
      langevin.setTemperature(305.0);
      langevin.setFriction(1.7);
      langevin.setRandomNumberSeed(9);
      assertEquals(305.0, langevin.getTemperature(), 0.0);
      assertEquals(1.7, langevin.getFriction(), 0.0);
      assertEquals(9, langevin.getRandomNumberSeed());
      assertEquals(310.0, alias.getTemperature(), 0.0);
      assertEquals(290.0, brownian.getTemperature(), 0.0);
      assertEquals(3.0, brownian.getFriction(), 0.0);
      brownian.setRandomNumberSeed(11);
      assertEquals(11, brownian.getRandomNumberSeed());
    }
  }

  @Test
  public void testVariableIntegratorParameters() {
    try (VariableVerletIntegrator verlet = new VariableVerletIntegrator(0.001);
         VariableLangevinIntegrator langevin =
             new VariableLangevinIntegrator(300.0, 1.0, 0.002)) {
      verlet.setMaximumStepSize(0.02);
      assertEquals(0.001, verlet.getErrorTolerance(), 0.0);
      assertEquals(0.02, verlet.getMaximumStepSize(), 0.0);
      verlet.setErrorTolerance(0.003);
      assertEquals(0.003, verlet.getErrorTolerance(), 0.0);

      langevin.setTemperature(315.0);
      langevin.setFriction(1.3);
      langevin.setErrorTolerance(0.004);
      langevin.setMaximumStepSize(0.03);
      langevin.setRandomNumberSeed(13);
      assertEquals(315.0, langevin.getTemperature(), 0.0);
      assertEquals(1.3, langevin.getFriction(), 0.0);
      assertEquals(0.004, langevin.getErrorTolerance(), 0.0);
      assertEquals(0.03, langevin.getMaximumStepSize(), 0.0);
      assertEquals(13, langevin.getRandomNumberSeed());
    }
  }
}
