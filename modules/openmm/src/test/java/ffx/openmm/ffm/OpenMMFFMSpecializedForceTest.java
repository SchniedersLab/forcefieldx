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

/** Native tests for specialized custom forces. */
public class OpenMMFFMSpecializedForceTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testCmapParameters() {
    try (CMAPTorsionForce force = new CMAPTorsionForce()) {
      assertEquals(0, force.addMap(2, new double[]{0.0, 1.0, 2.0, 3.0}));
      assertEquals(0, force.addTorsion(0, 0, 1, 2, 3, 4, 5, 6, 7));
      assertEquals(1, force.getNumMaps());
      assertEquals(1, force.getNumTorsions());
      CMAPTorsionForce.MapParameters map = force.getMapParameters(0);
      assertEquals(2, map.size());
      assertArrayEquals(new double[]{0.0, 1.0, 2.0, 3.0}, map.energy(), 0.0);
      assertEquals(
          new CMAPTorsionForce.TorsionParameters(0, 0, 1, 2, 3, 4, 5, 6, 7),
          force.getTorsionParameters(0));

      force.setMapParameters(0, 2, new double[]{3.0, 2.0, 1.0, 0.0});
      force.setTorsionParameters(0, 0, 7, 6, 5, 4, 3, 2, 1, 0);
      force.setUsesPeriodicBoundaryConditions(true);
      assertArrayEquals(new double[]{3.0, 2.0, 1.0, 0.0},
          force.getMapParameters(0).energy(), 0.0);
      assertEquals(
          new CMAPTorsionForce.TorsionParameters(0, 7, 6, 5, 4, 3, 2, 1, 0),
          force.getTorsionParameters(0));
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testGayBerneParameters() {
    try (GayBerneForce force = new GayBerneForce()) {
      force.addParticle(0.3, 1.2, 2, 3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9);
      assertEquals(1, force.getNumParticles());
      GayBerneForce.ParticleParameters particle = force.getParticleParameters(0);
      assertEquals(0.3, particle.sigma(), 0.0);
      assertEquals(1.2, particle.epsilon(), 0.0);
      assertEquals(2, particle.xparticle());
      assertEquals(3, particle.yparticle());
      assertEquals(0.4, particle.ex(), 0.0);
      assertEquals(0.5, particle.ey(), 0.0);
      assertEquals(0.6, particle.ez(), 0.0);
      assertEquals(0.7, particle.sx(), 0.0);
      assertEquals(0.8, particle.sy(), 0.0);
      assertEquals(0.9, particle.sz(), 0.0);

      force.setParticleParameters(0, 0.31, 1.3, 4, 5, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91);
      force.addException(0, 1, 0.25, 0.4, false);
      assertEquals(new GayBerneForce.ExceptionParameters(0, 1, 0.25, 0.4),
          force.getExceptionParameters(0));
      force.setExceptionParameters(0, 2, 3, 0.35, 0.45);
      assertEquals(new GayBerneForce.ExceptionParameters(2, 3, 0.35, 0.45),
          force.getExceptionParameters(0));
      force.setNonbondedMethod(1);
      force.setCutoffDistance(1.1);
      force.setUseSwitchingFunction(true);
      force.setSwitchingDistance(0.8);
      assertEquals(1, force.getNonbondedMethod());
      assertEquals(1.1, force.getCutoffDistance(), 0.0);
      assertTrue(force.isUseSwitchingFunction());
      assertEquals(0.8, force.getSwitchingDistance(), 0.0);
    }
  }

  @Test
  public void testRmsdParametersAndOwnership() {
    Vec3[] positions = {new Vec3(0.0, 0.1, 0.2), new Vec3(1.0, 1.1, 1.2)};
    try (RMSDForce force = new RMSDForce(positions, new int[]{1, 0})) {
      assertArrayEquals(new int[]{1, 0}, force.getParticles());
      assertArrayEquals(positions, force.getReferencePositions());
      force.setParticles(new int[]{0});
      force.setReferencePositions(new Vec3[]{new Vec3(2.0, 3.0, 4.0)});
      assertArrayEquals(new int[]{0}, force.getParticles());
      assertArrayEquals(new Vec3[]{new Vec3(2.0, 3.0, 4.0)}, force.getReferencePositions());
    }
  }

  @Test
  public void testAtmParametersAndNestedForceOwnership() {
    try (ATMForce force = new ATMForce("u1-u0");
         CustomVolumeForce nested = new CustomVolumeForce("volume")) {
      force.addGlobalParameter("lambda", 0.5);
      assertEquals(0, force.addParticle(new Vec3(1.0, 0.0, 0.0),
          new Vec3(-1.0, 0.0, 0.0)));
      assertEquals(0, force.addForce(nested));
      assertTrue(nested.isDestroyed());
      assertEquals(1, force.getNumForces());
      assertEquals(1, force.getNumParticles());
      assertEquals("lambda", force.getGlobalParameterName(0));
      assertEquals(0.5, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals(new ATMForce.ParticleParameters(
          new Vec3(1.0, 0.0, 0.0), new Vec3(-1.0, 0.0, 0.0)),
          force.getParticleParameters(0));

      force.setParticleParameters(0, new Vec3(0.5, 0.0, 0.0), new Vec3(-0.5, 0.0, 0.0));
      force.setEnergyFunction("2*(u1-u0)");
      force.setGlobalParameterName(0, "lambda2");
      force.setGlobalParameterDefaultValue(0, 0.6);
      assertEquals(new Vec3(0.5, 0.0, 0.0),
          force.getParticleParameters(0).displacement1());
      assertEquals("lambda2", force.getGlobalParameterName(0));
      assertEquals(0.6, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals("2*(u1-u0)", force.getEnergyFunction());
    }
  }
}
