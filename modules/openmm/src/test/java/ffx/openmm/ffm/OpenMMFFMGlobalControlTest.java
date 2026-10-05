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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native property tests for FFM global force controls. */
public class OpenMMFFMGlobalControlTest {
  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testControlProperties() {
    try (CMMotionRemover remover = new CMMotionRemover(5);
         AndersenThermostat thermostat = new AndersenThermostat(298.15, 1.0)) {
      assertEquals(5, remover.getFrequency());
      remover.setFrequency(9);
      assertEquals(9, remover.getFrequency());
      assertFalse(remover.usesPeriodicBoundaryConditions());

      thermostat.setDefaultTemperature(310.0);
      thermostat.setDefaultCollisionFrequency(2.5);
      thermostat.setRandomNumberSeed(71);
      assertEquals(310.0, thermostat.getDefaultTemperature(), 0.0);
      assertEquals(2.5, thermostat.getDefaultCollisionFrequency(), 0.0);
      assertEquals(71, thermostat.getRandomNumberSeed());
      assertFalse(thermostat.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testMonteCarloBarostatProperties() {
    try (MonteCarloBarostat barostat = new MonteCarloBarostat(1.0, 298.15, 25)) {
      barostat.setDefaultPressure(1.2);
      barostat.setDefaultTemperature(310.0);
      barostat.setFrequency(0);
      barostat.setRandomNumberSeed(19);
      assertEquals(1.2, barostat.getDefaultPressure(), 0.0);
      assertEquals(310.0, barostat.getDefaultTemperature(), 0.0);
      assertEquals(0, barostat.getFrequency());
      assertEquals(19, barostat.getRandomNumberSeed());
    }
  }

  @Test
  public void testMonteCarloAnisotropicBarostatProperties() {
    try (MonteCarloAnisotropicBarostat barostat =
        new MonteCarloAnisotropicBarostat(new Vec3(1.0, 2.0, 3.0), 298.15, true, false,
            true, 25)) {
      assertEquals(new Vec3(1.0, 2.0, 3.0), barostat.getDefaultPressure());
      assertTrue(barostat.getScaleX());
      assertFalse(barostat.getScaleY());
      assertTrue(barostat.getScaleZ());
      barostat.setDefaultPressure(new Vec3(4.0, 5.0, 6.0));
      assertEquals(new Vec3(4.0, 5.0, 6.0), barostat.getDefaultPressure());
    }
  }

  @Test
  public void testMonteCarloFlexibleBarostatProperties() {
    try (MonteCarloFlexibleBarostat barostat =
        new MonteCarloFlexibleBarostat(1.0, 298.15, 25, true)) {
      assertTrue(barostat.getScaleMoleculesAsRigid());
      barostat.setScaleMoleculesAsRigid(false);
      barostat.setDefaultPressure(1.5);
      barostat.setDefaultTemperature(310.0);
      barostat.setFrequency(0);
      assertFalse(barostat.getScaleMoleculesAsRigid());
      assertEquals(1.5, barostat.getDefaultPressure(), 0.0);
      assertEquals(310.0, barostat.getDefaultTemperature(), 0.0);
      assertEquals(0, barostat.getFrequency());
    }
  }

  @Test
  public void testMonteCarloMembraneBarostatModes() {
    try (MonteCarloMembraneBarostat barostat =
        new MonteCarloMembraneBarostat(1.0, 0.1, 298.15,
            MonteCarloMembraneBarostat.XYMode.ANISOTROPIC,
            MonteCarloMembraneBarostat.ZMode.CONSTANT_VOLUME, 25)) {
      assertEquals(MonteCarloMembraneBarostat.XYMode.ANISOTROPIC, barostat.getXYMode());
      assertEquals(MonteCarloMembraneBarostat.ZMode.CONSTANT_VOLUME, barostat.getZMode());
      barostat.setXYMode(MonteCarloMembraneBarostat.XYMode.ISOTROPIC);
      barostat.setZMode(MonteCarloMembraneBarostat.ZMode.FREE);
      assertEquals(MonteCarloMembraneBarostat.XYMode.ISOTROPIC, barostat.getXYMode());
      assertEquals(MonteCarloMembraneBarostat.ZMode.FREE, barostat.getZMode());
      barostat.setZMode(MonteCarloMembraneBarostat.ZMode.FIXED);
      assertEquals(MonteCarloMembraneBarostat.ZMode.FIXED, barostat.getZMode());
    }
  }
}
