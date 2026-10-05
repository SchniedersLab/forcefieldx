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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native parameter round-trip tests for the scalar custom-force façades. */
public class OpenMMFFMCustomScalarForceTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testCustomBondForce() {
    try (CustomBondForce force = new CustomBondForce("k*(r-r0)^2")) {
      assertEquals("k*(r-r0)^2", force.getEnergyFunction());
      assertEquals(0, force.addPerBondParameter("k"));
      assertEquals(1, force.addPerBondParameter("r0"));
      assertEquals(0, force.addGlobalParameter("scale", 2.0));
      force.addEnergyParameterDerivative("scale");
      assertEquals(0, force.addBond(2, 5, new double[]{3.0, 1.2}));
      assertEquals(1, force.getNumBonds());
      assertEquals(2, force.getNumPerBondParameters());
      assertEquals(1, force.getNumEnergyParameterDerivatives());
      CustomBondForce.BondParameters initial = force.getBondParameters(0);
      assertEquals(2, initial.particle1());
      assertEquals(5, initial.particle2());
      assertArrayEquals(new double[]{3.0, 1.2}, initial.parameters(), 0.0);
      force.setBondParameters(0, 1, 4, new double[]{6.0, 1.5});
      CustomBondForce.BondParameters updated = force.getBondParameters(0);
      assertEquals(1, updated.particle1());
      assertEquals(4, updated.particle2());
      assertArrayEquals(new double[]{6.0, 1.5}, updated.parameters(), 0.0);
      force.setEnergyFunction("2*k*(r-r0)^2");
      force.setGlobalParameterName(0, "scale2");
      force.setGlobalParameterDefaultValue(0, 3.0);
      force.setPerBondParameterName(0, "spring");
      assertEquals("2*k*(r-r0)^2", force.getEnergyFunction());
      assertEquals("scale2", force.getGlobalParameterName(0));
      assertEquals(3.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals("spring", force.getPerBondParameterName(0));
      force.setUsesPeriodicBoundaryConditions(true);
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomAngleAndTorsionForces() {
    try (CustomAngleForce angle = new CustomAngleForce("k*(theta-theta0)^2");
         CustomTorsionForce torsion = new CustomTorsionForce("k*(1-cos(theta-theta0))")) {
      angle.addPerAngleParameter("k");
      angle.addPerAngleParameter("theta0");
      angle.addGlobalParameter("scale", 1.0);
      angle.addEnergyParameterDerivative("scale");
      angle.addAngle(0, 1, 2, new double[]{2.0, 1.1});
      CustomAngleForce.AngleParameters angleParameters = angle.getAngleParameters(0);
      assertEquals(0, angleParameters.particle1());
      assertEquals(1, angleParameters.particle2());
      assertEquals(2, angleParameters.particle3());
      assertArrayEquals(new double[]{2.0, 1.1}, angleParameters.parameters(), 0.0);
      angle.setAngleParameters(0, 3, 4, 5, new double[]{3.0, 1.4});
      assertEquals(3, angle.getAngleParameters(0).particle1());
      assertEquals("theta0", angle.getPerAngleParameterName(1));
      assertEquals("scale", angle.getEnergyParameterDerivativeName(0));
      angle.setUsesPeriodicBoundaryConditions(true);
      assertTrue(angle.usesPeriodicBoundaryConditions());

      torsion.addPerTorsionParameter("k");
      torsion.addPerTorsionParameter("theta0");
      torsion.addGlobalParameter("scale", 1.0);
      torsion.addEnergyParameterDerivative("scale");
      torsion.addTorsion(0, 1, 2, 3, new double[]{4.0, 0.5});
      CustomTorsionForce.TorsionParameters torsionParameters = torsion.getTorsionParameters(0);
      assertEquals(0, torsionParameters.particle1());
      assertEquals(3, torsionParameters.particle4());
      assertArrayEquals(new double[]{4.0, 0.5}, torsionParameters.parameters(), 0.0);
      torsion.setTorsionParameters(0, 4, 5, 6, 7, new double[]{5.0, 0.8});
      assertEquals(4, torsion.getTorsionParameters(0).particle1());
      assertEquals("theta0", torsion.getPerTorsionParameterName(1));
      assertEquals("scale", torsion.getEnergyParameterDerivativeName(0));
      torsion.setUsesPeriodicBoundaryConditions(true);
      assertTrue(torsion.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomExternalForce() {
    try (CustomExternalForce force = new CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")) {
      assertEquals(0, force.addPerParticleParameter("x0"));
      assertEquals(1, force.addPerParticleParameter("y0"));
      assertEquals(2, force.addPerParticleParameter("z0"));
      assertEquals(0, force.addGlobalParameter("k", 10.0));
      force.addParticle(7, new double[]{1.0, 2.0, 3.0});
      CustomExternalForce.ParticleParameters initial = force.getParticleParameters(0);
      assertEquals(7, initial.particle());
      assertArrayEquals(new double[]{1.0, 2.0, 3.0}, initial.parameters(), 0.0);
      force.setParticleParameters(0, 8, new double[]{4.0, 5.0, 6.0});
      CustomExternalForce.ParticleParameters updated = force.getParticleParameters(0);
      assertEquals(8, updated.particle());
      assertArrayEquals(new double[]{4.0, 5.0, 6.0}, updated.parameters(), 0.0);
      force.setEnergyFunction("k*(x^2+y^2+z^2)");
      force.setGlobalParameterName(0, "spring");
      force.setGlobalParameterDefaultValue(0, 20.0);
      force.setPerParticleParameterName(2, "zTarget");
      assertEquals("k*(x^2+y^2+z^2)", force.getEnergyFunction());
      assertEquals("spring", force.getGlobalParameterName(0));
      assertEquals(20.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals("zTarget", force.getPerParticleParameterName(2));
      assertFalse(force.usesPeriodicBoundaryConditions());
    }
  }
}
