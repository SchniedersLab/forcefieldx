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

import ffx.openmm.ffm.bindings.OpenMMNative;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native Reference-platform tests for {@link HarmonicBondForce}. */
public class OpenMMFFMHarmonicBondTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testParametersEnergyAndContextUpdate() {
    double length = 0.1;
    double k = 1000.0;
    double displacement = 0.01;
    try (HarmonicBondForce force = new HarmonicBondForce();
         System system = new System();
         VerletIntegrator integrator = new VerletIntegrator(0.001);
         Platform platform = new Platform("Reference")) {
      int bond = force.addBond(0, 1, length, k);
      assertEquals(0, bond);
      assertEquals(1, force.getNumBonds());
      assertEquals(new HarmonicBondForce.BondParameters(0, 1, length, k),
          force.getBondParameters(bond));
      assertFalse(force.usesPeriodicBoundaryConditions());
      force.setUsesPeriodicBoundaryConditions(true);
      assertTrue(force.usesPeriodicBoundaryConditions());

      system.addParticle(12.0);
      system.addParticle(12.0);
      assertEquals(0, system.addForce(force));

      try (Context context = new Context(system, integrator, platform)) {
        context.setPositions(new double[]{0.0, 0.0, 0.0, length + displacement, 0.0, 0.0});
        assertEquals(0.5 * k * displacement * displacement, potentialEnergy(context), 1.0e-12);

        force.setBondParameters(bond, 0, 1, length + displacement, k);
        force.updateParametersInContext(context);
        assertEquals(0.0, potentialEnergy(context), 1.0e-12);
      }
    }
  }

  private static double potentialEnergy(Context context) {
    try (State state = context.getState(OpenMMNative.OpenMM_State_Energy(), false)) {
      return state.getPotentialEnergy();
    }
  }
}
