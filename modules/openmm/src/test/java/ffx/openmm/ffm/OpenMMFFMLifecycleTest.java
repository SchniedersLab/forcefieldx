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

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assume.assumeTrue;

/**
 * Native integration tests for the FFM OpenMM system, context, and state lifecycle.
 */
public class OpenMMFFMLifecycleTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testReferenceContextStateAndDestruction() {
    try (System system = new System();
         VerletIntegrator integrator = new VerletIntegrator(0.002);
         Platform platform = new Platform("Reference")) {
      system.addParticle(12.0);
      system.addParticle(16.0);

      Context context = new Context(system, integrator, platform);
      context.setPositions(new double[]{0.0, 0.0, 0.0, 0.1, 0.2, 0.3});
      context.setVelocities(new double[]{0.4, 0.5, 0.6, 0.7, 0.8, 0.9});
      context.setTime(1.25);
      context.setStepCount(17);

      int types = OpenMMNative.OpenMM_State_Positions()
          | OpenMMNative.OpenMM_State_Velocities();
      try (State state = context.getState(types, false)) {
        assertEquals(types, state.getDataTypes());
        assertArrayEquals(
            new double[]{0.0, 0.0, 0.0, 0.1, 0.2, 0.3}, state.getPositions(), 1.0e-12);
        assertArrayEquals(
            new double[]{0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, state.getVelocities(), 1.0e-12);
        assertEquals(1.25, state.getTime(), 0.0);
        assertEquals(17L, state.getStepCount());
      }

      context.destroy();
      context.destroy();
      assertNull(context.getIntegrator());
    }
  }
}
