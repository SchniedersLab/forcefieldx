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
import static org.junit.Assume.assumeTrue;

/** Native parameter round-trip tests for conventional FFM bonded forces. */
public class OpenMMFFMConventionalBondedForceTest {
  @BeforeClass public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test public void testParameterRoundTrips() {
    try (HarmonicAngleForce angle = new HarmonicAngleForce();
         PeriodicTorsionForce periodic = new PeriodicTorsionForce();
         RBTorsionForce rb = new RBTorsionForce()) {
      assertEquals(new HarmonicAngleForce.AngleParameters(0,1,2,1.5,4.0),
          angle.getAngleParameters(angle.addAngle(0,1,2,1.5,4.0)));
      assertEquals(new PeriodicTorsionForce.TorsionParameters(0,1,2,3,3,0.4,2.0),
          periodic.getTorsionParameters(periodic.addTorsion(0,1,2,3,3,0.4,2.0)));
      assertEquals(new RBTorsionForce.TorsionParameters(0,1,2,3,1,2,3,4,5,6),
          rb.getTorsionParameters(rb.addTorsion(0,1,2,3,1,2,3,4,5,6)));
    }
  }
}
