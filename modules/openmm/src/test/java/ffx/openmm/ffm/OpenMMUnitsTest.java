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

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Unit tests for OpenMM unit conversion constants.
 */
public class OpenMMUnitsTest {

  private static final double TOLERANCE = 1.0e-15;

  @Test
  public void testUnitConstants() {
    assertEquals(0.1, OpenMMUnits.NM_PER_ANGSTROM, TOLERANCE);
    assertEquals(10.0, OpenMMUnits.ANGSTROMS_PER_NM, TOLERANCE);
    assertEquals(0.001, OpenMMUnits.PS_PER_FS, TOLERANCE);
    assertEquals(1000.0, OpenMMUnits.FS_PER_PS, TOLERANCE);
    assertEquals(4.184, OpenMMUnits.KJ_PER_KCAL, TOLERANCE);
    assertEquals(1.0 / 4.184, OpenMMUnits.KCAL_PER_KJ, TOLERANCE);
    assertEquals(3.1415926535897932385 / 180.0, OpenMMUnits.RADIANS_PER_DEGREE, TOLERANCE);
    assertEquals(180.0 / 3.1415926535897932385, OpenMMUnits.DEGREES_PER_RADIAN, TOLERANCE);
    assertEquals(1.7817974362806786095, OpenMMUnits.SIGMA_PER_VDW_RADIUS, TOLERANCE);
    assertEquals(0.56123102415468649070, OpenMMUnits.VDW_RADIUS_PER_SIGMA, TOLERANCE);
  }

  @Test
  public void testReciprocalIdentities() {
    assertEquals(1.0, OpenMMUnits.NM_PER_ANGSTROM * OpenMMUnits.ANGSTROMS_PER_NM, TOLERANCE);
    assertEquals(1.0, OpenMMUnits.PS_PER_FS * OpenMMUnits.FS_PER_PS, TOLERANCE);
    assertEquals(1.0, OpenMMUnits.KJ_PER_KCAL * OpenMMUnits.KCAL_PER_KJ, TOLERANCE);
    assertEquals(1.0, OpenMMUnits.RADIANS_PER_DEGREE * OpenMMUnits.DEGREES_PER_RADIAN, TOLERANCE);
    assertEquals(1.0, OpenMMUnits.SIGMA_PER_VDW_RADIUS * OpenMMUnits.VDW_RADIUS_PER_SIGMA, 1.0e-14);
  }
}
