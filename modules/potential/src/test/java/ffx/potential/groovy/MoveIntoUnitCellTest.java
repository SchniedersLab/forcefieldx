// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.potential.groovy;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import ffx.potential.utils.PotentialTest;
import org.junit.Test;

/** Test the Cart2Frac script. */
public class MoveIntoUnitCellTest extends PotentialTest {

  @Test
  public void testMoveIntoUnitCell() {
    // Set-up the input arguments for the MoveIntoUnitCell script.
    String[] args = {"src/main/java/ffx/potential/structures/watertiny.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Contruct and evaluate the MoveIntoUnitCell script.
    MoveIntoUnitCell moveIntoUnitCell = new MoveIntoUnitCell(binding).run();
    potentialScript = moveIntoUnitCell;

    // Pull out the Cart2Frac results to check.
    double[][] origCoordinates = moveIntoUnitCell.origCoordinates;
    assertNotNull(origCoordinates);
    assertEquals(81, origCoordinates.length);
    double tolerance = 1.0e-6;
    assertEquals(-0.382446, origCoordinates[6][0], tolerance);
    assertEquals(1.447602, origCoordinates[6][1], tolerance);
    assertEquals(-0.456106, origCoordinates[6][2], tolerance);

    double[][] unitCellCoordinates = moveIntoUnitCell.unitCellCoordinates;
    assertNotNull(unitCellCoordinates);
    assertEquals(81, unitCellCoordinates.length);
    assertEquals(8.93905400, unitCellCoordinates[6][0], tolerance);
    assertEquals(1.44760200, unitCellCoordinates[6][1], tolerance);
    assertEquals(8.86539400, unitCellCoordinates[6][2], tolerance);
  }

  @Test
  public void testMoveIntoUnitCellHelp() {
    // Set-up the input arguments for the MoveIntoUnitCell script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Contruct and evaluate the MoveIntoUnitCell script.
    MoveIntoUnitCell moveIntoUnitCell = new MoveIntoUnitCell(binding).run();
    potentialScript = moveIntoUnitCell;
  }
}
