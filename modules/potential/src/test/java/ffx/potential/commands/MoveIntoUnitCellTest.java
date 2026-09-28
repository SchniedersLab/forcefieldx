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
package ffx.potential.commands;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import ffx.potential.utils.PotentialTest;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

/** Test the Cart2Frac script. */
public class MoveIntoUnitCellTest extends PotentialTest {

  @Test
  public void testMoveIntoUnitCell() {
    // Set-up the input arguments for the MoveIntoUnitCell script.
    String filepath = getResourcePath("watertiny.xyz");
    String[] args = {filepath};
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
  public void testMoveIntoUnitCellArc() throws IOException {
    // Set-up the input arguments for the MoveIntoUnitCell script.
    // Each snapshot of watertiny.arc is watertiny.xyz translated by a lattice vector:
    // (0, 0, 0), (+a, 0, 0) and (0, 0, -a).
    String filepath = getResourcePath("watertiny.arc");
    String[] args = {filepath};
    binding.setVariable("args", args);
    File baseDir = registerTemporaryDirectory().toFile();
    binding.setVariable("baseDir", baseDir);

    // Contruct and evaluate the MoveIntoUnitCell script.
    MoveIntoUnitCell moveIntoUnitCell = new MoveIntoUnitCell(binding).run();
    potentialScript = moveIntoUnitCell;

    // The coordinates are from the final (third) snapshot of the archive.
    double[][] origCoordinates = moveIntoUnitCell.origCoordinates;
    assertNotNull(origCoordinates);
    assertEquals(81, origCoordinates.length);
    double tolerance = 1.0e-6;
    assertEquals(-0.382446, origCoordinates[6][0], tolerance);
    assertEquals(1.447602, origCoordinates[6][1], tolerance);
    assertEquals(-9.777606, origCoordinates[6][2], tolerance);

    // Moving into the unit cell gives the same result as for watertiny.xyz.
    double[][] unitCellCoordinates = moveIntoUnitCell.unitCellCoordinates;
    assertNotNull(unitCellCoordinates);
    assertEquals(81, unitCellCoordinates.length);
    assertEquals(8.93905400, unitCellCoordinates[6][0], tolerance);
    assertEquals(1.44760200, unitCellCoordinates[6][1], tolerance);
    assertEquals(8.86539400, unitCellCoordinates[6][2], tolerance);

    // All three snapshots should be written to the output archive.
    File saveFile = new File(baseDir, "watertiny.arc");
    assertTrue(saveFile.exists());
    // Each snapshot is written with a line of unit cell parameters (a, b, c, alpha, beta, gamma).
    long nSnapshots = Files.readAllLines(saveFile.toPath()).stream()
        .map(line -> line.trim().split(" +"))
        .filter(tokens -> tokens.length == 6 && tokens[3].equals("90.00000000"))
        .count();
    assertEquals(3, nSnapshots);
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
