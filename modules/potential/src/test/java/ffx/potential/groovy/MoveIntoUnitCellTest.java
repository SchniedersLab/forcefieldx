/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.groovy;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import ffx.potential.groovy.MoveIntoUnitCell;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * Test the Cart2Frac script.
 */
public class MoveIntoUnitCellTest {

    Binding binding;
    MoveIntoUnitCell moveIntoUnitCell;

    @Before
    public void before() {
        binding = new Binding();
        moveIntoUnitCell = new MoveIntoUnitCell();
        moveIntoUnitCell.setBinding(binding);
    }

    @After
    public void after() {
        moveIntoUnitCell.destroyPotentials();
        System.gc();
    }

    @Test
    public void testMoveIntoUnitCellHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        moveIntoUnitCell.run();

        // Pull out the biotype results to check.
        Assert.assertNull(moveIntoUnitCell.origCoordinates);
        Assert.assertNull(moveIntoUnitCell.unitCellCoordinates);
    }

    @Test
    public void testMoveIntoUnitCell() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"src/main/java/ffx/potential/structures/watertiny.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("MoveIntoUnitCell");
            moveIntoUnitCell.setBaseDir(path.toFile());
        } catch (java.io.IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        moveIntoUnitCell.run();

        // Pull out the Cart2Frac results to check.
        double origCoordinates[][] = moveIntoUnitCell.origCoordinates;
        org.junit.Assert.assertNotNull(origCoordinates);
        org.junit.Assert.assertEquals(81, origCoordinates.length);
        // 7  O     -0.382446    1.447602   -0.456106     1     8     9
        org.junit.Assert.assertEquals(-0.382446, origCoordinates[6][0], 1.0e-6);
        org.junit.Assert.assertEquals(1.447602, origCoordinates[6][1], 1.0e-6);
        org.junit.Assert.assertEquals(-0.456106, origCoordinates[6][2], 1.0e-6);

        double unitCellCoordinates[][] = moveIntoUnitCell.unitCellCoordinates;
        org.junit.Assert.assertNotNull(unitCellCoordinates);
        org.junit.Assert.assertEquals(81, unitCellCoordinates.length);
        //  7   O    8.93905400    1.44760200    8.86539400     1       8       9
        org.junit.Assert.assertEquals(8.93905400, unitCellCoordinates[6][0], 1.0e-6);
        org.junit.Assert.assertEquals(1.44760200, unitCellCoordinates[6][1], 1.0e-6);
        org.junit.Assert.assertEquals(8.86539400, unitCellCoordinates[6][2], 1.0e-6);

        // Delate all created space grouop directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by Frac2Cart.");
        }
    }
}
