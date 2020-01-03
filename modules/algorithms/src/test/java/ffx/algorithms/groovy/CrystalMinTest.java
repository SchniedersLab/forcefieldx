//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
//******************************************************************************
package ffx.algorithms.groovy;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import ffx.algorithms.misc.PJDependentTest;
import ffx.algorithms.groovy.CrystalMin;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * Tests CrystalMin and CrystalMin.groovy scripts for -e -I -f and -c flags.
 *
 * @author Aaron J. Nessler
 * @author Mallory R. Tollefson
 */
public class CrystalMinTest extends PJDependentTest {

    Binding binding;
    CrystalMin xtalMin;

    @Before
    public void before() {
        binding = new Binding();
        xtalMin = new CrystalMin();
        xtalMin.setBinding(binding);
    }

    @After
    public void after() {
        xtalMin.destroyPotentials();
        System.gc();
    }

    @Test
    public void testCrystalMinHelp() {
        // Set-up the input arguments for the CrystalMin script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        xtalMin.run();
    }

    /**
     * Tests convergence criteria flag of the CrystalMin class.
     */
    @Test
    public void testCrystalMinConvergenceCriteria() {
        // Set-up the input arguments for the script.
        String[] args = {"-e", "0.25", "src/main/java/ffx/algorithms/structures/acetamide.xtal.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("CrystalMinTest");
            xtalMin.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        xtalMin.run();

        double expectedPotentialEnergy = -32.72657252765883;

        double actualPotentialEnergy = xtalMin.getPotentials().get(xtalMin.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(expectedPotentialEnergy, actualPotentialEnergy,1E-6);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by CrystalMinTest.");
        }
    }

    /**
     * Tests the iterations flag of the CrystalMin class.
     */
    @Test
    public void testCrystalMinIterations() {
        // Set-up the input arguments for the script.
        String[] args = {"-I", "1", "-e", "0.25", "src/main/java/ffx/algorithms/structures/acetamide.xtal.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("CrystalMinTest");
            xtalMin.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        xtalMin.run();

        double expectedPotentialEnergy = -32.63130589380454;

        double actualPotentialEnergy = xtalMin.getPotentials().get(xtalMin.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(expectedPotentialEnergy, actualPotentialEnergy, 1E-6);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by CrystalMinTest.");
        }
    }

    /**
     * Tests the fractional flag of the CrystalMin class.
     */
    @Test
    public void testCrystalMinFractional() {
        // Set-up the input arguments for the script.
        String[] args = {"-f", "ATOM", "-e", "0.24", "src/main/java/ffx/algorithms/structures/acetamide.xtal.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("CrystalMinTest");
            xtalMin.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        xtalMin.run();

        double expectedPotentialEnergy = -32.54442936288904;

        double actualPotentialEnergy = xtalMin.getPotentials().get(xtalMin.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(expectedPotentialEnergy, actualPotentialEnergy, 1E-6);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by CrystalMinTest.");
        }
    }

    /**
     * Tests the fractional flag of the CrystalMin class.
     */
    @Test
    public void testCrystalMinCoords() {
        // Set-up the input arguments for the script.
        String[] args = {"-c", "src/main/java/ffx/algorithms/structures/acetamide.xtal.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("CrystalMinTest");
            xtalMin.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        xtalMin.run();

        double expectedPotentialEnergy = -32.535694367226576;

        double actualPotentialEnergy = xtalMin.getPotentials().get(xtalMin.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(expectedPotentialEnergy, actualPotentialEnergy,1E-6);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by CrystalMinTest.");
        }
    }
}


