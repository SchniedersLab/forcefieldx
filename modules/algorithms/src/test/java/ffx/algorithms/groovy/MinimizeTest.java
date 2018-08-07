/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.algorithms.groovy;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import groovy.lang.Binding;

import ffx.utilities.DirectoryUtils;
import ffx.algorithms.PJDependentTest;
import ffx.algorithms.groovy.Minimizer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;

/**
 * Tests many body optimization and the many body groovy script under global, box and monte carlo parameter conditions.
 * @author Mallory R. Tollefson
 */
public class MinimizeTest extends PJDependentTest {

    Binding binding;
    Minimizer minimize;

    @Before
    public void before() {
        binding = new Binding();
        minimize = new Minimizer();
        minimize.setBinding(binding);
    }

    @After
    public void after() {
        minimize.destroyPotentials();
        System.gc();
    }

    @Test
    public void testMinimizeHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        minimize.run();
    }

    /**
     * Tests convergence criteria flag of the minimize class.
     */
    @Test
    public void testMinimizeConvergenceCriteria() {
        // Set-up the input arguments for the script.
        String[] args = {"-e", "2", "src/main/java/ffx/algorithms/structures/5awl.pdb"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("MinimizeTest");
            minimize.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        minimize.run();
        double expectedTotalPotential = -277.62004110;


        double actualTotalPotential = minimize.getPotentials().get(minimize.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(actualTotalPotential, expectedTotalPotential, 1E-8);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by MinimizeTest.");
        }
    }

    /**
     * Tests the iterations flag of the minimize class.
     */
    @Test
    public void testMinimizeIterations() {
        // Set-up the input arguments for the script.
        String[] args = {"-I", "5", "src/main/java/ffx/algorithms/structures/5awl.pdb"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("MinimizeTest");
            minimize.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        minimize.run();
        double expectedTotalPotential = -262.73461486;


        double actualTotalPotential = minimize.getPotentials().get(minimize.getPotentials().size() - 1).getTotalEnergy();
        Assert.assertEquals(actualTotalPotential, expectedTotalPotential, 1E-8);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by MinimizeTest.");
        }
    }
}

