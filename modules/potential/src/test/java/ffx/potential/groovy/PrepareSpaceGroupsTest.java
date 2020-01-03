/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2020.
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
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import ffx.potential.groovy.PrepareSpaceGroups;
import ffx.utilities.BaseFFXTest;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * Test the Energy script.
 */
public class PrepareSpaceGroupsTest extends BaseFFXTest {

    Binding binding;
    PrepareSpaceGroups prepareSpaceGroups;

    private static final Logger logger = Logger.getLogger(PrepareSpaceGroupsTest.class.getName());

    @Before
    public void before() {
        binding = new Binding();
        prepareSpaceGroups = new PrepareSpaceGroups();
        prepareSpaceGroups.setBinding(binding);
    }

    @After
    public void after() {
        prepareSpaceGroups.destroyPotentials();
        System.gc();
    }

    @Test
    public void testPrepareSpaceGroupHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        prepareSpaceGroups.run();

        // Pull out the biotype results to check.
        Assert.assertEquals(0, prepareSpaceGroups.numberCreated);
    }

    @Test
    public void testPrepareSpaceGroups() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"src/main/java/ffx/potential/structures/paracetamol.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("spacegroups");
            prepareSpaceGroups.baseDir = path.toFile();
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory");
        }

        // Evaluate the script.
        prepareSpaceGroups.run();

        // Pull out the Cart2Frac results to check.
        Assert.assertEquals(230, prepareSpaceGroups.numberCreated);

        // Delate all created space grouop directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by PrepareSpaceGroups.");
        }

    }
    @Test
    public void testPrepareSpaceGroupsChiral() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-c","src/main/java/ffx/potential/structures/paracetamol.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("spacegroups");
            prepareSpaceGroups.baseDir = path.toFile();
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory");
        }

        // Evaluate the script.
        prepareSpaceGroups.run();

        // Pull out the Cart2Frac results to check.
        Assert.assertEquals(65, prepareSpaceGroups.numberCreated);

        // Delete all created space group directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by PrepareSpaceGroups.");
        }

    }
    @Test
    public void testPrepareSpaceGroupsAChiral() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-a","src/main/java/ffx/potential/structures/paracetamol.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("spacegroups");
            prepareSpaceGroups.baseDir = path.toFile();
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory");
        }

        // Evaluate the script.
        prepareSpaceGroups.run();

        // Pull out the Cart2Frac results to check.
        Assert.assertEquals(165, prepareSpaceGroups.numberCreated);

        // Delete all created space group directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by PrepareSpaceGroups.");
        }

    }
}
