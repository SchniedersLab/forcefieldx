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

import edu.rit.pj.Comm;

import ffx.utilities.DirectoryUtils;
import ffx.algorithms.PJDependentTest;
import ffx.algorithms.groovy.ManyBody;
import ffx.potential.PotentialComponent;

import java.util.logging.Level;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;

/**
 * Tests many body optimization and the many body groovy script under global, box and monte carlo parameter conditions.
 * @author Mallory R. Tollefson
 */
public class ManyBodyTest extends PJDependentTest {

    Binding binding;
    ManyBody manyBody;

    @Before
    public void before() {
        binding = new Binding();
        manyBody = new ManyBody();
        manyBody.setBinding(binding);
    }

    @After
    public void after() {
        manyBody.destroyPotentials();
        System.gc();
    }

    @Test
    public void testManyBodyHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        manyBody.run();
    }

    @Test
    public void testManyBodyGlobal() {
        // Set-up the input arguments for the script.
        String[] args = {"-a", "2", "-L", "2", "--tC", "2",
            "src/main/java/ffx/algorithms/structures/5awl.pdb"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("ManyBodyTest");
            manyBody.setBaseDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        manyBody.run();
        double expectedTotalPotential = -219.8836543404126;
        
        
        double actualTotalPotential = manyBody.getPotential().getEnergyComponent(PotentialComponent.ForceFieldEnergy);
        Assert.assertEquals(actualTotalPotential, expectedTotalPotential, 1E-8);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by ManyBodyTest.");
        }

        manyBody.getManyBody().getRestartFile().delete();
    }

    /**
     * Tests ManyBody.groovy and RotamerOptimization.java by running a box optimization simulation on a small pdb file.
     */
    @Test
    public void testManyBodyBoxOptimization() {
        // Set-up the input arguments for the script.
        String[] args = {"-a", "5", "-L", "2", "--bL", "10", "--bB", "2", "--tC", "2", "--pr", "2",
            "src/main/java/ffx/algorithms/structures/5awl.pdb"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("ManyBodyTest");
            manyBody.setBaseDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        manyBody.run();
        
        double expectedTotalPotential = -219.8836543404126;
        double actualTotalPotential = manyBody.getPotential().getEnergyComponent(PotentialComponent.ForceFieldEnergy);
        Assert.assertEquals(actualTotalPotential, expectedTotalPotential, 1E-7);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by ManyBodyTest.");
        }

        manyBody.getManyBody().getPartial().delete();
    }

    /**
     * Tests ManyBody.groovy and RotamerOptimization.java by running a monte carlo optimization simulation on a small pdb file. Elimination criteria are not used during this test. A monte carlo search is done on the permuatations the protein experience.
     */
    @Test
    public void testManyBodyMonteCarlo(){
        // Initialize Parallel Java
        try {
            Comm.world();
        } catch (IllegalStateException ise) {
            try {
                String args[] = new String[0];
                Comm.init(args);
            } catch (Exception e) {
                String message = String.format(" Exception starting up the Parallel Java communication layer.");
                logger.log(Level.WARNING, message, e.toString());
                message = String.format(" Skipping rotamer optimization test.");
                logger.log(Level.WARNING, message, e.toString());
                return;
            }
        }

        System.setProperty("polarization", "direct");

        // Set-up the input arguments for the script.
        String[] args = {"-a", "2", "-L", "2", "--tC", "2", "--pr", "2", "--mC", "10000",
                "src/main/java/ffx/algorithms/structures/5awl.pdb"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("ManyBodyTest");
            manyBody.setBaseDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        manyBody.setTesting(true);
        manyBody.setMonteCarloTesting(true);

        // Evaluate the script.
        manyBody.run();

        double expectedTotalPotential = -203.72294133789995;
        double actualTotalPotential = manyBody.getPotential().getEnergyComponent(PotentialComponent.ForceFieldEnergy);
        Assert.assertEquals(actualTotalPotential, expectedTotalPotential, 1E-7);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by ManyBodyTest.");
        }

        // Clear properties and delete unneccesary files.
        manyBody.setMonteCarloTesting(false);
        manyBody.getManyBody().getRestartFile().delete();
        System.clearProperty("polarization");
    }

}
