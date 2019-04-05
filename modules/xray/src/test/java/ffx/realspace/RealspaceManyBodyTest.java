//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.realspace;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.testng.Assert;

import ffx.algorithms.misc.PJDependentTest;
import ffx.numerics.Potential;
import ffx.realspace.groovy.ManyBody;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * Tests realspace many body optimization and the realspace many body groovy script under varying parameters.
 *
 * @author Mallory R. Tollefson
 */
public class RealspaceManyBodyTest extends PJDependentTest {

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
        try {
            manyBody.run();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    //@Test
    public void testManyBodyGlobal() {
        // Set-up the input arguments for the script.
        String[] args = {"-a", "2", "-L", "2", "-s", "1", "--fi", "5",
                "src/main/java/ffx/realspace/structures/file.pdb", "src/main/java/ffx/realspace/structures/file.mtz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("RealspaceManyBodyTest");
            manyBody.setSaveDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        try {
            manyBody.run();
        } catch (AssertionError ex) {
            ex.printStackTrace();
            throw ex;
        }

        List<Potential> list = manyBody.getPotentials();
        double expectedPotential = 0;
        double actualPotential = list.get(0).getTotalEnergy();

        // The normalized difference between actual and expected potential should be very close to 0 no matter
        // the magnitude of the potentials.
        double differenceNorm = (expectedPotential - actualPotential) / actualPotential;
        Assert.assertEquals(differenceNorm, 0, 1E-8);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by RealspaceManyBodyTest.");
        }

        manyBody.getManyBody().getRestartFile().delete();
    }
}
