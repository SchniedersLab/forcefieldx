package ffx.xray;

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
import ffx.utilities.DirectoryUtils;
import ffx.xray.groovy.ManyBody;

import groovy.lang.Binding;

/**
 * Tests X-Ray many body optimization and the X-Ray many body groovy script under varying parameters.
 *
 * @author Mallory R. Tollefson
 */
public class XRayManyBodyTest extends PJDependentTest {

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

    @Test
    public void testManyBodyGlobal() {
        // Set-up the input arguments for the script.
        String[] args = {"-a", "2", "-L", "2", "-s", "1", "--fi", "5",
                "src/main/java/ffx/xray/structures/5awl.pdb", "src/main/java/ffx/xray/structures/5awl.mtz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("ManyBodyTest");
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
        double expectedPotential = 2.8613333175581248E16;
        double actualPotential = list.get(0).getTotalEnergy();

        /**
         * The normalized difference between actual and expected potential should be very close to 0 no matter
         * the magnitude of the potentials.
         */
        double differenceNorm = (expectedPotential - actualPotential) / actualPotential;
        Assert.assertEquals(differenceNorm, 0, 1E-8);

        // Delete all created directories and files.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by ManyBodyTest.");
        }

        manyBody.getManyBody().getRestartFile().delete();
    }
}
