package ffx.potential.groovy;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import ffx.potential.groovy.Volume;

import groovy.lang.Binding;

/**
 * JUnit Tests for the Volume Script
 */
public class VolumeTest {

    private Volume volume;
    private Binding binding;

    @Before
    public void before() {
        binding = new Binding();
        System.clearProperty("platform");

        volume = new Volume();
        volume.setBinding(binding);
    }

    @After
    public void after() {
        volume.destroyPotentials();
        System.gc();
    }

    @Test
    public void testVolumeButane() {
        // Set-up the input arguments for the Volume script using butane as a test case
        // *Note: updating this test case will require updating butane.xyz, butane.properties, and the
        // parameter file referenced in butane.properties (in this case, the amoeba09.prm file.
        String[] args = {"-p", "1.4", "-o", "0.4", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script.
        volume.run();

        org.junit.Assert.assertEquals(383.2278278467365, volume.totalVolume, 0.01);
        org.junit.Assert.assertEquals(134.25693852184395, volume.totalSurfaceArea, 0.01);

    }

    @Test
    public void testVolumeEthylbenzene() {
        // Set-up the input arguments for the Volume script using butane as a test case
        // *Note: updating this test case will require updating ethylbenzene.xyz, ethylbenzene.properties, and the
        // parameter file referenced in ethylbenzene.properties (in this case, the amoeba09.prm file.
        String[] args = {"-p", "1.4", "-o", "0.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();

        org.junit.Assert.assertEquals(521.7024513160144, volume.totalVolume, 0.001);
        org.junit.Assert.assertEquals(193.15447011592823, volume.totalSurfaceArea, 0.001);

    }

}
