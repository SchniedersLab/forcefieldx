package ffx.potential.groovy;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

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

    /**
     *  Test GaussVol without hydrogen and a 0.4 A radii offset.
     */
    @Test
    public void testGaussVolButane() {
        String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script.
        volume.run();
        assertEquals(125.5120767378517, volume.totalVolume, 0.01);
        assertEquals(134.25693852184395, volume.totalSurfaceArea, 0.01);
    }

    /**
     *  Test GaussVol without hydrogen and a 0.4 A radii offset.
     */
    @Test
    public void testGaussVolEthylbenzene() {
        String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();
        assertEquals(194.44960348422916, volume.totalVolume, 0.001);
        assertEquals(193.15447011592823, volume.totalSurfaceArea, 0.001);
    }

    /**
     * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollyMolecularButane() {
        String[] args = {"-c", "-p", "1.4", "-m", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        volume.run();
        assertEquals(79.85081214917432, volume.totalVolume, 0.01);
        assertEquals(98.14422780202705, volume.totalSurfaceArea, 0.01);
    }

    /**
     * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollyMolecularHydrogenButane() {
        String[] args = {"-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        volume.run();
        assertEquals(119.0512090546068, volume.totalVolume, 0.01);
        assertEquals(132.598010667639, volume.totalSurfaceArea, 0.01);
    }

    /**
     * Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollySEVButane() {
        String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        volume.run();
        assertEquals(302.0524178356785, volume.totalVolume, 0.01);
        assertEquals(227.49657650050898, volume.totalSurfaceArea, 0.01);
    }

    /**
     * Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollySEVHydrogenButane() {
        String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/butane.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        volume.run();
        assertEquals(402.708673112844, volume.totalVolume, 0.01);
        assertEquals(280.836964340470, volume.totalSurfaceArea, 0.01);
    }

    /**
     * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollyMolecularEthylbenzene() {
        String[] args = {"-c", "-p", "1.4", "-m", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();
        assertEquals(126.96932864947797, volume.totalVolume, 0.001);
        assertEquals(137.7753894764906, volume.totalSurfaceArea, 0.001);
    }

    /**
     * Test Connolly molecular surface area and volume with hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollyMolecularHydrogenEthylbenzene() {
        String[] args = {"-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();
        assertEquals(165.83018235892447, volume.totalVolume, 0.001);
        assertEquals(171.91980284282084, volume.totalSurfaceArea, 0.001);
    }

    /**
     * Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollySEVEthylbenzene() {
        String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();
        assertEquals(418.96380636214514, volume.totalVolume, 0.001);
        assertEquals(287.56961427851286, volume.totalSurfaceArea, 0.001);
    }

    /**
     * Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius.
     */
    @Test
    public void testConnollySEVHydrogenEthylbenzene() {
        String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script
        volume.run();
        assertEquals(518.612603965319, volume.totalVolume, 0.001);
        assertEquals(340.264998320387, volume.totalSurfaceArea, 0.001);
    }


}
