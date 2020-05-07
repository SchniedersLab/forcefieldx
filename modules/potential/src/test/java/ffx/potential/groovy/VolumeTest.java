package ffx.potential.groovy;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

import ffx.potential.groovy.test.Gradient;
import groovy.lang.Binding;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/** JUnit Tests for the Volume Script */
public class VolumeTest {

  private Volume volume;
  private Binding binding;
  private double tolerance = 0.001;

  @After
  public void after() {
    volume.destroyPotentials();
    System.gc();
  }

  @Before
  public void before() {
    binding = new Binding();
    volume = new Volume();
    volume.setBinding(binding);
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
    assertEquals(79.85081214917432, volume.totalVolume, tolerance);
    assertEquals(98.14422780202705, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius. */
  @Test
  public void testConnollyMolecularCrambin() {
    String[] args = {"-c", "-m", "-p", "1.4", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(5222.628196815338, volume.totalVolume, tolerance);
    assertEquals(2326.375086471378, volume.totalSurfaceArea, tolerance);
  }

  /**
   * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
   */
  @Test
  public void testConnollyMolecularEthylbenzene() {
    String[] args = {
      "-c", "-p", "1.4", "-m", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(126.96932864947797, volume.totalVolume, tolerance);
    assertEquals(137.7753894764906, volume.totalSurfaceArea, tolerance);
  }

  /**
   * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
   */
  @Test
  public void testConnollyMolecularHydrogenButane() {
    String[] args = {
      "-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/butane.xyz"
    };
    binding.setVariable("args", args);
    // Evaluate the script.
    volume.run();
    assertEquals(119.0512090546068, volume.totalVolume, tolerance);
    assertEquals(132.598010667639, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius and hydrogen atoms. */
  @Test
  public void testConnollyMolecularHydrogenCrambin() {
    String[] args = {
      "-c", "-m", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/crambin.xyz"
    };
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(5890.028072142157, volume.totalVolume, tolerance);
    assertEquals(2456.835858765312, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly molecular surface area and volume with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollyMolecularHydrogenEthylbenzene() {
    String[] args = {
      "-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(165.83018235892447, volume.totalVolume, tolerance);
    assertEquals(171.91980284282084, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVButane() {
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);
    // Evaluate the script.
    volume.run();
    assertEquals(302.0524178356785, volume.totalVolume, tolerance);
    assertEquals(227.49657650050898, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVCrambin() {
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(8956.620463626994, volume.totalVolume, tolerance);
    assertEquals(3015.7687533888334, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVEthylbenzene() {
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(418.96380636214514, volume.totalVolume, tolerance);
    assertEquals(287.56961427851286, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenButane() {
    String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);
    // Evaluate the script.
    volume.run();
    assertEquals(402.708673112844, volume.totalVolume, tolerance);
    assertEquals(280.836964340470, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius and hydrogen atoms. */
  @Test
  public void testConnollySEVHydrogenCrambin() {
    String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(9804.262055253388, volume.totalVolume, tolerance);
    assertEquals(3142.106787870207, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenEthylbenzene() {
    String[] args = {
      "-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(518.612603965319, volume.totalVolume, tolerance);
    assertEquals(340.264998320387, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenEthylbenzeneDerivatives() {
    Gradient gradient = new Gradient();
    binding = new Binding();
    gradient.setBinding(binding);

    String[] args = {"--tol", "5.0e-2", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    gradient.run();

    assertEquals("Ethylbenzene gradient failures: ", 0, gradient.nFailures);
    gradient.destroyPotentials();
    System.gc();
  }

  /** Test Connolly vdW volume and surface area (probe radius = 0.0). */
  @Test
  public void testConnollyVDWCrambin() {
    String[] args = {"-c", "-p", "0.0", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(4418.303482563956, volume.totalVolume, tolerance);
    assertEquals(4168.547763834282, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol without hydrogen and a 0.4 A radii offset. */
  @Test
  public void testGaussVolButane() {
    String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script.
    volume.run();
    assertEquals(125.5120767378517, volume.totalVolume, tolerance);
    assertEquals(134.50524040566165, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol without hydrogen and a 0.0 A radii offset. */
  @Test
  public void testGaussVolCrambin() {
    String[] args = {"-o", "0.0", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(4371.667466648112, volume.totalVolume, tolerance);
    assertEquals(3971.0619085859435, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol derivatives. */
  @Test
  public void testGaussVolCrambinDerivatives() {
    Gradient gradient = new Gradient();
    binding = new Binding();
    gradient.setBinding(binding);

    System.setProperty("gkterm", "true");
    System.setProperty("cavmodel", "gauss-disp");

    // Choose a random atom to test.
    int atomID = (int) floor(random() * 642) + 1;

    String[] args = {
      "--ga", Integer.toString(atomID), "src/main/java/ffx/potential/structures/crambin.xyz"
    };
    binding.setVariable("args", args);

    // Evaluate the script
    gradient.run();

    System.clearProperty("gkterm");
    System.clearProperty("cavmodel");
    assertEquals("Crambin gradient failures: ", 0, gradient.nFailures);
    gradient.destroyPotentials();
    System.gc();
  }

  /** Test GaussVol without hydrogen and a 0.4 A radii offset. */
  @Test
  public void testGaussVolEthylbenzene() {
    String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    volume.run();
    assertEquals(194.44960348422916, volume.totalVolume, tolerance);
    assertEquals(193.55314214639066, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol derivatives. */
  @Test
  public void testGaussVolHydrogenEthylbenzeneDerivatives() {
    Gradient gradient = new Gradient();
    binding = new Binding();
    gradient.setBinding(binding);

    System.setProperty("cavmodel", "gauss-disp");

    String[] args = {"src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Evaluate the script
    gradient.run();

    System.clearProperty("cavmodel");
    assertEquals("Ethylbenzene gradient failures: ", 0, gradient.nFailures);
    gradient.destroyPotentials();
    System.gc();
  }
}
