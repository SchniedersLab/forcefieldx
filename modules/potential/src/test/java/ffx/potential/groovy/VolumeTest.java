// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.potential.groovy;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

import ffx.potential.groovy.test.Gradient;
import ffx.potential.utils.PotentialTest;
import org.junit.Test;

/** JUnit Tests for the Volume Script */
public class VolumeTest extends PotentialTest {

  private final double tolerance = 0.001;

  /**
   * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
   */
  @Test
  public void testConnollyMolecularButane() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "-m", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);
    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(79.85081214917432, volume.totalVolume, tolerance);
    assertEquals(98.14422780202705, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius. */
  @Test
  public void testConnollyMolecularCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-m", "-p", "1.4", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(5222.628196815338, volume.totalVolume, tolerance);
    assertEquals(2326.375086471378, volume.totalSurfaceArea, tolerance);
  }

  /**
   * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
   */
  @Test
  public void testConnollyMolecularEthylbenzene() {
    // Configure input arguments for the Volume script.
    String[] args = {
        "-c", "-p", "1.4", "-m", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(126.96932864947797, volume.totalVolume, tolerance);
    assertEquals(137.7753894764906, volume.totalSurfaceArea, tolerance);
  }

  /**
   * Test Connolly molecular surface area and volume without hydrogen and a 1.4 A exclude radius.
   */
  @Test
  public void testConnollyMolecularHydrogenButane() {
    // Configure input arguments for the Volume script.
    String[] args = {
        "-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/butane.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(119.0512090546068, volume.totalVolume, tolerance);
    assertEquals(132.598010667639, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius and hydrogen atoms. */
  @Test
  public void testConnollyMolecularHydrogenCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {
        "-c", "-m", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/crambin.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(5890.028072142157, volume.totalVolume, tolerance);
    assertEquals(2456.835858765312, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly molecular surface area and volume with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollyMolecularHydrogenEthylbenzene() {
    // Configure input arguments for the Volume script.
    String[] args = {
        "-c", "-p", "1.4", "-m", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(165.83018235892447, volume.totalVolume, tolerance);
    assertEquals(171.91980284282084, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVButane() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(302.0524178356785, volume.totalVolume, tolerance);
    assertEquals(227.49657650050898, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(8956.620463626994, volume.totalVolume, tolerance);
    assertEquals(3015.7687533888334, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA without hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVEthylbenzene() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(418.96380636214514, volume.totalVolume, tolerance);
    assertEquals(287.56961427851286, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenButane() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(402.708673112844, volume.totalVolume, tolerance);
    assertEquals(280.836964340470, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with a 1.4 A exclude radius and hydrogen atoms. */
  @Test
  public void testConnollySEVHydrogenCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(9804.262055253388, volume.totalVolume, tolerance);
    assertEquals(3142.106787870207, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenEthylbenzene() {
    // Configure input arguments for the Volume script.
    String[] args = {
        "-c", "-p", "1.4", "-y", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(518.612603965319, volume.totalVolume, tolerance);
    assertEquals(340.264998320387, volume.totalSurfaceArea, tolerance);
  }

  /** Test Connolly SEV & SASA with hydrogen and a 1.4 A exclude radius. */
  @Test
  public void testConnollySEVHydrogenEthylbenzeneDerivatives() {
    // Configure input arguments for the Gradient script.
    String[] args = {"--tol", "5.0e-2", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Gradient script
    Gradient gradient = new Gradient(binding).run();
    potentialScript = gradient;
    assertEquals("Ethylbenzene gradient failures: ", 0, gradient.nFailures);
  }

  /** Test Connolly vdW volume and surface area (probe radius = 0.0). */
  @Test
  public void testConnollyVDWCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {"-c", "-p", "0.0", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(4418.303482563956, volume.totalVolume, tolerance);
    assertEquals(4168.547763834282, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol without hydrogen and a 0.4 A radii offset. */
  @Test
  public void testGaussVolButane() {
    // Configure input arguments for the Volume script.
    String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/butane.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(125.5120767378517, volume.totalVolume, tolerance);
    assertEquals(134.50524040566165, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol without hydrogen and a 0.0 A radii offset. */
  @Test
  public void testGaussVolCrambin() {
    // Configure input arguments for the Volume script.
    String[] args = {"-o", "0.0", "src/main/java/ffx/potential/structures/crambin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(4371.667466648112, volume.totalVolume, tolerance);
    assertEquals(3971.0619085859435, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol derivatives. */
  @Test
  public void testGaussVolCrambinDerivatives() {
    // Configure input arguments for the Gradient script.
    System.setProperty("gkterm", "true");
    System.setProperty("cavmodel", "gauss-disp");
    // Choose a random atom to test.
    int atomID = (int) floor(random() * 642) + 1;
    String[] args = {
        "--ga", Integer.toString(atomID), "src/main/java/ffx/potential/structures/crambin.xyz"
    };
    binding.setVariable("args", args);

    // Construct and evaluate the Gradient script
    Gradient gradient = new Gradient(binding).run();
    potentialScript = gradient;
    assertEquals("Crambin gradient failures: ", 0, gradient.nFailures);
  }

  /** Test GaussVol without hydrogen and a 0.4 A radii offset. */
  @Test
  public void testGaussVolEthylbenzene() {
    // Configure input arguments for the Volume script.
    String[] args = {"-o", "0.4", "src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Volume script.
    Volume volume = new Volume(binding).run();
    potentialScript = volume;
    assertEquals(194.44960348422916, volume.totalVolume, tolerance);
    assertEquals(193.55314214639066, volume.totalSurfaceArea, tolerance);
  }

  /** Test GaussVol derivatives. */
  @Test
  public void testGaussVolHydrogenEthylbenzeneDerivatives() {
    // Configure input arguments for the Gradient script.
    System.setProperty("gkterm", "true");
    System.setProperty("cavmodel", "gauss-disp");
    String[] args = {"src/main/java/ffx/potential/structures/ethylbenzene.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Gradient script
    Gradient gradient = new Gradient(binding).run();
    potentialScript = gradient;
    assertEquals("Ethylbenzene gradient failures: ", 0, gradient.nFailures);
  }
}
