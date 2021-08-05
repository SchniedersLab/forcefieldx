/**
 * Title: Force Field X.
 *
 * <p>Description: Force Field X - Software for Molecular Biophysics.
 *
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * <p>This file is part of Force Field X.
 *
 * <p>Force Field X is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * <p>Force Field X is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * <p>You should have received a copy of the GNU General Public License along with Force Field X; if
 * not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 * <p>Linking this library statically or dynamically with other modules is making a combined work
 * based on this library. Thus, the terms and conditions of the GNU General Public License cover the
 * whole combination.
 *
 * <p>As a special exception, the copyright holders of this library give you permission to link this
 * library with independent modules to produce an executable, regardless of the license terms of
 * these independent modules, and to copy and distribute the resulting executable under terms of
 * your choice, provided that you also meet, for each linked independent module, the terms and
 * conditions of the license of that module. An independent module is a module which is not derived
 * from or based on this library. If you modify this library, you may extend this exception to your
 * version of the library, but you are not obligated to do so. If you do not wish to do so, delete
 * this exception statement from your version.
 */
package ffx.potential.bonded;

import static org.junit.Assert.assertArrayEquals;

import ffx.potential.MolecularAssembly;
import ffx.potential.utils.Loop;
import ffx.potential.utils.PotentialsUtils;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import org.junit.After;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/** @author Mallory R. Tollefson */
@RunWith(Parameterized.class)
public class LoopClosureTest {

  private final MolecularAssembly molecularAssembly;
  private final Loop loop;

  private double[][] xyzNTest;
  private double[][] xyzCTest;
  private double[][] xyzATest;

  public LoopClosureTest(
      double[][] xyzNTest, double[][] xyzATest, double[][] xyzCTest, double[][] xyzOTest) {
    int startResidue = 2;
    int endResidue = 4;
    ClassLoader classLoader = getClass().getClassLoader();
    File structure =
        new File(classLoader.getResource("ffx/potential/structures/LoopClosureTest.pdb").getPath());
    PotentialsUtils potentialsUtils = new PotentialsUtils();
    molecularAssembly = potentialsUtils.open(structure);
    loop = new Loop(molecularAssembly, startResidue, endResidue);

    this.xyzNTest = xyzNTest;
    this.xyzATest = xyzATest;
    this.xyzCTest = xyzCTest;
  }

  @Parameters
  public static Collection<Object[]> data() {

    double[][] xyzNTest = new double[3][3];
    double[][] xyzATest = new double[3][3];
    double[][] xyzCTest = new double[3][3];
    double[][] xyzOTest = new double[3][3];

    // Residue 1
    xyzNTest[0][0] = 7.773;
    xyzNTest[0][1] = -9.71;
    xyzNTest[0][2] = -7.32;
    xyzATest[0][0] = 6.331;
    xyzATest[0][1] = -9.839;
    xyzATest[0][2] = -7.259;
    xyzCTest[0][0] = 5.886372894231285;
    xyzCTest[0][1] = -10.55641925089512;
    xyzCTest[0][2] = -5.994873283542817;
    xyzOTest[0][0] = 4.7066623518635335;
    xyzOTest[0][1] = -10.772063009151791;
    xyzOTest[0][2] = -5.7426213147669065;
    // Residue 2
    xyzNTest[1][0] = 6.267265566616004;
    xyzNTest[1][1] = -11.821304411459156;
    xyzNTest[1][2] = -5.840321341761048;
    xyzATest[1][0] = 5.873570174757412;
    xyzATest[1][1] = -12.55730694668949;
    xyzATest[1][2] = -4.654655197113309;
    xyzCTest[1][0] = 4.522327673119161;
    xyzCTest[1][1] = -13.229312851952344;
    xyzCTest[1][2] = -4.836181407502477;
    xyzOTest[1][0] = 4.000608042124467;
    xyzOTest[1][1] = -13.903386108027433;
    xyzOTest[1][2] = -3.955679208709603;
    // Residue 3
    xyzNTest[2][0] = 3.9321584783416297;
    xyzNTest[2][1] = -13.724556941299172;
    xyzNTest[2][2] = -3.7520533645561343;
    xyzATest[2][0] = 2.642;
    xyzATest[2][1] = -14.377;
    xyzATest[2][2] = -3.863;
    xyzCTest[2][0] = 1.658;
    xyzCTest[2][1] = -13.856;
    xyzCTest[2][2] = -2.821;
    xyzOTest[2][0] = 0.5084362754396345;
    xyzOTest[2][1] = -14.272368583539699;
    xyzOTest[2][2] = -2.7373896189696216;

    return Arrays.asList(
        new Object[][] {
          {xyzNTest, xyzATest, xyzCTest, xyzOTest}, // constructor arguments for test set 1
        });
  }

  @After
  public void after() {
    molecularAssembly.getPotentialEnergy().destroy();
    System.gc();
  }

  @Test
  public void loopTest() {
    double[][] rA;
    double[][] rC;
    double[][] rN;

    rA = loop.getRA();
    rC = loop.getRC();
    rN = loop.getRN();

    for (int j = 0; j < 3; j++) {
      assertArrayEquals(rA[j], xyzATest[j], 1e-8);
    }

    for (int j = 0; j < 3; j++) {
      assertArrayEquals(rC[j], xyzCTest[j], 1e-8);
    }

    for (int j = 0; j < 3; j++) {
      assertArrayEquals(rN[j], xyzNTest[j], 1e-8);
    }
  }
}
