// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.ommffm;

import ffx.openmm.ffm.OpenMMRuntime;
import ffx.potential.MolecularAssembly;
import ffx.potential.Platform;
import ffx.potential.utils.PotentialTest;
import ffx.potential.utils.PotentialsUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.net.URL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/**
 * Targeted tests for OpenMM FFM bond energy and gradient on the Reference platform.
 */
public class OpenMMFFMBondEnergyTest extends PotentialTest {

  private static final double TOLERANCE = 1.0e-8;

  @BeforeClass
  public static void setUp() {
    assumeTrue(System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Before
  public void isolateBondTerm() {
    System.setProperty("bondterm", "true");
    System.setProperty("angleterm", "false");
    System.setProperty("strbndterm", "false");
    System.setProperty("ureyterm", "false");
    System.setProperty("opbendterm", "false");
    System.setProperty("improperterm", "false");
    System.setProperty("torsionterm", "false");
    System.setProperty("pitorsterm", "false");
    System.setProperty("imptorterm", "false");
    System.setProperty("tortorterm", "false");
    System.setProperty("strtorterm", "false");
    System.setProperty("angtorterm", "false");
    System.setProperty("pitorterm", "false");
    System.setProperty("vdwterm", "false");
    System.setProperty("mpoleterm", "false");
    System.setProperty("polarizeterm", "false");
    System.setProperty("gkterm", "false");
  }

  @After
  public void clearProperties() {
    System.clearProperty("bondterm");
    System.clearProperty("angleterm");
    System.clearProperty("strbndterm");
    System.clearProperty("ureyterm");
    System.clearProperty("opbendterm");
    System.clearProperty("torsionterm");
    System.clearProperty("pitorsterm");
    System.clearProperty("improperterm");
    System.clearProperty("imptorterm");
    System.clearProperty("tortorterm");
    System.clearProperty("strtorterm");
    System.clearProperty("angtorterm");
    System.clearProperty("pitorterm");
    System.clearProperty("vdwterm");
    System.clearProperty("mpoleterm");
    System.clearProperty("polarizeterm");
    System.clearProperty("gkterm");
  }

  @Test
  public void testOrganicMoleculeBondEnergyAndGradient() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("ethylbenzene.xyz");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);

    OpenMMEnergy openMMEnergy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);

    int nVar = openMMEnergy.getNumberOfVariables();
    double[] x = new double[nVar];
    openMMEnergy.getCoordinates(x);

    double eFFX = openMMEnergy.energyFFX(x, true);
    double[] gFFX = new double[nVar];
    openMMEnergy.energyAndGradientFFX(x, gFFX);

    double[] gOMM = new double[nVar];
    double eOMM = openMMEnergy.energyAndGradient(x, gOMM);

    assertEquals(eFFX, eOMM, TOLERANCE);
    for (int i = 0; i < nVar; i++) {
      assertEquals(gFFX[i], gOMM[i], TOLERANCE);
    }

    // Displace several coordinates
    for (int i = 0; i < nVar; i += 3) {
      x[i] += 0.02 * (i % 5 - 2);
    }

    double eFFXDisplaced = openMMEnergy.energyFFX(x);
    openMMEnergy.energyAndGradientFFX(x, gFFX);
    double eOMMDisplaced = openMMEnergy.energyAndGradient(x, gOMM);

    assertTrue(eOMMDisplaced > 0.0);
    assertEquals(eFFXDisplaced, eOMMDisplaced, TOLERANCE);
    for (int i = 0; i < nVar; i++) {
      assertEquals(gFFX[i], gOMM[i], TOLERANCE);
    }

    openMMEnergy.destroy();
  }

  @Test
  public void testLifecycleLeakFree() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("ethylbenzene.xyz");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);

    for (int iter = 0; iter < 100; iter++) {
      OpenMMEnergy energy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);
      int nVar = energy.getNumberOfVariables();
      double[] x = new double[nVar];
      energy.getCoordinates(x);
      double e = energy.energy(x);
      assertTrue(Double.isFinite(e));
      energy.destroy();
    }
  }
}
