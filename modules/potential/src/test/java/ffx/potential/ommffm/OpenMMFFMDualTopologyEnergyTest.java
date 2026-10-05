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

import ffx.numerics.switching.PowerSwitch;
import ffx.numerics.switching.UnivariateSwitchingFunction;
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
import static org.junit.Assume.assumeTrue;

/**
 * Validation tests for OpenMM FFM Dual Topology energy evaluation across lambda paths.
 */
public class OpenMMFFMDualTopologyEnergyTest extends PotentialTest {

  private static final double TOLERANCE = 1.0e-4;

  @BeforeClass
  public static void setUp() {
    assumeTrue(System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testDualTopologyPeptide() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly topology1 = potentialsUtils.open(structure);
    MolecularAssembly topology2 = potentialsUtils.open(structure);

    UnivariateSwitchingFunction switchFunction = new PowerSwitch();
    OpenMMDualTopologyEnergy dualEnergy = new OpenMMDualTopologyEnergy(topology1, topology2, switchFunction, Platform.OMM_REF);

    int nVar = dualEnergy.getNumberOfVariables();
    double[] x = new double[nVar];
    dualEnergy.getCoordinates(x);

    double[] lambdas = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double lambda : lambdas) {
      dualEnergy.setLambda(lambda);
      double eFFX = dualEnergy.energyFFX(x, true);
      double eOMM = dualEnergy.energy(x, true);

      String message = String.format("Dual topology energy mismatch at lambda=%.2f: FFX=%12.6f, OMM=%12.6f",
          lambda, eFFX, eOMM);
      assertEquals(message, eFFX, eOMM, TOLERANCE);
    }

    dualEnergy.destroy();
  }

  @Test
  public void testDualTopologyDisplacedCoordinates() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly topology1 = potentialsUtils.open(structure);
    MolecularAssembly topology2 = potentialsUtils.open(structure);

    UnivariateSwitchingFunction switchFunction = new PowerSwitch();
    OpenMMDualTopologyEnergy dualEnergy = new OpenMMDualTopologyEnergy(topology1, topology2, switchFunction, Platform.OMM_REF);

    int nVar = dualEnergy.getNumberOfVariables();
    double[] x = new double[nVar];
    dualEnergy.getCoordinates(x);

    // Displace coordinates
    for (int i = 0; i < nVar; i += 3) {
      x[i] += 0.02 * (i % 5 - 2);
    }

    double[] lambdas = {0.0, 0.5, 1.0};
    for (double lambda : lambdas) {
      dualEnergy.setLambda(lambda);
      double eFFX = dualEnergy.energyFFX(x, true);
      double eOMM = dualEnergy.energy(x, true);

      String message = String.format("Dual topology displaced energy mismatch at lambda=%.2f: FFX=%12.6f, OMM=%12.6f",
          lambda, eFFX, eOMM);
      assertEquals(message, eFFX, eOMM, TOLERANCE);
    }

    dualEnergy.destroy();
  }
}
