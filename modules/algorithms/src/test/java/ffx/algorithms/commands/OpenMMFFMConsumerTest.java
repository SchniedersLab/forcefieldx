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
package ffx.algorithms.commands;

import ffx.algorithms.dynamics.MolecularDynamicsOpenMM;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.misc.AlgorithmsTest;
import ffx.algorithms.optimize.MinimizeOpenMM;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Platform;
import ffx.potential.utils.PotentialsUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.net.URL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/**
 * End-to-end integration tests for high-level consumer algorithms running over the OpenMM FFM backend.
 */
public class OpenMMFFMConsumerTest extends AlgorithmsTest {

  @BeforeClass
  public static void setUp() {
    assumeTrue(System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
    System.setProperty("ffx.openmm.backend", "ffm");
  }

  @AfterClass
  public static void tearDown() {
    System.clearProperty("ffx.openmm.backend");
  }

  @Test
  public void testMinimizeOpenMMFFM() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("trpcage.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly assembly = potentialsUtils.open(structure);
    assertNotNull(assembly);

    assembly.getForceField().addProperty("platform", "OMM_REF");
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(assembly);
    assertTrue("Energy should be FFM OpenMMEnergy", energy instanceof ffx.potential.ommffm.OpenMMEnergy);

    double eInitial = energy.energy();

    MinimizeOpenMM minimizer = new MinimizeOpenMM(assembly, energy);
    minimizer.minimize(1.0, 50);

    double eFinal = energy.energy();
    assertTrue("Final energy should be lower than initial energy after minimization", eFinal < eInitial);
  }

  @Test
  public void testMolecularDynamicsOpenMMFFMVerlet() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("trpcage.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly assembly = potentialsUtils.open(structure);
    assertNotNull(assembly);

    assembly.getForceField().addProperty("platform", "OMM_REF");
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(assembly);
    assertTrue("Energy should be FFM OpenMMEnergy", energy instanceof ffx.potential.ommffm.OpenMMEnergy);

    MolecularDynamicsOpenMM md = new MolecularDynamicsOpenMM(assembly, energy, null,
        ThermostatEnum.ADIABATIC, IntegratorEnum.VERLET);

    // Run 50 steps of 1 fs
    md.dynamic(50, 1.0, 10.0, 10.0, 298.15, true, null);

    double totalEnergy = md.getTotalEnergy();
    assertTrue("Total energy should be finite", !Double.isNaN(totalEnergy) && !Double.isInfinite(totalEnergy));
  }

  @Test
  public void testMolecularDynamicsOpenMMFFMLangevin() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("trpcage.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly assembly = potentialsUtils.open(structure);
    assertNotNull(assembly);

    assembly.getForceField().addProperty("platform", "OMM_REF");
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(assembly);
    assertTrue("Energy should be FFM OpenMMEnergy", energy instanceof ffx.potential.ommffm.OpenMMEnergy);

    MolecularDynamicsOpenMM md = new MolecularDynamicsOpenMM(assembly, energy, null,
        ThermostatEnum.BUSSI, IntegratorEnum.LANGEVIN);

    // Run 50 steps of 1 fs
    md.dynamic(50, 1.0, 10.0, 10.0, 298.15, true, null);

    double temp = md.getTemperature();
    assertTrue("Temperature should be positive and finite", temp > 0.0 && !Double.isNaN(temp) && !Double.isInfinite(temp));
  }
}
