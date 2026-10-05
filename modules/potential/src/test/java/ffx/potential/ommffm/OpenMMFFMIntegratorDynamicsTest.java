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
import ffx.openmm.ffm.bindings.OpenMMNative;
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
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/**
 * Validates OpenMM Integrator execution and state dynamics using the FFM backend.
 */
public class OpenMMFFMIntegratorDynamicsTest extends PotentialTest {

  @BeforeClass
  public static void setUpClass() {
    assumeTrue(System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Before
  public void setUp() {
    System.setProperty("polarization", "direct");
    System.setProperty("precision", "double");
    System.setProperty("ffx.openmm.backend", "ffm");
  }

  @After
  public void tearDown() {
    System.clearProperty("polarization");
    System.clearProperty("precision");
    System.clearProperty("ffx.openmm.backend");
  }

  @Test
  public void testVerletIntegratorNVE() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);
    assertNotNull(molecularAssembly);

    OpenMMEnergy openMMEnergy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);
    OpenMMContext context = openMMEnergy.getContext();
    assertNotNull(context);

    // Initialize Verlet dynamics with 0.1 fs timestep.
    context.update("VERLET", 0.0001, 298.15, true);
    context.setVelocitiesToTemperature(298.15, 42);

    int energyMask = OpenMMNative.OpenMM_State_Energy()
        | OpenMMNative.OpenMM_State_Positions()
        | OpenMMNative.OpenMM_State_Velocities();

    OpenMMState state0 = context.getOpenMMState(energyMask);
    double eTotal0 = state0.getTotalEnergy();
    double ePot0 = state0.getPotentialEnergy();
    double eKin0 = state0.getKineticEnergy();
    assertTrue(Double.isFinite(eTotal0));
    assertTrue(eKin0 > 0.0);
    state0.destroy();

    // Integrate 100 steps (0.1 ps).
    context.integrate(100);

    OpenMMState state1 = context.getOpenMMState(energyMask);
    double eTotal1 = state1.getTotalEnergy();
    double time1 = state1.getTime();
    assertTrue(Double.isFinite(eTotal1));
    assertEquals(0.0, time1, 1e-1);

    // Check NVE total energy conservation within numerical tolerance.
    double deltaE = Math.abs(eTotal1 - eTotal0);
    assertTrue("Total energy should be conserved in NVE Verlet integration: deltaE=" + deltaE, deltaE < 0.2);

    state1.destroy();
    openMMEnergy.destroy();
  }

  @Test
  public void testCustomMTSIntegratorNVE() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);
    assertNotNull(molecularAssembly);

    OpenMMEnergy openMMEnergy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);
    OpenMMContext context = openMMEnergy.getContext();
    assertNotNull(context);

    // Initialize Custom MTS integrator with 0.5 fs outer timestep.
    context.update("MTS", 0.0005, 298.15, true);
    context.setVelocitiesToTemperature(298.15, 101);

    int energyMask = OpenMMNative.OpenMM_State_Energy()
        | OpenMMNative.OpenMM_State_Positions()
        | OpenMMNative.OpenMM_State_Velocities();

    OpenMMState state0 = context.getOpenMMState(energyMask);
    double eTotal0 = state0.getTotalEnergy();
    double eKin0 = state0.getKineticEnergy();
    assertTrue(Double.isFinite(eTotal0));
    assertTrue(eKin0 > 0.0);
    state0.destroy();

    // Integrate 50 steps (0.1 ps).
    context.integrate(50);

    OpenMMState state1 = context.getOpenMMState(energyMask);
    double eTotal1 = state1.getTotalEnergy();
    double time1 = state1.getTime();
    assertTrue(Double.isFinite(eTotal1));
    assertEquals(0.0, time1, 1e-1);

    double deltaE = Math.abs(eTotal1 - eTotal0);
    assertTrue("Total energy should be conserved in Custom MTS integration: deltaE=" + deltaE, deltaE < 0.25);

    state1.destroy();
    openMMEnergy.destroy();
  }

  @Test
  public void testLangevinIntegratorsNVT() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);
    assertNotNull(molecularAssembly);

    OpenMMEnergy openMMEnergy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);
    OpenMMContext context = openMMEnergy.getContext();
    assertNotNull(context);

    int energyMask = OpenMMNative.OpenMM_State_Energy()
        | OpenMMNative.OpenMM_State_Positions()
        | OpenMMNative.OpenMM_State_Velocities();

    // 1. Standard Langevin Integrator
    context.update("LANGEVIN", 0.001, 298.15, true);
    context.integrate(50);
    OpenMMState langevinState = context.getOpenMMState(energyMask);
    assertTrue(Double.isFinite(langevinState.getPotentialEnergy()));
    assertTrue(langevinState.getKineticEnergy() > 0.0);
    assertEquals(0.05, langevinState.getTime(), 1e-5);
    langevinState.destroy();

    // 2. Custom MTS Langevin Integrator
    context.update("LANGEVIN-MTS", 0.002, 298.15, true);
    context.integrate(50);
    OpenMMState mtsLangevinState = context.getOpenMMState(energyMask);
    assertTrue(Double.isFinite(mtsLangevinState.getPotentialEnergy()));
    assertTrue(mtsLangevinState.getKineticEnergy() > 0.0);
    assertEquals(0.1, mtsLangevinState.getTime(), 1e-5);
    mtsLangevinState.destroy();

    openMMEnergy.destroy();
  }

  @Test
  public void testLocalEnergyMinimizer() {
    ClassLoader classLoader = getClass().getClassLoader();
    URL url = classLoader.getResource("peptide.pdb");
    assumeTrue(url != null);
    File structure = new File(url.getPath());

    PotentialsUtils potentialsUtils = new PotentialsUtils();
    MolecularAssembly molecularAssembly = potentialsUtils.open(structure);
    assertNotNull(molecularAssembly);

    OpenMMEnergy openMMEnergy = new OpenMMEnergy(molecularAssembly, Platform.OMM_REF, 1);
    OpenMMContext context = openMMEnergy.getContext();
    assertNotNull(context);

    int nVar = openMMEnergy.getNumberOfVariables();
    double[] x = new double[nVar];
    openMMEnergy.getCoordinates(x);

    // Displace coordinates to move away from minimum
    for (int i = 0; i < nVar; i += 3) {
      x[i] += 0.05 * (i % 3 + 1);
    }
    context.setPositions(x);

    OpenMMState stateUnmin = context.getOpenMMState(OpenMMNative.OpenMM_State_Energy());
    double eInitial = stateUnmin.getPotentialEnergy();
    stateUnmin.destroy();

    // Optimize
    context.optimize(1.0, 50);

    OpenMMState stateMin = context.getOpenMMState(OpenMMNative.OpenMM_State_Energy());
    double eFinal = stateMin.getPotentialEnergy();
    stateMin.destroy();

    assertTrue("Minimization should reduce potential energy: initial=" + eInitial + ", final=" + eFinal, eFinal < eInitial);

    openMMEnergy.destroy();
  }
}
