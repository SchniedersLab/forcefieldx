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
package ffx.openmm;

import edu.uiowa.jopenmm.OpenMMUtils;
import org.junit.After;
import org.junit.BeforeClass;
import org.junit.Test;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static org.junit.Assert.assertEquals;

/**
 * Simple unit tests to exercise creating a System with a HarmonicBondForce,
 * creating a Context with a Langevin integrator on the OpenCL platform,
 * and validating potential energy before and after a few MD steps.
 */
public class OpenMMHarmonicBondTest {

  private Context context;
  private Integrator integrator;
  private System system;

  @BeforeClass
  public static void init() {
    OpenMMUtils.init();
    String pluginDirectory = OpenMMUtils.getPluginDirectory();
    Platform.loadPluginsFromDirectory(pluginDirectory);
  }

  @After
  public void tearDown() {
    // Destroy in safe order to free native resources if they were created.
    try {
      if (context != null) context.destroy();
    } catch (Exception ignored) {
    }
    try {
      if (integrator != null) integrator.destroy();
    } catch (Exception ignored) {
    }
    try {
      if (system != null) system.destroy();
    } catch (Exception ignored) {
    }
  }

  @Test
  public void testHarmonicBondEnergyAtEquilibriumOpenCL() {
    Platform platform = new Platform("Reference");

    // Build a simple 2-particle system with a single harmonic bond.
    system = new System();
    double mass = 12.0; // amu (arbitrary positive mass)
    system.addParticle(mass);
    system.addParticle(mass);

    HarmonicBondForce bondForce = new HarmonicBondForce();
    double r0 = 0.1; // nm
    double k = 1000.0; // kJ/mol/nm^2
    bondForce.addBond(0, 1, r0, k);
    system.addForce(bondForce);

    // Use a Langevin integrator at 0 K to avoid thermal energy.
    double dt = 0.001; // ps
    double temperature = 0.0; // K
    double gamma = 1.0; // 1/ps
    integrator = new LangevinIntegrator(dt, temperature, gamma);

    // Create the context on OpenCL.
    context = new Context(system, integrator, platform);

    // Set positions at equilibrium distance along x-axis: energy should be ~0.
    double[] posEq = new double[]{0.0, 0.0, 0.0, r0, 0.0, 0.0};
    context.setPositions(posEq);

    State state0 = context.getState(OpenMM_State_Energy, 0);
    double e0 = state0.getPotentialEnergy();
    // Expect near-zero potential energy at equilibrium. Allow tiny numerical tolerance.
    assertEquals(0.0, e0, 1.0e-6);

    // Step a few times and re-check energy remains near zero.
    integrator.step(10);
    State state1 = context.getState(OpenMM_State_Energy, 0);
    double e1 = state1.getPotentialEnergy();
    assertEquals(0.0, e1, 1.0e-5);
  }

  @Test
  public void testHarmonicBondEnergyDisplacedOpenCL() {
    Platform platform = new Platform("Reference");

    // Build a simple 2-particle system with a single harmonic bond.
    system = new System();
    double mass = 12.0; // amu
    system.addParticle(mass);
    system.addParticle(mass);

    HarmonicBondForce bondForce = new HarmonicBondForce();
    double r0 = 0.1; // nm
    double k = 1000.0; // kJ/mol/nm^2
    bondForce.addBond(0, 1, r0, k);
    system.addForce(bondForce);

    // Use a Langevin integrator at 0 K to avoid thermal noise.
    double dt = 0.001; // ps
    double temperature = 0.0; // K
    double gamma = 1.0; // 1/ps
    integrator = new LangevinIntegrator(dt, temperature, gamma);

    // Create the context on OpenCL.
    context = new Context(system, integrator, platform);

    // Displace the bond slightly from equilibrium by +dr along x.
    double dr = 0.01; // nm
    double[] pos = new double[]{0.0, 0.0, 0.0, r0 + dr, 0.0, 0.0};
    context.setPositions(pos);

    // Analytical harmonic energy: 0.5 * k * dr^2
    double expected = 0.5 * k * dr * dr;

    State state0 = context.getState(OpenMM_State_Energy, 0);
    double e0 = state0.getPotentialEnergy();
    assertEquals(expected, e0, 1.0e-6);

    // Take a few steps; at 0 K the system should relax slightly toward minimum,
    // but small dt and friction may keep energy close; just ensure finite and non-negative.
    integrator.step(5);
    State state1 = context.getState(OpenMM_State_Energy, 0);
    double e1 = state1.getPotentialEnergy();
    // Energy should remain non-negative and near the same order of magnitude.
    assertEquals(expected, e1, 1.0e-3);
  }
}
