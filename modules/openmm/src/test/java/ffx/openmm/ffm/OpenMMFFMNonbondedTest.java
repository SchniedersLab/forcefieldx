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
package ffx.openmm.ffm;

import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native parameter and context-update tests for nonbonded FFM forces. */
public class OpenMMFFMNonbondedTest {
  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testNonbondedParticleExceptionAndOptions() {
    try (NonbondedForce force = new NonbondedForce()) {
      assertEquals(0, force.addParticle(0.5, 0.3, 0.8));
      assertEquals(new NonbondedForce.ParticleParameters(0.5, 0.3, 0.8),
          force.getParticleParameters(0));
      force.setParticleParameters(0, -0.5, 0.31, 0.9);
      assertEquals(new NonbondedForce.ParticleParameters(-0.5, 0.31, 0.9),
          force.getParticleParameters(0));

      int exception = force.addException(0, 1, -0.25, 0.2, 0.4, false);
      assertEquals(0, exception);
      assertEquals(new NonbondedForce.ExceptionParameters(0, 1, -0.25, 0.2, 0.4),
          force.getExceptionParameters(exception));
      force.setExceptionParameters(exception, 1, 2, 0.1, 0.25, 0.5);
      assertEquals(new NonbondedForce.ExceptionParameters(1, 2, 0.1, 0.25, 0.5),
          force.getExceptionParameters(exception));
      assertEquals(1, force.getNumParticles());
      assertEquals(1, force.getNumExceptions());

      force.setNonbondedMethod(NonbondedForce.NonbondedMethod.CUTOFF_PERIODIC);
      force.setCutoffDistance(1.2);
      force.setUseSwitchingFunction(true);
      force.setSwitchingDistance(0.9);
      force.setUseDispersionCorrection(false);
      force.setReactionFieldDielectric(1.0);
      force.setEwaldErrorTolerance(1.0e-5);
      force.setPMEParameters(2.5, 32, 34, 36);
      assertEquals(NonbondedForce.NonbondedMethod.CUTOFF_PERIODIC, force.getNonbondedMethod());
      assertEquals(1.2, force.getCutoffDistance(), 0.0);
      assertTrue(force.getUseSwitchingFunction());
      assertEquals(0.9, force.getSwitchingDistance(), 0.0);
      assertFalse(force.getUseDispersionCorrection());
      assertEquals(1.0, force.getReactionFieldDielectric(), 0.0);
      assertEquals(1.0e-5, force.getEwaldErrorTolerance(), 0.0);
      assertEquals(new NonbondedForce.PMEParameters(2.5, 32, 34, 36), force.getPMEParameters());
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testGBSAOBCParameters() {
    try (GBSAOBCForce force = new GBSAOBCForce()) {
      assertEquals(0, force.addParticle(-0.4, 0.15, 0.8));
      assertEquals(new GBSAOBCForce.ParticleParameters(-0.4, 0.15, 0.8),
          force.getParticleParameters(0));
      force.setParticleParameters(0, 0.4, 0.16, 0.7);
      force.setSoluteDielectric(2.0);
      force.setSolventDielectric(78.5);
      force.setSurfaceAreaEnergy(2.25936);
      force.setNonbondedMethod(GBSAOBCForce.NonbondedMethod.CUTOFF_PERIODIC);
      force.setCutoffDistance(1.1);
      assertEquals(new GBSAOBCForce.ParticleParameters(0.4, 0.16, 0.7),
          force.getParticleParameters(0));
      assertEquals(2.0, force.getSoluteDielectric(), 0.0);
      assertEquals(78.5, force.getSolventDielectric(), 0.0);
      assertEquals(2.25936, force.getSurfaceAreaEnergy(), 0.0);
      assertEquals(GBSAOBCForce.NonbondedMethod.CUTOFF_PERIODIC, force.getNonbondedMethod());
      assertEquals(1.1, force.getCutoffDistance(), 0.0);
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testLennardJonesEnergyAndLiveParameterUpdate() {
    double sigma = 0.3;
    try (NonbondedForce force = new NonbondedForce();
         System system = new System();
         VerletIntegrator integrator = new VerletIntegrator(0.001);
         Platform platform = new Platform("Reference")) {
      force.setNonbondedMethod(NonbondedForce.NonbondedMethod.NO_CUTOFF);
      force.addParticle(0.0, sigma, 0.8);
      force.addParticle(0.0, sigma, 0.8);
      system.addParticle(12.0);
      system.addParticle(12.0);
      system.addForce(force);
      try (Context context = new Context(system, integrator, platform)) {
        double minimumDistance = Math.pow(2.0, 1.0 / 6.0) * sigma;
        context.setPositions(new double[]{0.0, 0.0, 0.0, minimumDistance, 0.0, 0.0});
        assertEquals(-0.8, potentialEnergy(context), 1.0e-10);
        force.setParticleParameters(1, 0.0, sigma, 1.8);
        force.updateParametersInContext(context);
        assertEquals(-1.2, potentialEnergy(context), 1.0e-10);
      }
    }
  }

  @Test
  public void testGBSAParameterUpdateInContext() {
    try (NonbondedForce nonbonded = new NonbondedForce();
         GBSAOBCForce gbsa = new GBSAOBCForce();
         System system = new System();
         VerletIntegrator integrator = new VerletIntegrator(0.001);
         Platform platform = new Platform("Reference")) {
      for (int particle = 0; particle < 2; particle++) {
        nonbonded.addParticle(0.0, 0.3, 0.0);
        gbsa.addParticle(0.0, 0.15, 0.8);
        system.addParticle(12.0);
      }
      system.addForce(nonbonded);
      system.addForce(gbsa);
      try (Context context = new Context(system, integrator, platform)) {
        context.setPositions(new double[]{0.0, 0.0, 0.0, 0.3, 0.0, 0.0});
        gbsa.setParticleParameters(0, 0.0, 0.16, 0.75);
        gbsa.updateParametersInContext(context);
        assertTrue(Double.isFinite(potentialEnergy(context)));
      }
    }
  }

  private static double potentialEnergy(Context context) {
    try (State state = context.getState(
        ffx.openmm.ffm.bindings.OpenMMNative.OpenMM_State_Energy(), false)) {
      return state.getPotentialEnergy();
    }
  }
}
