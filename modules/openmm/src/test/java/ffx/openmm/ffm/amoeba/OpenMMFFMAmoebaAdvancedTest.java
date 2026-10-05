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
package ffx.openmm.ffm.amoeba;

import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.OpenMMRuntime;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assume.assumeTrue;

/** Native tests for AMOEBA solvent, multipole, and HIPPO FFM façades. */
public class OpenMMFFMAmoebaAdvancedTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testGeneralizedKirkwoodParameters() {
    try (GeneralizedKirkwoodForce force = new GeneralizedKirkwoodForce()) {
      assertEquals(0, force.addParticle(-0.2, 0.15, 0.8));
      assertEquals(1, force.addParticle(0.2, 0.16, 0.9, 0.17, 0.12));
      assertEquals(2, force.getNumParticles());
      assertEquals(new GeneralizedKirkwoodForce.ParticleParameters(0.2, 0.16, 0.9, 0.17, 0.12),
          force.getParticleParameters(1));

      force.setParticleParameters(1, new GeneralizedKirkwoodForce.ParticleParameters(
          0.3, 0.18, 0.7, 0.19, 0.11));
      assertEquals(new GeneralizedKirkwoodForce.ParticleParameters(0.3, 0.18, 0.7, 0.19, 0.11),
          force.getParticleParameters(1));

      force.setSoluteDielectric(1.2);
      force.setSolventDielectric(78.3);
      force.setProbeRadius(0.14);
      force.setDescreenOffset(0.01);
      force.setDielectricOffset(0.02);
      force.setSurfaceAreaFactor(0.03);
      force.setIncludeCavityTerm(1);
      force.setTanhParameters(0.8, 0.9, 1.1);
      assertEquals(1.2, force.getSoluteDielectric(), 0.0);
      assertEquals(78.3, force.getSolventDielectric(), 0.0);
      assertEquals(0.14, force.getProbeRadius(), 0.0);
      assertEquals(0.01, force.getDescreenOffset(), 0.0);
      assertEquals(0.02, force.getDielectricOffset(), 0.0);
      assertEquals(0.03, force.getSurfaceAreaFactor(), 0.0);
      assertEquals(1, force.getIncludeCavityTerm());
      assertEquals(new GeneralizedKirkwoodForce.TanhParameters(0.8, 0.9, 1.1),
          force.getTanhParameters());
      assertFalse(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testMultipoleParameters() {
    try (MultipoleForce force = new MultipoleForce();
         DoubleArray dipole = values(0.1, 0.2, 0.3);
         DoubleArray quadrupole = values(0.0, 0.1, 0.2, 0.1, 0.0, 0.3, 0.2, 0.3, 0.0)) {
      assertEquals(0, force.addMultipole(0.25, dipole, quadrupole, 5, -1, -1, -1, 0.4, 0.5, 0.6));
      assertEquals(1, force.getNumMultipoles());
      MultipoleForce.MultipoleParameters actual = force.getMultipoleParameters(0);
      assertEquals(0.25, actual.charge(), 0.0);
      assertArrayEquals(new double[]{0.1, 0.2, 0.3}, actual.molecularDipole(), 0.0);
      assertArrayEquals(new double[]{0.0, 0.1, 0.2, 0.1, 0.0, 0.3, 0.2, 0.3, 0.0},
          actual.molecularQuadrupole(), 0.0);
      assertEquals(5, actual.axisType());
      assertEquals(-1, actual.multipoleAtomZ());
      assertEquals(-1, actual.multipoleAtomX());
      assertEquals(-1, actual.multipoleAtomY());
      assertEquals(0.4, actual.thole(), 0.0);
      assertEquals(0.5, actual.dampingFactor(), 0.0);
      assertEquals(0.6, actual.polarity(), 0.0);
      force.setPMEParameters(0.0, 32, 32, 32);
      assertEquals(new MultipoleForce.PmeParameters(0.0, 32, 32, 32), force.getPMEParameters());
    }
  }

  @Test
  public void testHippoParameters() {
    try (HippoNonbondedForce force = new HippoNonbondedForce()) {
      HippoNonbondedForce.ParticleParameters particle = new HippoNonbondedForce.ParticleParameters(
          0.1, new double[]{0.1, 0.0, 0.0}, new double[9], 0.2, 0.3, 0.4, 0.5, 0.6,
          0.7, 0.8, 0.9, 1.0, 5, -1, -1, -1);
      assertEquals(0, force.addParticle(particle));
      assertEquals(1, force.getNumParticles());
      HippoNonbondedForce.ParticleParameters actual = force.getParticleParameters(0);
      assertEquals(particle.charge(), actual.charge(), 0.0);
      assertArrayEquals(particle.dipole(), actual.dipole(), 0.0);
      assertArrayEquals(particle.quadrupole(), actual.quadrupole(), 0.0);
      assertEquals(particle.coreCharge(), actual.coreCharge(), 0.0);
      assertEquals(particle.alpha(), actual.alpha(), 0.0);
      assertEquals(particle.epsilon(), actual.epsilon(), 0.0);
      assertEquals(particle.damping(), actual.damping(), 0.0);
      assertEquals(particle.c6(), actual.c6(), 0.0);
      assertEquals(particle.pauliK(), actual.pauliK(), 0.0);
      assertEquals(particle.pauliQ(), actual.pauliQ(), 0.0);
      assertEquals(particle.pauliAlpha(), actual.pauliAlpha(), 0.0);
      assertEquals(particle.polarizability(), actual.polarizability(), 0.0);
      assertEquals(particle.axisType(), actual.axisType());
      assertEquals(particle.multipoleAtomZ(), actual.multipoleAtomZ());
      assertEquals(particle.multipoleAtomX(), actual.multipoleAtomX());
      assertEquals(particle.multipoleAtomY(), actual.multipoleAtomY());

      HippoNonbondedForce.ExceptionParameters exception = new HippoNonbondedForce.ExceptionParameters(
          0, 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
      assertEquals(0, force.addException(0, 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, true));
      assertEquals(1, force.getNumExceptions());
      assertEquals(exception, force.getExceptionParameters(0));
      force.setExceptionParameters(0, new HippoNonbondedForce.ExceptionParameters(
          0, 1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7));
      assertEquals(0.7, force.getExceptionParameters(0).chargeTransferScale(), 0.0);

      force.setPMEParameters(0.0, 32, 32, 32);
      force.setDPMEParameters(0.0, 16, 16, 16);
      assertEquals(new HippoNonbondedForce.PmeParameters(0.0, 32, 32, 32),
          force.getPMEParameters());
      assertEquals(new HippoNonbondedForce.PmeParameters(0.0, 16, 16, 16),
          force.getDPMEParameters());
    }
  }

  private static DoubleArray values(double... values) {
    DoubleArray result = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) {
      result.set(i, values[i]);
    }
    return result;
  }
}
