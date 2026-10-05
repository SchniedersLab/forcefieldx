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
import ffx.openmm.ffm.IntArray;
import ffx.openmm.ffm.OpenMMRuntime;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native tests for the basic AMOEBA FFM façades. */
public class OpenMMFFMAmoebaBasicTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testDoubleArray3D() {
    try (DoubleArray3D array = new DoubleArray3D(1, 1, 2);
         DoubleArray values = new DoubleArray(0)) {
      values.append(1.25);
      values.append(-2.5);
      array.set(0, 0, values);
    }
  }

  @Test
  public void testVdwParticleParametersAndControls() {
    try (VdwForce force = new VdwForce();
         IntArray exclusions = new IntArray(0)) {
      assertEquals(0, force.addParticle(-1, 0.31, 0.42, 0.8, true, 0.75));
      assertEquals(1, force.addParticle(-1, 0.29, 0.33, 1.0, false, 1.0));
      assertEquals(2, force.getNumParticles());
      assertEquals(new VdwForce.ParticleParameters(-1, 0.31, 0.42, 0.8, true, -1, 0.75),
          force.getParticleParameters(0));

      force.setParticleParameters(0, 1, 0.32, 0.5, 0.7, false, -1, 0.9);
      assertEquals(new VdwForce.ParticleParameters(1, 0.32, 0.5, 0.7, false, -1, 0.9),
          force.getParticleParameters(0));

      force.setLambdaName("testVdwLambda");
      force.setPotentialFunction(1);
      force.setAlchemicalMethod(2);
      force.setNonbondedMethod(1);
      force.setCutoffDistance(1.4);
      force.setSoftcoreAlpha(0.6);
      force.setSoftcorePower(4);
      force.setUseDispersionCorrection(false);
      force.setSigmaCombiningRule("GEOMETRIC");
      force.setEpsilonCombiningRule("HARMONIC");
      assertEquals("testVdwLambda", force.getLambda());
      assertEquals(1, force.getPotentialFunction());
      assertEquals(2, force.getAlchemicalMethod());
      assertEquals(1, force.getNonbondedMethod());
      assertEquals(1.4, force.getCutoffDistance(), 0.0);
      assertEquals(0.6, force.getSoftcoreAlpha(), 0.0);
      assertEquals(4, force.getSoftcorePower());
      assertFalse(force.getUseDispersionCorrection());
      assertEquals("GEOMETRIC", force.getSigmaCombiningRule());
      assertEquals("HARMONIC", force.getEpsilonCombiningRule());
      assertTrue(force.usesPeriodicBoundaryConditions());

      exclusions.append(1);
      force.setParticleExclusions(0, exclusions);
      try (IntArray actual = force.getParticleExclusions(0)) {
        assertEquals(1, actual.getSize());
        assertEquals(1, actual.get(0));
      }
    }
  }

  @Test
  public void testVdwParticleTypesAndTypePair() {
    try (VdwForce force = new VdwForce()) {
      int firstType = force.addParticleType(0.3, 0.4);
      int secondType = force.addParticleType(0.5, 0.6);
      assertEquals(0, firstType);
      assertEquals(1, secondType);
      assertEquals(new VdwForce.TypeParameters(0.5, 0.6),
          force.getParticleTypeParameters(secondType));

      assertEquals(0, force.addParticle(-1, firstType, 1.0, false, 1.0));
      assertEquals(1, force.addParticle(-1, secondType, 1.0, true, 0.8));
      assertTrue(force.getUseParticleTypes());
      assertEquals(1, force.getParticleParameters(1).typeIndex());

      int pairIndex = force.addTypePair(firstType, secondType, 0.37, 0.45);
      assertEquals(0, pairIndex);
      assertEquals(new VdwForce.TypePairParameters(firstType, secondType, 0.37, 0.45),
          force.getTypePairParameters(pairIndex));
      force.setTypePairParameters(pairIndex, secondType, firstType, 0.38, 0.46);
      assertEquals(new VdwForce.TypePairParameters(secondType, firstType, 0.38, 0.46),
          force.getTypePairParameters(pairIndex));
    }
  }

  @Test
  public void testWcaDispersionParameters() {
    try (WcaDispersionForce force = new WcaDispersionForce()) {
      assertEquals(0, force.addParticle(0.14, 0.8));
      assertEquals(1, force.getNumParticles());
      assertEquals(new WcaDispersionForce.ParticleParameters(0.14, 0.8),
          force.getParticleParameters(0));
      force.setParticleParameters(0, 0.16, 0.9);
      assertEquals(new WcaDispersionForce.ParticleParameters(0.16, 0.9),
          force.getParticleParameters(0));

      force.setEpso(0.17);
      force.setEpsh(0.02);
      force.setRmino(0.17);
      force.setRminh(0.12);
      force.setAwater(0.0334);
      force.setShctd(0.82);
      force.setDispoff(0.18);
      force.setSlevy(1.2);
      assertEquals(0.17, force.getEpso(), 0.0);
      assertEquals(0.02, force.getEpsh(), 0.0);
      assertEquals(0.17, force.getRmino(), 0.0);
      assertEquals(0.12, force.getRminh(), 0.0);
      assertEquals(0.0334, force.getAwater(), 0.0);
      assertEquals(0.82, force.getShctd(), 0.0);
      assertEquals(0.18, force.getDispoff(), 0.0);
      assertEquals(1.2, force.getSlevy(), 0.0);
      assertFalse(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test(expected = UnsupportedOperationException.class)
  public void testGkCavitationReportsMissingNativeApi() {
    new GKCavitationForce();
  }
}
