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

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native tests for multi-particle custom force APIs. */
public class OpenMMFFMCustomMultiParticleForceTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testCustomNonbondedParameters() {
    try (CustomNonbondedForce force = new CustomNonbondedForce("scale*a1*a2/r");
         IntSet group1 = new IntSet();
         IntSet group2 = new IntSet()) {
      assertEquals("scale*a1*a2/r", force.getEnergyFunction());
      force.addPerParticleParameter("a");
      force.addGlobalParameter("scale", 2.0);
      force.addEnergyParameterDerivative("scale");
      assertEquals(0, force.addParticle(new double[]{1.0}));
      force.addParticle(new double[]{2.0});
      force.addParticle(new double[]{3.0});
      force.addExclusion(0, 1);
      assertArrayEquals(new double[]{2.0}, force.getParticleParameters(1), 0.0);
      assertEquals(new CustomNonbondedForce.Exclusion(0, 1), force.getExclusionParticles(0));

      group1.insert(0);
      group1.insert(1);
      group2.insert(2);
      force.addInteractionGroup(group1, group2);
      assertEquals(1, force.getNumInteractionGroups());
      try (IntSet actual1 = new IntSet(); IntSet actual2 = new IntSet()) {
        force.getInteractionGroupParameters(0, actual1, actual2);
        assertEquals(2, actual1.getSize());
        assertEquals(1, actual2.getSize());
        force.setInteractionGroupParameters(0, group2, group1);
      }

      force.setParticleParameters(1, new double[]{3.0});
      assertArrayEquals(new double[]{3.0}, force.getParticleParameters(1), 0.0);
      force.setGlobalParameterName(0, "scale2");
      force.setGlobalParameterDefaultValue(0, 4.0);
      force.setPerParticleParameterName(0, "amplitude");
      assertEquals("scale2", force.getGlobalParameterName(0));
      assertEquals(4.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals("amplitude", force.getPerParticleParameterName(0));
      assertEquals("scale2", force.getEnergyParameterDerivativeName(0));

      force.setNonbondedMethod(1);
      force.setCutoffDistance(1.2);
      force.setUseSwitchingFunction(true);
      force.setSwitchingDistance(0.8);
      force.setUseLongRangeCorrection(true);
      assertEquals(1, force.getNonbondedMethod());
      assertEquals(1.2, force.getCutoffDistance(), 0.0);
      assertTrue(force.getUseSwitchingFunction());
      assertEquals(0.8, force.getSwitchingDistance(), 0.0);
      assertTrue(force.getUseLongRangeCorrection());
    }
  }

  @Test
  public void testCustomGBParameters() {
    try (CustomGBForce force = new CustomGBForce()) {
      force.addPerParticleParameter("charge");
      force.addGlobalParameter("scale", 1.0);
      force.addEnergyParameterDerivative("scale");
      force.addComputedValue("q2", "charge*charge", 0);
      force.addEnergyTerm("scale*q2", 0);
      force.addParticle(new double[]{-0.5});
      force.addParticle(new double[]{0.5});
      force.addExclusion(0, 1);
      assertEquals(2, force.getNumParticles());
      assertEquals(1, force.getNumComputedValues());
      assertEquals(1, force.getNumEnergyTerms());
      assertEquals(1, force.getNumExclusions());
      assertEquals("charge", force.getPerParticleParameterName(0));
      assertEquals("scale", force.getEnergyParameterDerivativeName(0));
      assertArrayEquals(new double[]{-0.5}, force.getParticleParameters(0), 0.0);
      assertEquals(new CustomGBForce.Exclusion(0, 1), force.getExclusionParticles(0));
      force.setComputedValueParameters(0, "qSquared", "charge*charge", 0);
      force.setEnergyTermParameters(0, "scale*qSquared", 0);
      force.setParticleParameters(0, new double[]{-0.7});
      assertArrayEquals(new double[]{-0.7}, force.getParticleParameters(0), 0.0);
      force.setGlobalParameterName(0, "scale2");
      force.setGlobalParameterDefaultValue(0, 2.0);
      assertEquals("scale2", force.getGlobalParameterName(0));
      assertEquals(2.0, force.getGlobalParameterDefaultValue(0), 0.0);
      force.setNonbondedMethod(1);
      force.setCutoffDistance(1.1);
      assertEquals(1, force.getNonbondedMethod());
      assertEquals(1.1, force.getCutoffDistance(), 0.0);
    }
  }

  @Test
  public void testCustomHbondParameters() {
    try (CustomHbondForce force = new CustomHbondForce("k*distance(a1,d1)^2")) {
      force.addPerDonorParameter("donorWeight");
      force.addPerAcceptorParameter("acceptorWeight");
      force.addGlobalParameter("k", 2.0);
      force.addDonor(0, 1, -1, new double[]{0.5});
      force.addAcceptor(2, 3, -1, new double[]{1.5});
      force.addExclusion(0, 0);
      assertEquals(1, force.getNumDonors());
      assertEquals(1, force.getNumAcceptors());
      assertEquals(1, force.getNumExclusions());
      CustomHbondForce.GroupParameters donor = force.getDonorParameters(0);
      assertEquals(0, donor.particle1());
      assertEquals(1, donor.particle2());
      assertEquals(-1, donor.particle3());
      assertArrayEquals(new double[]{0.5}, donor.parameters(), 0.0);
      CustomHbondForce.GroupParameters acceptor = force.getAcceptorParameters(0);
      assertEquals(2, acceptor.particle1());
      assertEquals(3, acceptor.particle2());
      assertArrayEquals(new double[]{1.5}, acceptor.parameters(), 0.0);
      assertEquals(new CustomHbondForce.Exclusion(0, 0), force.getExclusionParticles(0));
      force.setDonorParameters(0, 4, 5, 6, new double[]{0.75});
      force.setAcceptorParameters(0, 7, 8, 9, new double[]{1.75});
      assertEquals(4, force.getDonorParameters(0).particle1());
      assertEquals(7, force.getAcceptorParameters(0).particle1());
      force.setGlobalParameterName(0, "strength");
      force.setGlobalParameterDefaultValue(0, 3.0);
      force.setPerDonorParameterName(0, "dWeight");
      force.setPerAcceptorParameterName(0, "aWeight");
      force.setNonbondedMethod(1);
      force.setCutoffDistance(0.9);
      assertEquals("strength", force.getGlobalParameterName(0));
      assertEquals("dWeight", force.getPerDonorParameterName(0));
      assertEquals("aWeight", force.getPerAcceptorParameterName(0));
      assertEquals(3.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals(1, force.getNonbondedMethod());
      assertEquals(0.9, force.getCutoffDistance(), 0.0);
      assertFalse(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomManyParticleParametersAndOwnership() {
    try (CustomManyParticleForce force = new CustomManyParticleForce(3, "scale*a1*a2*a3")) {
      force.addPerParticleParameter("a");
      force.addGlobalParameter("scale", 1.0);
      force.addParticle(new double[]{1.0}, 2);
      force.addParticle(new double[]{2.0}, 3);
      force.addParticle(new double[]{3.0}, 4);
      assertEquals(3, force.getNumParticles());
      assertEquals(3, force.getNumParticlesPerSet());
      CustomManyParticleForce.ParticleParameters initial = force.getParticleParameters(1);
      assertEquals(3, initial.type());
      assertArrayEquals(new double[]{2.0}, initial.parameters(), 0.0);
      force.setParticleParameters(1, new double[]{4.0}, 8);
      assertEquals(8, force.getParticleParameters(1).type());
      assertArrayEquals(new double[]{4.0}, force.getParticleParameters(1).parameters(), 0.0);
      force.addExclusion(0, 2);
      assertEquals(new CustomManyParticleForce.Exclusion(0, 2), force.getExclusionParticles(0));
      force.setPermutationMode(1);
      force.setNonbondedMethod(2);
      force.setCutoffDistance(1.0);
      assertEquals(1, force.getPermutationMode());
      assertEquals(2, force.getNonbondedMethod());
      assertEquals(1.0, force.getCutoffDistance(), 0.0);

      try (IntSet types = new IntSet()) {
        types.insert(2);
        types.insert(4);
        force.setTypeFilter(0, types);
        try (IntSet actualTypes = new IntSet()) {
          force.getTypeFilter(0, actualTypes);
          assertEquals(2, actualTypes.getSize());
        }
      }
    }
  }
}
