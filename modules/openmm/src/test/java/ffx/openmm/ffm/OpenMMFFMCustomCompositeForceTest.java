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
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeTrue;

/** Native tests for nested and composite custom forces. */
public class OpenMMFFMCustomCompositeForceTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testCustomCentroidBondParameters() {
    try (CustomCentroidBondForce force =
             new CustomCentroidBondForce(2, "k*distance(g1,g2)^2")) {
      force.addPerBondParameter("k");
      force.addGlobalParameter("scale", 2.0);
      assertEquals(0, force.addGroup(new int[]{0, 1}, new double[]{1.0, 3.0}));
      assertEquals(1, force.addGroup(new int[]{2}, new double[0]));
      assertEquals(0, force.addBond(new int[]{0, 1}, new double[]{5.0}));
      assertEquals(2, force.getNumGroups());
      assertEquals(2, force.getNumGroupsPerBond());
      assertEquals(1, force.getNumBonds());
      assertEquals("k", force.getPerBondParameterName(0));
      assertEquals("scale", force.getGlobalParameterName(0));
      assertEquals("k*distance(g1,g2)^2", force.getEnergyFunction());

      CustomCentroidBondForce.GroupParameters group = force.getGroupParameters(0);
      assertArrayEquals(new int[]{0, 1}, group.particles());
      assertArrayEquals(new double[]{1.0, 3.0}, group.weights(), 0.0);
      CustomCentroidBondForce.BondParameters bond = force.getBondParameters(0);
      assertArrayEquals(new int[]{0, 1}, bond.groups());
      assertArrayEquals(new double[]{5.0}, bond.parameters(), 0.0);

      try (IntArray particles = toIntArray(1, 2);
           DoubleArray weights = toDoubleArray(2.0, 4.0);
           IntArray groups = toIntArray(1, 0);
           DoubleArray parameters = toDoubleArray(7.0)) {
        force.setGroupParameters(0, particles, weights);
        force.setBondParameters(0, groups, parameters);
      }
      force.setEnergyFunction("scale*k*distance(g1,g2)^2");
      force.setGlobalParameterName(0, "scale2");
      force.setGlobalParameterDefaultValue(0, 4.0);
      force.setPerBondParameterName(0, "strength");
      force.setUsesPeriodicBoundaryConditions(true);
      assertArrayEquals(new int[]{1, 2}, force.getGroupParameters(0).particles());
      assertArrayEquals(new int[]{1, 0}, force.getBondParameters(0).groups());
      assertEquals("scale2", force.getGlobalParameterName(0));
      assertEquals("strength", force.getPerBondParameterName(0));
      assertEquals(4.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomCompoundBondParameters() {
    try (CustomCompoundBondForce force =
             new CustomCompoundBondForce(2, "k*(distance(p1,p2)-r0)^2")) {
      force.addGlobalParameter("k", 10.0);
      force.addPerBondParameter("r0");
      assertEquals(0, force.addBond(new int[]{3, 4}, new double[]{0.25}));
      assertEquals(1, force.getNumBonds());
      assertEquals(2, force.getNumParticlesPerBond());
      assertEquals(1, force.getNumGlobalParameters());
      assertEquals(1, force.getNumPerBondParameters());
      assertEquals("k", force.getGlobalParameterName(0));
      assertEquals("r0", force.getPerBondParameterName(0));
      assertEquals("k*(distance(p1,p2)-r0)^2", force.getEnergyFunction());

      CustomCompoundBondForce.BondParameters bond = force.getBondParameters(0);
      assertArrayEquals(new int[]{3, 4}, bond.particles());
      assertArrayEquals(new double[]{0.25}, bond.parameters(), 0.0);
      try (IntArray particles = toIntArray(1, 2);
           DoubleArray parameters = toDoubleArray(0.5)) {
        force.setBondParameters(0, particles, parameters);
      }
      force.setEnergyFunction("2*k*(distance(p1,p2)-r0)^2");
      force.setGlobalParameterName(0, "spring");
      force.setGlobalParameterDefaultValue(0, 20.0);
      force.setPerBondParameterName(0, "length");
      force.setUsesPeriodicBoundaryConditions(true);
      assertArrayEquals(new int[]{1, 2}, force.getBondParameters(0).particles());
      assertArrayEquals(new double[]{0.5}, force.getBondParameters(0).parameters(), 0.0);
      assertEquals("spring", force.getGlobalParameterName(0));
      assertEquals(20.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertEquals("length", force.getPerBondParameterName(0));
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomCVForceOwnershipAndParameters() {
    try (CustomCVForce force = new CustomCVForce("scale*volume^2");
         CustomVolumeForce variable = new CustomVolumeForce("volume")) {
      force.addGlobalParameter("scale", 1.5);
      assertEquals(0, force.addCollectiveVariable("volume", variable));
      assertTrue(variable.isDestroyed());
      assertEquals(1, force.getNumCollectiveVariables());
      assertEquals("volume", force.getCollectiveVariableName(0));
      assertTrue(force.getCollectiveVariable(0).address() != 0L);
      assertEquals("scale*volume^2", force.getEnergyFunction());
      assertEquals("scale", force.getGlobalParameterName(0));
      assertEquals(1.5, force.getGlobalParameterDefaultValue(0), 0.0);

      force.setEnergyFunction("2*scale*volume^2");
      force.setGlobalParameterName(0, "coefficient");
      force.setGlobalParameterDefaultValue(0, 2.5);
      assertEquals("2*scale*volume^2", force.getEnergyFunction());
      assertEquals("coefficient", force.getGlobalParameterName(0));
      assertEquals(2.5, force.getGlobalParameterDefaultValue(0), 0.0);
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  @Test
  public void testCustomVolumeForceParameters() {
    try (CustomVolumeForce force = new CustomVolumeForce("pressure*volume")) {
      force.addGlobalParameter("pressure", 1.0);
      assertEquals("pressure*volume", force.getEnergyFunction());
      assertEquals(1, force.getNumGlobalParameters());
      assertEquals("pressure", force.getGlobalParameterName(0));
      assertEquals(1.0, force.getGlobalParameterDefaultValue(0), 0.0);
      force.setEnergyFunction("2*pressure*volume");
      force.setGlobalParameterName(0, "p");
      force.setGlobalParameterDefaultValue(0, 2.0);
      assertEquals("2*pressure*volume", force.getEnergyFunction());
      assertEquals("p", force.getGlobalParameterName(0));
      assertEquals(2.0, force.getGlobalParameterDefaultValue(0), 0.0);
      assertTrue(force.usesPeriodicBoundaryConditions());
    }
  }

  private static IntArray toIntArray(int... values) {
    IntArray array = new IntArray(values.length);
    for (int i = 0; i < values.length; i++) {
      array.set(i, values[i]);
    }
    return array;
  }

  private static DoubleArray toDoubleArray(double... values) {
    DoubleArray array = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) {
      array.set(i, values[i]);
    }
    return array;
  }
}
