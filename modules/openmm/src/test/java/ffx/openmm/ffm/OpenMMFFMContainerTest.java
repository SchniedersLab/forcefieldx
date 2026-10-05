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

/**
 * Integration tests for FFM-backed OpenMM foundation containers.
 */
public class OpenMMFFMContainerTest {

  @BeforeClass
  public static void initializeRuntime() {
    assumeTrue(java.lang.System.getenv(OpenMMRuntime.LIBRARY_DIRECTORY_ENVIRONMENT) != null);
    OpenMMRuntime.initialize();
  }

  @Test
  public void testBondArray() {
    try (BondArray bonds = new BondArray(1)) {
      bonds.set(0, 2, 3);
      bonds.append(5, 7);

      assertEquals(2, bonds.getSize());
      assertEquals(new BondArray.Bond(2, 3), bonds.get(0));
      assertEquals(new BondArray.Bond(5, 7), bonds.get(1));

      bonds.resize(1);
      assertEquals(1, bonds.getSize());
    }
  }

  @Test
  public void testDoubleAndIntArrays() {
    try (DoubleArray doubles = new DoubleArray(1); IntArray ints = new IntArray(1)) {
      doubles.set(0, 1.25);
      doubles.append(2.5);
      ints.set(0, 4);
      ints.append(8);

      assertEquals(2, doubles.getSize());
      assertEquals(1.25, doubles.get(0), 0.0);
      assertEquals(2.5, doubles.get(1), 0.0);
      assertEquals(2, ints.getSize());
      assertEquals(4, ints.get(0));
      assertEquals(8, ints.get(1));
    }
  }

  @Test
  public void testIntSet() {
    try (IntSet values = new IntSet()) {
      values.insert(3);
      values.insert(3);
      values.insert(5);

      assertEquals(2, values.getSize());
    }
  }

  @Test
  public void testStringArray() {
    try (StringArray strings = new StringArray(1)) {
      strings.set(0, "alpha");
      strings.append("beta-\u03b2");

      assertEquals(2, strings.getSize());
      assertEquals("alpha", strings.get(0));
      assertEquals("beta-\u03b2", strings.get(1));
      assertEquals(null, strings.get(-1));
      assertEquals(null, strings.get(2));
    }
  }

  @Test
  public void testVec3Array() {
    try (Vec3Array vectors = new Vec3Array(1)) {
      vectors.set(0, new Vec3(1.0, 2.0, 3.0));
      vectors.append(new Vec3(4.0, 5.0, 6.0));

      assertEquals(new Vec3(1.0, 2.0, 3.0), vectors.get(0));
      assertEquals(new Vec3(4.0, 5.0, 6.0), vectors.get(1));
      assertArrayEquals(new double[]{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, vectors.getArray(), 0.0);
    }

    try (Vec3Array vectors = Vec3Array.toVec3Array(new double[]{7.0, 8.0, 9.0})) {
      assertEquals(new Vec3(7.0, 8.0, 9.0), vectors.get(0));
    }
  }

  @Test
  public void testBooleanConversion() {
    assertTrue(OpenMMBooleans.fromNative(1));
    assertFalse(OpenMMBooleans.fromNative(0));
    assertEquals(1, OpenMMBooleans.toNative(true));
    assertEquals(0, OpenMMBooleans.toNative(false));
  }
}
