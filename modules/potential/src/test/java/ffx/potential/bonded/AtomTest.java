//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
//******************************************************************************
package ffx.potential.bonded;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import ffx.potential.parameters.AtomType;
import org.jogamp.vecmath.Vector3d;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/** Unit tests for the Atom class. */
public class AtomTest {

  private final double d[] = {0.0, 0.0, 0.0};
  private Atom atom = null;
  private Atom atom2 = null;
  private Atom atom3 = null;
  private Bond b = null;

  @Test(timeout = 500)
  public void Atom_addTrajectoryCoords() {
    int position = 1;
    Vector3d v3d = new Vector3d(1.0, 2.0, 3.0);
    atom.addTrajectoryCoords(v3d, position);
    atom.setCurrentCycle(2);
    atom.getV3D(v3d);
    assertEquals(1.0, v3d.x, 0.0);
    assertEquals(2.0, v3d.y, 0.0);
    assertEquals(3.0, v3d.z, 0.0);
  }

  @Test(timeout = 500)
  public void Atom_constructor() {
    String s = "Carbon";
    Atom actualReturn = new Atom(s);
    assertEquals("return value", s, actualReturn.getName());
  }

  @Test(timeout = 500)
  public void Atom_constructor3() {
    int index = 0;
    String i = "N";
    AtomType mmdata = null;
    Atom actualReturn = new Atom(index, i, mmdata, d);
    assertEquals("return value", 0, actualReturn.getIndex());
  }

  @Test(timeout = 500)
  public void Atom_constructor4() {
    int xyznum = 0;
    String id = "O";
    Character altLoc = 'A';
    String r = "GLY";
    int n = 1;
    Character c = 'A';
    Atom actualReturn = new Atom(xyznum, id, altLoc, d, r, n, c, 0.0, 0.0, "A");
    assertEquals("return value", "GLY", actualReturn.getResidueName());
  }

  @Test(timeout = 500)
  public void Atom_destroy() {
    boolean expectedReturn = true;
    boolean actualReturn = atom.destroy();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  @Test(timeout = 500)
  public void Atom_equals() {
    boolean expectedReturn = true;
    boolean actualReturn = atom.equals(atom);
    assertEquals("return value", expectedReturn, actualReturn);
    expectedReturn = false;
    actualReturn = atom.equals(null);
    assertEquals("return value", expectedReturn, actualReturn);
    Atom atom4 = new Atom(atom.getName());
    expectedReturn = false;
    actualReturn = atom.equals(atom4);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  @Test(timeout = 500)
  public void Atom_getBond() {
    Bond expectedReturn = b;
    Bond actualReturn = atom.getBond(atom2);
    assertEquals("return value", expectedReturn, actualReturn);
    actualReturn = atom.getBond(null);
    assertNull(actualReturn);
    actualReturn = atom.getBond(atom3);
    assertNull(actualReturn);
  }

  @Test(timeout = 500)
  public void Atom_getNumBonds() {
    int expectedReturn = 1;
    int actualReturn = atom.getNumBonds();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  @Test(timeout = 500)
  public void Atom_isBonded() {
    boolean expectedReturn = true;
    boolean actualReturn = atom.isBonded(atom2);
    assertEquals("return value", expectedReturn, actualReturn);
    expectedReturn = false;
    actualReturn = atom.isBonded(atom3);
    assertEquals("return value", expectedReturn, actualReturn);
    expectedReturn = false;
    actualReturn = atom.isBonded(null);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  @Before
  public void setUp() {
    atom = new Atom(1, "C", 'A', d, "GLY", 1, 'A', 0.0, 0.0, "A");
    atom2 = new Atom("C");
    atom3 = new Atom("O");
    b = new Bond(atom, atom2);
  }

  @After
  public void tearDown() {
    assertTrue(atom.destroy());
    assertTrue(atom2.destroy());
    assertTrue(atom3.destroy());
    assertTrue(b.destroy());
    atom = null;
    atom2 = null;
    atom3 = null;
    b = null;
  }
}
