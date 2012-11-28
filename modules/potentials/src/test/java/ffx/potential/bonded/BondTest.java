/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential.bonded;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Unit tests for the Bond class.
 */
public class BondTest {

    private Bond bond = null;
    private Bond bond2 = null;
    private Bond bond3 = null;
    private final double[] a1d = {0.0, 0.0, 0.0};
    private final double[] a2d = {1.0, 0.0, 0.0};
    private final double[] a3d = {1.0, 1.0, 0.0};
    private final double[] a4d = {0.0, 1.0, 1.0};
    private Atom a1;
    private Atom a2;
    private Atom a3;
    private Atom a4;
    private Angle an1;

    @Test(timeout = 500)
    public void Bond_constructor() {
        bond = new Bond("Empty Bond");
        assertNotNull(bond);
    }

    @Test(timeout = 500)
    public void Bond_formsAngleWith() {
        boolean expectedReturn = true;
        boolean actualReturn = bond.formsAngleWith(bond2);
        assertEquals("return value", expectedReturn, actualReturn);
        expectedReturn = false;
        actualReturn = bond.formsAngleWith(bond3);
        assertEquals("return value", expectedReturn, actualReturn);
    }

    @Test(timeout = 500)
    public void Bond_getCommonAtom() {
        Atom expectedReturn = a2;
        Atom actualReturn = bond.getCommonAtom(bond2);
        assertEquals("return value", expectedReturn, actualReturn);
        actualReturn = bond.getCommonAtom(bond3);
        assertNull(actualReturn);
        actualReturn = bond.getCommonAtom(null);
        assertNull(actualReturn);
        actualReturn = bond.getCommonAtom(bond);
        assertNull(actualReturn);
    }

    @Test(timeout = 500)
    public void Bond_getOtherAtom() {
        Atom expectedReturn = a2;
        Atom actualReturn = bond.get1_2(a1);
        assertEquals("return value", expectedReturn, actualReturn);
        expectedReturn = a1;
        actualReturn = bond.get1_2(a2);
        assertEquals("return value", expectedReturn, actualReturn);
        actualReturn = bond.get1_2(a3);
        assertNull(actualReturn);
        Atom nullAtom = null;
        actualReturn = bond.get1_2(nullAtom);
        assertNull(actualReturn);
    }

    @Test(timeout = 500)
    public void Bond_getOtherAtom2() {
        Atom expectedReturn = a1;
        Atom actualReturn = bond.getOtherAtom(bond2);
        assertEquals("return value", expectedReturn, actualReturn);
        actualReturn = bond.getOtherAtom(bond3);
        assertNull(actualReturn);
        Bond nullBond = null;
        actualReturn = bond.getOtherAtom(nullBond);
        assertNull(actualReturn);
    }

    @Test(timeout = 500)
    public void Bond_update() {
        // Atom 2 position == Atom 1 position
        a2.moveTo(a1.getX(), a1.getY(), a1.getZ());
        bond.update();
        assertEquals(0.0, bond.getValue(), 0.0);
        // PI
        a1.moveTo(0.0, 0.0, 0.0);
        a2.moveTo(Math.PI, 0.0, 0.0);
        bond.update();
        assertEquals(Math.PI, bond.getValue(), 0.0);
    }

    @Before
    public void setUp() {
        a1 = new Atom("A1");
        a1.setXYZ(a1d);
        a2 = new Atom("A2");
        a2.setXYZ(a2d);
        a3 = new Atom("A3");
        a3.setXYZ(a3d);
        a4 = new Atom("A4");
        a4.setXYZ(a4d);
        assertNotNull(a1);
        assertNotNull(a2);
        assertNotNull(a3);
        assertNotNull(a4);
        bond = new Bond(a1, a2);
        bond2 = new Bond(a2, a3);
        bond3 = new Bond(a3, a4);
        assertNotNull(bond);
        assertNotNull(bond2);
        assertNotNull(bond3);
        an1 = new Angle(bond, bond2);
        assertNotNull(an1);
    }

    @After
    public void tearDown() {
        assertTrue(a1.destroy());
        assertTrue(a2.destroy());
        assertTrue(a3.destroy());
        assertTrue(a4.destroy());
        assertTrue(bond.destroy());
        assertTrue(bond2.destroy());
        assertTrue(bond3.destroy());
        assertTrue(an1.destroy());
        a1 = null;
        a2 = null;
        a3 = null;
        a4 = null;
        bond = null;
        bond2 = null;
        bond3 = null;
    }
}
