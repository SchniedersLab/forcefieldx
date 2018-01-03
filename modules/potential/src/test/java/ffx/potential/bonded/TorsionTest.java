/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.bonded;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Unit tests for the Torsion class.
 */
public class TorsionTest {

    private final double[] a1d = {0.0, 0.0, 0.0};
    private final double[] a2d = {1.0, 0.0, 0.0};
    private final double[] a3d = {1.0, 1.0, 0.0};
    private final double[] a4d = {1.0, 1.0, 1.0};
    private Atom a1;
    private Atom a2;
    private Atom a3;
    private Atom a4;
    private Bond b1;
    private Bond b2;
    private Bond b3;
    private Angle an1;
    private Angle an2;
    private Torsion torsion = null;

    @Test(timeout = 500)
    public void Torsion_constructor() {
        torsion = new Torsion(an1, b3);
        assertNotNull(torsion);
    }

    @Test(timeout = 500)
    public void Torsion_constructor2() {
        String n = "Empty Dihedral";
        torsion = new Torsion(n);
        assertNotNull(torsion);
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
        b1 = new Bond(a1, a2);
        b2 = new Bond(a2, a3);
        b3 = new Bond(a3, a4);
        assertNotNull(b1);
        assertNotNull(b2);
        assertNotNull(b3);
        an1 = new Angle(b1, b2);
        an2 = new Angle(b2, b3);
        assertNotNull(an1);
        assertNotNull(an2);
        torsion = new Torsion(an1, an2);
        assertNotNull(torsion);
    }

    @After
    public void tearDown() {
        assertTrue(a1.destroy());
        assertTrue(a2.destroy());
        assertTrue(a3.destroy());
        assertTrue(a4.destroy());
        assertTrue(b1.destroy());
        assertTrue(b2.destroy());
        assertTrue(b3.destroy());
        assertTrue(an1.destroy());
        assertTrue(an2.destroy());
        assertTrue(torsion.destroy());
        a1 = null;
        a2 = null;
        a3 = null;
        a4 = null;
        b1 = null;
        b2 = null;
        b3 = null;
        an1 = null;
        an2 = null;
        torsion = null;
    }
}
