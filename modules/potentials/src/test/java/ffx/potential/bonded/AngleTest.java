/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Biophysics Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */

package ffx.potential.bonded;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Unit tests for the Angle class.
 */
public class AngleTest {
	private final double[] a1d = { 0.0, 0.0, 0.0 };

	private final double[] a2d = { 1.0, 0.0, 0.0 };

	private final double[] a3d = { 1.0, 1.0, 0.0 };

	private final double[] a4d = { 1.0, 1.0, 1.0 };

	private final double[] a5d = { 1.0, 1.0, 2.0 };

	private Atom a1;

	private Atom a2;

	private Atom a3;

	private Atom a4;

	private Atom a5;

	private Bond b1;

	private Bond b2;

	private Bond b3;

	private Bond b4;

	private Angle angle;

	private Angle angle2;

	private Angle angle3;

	@Test(timeout = 500)
	public void Angle_getCommonBond() {
		Bond expectedReturn = b2;
		Bond actualReturn = angle.getCommonBond(angle2);
		assertEquals("return value", expectedReturn, actualReturn);
		actualReturn = angle.getCommonBond(null);
		assertNull(actualReturn);
		actualReturn = angle.getCommonBond(angle3);
		assertNull(actualReturn);
	}

	@Test(timeout = 500)
	public void Angle_getOtherBond() {
		Bond expectedReturn = b2;
		Bond actualReturn = angle.getOtherBond(b1);
		assertEquals("return value", expectedReturn, actualReturn);
		expectedReturn = b1;
		actualReturn = angle.getOtherBond(b2);
		assertEquals("return value", expectedReturn, actualReturn);
		actualReturn = angle.getOtherBond(null);
		assertNull("return value", actualReturn);
		actualReturn = angle.getOtherBond(b3);
		assertNull("return value", actualReturn);
	}

	@Before
	public void setUp() {
		a1 = new Atom("A1"); a1.setXYZ(a1d);
		a2 = new Atom("A2"); a2.setXYZ(a2d);
		a3 = new Atom("A3"); a3.setXYZ(a3d);
		a4 = new Atom("A4"); a4.setXYZ(a4d);
		a5 = new Atom("A5"); a5.setXYZ(a5d);
		b1 = new Bond(a1, a2);
		b2 = new Bond(a2, a3);
		b3 = new Bond(a3, a4);
		b4 = new Bond(a4, a5);
		angle = new Angle(b1, b2);
		assertNotNull(angle);
		angle2 = new Angle(b2, b3);
		assertNotNull(angle2);
		angle3 = new Angle(b3, b4);
		assertNotNull(angle3);
	}

	@After
	public void tearDown() {
		assertTrue(a1.destroy());
		assertTrue(a2.destroy());
		assertTrue(a3.destroy());
		assertTrue(a4.destroy());
		assertTrue(a5.destroy());
		assertTrue(b1.destroy());
		assertTrue(b2.destroy());
		assertTrue(b3.destroy());
		assertTrue(b4.destroy());
		a1 = null;
		a2 = null;
		a3 = null;
		a4 = null;
		a5 = null;
		b1 = null;
		b2 = null;
		b3 = null;
		assertTrue(angle.destroy());
		assertTrue(angle2.destroy());
		assertTrue(angle3.destroy());
		angle = null;
		angle2 = null;
		angle3 = null;
	}
}
