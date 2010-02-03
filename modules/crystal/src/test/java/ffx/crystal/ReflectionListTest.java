/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.crystal;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fennt
 */
public class ReflectionListTest {

    Crystal crystal1 = new Crystal(86.031, 92.854, 98.312, 90.00, 90.00, 90.00, "I222");
    Resolution resolution1 = new Resolution(1.0);
    ReflectionList i222list = new ReflectionList(crystal1, resolution1);
    Crystal crystal2 = new Crystal(25.015, 29.415, 52.761, 89.54, 86.1, 82.39, "P1");
    Resolution resolution2 = new Resolution(1.0);
    ReflectionList p1list = new ReflectionList(crystal2, resolution2);

    public ReflectionListTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testP1sizes() {
        assertEquals("hash map and arraylist should have equal length",
                p1list.hklmap.size(), p1list.hkllist.size());
        assertEquals("P1 reflection list should have 80373 elements",
                80373, p1list.hkllist.size());
    }

    @Test
    public void testI222sizes() {
        assertEquals("hash map and arraylist should have equal length",
                i222list.hklmap.size(), i222list.hkllist.size());
        assertEquals("I222 reflection list should have 210663 elements",
                210663, i222list.hkllist.size());
    }

    @Test
    public void testP1reflections() {
        assertNotNull("P1 list should have -24 -10 1 reflection",
                p1list.getHKL(-24, -10, 1));
        assertNull("P1 list should NOT have -24 -10 8 reflection",
                p1list.getHKL(-24, -10, 8));

        HKL hkl = p1list.getHKL(0, 0, 0);
        assertEquals("P1 list 0 0 0 reflection should have epsilon of 1",
                1, hkl.epsilon());
        assertEquals("P1 list 0 0 0 reflection should NOT be allowed",
                0, hkl.allowed);
    }

    @Test
    public void testI222reflections() {
        assertNotNull("I222 list should have 0 0 4 reflection",
                i222list.getHKL(0, 0, 4));
        assertNull("I222 list should NOT have 0 0 3 reflection",
                i222list.getHKL(0, 0, 3));

        HKL hkl = i222list.getHKL(0, 0, 0);
        assertEquals("I222 list 0 0 0 reflection should have epsilon of 8",
                8, hkl.epsilon());
        assertEquals("I222 list 0 0 0 reflection should NOT be allowed",
                0, hkl.allowed);
    }

    /*
    @Test
    public void listreflections() {
    int nhkl = 0;
    int nbin[] = new int[10];
    for (HKL i : i222list.hkllist) {
    System.out.println(nhkl + ": " + i.h() + " " + i.k() + " " + i.l()
    + " eps: " + i.epsilon()
    + " allowed: " + i.allowed()
    + " bin: " + i.bin
    + " res: " + Crystal.res(i222list.crystal, i));
    nhkl++;
    nbin[i.bin]++;
    }
    System.out.println("bin  # HKL");
    for (int i = 0; i < 10; i++) {
    System.out.println(i + "  " + nbin[i]);
    }
    }
     */

    /*
    for (Iterator i = hkls.hkls.entrySet().iterator(); i.hasNext(); nhkl++) {
    Map.Entry ei = (Map.Entry) i.next();
    Object key = ei.getKey();
    HKL ih = (HKL) ei.getValue();

    System.out.println(nhkl + ": " + key + " " + ih.h() + " " + ih.k() + " " + ih.l()
    + " eps: " + ih.epsilon() + " res: " + Crystal.res(crystal, ih));
    }
     */
}
