/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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

import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.*;

/**
 *
 * @author fennt
 */
@RunWith(Parameterized.class)
public class ReflectionListTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                    {false,
                        "P6 test",
                        115.996, 115.996, 44.13, 90.0, 90.0, 120.0, "P6",
                        1.89631,
                        27200,
                        new HKL(0, 0, 4),
                        new HKL(-1, 0, 0),
                        6, 0},
                    {true,
                        "I222 test",
                        86.031, 92.854, 98.312, 90.0, 90.0, 90.0, "I222",
                        1.0,
                        210663,
                        new HKL(0, 0, 4),
                        new HKL(0, 0, 3),
                        8, 0},
                    {true,
                        "P1 test",
                        25.015, 29.415, 52.761, 89.54, 86.1, 82.39, "P1",
                        1.0,
                        80373,
                        new HKL(-24, -10, 1),
                        new HKL(-24, -10, 8),
                        1, 0}
                });
    }
    private final String info;
    private final int size;
    private final HKL valid;
    private final HKL invalid;
    private final int epsilon;
    private final int allowed;
    private final boolean ci;
    private final boolean ciOnly;
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;

    public ReflectionListTest(boolean ciOnly,
            String info, double a, double b, double c,
            double alpha, double beta, double gamma, String sg,
            double resolution,
            int size, HKL valid, HKL invalid, int epsilon, int allowed) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.size = size;
        this.valid = valid;
        this.invalid = invalid;
        this.epsilon = epsilon;
        this.allowed = allowed;

        String ffxCi = System.getProperty("ffx.ci");
        if (ffxCi != null && ffxCi.equalsIgnoreCase("true")) {
            ci = true;
        } else {
            ci = false;
        }

        if (!ci && ciOnly) {
            this.crystal = null;
            this.resolution = null;
            this.reflectionlist = null;
            return;
        }

        this.crystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
        this.resolution = new Resolution(resolution);
        this.reflectionlist = new ReflectionList(this.crystal, this.resolution);
    }

    @Test
    public void testsize() {
        if (!ci && ciOnly) {
            return;
        }

        assertEquals(info + " hash map and arraylist should have equal length",
                reflectionlist.hklmap.size(), reflectionlist.hkllist.size());
        assertEquals(info + " reflection list should have correct size",
                size, reflectionlist.hkllist.size());
    }

    @Test
    public void testreflections() {
        if (!ci && ciOnly) {
            return;
        }

        assertNotNull(info + " list should have valid reflection",
                reflectionlist.getHKL(valid));
        assertNull(info + " list should NOT have invalid reflection",
                reflectionlist.getHKL(invalid));

        HKL hkl = reflectionlist.getHKL(0, 0, 0);
        assertEquals(info + " list 0 0 0 reflection should have correct epsilon",
                epsilon, hkl.epsilon());
        assertEquals(info + " list 0 0 0 reflection should have correct allowance",
                allowed, hkl.allowed);
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
