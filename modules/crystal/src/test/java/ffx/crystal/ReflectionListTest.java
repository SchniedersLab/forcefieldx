//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.crystal;

import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import ffx.utilities.BaseFFXTest;

/**
 * Test the ReflectionList class.
 *
 * @author Timothy D. Fenn
 */
@RunWith(Parameterized.class)
public class ReflectionListTest extends BaseFFXTest {

    private final String info;
    private final int size;
    private final HKL valid;
    private final HKL invalid;
    private final int epsilon;
    private final int allowed;
    private final boolean ciOnly;
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

        if (!ffxCI && ciOnly) {
            this.reflectionlist = null;
            return;
        }

        Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
        Resolution res = new Resolution(resolution);
        this.reflectionlist = new ReflectionList(crystal, res);
    }

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

    @Test
    public void testreflections() {
        if (!ffxCI && ciOnly) {
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

    @Test
    public void testsize() {
        if (!ffxCI && ciOnly) {
            return;
        }

        assertEquals(info + " hash map and arraylist should have equal length",
                reflectionlist.hklmap.size(), reflectionlist.hkllist.size());
        assertEquals(info + " reflection list should have correct size",
                size, reflectionlist.hkllist.size());
    }
}
