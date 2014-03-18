/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.xray;

import java.io.File;

import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.*;

import static org.junit.Assert.*;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.utilities.Keyword;

/**
 *
 * @author fenn
 */
public class CIFFilterTest {

    public CIFFilterTest() {
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
    public void testCIFFilter3DYC() {
        String filename = "ffx/xray/structures/3DYC.ent";
        ClassLoader cl = this.getClass().getClassLoader();
        File ciffile = new File(cl.getResource(filename).getPath());
        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(ciffile);

        CIFFilter ciffilter = new CIFFilter();
        ReflectionList reflectionlist = ciffilter.getReflectionList(ciffile);
        assertNotNull("reflection list should be nonnull", reflectionlist);
        DiffractionRefinementData refinementdata = new DiffractionRefinementData(properties,
                reflectionlist);

        assertTrue("CIF data shoud be read in correctly",
                ciffilter.readFile(ciffile, reflectionlist, refinementdata,
                properties));

        HKL hkl = reflectionlist.getHKL(58, 0, 13);
        assertEquals("58 0 13 F should be correct",
                99.7, refinementdata.get_f(hkl.index()), 0.01);
        assertEquals("58 0 13 sigF should be correct",
                69.7, refinementdata.get_sigf(hkl.index()), 0.01);
        assertEquals("58 0 13 freeR value should be correct",
                1, refinementdata.freer[hkl.index()]);

        hkl = reflectionlist.getHKL(28, 20, 5);
        assertEquals("28 20 5 F should be correct",
                428.1, refinementdata.get_f(hkl.index()), 0.01);
        assertEquals("28 20 5 sigF should be correct",
                10.1, refinementdata.get_sigf(hkl.index()), 0.01);
        assertEquals("28 20 5 freeR value should be correct",
                0, refinementdata.freer[hkl.index()]);
    }

    @Test
    public void testCIFFilter2DRM() {
        String filename = "ffx/xray/structures/2DRM.cif";
        ClassLoader cl = this.getClass().getClassLoader();
        File ciffile = new File(cl.getResource(filename).getPath());
        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(ciffile);

        CIFFilter ciffilter = new CIFFilter();
        ReflectionList reflectionlist = ciffilter.getReflectionList(ciffile);
        assertNull("reflection list should be null", reflectionlist);

        Crystal crystal = new Crystal(29.969, 37.861, 44.506,
                90.28, 90.11, 90.64, "P1");
        Resolution resolution = new Resolution(1.30);
        reflectionlist = new ReflectionList(crystal, resolution);
        DiffractionRefinementData refinementdata = new DiffractionRefinementData(properties,
                reflectionlist);

        assertTrue("CIF data shoud be read in correctly",
                ciffilter.readFile(ciffile, reflectionlist, refinementdata,
                properties));

        HKL hkl = reflectionlist.getHKL(-21, -6, 7);
        assertEquals("-21 -6 7 F should be correct",
                18.6, refinementdata.get_f(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF should be correct",
                3.6, refinementdata.get_sigf(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value should be correct",
                0, refinementdata.freer[hkl.index()]);

        hkl = reflectionlist.getHKL(-21, -6, 8);
        assertEquals("-21 -6 7 F should be correct",
                20.2, refinementdata.get_f(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF should be correct",
                5.0, refinementdata.get_sigf(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value should be correct",
                1, refinementdata.freer[hkl.index()]);
    }
}
