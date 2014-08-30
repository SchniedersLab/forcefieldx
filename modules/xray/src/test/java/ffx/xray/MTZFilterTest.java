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
 * @author Timothy D. Fenn
 */
public class MTZFilterTest {

    String filename = "ffx/xray/structures/2DRM.mtz";
    ClassLoader cl = this.getClass().getClassLoader();
    File mtzfile = new File(cl.getResource(filename).getPath());
    // load any properties associated with it
    CompositeConfiguration properties = Keyword.loadProperties(mtzfile);
    // set up the crystal data
    Crystal crystal
            = new Crystal(29.97, 37.86, 44.51, 90.28, 90.11, 90.64, "P1");
    Resolution resolution = new Resolution(1.30);
    ReflectionList reflectionlist = new ReflectionList(crystal, resolution);
    DiffractionRefinementData refinementdata = new DiffractionRefinementData(properties, reflectionlist);

    public MTZFilterTest() {
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
    public void testMTZReadFile() {
        MTZFilter mtzfilter = new MTZFilter();
        assertTrue("mtz file should be read in without errors",
                mtzfilter.readFile(mtzfile, reflectionlist, refinementdata, null));
    }

    @Test
    public void testMTZParams() {
        MTZFilter mtzfilter = new MTZFilter();
        assertTrue("mtz file should be read in without errors",
                mtzfilter.readFile(mtzfile, reflectionlist, refinementdata, null));

        assertEquals("mtz file should have correct number of columns",
                6, mtzfilter.nColumns);

        assertEquals("mtz file should have correct number of reflections",
                48115, mtzfilter.nReflections);

        assertEquals("mtz file should have correct space group number",
                1, mtzfilter.sgnum);
    }

    @Test
    public void testMTZHKL() {
        MTZFilter mtzfilter = new MTZFilter();
        assertTrue("mtz file should be read in without errors",
                mtzfilter.readFile(mtzfile, reflectionlist, refinementdata, null));

        HKL hkl = reflectionlist.getHKL(-10, 1, 1);
        assertEquals("-10 1 1 FP value should be correct",
                229.90, refinementdata.get_f(hkl.index()), 0.02);
        assertEquals("-10 1 1 SIGFP value should be correct",
                2.50, refinementdata.get_sigf(hkl.index()), 0.02);
        assertEquals("-10 1 1 FREE value should be correct",
                1, refinementdata.get_freer(hkl.index()));

        hkl = reflectionlist.getHKL(-10, 1, 10);
        assertEquals("-10 1 10 FP value should be NaN",
                Double.NaN, refinementdata.get_f(hkl.index()), 0.1);
        assertEquals("-10 1 10 SIGFP value should be NaN",
                Double.NaN, refinementdata.get_sigf(hkl.index()), 0.1);
    }

    @Test
    public void testMTZReflectionList() {
        MTZFilter mtzfilter = new MTZFilter();

        ReflectionList tmp = mtzfilter.getReflectionList(mtzfile);
        assertNotNull("mtz file should return a reflection list", tmp);

        assertEquals("reflection list should have correct space group",
                1, tmp.spaceGroup.number);

        assertEquals("reflection list should have correct a cell dimension",
                29.969, tmp.crystal.a, 0.01);

        assertEquals("reflection list should have correct resolution",
                1.30, tmp.resolution.res_limit(), 0.01);
    }
}
