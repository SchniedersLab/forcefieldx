//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.xray.parsers;

import java.io.File;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.utilities.Keyword;
import ffx.xray.DiffractionRefinementData;

/**
 * Test the MTZ Filter.
 *
 * @author Timothy D. Fenn
 */
public class MTZFilterTest {

    private String filename = "ffx/xray/structures/2DRM.mtz";
    private ClassLoader cl = this.getClass().getClassLoader();
    private File mtzFile = new File(cl.getResource(filename).getPath());
    // Load any properties associated with it
    private CompositeConfiguration properties = Keyword.loadProperties(mtzFile);
    // set up the crystal data
    private Crystal crystal = new Crystal(29.97, 37.86, 44.51, 90.28, 90.11, 90.64, "P1");
    private Resolution resolution = new Resolution(1.30);
    private ReflectionList reflectionList = new ReflectionList(crystal, resolution);
    private DiffractionRefinementData refinementData = new DiffractionRefinementData(properties, reflectionList);

    @Test
    public void testMTZReadFile() {
        MTZFilter mtzFilter = new MTZFilter();
        assertTrue("mtz file should be read in without errors",
                mtzFilter.readFile(mtzFile, reflectionList, refinementData, null));
    }

    @Test
    public void testMTZParams() {
        MTZFilter mtzFilter = new MTZFilter();
        assertTrue("mtz file should be read in without errors",
                mtzFilter.readFile(mtzFile, reflectionList, refinementData, null));
        assertEquals("mtz file number of columns",
                6, mtzFilter.getnColumns());
        assertEquals("mtz file number of reflections",
                48115, mtzFilter.getnReflections());
        assertEquals("mtz file space group number",
                1, mtzFilter.getSpaceGroupNum());
    }

    @Test
    public void testMTZHKL() {
        MTZFilter mtzFilter = new MTZFilter();
        assertTrue("mtz file errors",
                mtzFilter.readFile(mtzFile, reflectionList, refinementData, null));
        HKL hkl = reflectionList.getHKL(-10, 1, 1);
        assertEquals("-10 1 1 FP value",
                229.90, refinementData.getF(hkl.index()), 0.02);
        assertEquals("-10 1 1 SIGFP value",
                2.50, refinementData.getSigF(hkl.index()), 0.02);
        assertEquals("-10 1 1 FREE value",
                1, refinementData.getFreeR(hkl.index()));
        hkl = reflectionList.getHKL(-10, 1, 10);
        assertEquals("-10 1 10 FP value should be NaN",
                Double.NaN, refinementData.getF(hkl.index()), 0.1);
        assertEquals("-10 1 10 SIGFP value should be NaN",
                Double.NaN, refinementData.getSigF(hkl.index()), 0.1);
    }

    @Test
    public void testMTZReflectionList() {
        MTZFilter mtzFilter = new MTZFilter();
        ReflectionList tmp = mtzFilter.getReflectionList(mtzFile);
        assertNotNull("mtz file should return a reflection list", tmp);
        assertEquals("reflection list space group",
                1, tmp.spaceGroup.number);
        assertEquals("reflection list cell dimension",
                29.969, tmp.crystal.a, 0.01);
        assertEquals("reflection list resolution",
                1.30, tmp.resolution.resolutionLimit(), 0.01);
    }
}
