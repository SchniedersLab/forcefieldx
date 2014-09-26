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
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

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

    @Test
    public void testCIFFilter3DYC() {
        String filename = "ffx/xray/structures/3DYC.ent";
        ClassLoader cl = this.getClass().getClassLoader();
        File cifFile = new File(cl.getResource(filename).getPath());
        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(cifFile);

        CIFFilter cifFilter = new CIFFilter();
        ReflectionList reflectionList = cifFilter.getReflectionList(cifFile);
        assertNotNull(" Reflection list should not be null", reflectionList);
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties,
                reflectionList);

        assertTrue(" CIF data not read in correctly",
                cifFilter.readFile(cifFile, reflectionList, refinementData, properties));

        HKL hkl = reflectionList.getHKL(58, 0, 13);
        assertEquals("58 0 13 F",
                99.7, refinementData.getF(hkl.index()), 0.01);
        assertEquals("58 0 13 sigF",
                69.7, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("58 0 13 freeR value",
                1, refinementData.freer[hkl.index()]);

        hkl = reflectionList.getHKL(28, 20, 5);
        assertEquals("28 20 5 F",
                428.1, refinementData.getF(hkl.index()), 0.01);
        assertEquals("28 20 5 sigF",
                10.1, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("28 20 5 freeR value",
                0, refinementData.freer[hkl.index()]);
    }

    @Test
    public void testCIFFilter2DRM() {
        String filename = "ffx/xray/structures/2DRM.cif";
        ClassLoader cl = this.getClass().getClassLoader();
        File cifFile = new File(cl.getResource(filename).getPath());
        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(cifFile);

        CIFFilter cifFilter = new CIFFilter();
        ReflectionList reflectionList = cifFilter.getReflectionList(cifFile);
        assertNull(" Reflection list should be null", reflectionList);

        Crystal crystal = new Crystal(29.969, 37.861, 44.506,
                90.28, 90.11, 90.64, "P1");
        Resolution resolution = new Resolution(1.30);
        reflectionList = new ReflectionList(crystal, resolution);
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties,
                reflectionList);

        assertTrue(" CIF data not read correctly",
                cifFilter.readFile(cifFile, reflectionList, refinementData,
                        properties));

        HKL hkl = reflectionList.getHKL(-21, -6, 7);
        assertEquals("-21 -6 7 F",
                18.6, refinementData.getF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF",
                3.6, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value",
                0, refinementData.freer[hkl.index()]);

        hkl = reflectionList.getHKL(-21, -6, 8);
        assertEquals("-21 -6 7 F",
                20.2, refinementData.getF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF",
                5.0, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value",
                1, refinementData.freer[hkl.index()]);
    }
}
