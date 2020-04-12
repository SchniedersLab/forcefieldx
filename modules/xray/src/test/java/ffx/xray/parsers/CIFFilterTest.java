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
package ffx.xray.parsers;

import java.io.File;

import org.apache.commons.configuration2.CompositeConfiguration;
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
import ffx.xray.DiffractionRefinementData;

/**
 * Test the CIF Filter.
 *
 * @author Timoth D. Fenn
 */
public class CIFFilterTest {

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
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties, reflectionList);

        assertTrue(" CIF data not read correctly",
                cifFilter.readFile(cifFile, reflectionList, refinementData, properties));

        HKL hkl = reflectionList.getHKL(-21, -6, 7);
        assertEquals("-21 -6 7 F", 18.6, refinementData.getF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF", 3.6, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value", 0, refinementData.freeR[hkl.index()]);

        hkl = reflectionList.getHKL(-21, -6, 8);
        assertEquals("-21 -6 7 F", 20.2, refinementData.getF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 sigF", 5.0, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("-21 -6 7 freeR value", 1, refinementData.freeR[hkl.index()]);
    }

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
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties, reflectionList);

        assertTrue(" CIF data not read in correctly",
                cifFilter.readFile(cifFile, reflectionList, refinementData, properties));

        HKL hkl = reflectionList.getHKL(58, 0, 13);
        assertEquals("58 0 13 F", 99.7, refinementData.getF(hkl.index()), 0.01);
        assertEquals("58 0 13 sigF", 69.7, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("58 0 13 freeR value", 1, refinementData.freeR[hkl.index()]);

        hkl = reflectionList.getHKL(28, 20, 5);
        assertEquals("28 20 5 F", 428.1, refinementData.getF(hkl.index()), 0.01);
        assertEquals("28 20 5 sigF", 10.1, refinementData.getSigF(hkl.index()), 0.01);
        assertEquals("28 20 5 freeR value", 0, refinementData.freeR[hkl.index()]);
    }
}
