/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.xray;

import ffx.crystal.HKL;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AtomType;

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
public class FormFactorTest {

    private Atom carbon;
    private FormFactor carbonff;

    public FormFactorTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        double d[] = new double[3];
        d[0] = d[1] = d[2] = 0.0;
        double anisou[] = new double[6];
        anisou[0] = anisou[1] = anisou[2] = 1.0;
        anisou[3] = anisou[4] = anisou[5] = 0.0;
        carbon = new Atom(1, "C", 'A', d, "ALA", 1, "A", 1.0, 20.0);
        AtomType atomType = new AtomType(1, 1, "C", null, 6, 12.01, 1);
        carbon.setAtomType(atomType);
        carbon.setAltLoc('A');
        carbon.setAnisou(anisou);
        carbonff = new FormFactor(carbon);
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testCarbonFF() {
        double ff[][] = new double[2][6];

        assertNotNull("carbon form factors should exist",
                FormFactor.getFormFactor("6"));

        ff = FormFactor.getFormFactor("6");

        assertEquals("carbon form factors should be correct",
                2.09921, ff[0][0], 0.0001);
        assertEquals("carbon form factors should be correct",
                13.18997, ff[1][0], 0.0001);
    }

    @Test
    public void testCarbonfrho() {
        HKL hkl = new HKL(1, 1, 1);
        assertEquals("carbon (1 1 1) structure factor should be correct",
                2.39874e-26, carbonff.f(hkl), 1e-30);

        double xyz[] = {1.0, 1.0, 1.0};
        assertEquals("carbon (1 1 1) electron density should be correct",
                0.0819388, carbonff.rho(xyz), 0.000001);
    }
}
