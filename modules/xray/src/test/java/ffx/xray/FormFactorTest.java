/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import ffx.crystal.HKL;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AtomType;

/**
 * @author Timothy D. Fenn
 */
public class FormFactorTest {

    private Atom carbon;
    private XRayFormFactor carbonFormFactor;

    @Before
    public void setUp() {
        double d[] = new double[3];
        d[0] = d[1] = d[2] = 0.0;
        double anisou[] = new double[6];
        anisou[0] = anisou[1] = anisou[2] = 1.0;
        anisou[3] = anisou[4] = anisou[5] = 0.0;
        carbon = new Atom(1, "C", 'A', d, "ALA", 1, 'A', 1.0, 20.0, "A");
        AtomType atomType = new AtomType(1, 1, "C", null, 6, 12.01, 1);
        carbon.setAtomType(atomType);
        carbon.setAltLoc('A');
        carbon.setAnisou(anisou);
        carbonFormFactor = new XRayFormFactor(carbon, false);
    }

    @Test
    public void testCarbonFF() {
        double formFactor[][] = new double[2][6];

        assertNotNull(" Carbon form factors should exist",
                XRayFormFactor.getFormFactor("6"));

        formFactor = XRayFormFactor.getFormFactor("6");

        assertEquals(" Carbon form factors",
                5, (int) formFactor[0][0]);
        assertEquals(" Carbon form factors",
                2.09921, formFactor[1][0], 0.0001);
        assertEquals(" Carbon form factors",
                13.18997, formFactor[2][0], 0.0001);
    }

    @Test
    public void testCarbonfrho() {
        HKL hkl = new HKL(1, 1, 1);
        assertEquals("carbon (1 1 1) structure factor should be correct",
                2.3986e-26, carbonFormFactor.f(hkl), 1e-30);

        double xyz[] = {1.0, 1.0, 1.0};
        assertEquals("carbon (1 1 1) electron density should be correct",
                0.081937, carbonFormFactor.rho(0.0, 1.0, xyz), 0.000001);
    }
}
