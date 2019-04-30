/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.potential.bonded;

import java.util.ArrayList;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Unit tests for the MSNode class.
 */
public class MSNodeTest {

    private MSNode dataNode = null;

    @Test(timeout = 500)
    public void MSNode_constructor() {
        String n = "Test";
        assertEquals("MSNode", n, dataNode.getName());
    }

    @Test(timeout = 500)
    public void MSNode_destroy() {
        dataNode.setSelected(true);
        boolean expectedReturn = true;
        boolean actualReturn = dataNode.destroy();
        assertEquals("return value", expectedReturn, actualReturn);
        assertTrue(!dataNode.isSelected());
        assertNull(dataNode.getName());
        assertNull(dataNode.getParent());
        assertNotNull(dataNode.getAtomList());
        assertNotNull(dataNode.getBondList());
        assertNotNull(dataNode.getList(BondedTerm.class, new ArrayList<>()));
        assertNotNull(dataNode.getChildList());
    }

    @Test(timeout = 500)
    public void MSNode_equals() {
        Object object = "";
        boolean expectedReturn = false;
        boolean actualReturn = dataNode.equals(object);
        assertEquals("return value", expectedReturn, actualReturn);
        actualReturn = dataNode.equals(null);
        assertEquals("return value", expectedReturn, actualReturn);
        object = new MSNode("Test");
        expectedReturn = true;
        actualReturn = dataNode.equals(object);
        assertEquals("return value", expectedReturn, actualReturn);
        expectedReturn = true;
        actualReturn = dataNode.equals(dataNode);
        assertEquals("return value", expectedReturn, actualReturn);
    }

    @Before
    public void setUp() {
        dataNode = new MSNode("Test");
    }

    @After
    public void tearDown() {
        dataNode = null;
    }
}
