/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.potential.parsers;

import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.Test;

import ffx.potential.parameters.ForceField;
import ffx.utilities.Keyword;

/**
 * Test creation of a Force Field object from a CompositeConfiguration.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class ForceFieldFilterTest {

    public ForceFieldFilterTest() {
    }

    /**
     * Test of parse method, of class ForceFieldFilter.
     */
    @Test
    public void testParse() {
        CompositeConfiguration properties = Keyword.loadProperties(null);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();
    }
}
