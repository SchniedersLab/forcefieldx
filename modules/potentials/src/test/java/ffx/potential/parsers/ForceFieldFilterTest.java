/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
