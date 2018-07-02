/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential;

import java.io.File;
import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import groovy.lang.Binding;
import groovy.lang.GroovyShell;

/**
 * Test the Biotype script.
 */
public class BiotypeTest {
    private GroovyShell shell;
    private Binding binding;

    @Before
    public void setUp() {
        binding = new Binding();
        shell = new GroovyShell(binding);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testBiotype() {
        try {
            // Set-up the input arguments for the Biotype script.
            String args[] = {"src/main/java/ffx/potential/structures/acetanilide.xyz"};
            binding.setVariable("args", args);

            // Evaluate the script.
            shell.evaluate(new File("src/main/groovy/ffx/potential/Biotype.groovy"));

            // Pull out the biotype results to check.
            List<String> biotypes = (List<String>) binding.getVariable("biotypes");
            Assert.assertNotNull(biotypes);
            Assert.assertEquals(19, biotypes.size());
            Assert.assertTrue(" Check the value of the first Biotype.",
                    biotypes.get(0).trim().equalsIgnoreCase("biotype   1    C \"ace\" 405    C    C    N"));
        } catch (Exception e) {
            System.out.println(e.toString());
            Assert.fail(" Failed to evaluate the Biotype script.");
        }
    }
}
