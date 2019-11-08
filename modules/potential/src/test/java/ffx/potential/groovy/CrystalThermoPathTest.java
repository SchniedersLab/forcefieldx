/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.groovy;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import ffx.potential.cli.PotentialScript;
import ffx.potential.groovy.test.LambdaGradient;

import groovy.lang.Binding;

/**
 * Tests test.LambdaGradient command to determine that the end state potentials and derivatives are correct.
 *
 * @author Aaron J. Nessler
 * @author Mallory R. Tollefson
 */
public class CrystalThermoPathTest extends PotentialScript {

    Binding binding;
    LambdaGradient lambdaGradient;

    @Before
    public void before() {
        binding = new Binding();
        lambdaGradient = new ffx.potential.groovy.test.LambdaGradient();
        lambdaGradient.setBinding(binding);
    }

    @After
    public void after() {
        lambdaGradient.destroyPotentials();
        System.gc();
    }

    @Test
    public void testLambdaGradientHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        lambdaGradient.run();
    }

    /**
     * Tests the End States of the LambdaGradient class.
     */
    @Test
    public void testLambdaGradient() {

        String[] args = {"--s1", "1",
                "--f1", "26",
                "src/main/java/ffx/potential/structures/phenacetin.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script.
        lambdaGradient.run();

        System.clearProperty("lambdaterm");

        double expectedPotentialEnergyVac =   -7.70864934;
        double expectedPotentialEnergyXtal = -36.22707541638729;
        double actualPotentialEnergyVac = lambdaGradient.e0;
        double actualPotentialEnergyXtal = lambdaGradient.e1;

        assertEquals(expectedPotentialEnergyVac, actualPotentialEnergyVac, 1E-6);
        assertEquals(expectedPotentialEnergyXtal, actualPotentialEnergyXtal, 1E-6);
        assertEquals(0, lambdaGradient.ndEdLFailures);
        assertEquals(0, lambdaGradient.nd2EdL2Failures);
        assertEquals(0, lambdaGradient.ndEdXdLFailures);
        assertEquals(0, lambdaGradient.ndEdXFailures);
    }

    /**
     * Tests the End States of the LambdaGradient class when softcore is active.
     */
    @Test
    public void testLambdaGradientIntermolecularSoftcore() {

        String[] args = {"--s1", "1",
                "--f1", "44",
                "src/main/java/ffx/potential/structures/ethylparaben.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script.
        lambdaGradient.run();

        System.clearProperty("lambdaterm");

        double expectedPotentialEnergyVac =   -6.67890842;
        double expectedPotentialEnergyXtal = -57.410945848200896;
        double actualPotentialEnergyVac = lambdaGradient.e0;
        double actualPotentialEnergyXtal = lambdaGradient.e1;

        assertEquals(expectedPotentialEnergyVac, actualPotentialEnergyVac, 1E-6);
        assertEquals(expectedPotentialEnergyXtal, actualPotentialEnergyXtal, 1E-6);
        assertEquals(0, lambdaGradient.ndEdLFailures);
        assertEquals(0, lambdaGradient.nd2EdL2Failures);
        assertEquals(0, lambdaGradient.ndEdXdLFailures);
        assertEquals(0, lambdaGradient.ndEdXFailures);
    }
}


