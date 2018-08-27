/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.algorithms.groovy;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;

import ffx.algorithms.groovy.Dynamics;
import ffx.algorithms.MolecularDynamics;

import groovy.lang.Binding;

/**
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsStochasticTest extends ffx.algorithms.PJDependentTest {

    private String info;
    private String filename;
    private double endKineticEnergy;
    private double endPotentialEnergy;
    private double endTotalEnergy;
    private boolean testSeed;
    private boolean testFriction00;
    private boolean testFriction01;
    private double tolerance = 0.1;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsStochasticTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide Peptide Restart and Stochastic Random Seed", // info
                "ffx/algorithms/structures/acetamide_res_stoch.xyz", // filename
                6.8546, // endKineticEnergy
                -26.9921, // endPotentialEnergy
                -20.1375, // endTotalEnergy
                true, // testSeed
                false, // testFriction00
                false, // testFriction01
            },
            {
                "Acetamide Peptide Restart, Stochastic Random Seed and Friction 0.0", // info
                "ffx/algorithms/structures/acetamide_res_stoch.xyz", //filename
                4.5625, // endKineticEnergy
                -29.8043, // endPotentialEnergy
                -25.2418, // endTotalEnergy
                false, // testSeed
                true, // testFriction00
                false // testFriction01
            },
            {
                "Acetamide Peptide Restart, Stochastic Random Seed and Friction 0.1", // info
                "ffx/algorithms/structures/acetamide_res_stoch.xyz", // filename
                4.5743, // endKineticEnergy
                -29.7373, // endPotentialEnergy
                -25.1630, // endTotalEnergy
                false, // testSeed
                false, // testFriction00
                true // testFriction01
            }
        });
    }

    public DynamicsStochasticTest(String info, String filename, double endKineticEnergy,
            double endPotentialEnergy, double endTotalEnergy, boolean testSeed, boolean testFriction00, boolean testFriction01) {

        this.info = info;
        this.filename = filename;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
        this.testSeed = testSeed;
        this.testFriction00 = testFriction00;
        this.testFriction01 = testFriction01;
    }

    @Before
    public void before() {
        binding = new Binding();
        dynamics = new Dynamics();
        dynamics.setBinding(binding);
    }

    @After
    public void after() {
        dynamics.destroyPotentials();
        System.gc();
    }

    @Test
    public void testDynamicsStochasticRandomSeed() {
        
        if (!testSeed){
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        //Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that the end energies are the same meaning that the Stochastic integrator works as intended.
        assertEquals(info + " Final kinetic energy", endKineticEnergy, molDyn.getKineticEnergy(), tolerance);
        assertEquals(info + " Final potential energy", endPotentialEnergy, molDyn.getPotentialEnergy(), tolerance);
        assertEquals(info + " Final total energy", endTotalEnergy, molDyn.getTotalEnergy(), tolerance);
        
        //dynamics.destroyPotentials();
    }

    @Test
    public void testDynamicsStochasticRandomSeedFriction0() {
        
        if (!testFriction00){
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        System.setProperty("friction", "0.0");

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that the end energies are the same meaning that the Stochastic integrator works as intended with friciton 0.0 collisions/picosecond.
        assertEquals(info + " Final kinetic energy", endKineticEnergy, molDyn.getKineticEnergy(), tolerance);
        assertEquals(info + " Final potential energy", endPotentialEnergy, molDyn.getPotentialEnergy(), tolerance);
        assertEquals(info + " Final total energy", endTotalEnergy, molDyn.getTotalEnergy(), tolerance);

        System.clearProperty("friction");
        
        //dynamics.destroyPotentials();
    }

    @Test
    public void testDynamicsStochasticRandomSeedFriction01() {
        
        if (!testFriction01){
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        System.setProperty("friction", "0.1");

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that the end energies are the same meaning that the Stochastic integrator works as intended with friciton 0.0 collisions/picosecond.
        assertEquals(info + " Final kinetic energy", endKineticEnergy, molDyn.getKineticEnergy(), tolerance);
        assertEquals(info + " Final potential energy", endPotentialEnergy, molDyn.getPotentialEnergy(), tolerance);
        assertEquals(info + " Final total energy", endTotalEnergy, molDyn.getTotalEnergy(), tolerance);

        System.clearProperty("friction");
    }

}
