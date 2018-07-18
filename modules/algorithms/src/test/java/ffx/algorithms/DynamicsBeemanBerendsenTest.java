/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;

import edu.rit.pj.Comm;

import ffx.algorithms.groovy.Dynamics;

import groovy.lang.Binding;

/**
 *
 * @author hbernabe
 */
@RunWith(Parameterized.class)
public class DynamicsBeemanBerendsenTest {

    private String info;
    private String filename;
    private String restartFile;
    private double endKineticEnergy;
    private double kineticEnergyTolerance = 0.5;
    private double endPotentialEnergy;
    private double potentialEnergyTolerance = 0.5;
    private double endTotalEnergy;
    private double totalEnergyTolerance = 0.5;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsBeemanBerendsenTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide Beeman Berendsen", // info
                "ffx/algorithms/structures/acetamide_NVE.xyz", // filename
                "ffx/algorithms/structures/acetamide_NVE.dyn", // restartFile
                4.8087, // endKineticEnergy
                -29.7186, // endPotentialEnergy
                -24.9099 // endTotalEnergy
            }

        });
    }

    public DynamicsBeemanBerendsenTest(String info, String filename, String restartFile, double endKineticEnergy, double endPotentialEnergy,
            double endTotalEnergy) {

        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;

    }

    @Before
    public void before() {

        binding = new Binding();
        dynamics = new Dynamics();
        dynamics.setBinding(binding);

        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
        }
    }

    @Test
    public void testDynamicsBeemanBerendsen() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Beeman", "-b", "Berendsen", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that end energies are within the tolerance for the dynamics trajectory
        assertEquals(info + "End kinetic energy for Beeman integrator and Berendsen thermostat", endKineticEnergy, molDyn.getKineticEnergy(), kineticEnergyTolerance);

        assertEquals(info + "End potential energy for Beeman integrator and Berendsen thermostat", endPotentialEnergy, molDyn.getPotentialEnergy(), potentialEnergyTolerance);

        assertEquals(info + "End total energy for Beeman integrator and Berendsen thermostat", endTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }

}
