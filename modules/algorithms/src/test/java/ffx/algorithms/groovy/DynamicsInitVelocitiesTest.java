/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms.groovy;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import edu.rit.pj.Comm;

import ffx.algorithms.groovy.Dynamics;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.PJDependentTest;

import groovy.lang.Binding;

/**
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsInitVelocitiesTest extends PJDependentTest {

    private String info;
    private String filename;
    private double endKineticEnergy;
    private double kineticEnergyTolerance = 5.0;
    private double endPotentialEnergy;
    private double potentialEnergyTolerance = 5.0;
    private double endTotalEnergy;
    private double totalEnergyTolerance = 5.0;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsInitVelocitiesTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "Water Tiny Initialize Velocities",  //info
                        "ffx/algorithms/structures/watertiny.xyz",  // filename
                        64.8511,  // endKineticEnergy
                        -218.0760,  // endPotentialEnergy
                        -153.2249  // endTotalEnergy
                }

        });
    }

    public DynamicsInitVelocitiesTest(String info, String filename, double endKineticEnergy, double endPotentialEnergy, double endTotalEnergy) {

        this.info = info;
        this.filename = filename;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
    }

    @Before
    public void before() {
        binding = new Binding();
        dynamics = new Dynamics();
        dynamics.setBinding(binding);
    }

    @Test
    public void testDynamicsInitVelocities() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "VelocityVerlet", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that final energies are within the tolerance for the dynamics trajectory.
        assertEquals(info + "End kinetic energy for initializing velocities test", endKineticEnergy, molDyn.getKineticEnergy(), kineticEnergyTolerance);

        assertEquals(info + "End potential energy for initializing velocities test", endPotentialEnergy, molDyn.getPotentialEnergy(), potentialEnergyTolerance);

        assertEquals(info + "End total energy for initializing velocities test", endTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }

}
