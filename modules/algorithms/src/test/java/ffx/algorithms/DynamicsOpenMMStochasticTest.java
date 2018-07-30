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

import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import edu.rit.pj.Comm;

import ffx.algorithms.groovy.DynamicsOpenMM;

import groovy.lang.Binding;

/**
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsOpenMMStochasticTest {

    private String info;
    private String filename;
    private double endKineticEnergy;
    private double endPotentialEnergy;
    private double endTotalEnergy;
    private double tolerance = 5.0;
    private boolean ffxOpenMM;

    private Binding binding;
    private DynamicsOpenMM dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsOpenMMStochasticTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "System OpenMM Stochastic", // info
                        "ffx/algorithms/structures/waterbox_eq.xyz", // filename
                        // 11796.9508,  endKineticEnergy
                        // -36782.1559, endPotentialEnergy
                        // -24985.2051  endTotalEnergy
                        11017.6784,  -36276.1969, -25258.5186
                }
        });
    }

    public DynamicsOpenMMStochasticTest(String info, String filename, double endKineticEnergy, double endPotentialEnergy,
                                        double endTotalEnergy) {

        this.info = info;
        this.filename = filename;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;

        ffxOpenMM = System.getProperty("ffx.openMM", "false").equalsIgnoreCase("true");
    }

    @Before
    public void before() {
        binding = new Binding();
        dynamics = new DynamicsOpenMM();
        dynamics.setBinding(binding);
    }

    @BeforeClass
    public static void beforeClass() {
        // Initialize Parallel Java

        System.setProperty("platform", "omm");
        try {
            Comm.world();
        } catch (IllegalStateException ise) {
            try {
                String args[] = new String[0];
                Comm.init(args);
            } catch (Exception e) {
                String message = String.format(" Exception starting up the Parallel Java communication layer.");
                logger.log(Level.WARNING, message, e.toString());
                message = String.format(" Skipping Beeman/Berendsen NVT dynamics test.");
                logger.log(Level.WARNING, message, e.toString());
                fail();
            }
        }
    }

    @AfterClass
    public static void afterClass() {
        System.clearProperty("platform");
    }

    @Test
    public void testDynamicsOpenMMStochastic() {

        if (!ffxOpenMM) {
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-z", "1", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate script
        dynamics.run();

        MolecularDynamicsOpenMM molDynOpenMM = dynamics.getMolecularDynamics();

        // Assert that the end energies are within the threshold for the dynamics trajectory.
        assertEquals(info + "End kinetic energy for OpenMM Langevin(Stochastic) integrator",
                endKineticEnergy, molDynOpenMM.getKineticEnergy(), tolerance);
        assertEquals(info + "End potential energy for OpenMM Langevin(Stochastic) integrator",
                endPotentialEnergy, molDynOpenMM.getPotentialEnergy(), tolerance);
        assertEquals(info + "End total energy for OpenMM Langevin(Stochastic) integrator",
                endTotalEnergy, molDynOpenMM.getTotalEnergy(), tolerance);
    }

}
