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

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;

import edu.rit.pj.Comm;

import ffx.algorithms.groovy.DynamicsOpenMM;

import groovy.lang.Binding;

/**
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsOpenMMRESPANVETest {

    private String info;
    private String filename;
    private double startingTotalEnergy;
    private double totalEnergyTolerance = 5.0;
    private boolean ffxOpenMM;

    private Binding binding;
    private DynamicsOpenMM dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsOpenMMRESPANVETest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "System OpenMM RESPA NVE", // info
                        "ffx/algorithms/structures/waterbox_eq.xyz", // filename
                        -25240.5415 // startingTotalEnergy
                }

        });

    }

    public DynamicsOpenMMRESPANVETest(String info, String filename, double startingTotalEnergy) {

        this.info = info;
        this.filename = filename;
        this.startingTotalEnergy = startingTotalEnergy;

        ffxOpenMM = System.getProperty("ffx.openMM", "false").equalsIgnoreCase("true");
    }

    @Before
    public void before() {

        binding = new Binding();
        dynamics = new DynamicsOpenMM();
        dynamics.setBinding(binding);

        System.setProperty("platform", "omm");

        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
        }
    }

    @After
    public void after() {
        System.clearProperty("platform");
    }

    @Test
    public void testDynamicsOpenMMRESPANVE() {

        if (!ffxOpenMM) {
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-z", "1", "-t", "298.15", "-i", "RESPA", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate script.
        dynamics.run();

        MolecularDynamicsOpenMM molDyn = dynamics.getMolecularDynamics();

        // Assert that the end total energy is within the threshold for the dynamics trajectory.
        assertEquals(info + "End total energy for OpenMM RESPA integrator under the NVE ensemble", startingTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }

}
