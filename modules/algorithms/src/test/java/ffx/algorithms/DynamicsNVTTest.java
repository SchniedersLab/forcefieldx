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
 * @author hbernabe
 */
@RunWith(Parameterized.class)
public class DynamicsNVTTest {

    private String info;
    private String filename;
    private String restartFile;
    private double startingTemp;
    private double tempTolerance = 6.0;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsNVTTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "Water Box NVT",                                                // info
                        "ffx/algorithms/structures/waterbox_eq.xyz",                    // filename
                        "ffx/algorithms/structures/waterbox_eq.dyn",                    // restartFile
                        298.15                                                          // startingTemp
                }
        });
    }

    public DynamicsNVTTest(String info, String filename, String restartFile, double startingTemp) {
        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.startingTemp = startingTemp;
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
    public void testDynamicsNVT() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "VelocityVerlet", "-b", "Bussi", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that temperature is within tolerance at the end of the dynamics trajectory.
        assertEquals(info + " End temperature for NVT test", startingTemp, molDyn.getTemperature(), tempTolerance);
    }
}
