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
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsRESPANVETest {

    private String info;
    private String filename;
    private String restartFile;
    private double startingTotalEnergy;
    private double totalEnergyTolerance = 0.5;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsRESPANVETest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide RESPA NVE", // info
                "ffx/algorithms/structures/acetamide_NVE.xyz", // filename
                "ffx/algorithms/structures/acetamide_NVE.dyn", // restartFile
                -25.1958 // startingTotalEnergy
            }

        });

    }

    public DynamicsRESPANVETest(String info, String filename, String restartFile, double startingTotalEnergy) {

        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.startingTotalEnergy = startingTotalEnergy;
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
    public void testRESPANVE() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "RESPA", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that the final total energy is within the tolerance for the molecular dynamics trajectory
        assertEquals(info + "End total energy for RESPA integrator NVE", startingTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }

}
