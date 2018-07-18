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
public class DynamicsNVETest {

    private String info;
    private String filename;
    private String restartFile;
    private double startingTotalEnergy;
    private double energyTolerance = 0.5;
    private double startingTemp;
    private double tempTolerance = 6.0;
    private double endKineticEnergy;
    private double kineticEnergyTolerance = 2.0;
    private double endPotentialEnergy;
    private double potentialEnergyTolerance = 2.0;
    private double endTotalEnergy;
    private double totalEnergyTolerance = 2.0;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsNVETest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide Peptide NVE", // info
                "ffx/algorithms/structures/acetamide_NVE.xyz", // filename
                "ffx/algorithms/structures/acetamide_NVE.dyn", // restartFile
                -25.1958, // startingTotalEnergy
                298.15, // startingTemp
                4.5625, // endKineticEnergy
                -29.8043, // endPotentialEnergy
                -25.2418 // endTotalEnergy
            }
        });
    }

    public DynamicsNVETest(String info, String filename, String restartFile, double startingTotalEnergy, double startingTemp, double endKineticEnergy,
            double endPotentialEnergy, double endTotalEnergy) {
        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.startingTotalEnergy = startingTotalEnergy;
        this.startingTemp = startingTemp;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
    }

    @Before
    public void before() {
        binding = new Binding();
        dynamics = new Dynamics();
        dynamics.setBinding(binding);

        // Initialize Parallel Java
        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
        }

    }

    @Test
    public void testDynamicsHelp() {
        // Set-up the input arguments for the script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();
    }

    @Test
    public void testDynamicsNVE() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "VelocityVerlet", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that energy is conserved at the end of the dynamics trajectory.
        assertEquals(info + " End total energy for NVE test", startingTotalEnergy, molDyn.getTotalEnergy(), energyTolerance);
    }

}
