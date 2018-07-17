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

import edu.rit.pj.Comm;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assume.assumeTrue;

import groovy.lang.Binding;

import ffx.algorithms.groovy.Dynamics;

/**
 *
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
    private boolean testNVE;
    private boolean testNVT;
    private boolean testRestart;
    private boolean testStochasticRandomSeed;
    
    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsNVETest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide Peptide NVE",                                        // info
                "ffx/algorithms/structures/acetamide_NVE.xyz",                  // filename
                "ffx/algorithms/structures/acetamide_NVE.dyn",                  // restartFile
                6.8672,                                                         // startingTotalEnergy
                298.15,                                                         // startingTemp
                5.2671,                                                         // endKineticEnergy
                1.4667,                                                         // endPotentialEnergy
                6.7338,                                                         // endTotalEnergy
                true,                                                           // testNVE
                false,                                                          // testNVT
                false,                                                          // testRestart
                false                                                           // testStochasticRandomSeed
            }
        });
    }

    public DynamicsNVETest(String info, String filename, String restartFile, double startingTotalEnergy, double startingTemp, double endKineticEnergy,
            double endPotentialEnergy, double endTotalEnergy, boolean testNVE, boolean testNVT, boolean testRestart, boolean testStochasticRandomSeed) {
        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.startingTotalEnergy = startingTotalEnergy;
        this.startingTemp = startingTemp;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
        this.testNVE = testNVE;
        this.testNVT = testNVT;
        this.testRestart = testRestart;
        this.testStochasticRandomSeed = testStochasticRandomSeed;
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

    /*
    //@Test
    public void testDynamicsHelp() {
                // Set-up the input arguments for the script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        dynamics.run();
    }*/

    @Test
    public void testDynamicsNVE() {
        
        //DynamicsNVETest.class.getClassLoader().setClassAssertionStatus(DynamicsNVETest.class.getName(), false);
        
        //assumeTrue(assertionsDisabled());
        
        logger.info("In the NVE test, before passing in arguments");
        
        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "VelocityVerlet", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);
        
        logger.info("Arguments passed in to NVE test dynamics object");

        // Evaluate the script.
        dynamics.run();
        
        logger.info("Dynamics run complete");
        
        MolecularDynamics molDyn = dynamics.getMolecularDynamics();
        
        logger.info(String.format(" Starting total energy %f", startingTotalEnergy));
        
        logger.info(String.format(" End total energy %f", molDyn.getTotalEnergy()));
        
        // Assert that energy is conserved at the end of the dynamics trajectory.
        assertEquals(info + " End total energy for NVE test", startingTotalEnergy, molDyn.getTotalEnergy(), energyTolerance);
    }
    
    public static boolean assertionsDisabled(){
        
        return !DynamicsNVETest.class.desiredAssertionStatus();
    }
    
    /*
    @Test
    public void testDynamicsRestart(){
        
        if(!testRestart){
            return;
        }
        
        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Stochastic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);
        
        // Evaluate the script.
        dynamics.run();
        
        MolecularDynamics molDyn = dynamics.getMolecularDynamics();
        
        // Asssert that the end kinetic, potential and total energies are the same meaning that the simulation started from the appropriate .dyn file.        
        assertEquals(info + " End kinetic energy for restart test", endKineticEnergy, molDyn.getKineticEnergy(), kineticEnergyTolerance);
        
        assertEquals(info + " End potential energy for restart test", endPotentialEnergy, molDyn.getPotentialEnergy(), potentialEnergyTolerance);
        
        assertEquals(info + " End total energy for restart test", endTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }*/
}
