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
public class DynamicsStochasticTest {

    private String info;
    private String filename;
    private String restartFile;
    private double endKineticEnergy;
    private double kineticEnergyTolerance = 2.0;
    private double endPotentialEnergy;
    private double potentialEnergyTolerance = 2.0;
    private double endTotalEnergy;
    private double totalEnergyTolerance = 2.0;

    private Binding binding;
    private Dynamics dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsStochasticTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "Acetamide Peptide Restart and Stochastic Random Seed",  // info
                "ffx/algorithms/structures/acetamide_res_stoch.xyz",     // filename
                "ffx/algorithms/structures/acetamide_res_stoch.dyn",     // restartFile
                6.7599,                                                  // endKineticEnergy
                4.2177,                                                  // endPotentialEnergy
                10.9776                                                  // endTotalEnergy
            }
        });
    }

    public DynamicsStochasticTest(String info, String filename, String restartFile, double endKineticEnergy,
            double endPotentialEnergy, double endTotalEnergy) {

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
    public void testDynamicsStochasticRandomSeed() {
        
        //DynamicsStochasticTest.class.getClassLoader().setClassAssertionStatus(DynamicsStochasticTest.class.getName(), false);
        
        //assumeTrue(assertionsDisabled());

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);

        //Evaluate the script.
        dynamics.run();

        MolecularDynamics molDyn = dynamics.getMolecularDynamics();

        // Assert that the end energies are the same meaning that the Stochastic integrator works as intended.
        assertEquals(info + " End total energy for stochastic random seed test", endTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }
    
    public static boolean assertionsDisabled(){
        
        return !DynamicsStochasticTest.class.desiredAssertionStatus();
    }

}
