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

import ffx.algorithms.groovy.DynamicsOpenMM;

import groovy.lang.Binding;

/**
 *
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsOpenMMStochasticTest {

    private String info;
    private String filename;
    private String restartFile;
    private double endKineticEnergy;
    private double kineticEnergyTolerance = 5.0;
    private double endPotentialEnergy;
    private double potentialEnergyTolerance = 5.0;
    private double endTotalEnergy;
    private double totalEnergyTolerance = 5.0;
    private boolean ffxOpenMM;

    private Binding binding;
    private DynamicsOpenMM dynamics;

    private static final Logger logger = Logger.getLogger(DynamicsOpenMMStochasticTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "System OpenMM Stochastic", // info
                "ffx/algorithms/structues/waterbox_eq.xyz", // filename
                "ffx/algorithms/structues/waterbox_eq.dyn", // restartFile
                11796.9508, // endKineticEnergy
                -36782.1559, // endPotentialEnergy
                -24985.2051 // endTotalEnergy
            }
        });
    }
    
    public DynamicsOpenMMStochasticTest(String info, String filename, String restartFile, double endKineticEnergy, double endPotentialEnergy,
            double endTotalEnergy){
        
        this.info = info;
        this.filename = filename;
        this.restartFile = restartFile;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
        
        ffxOpenMM = System.getProperty("ffx.openMM","false").equalsIgnoreCase("true");
    }
    
    @Before
    public void before(){
        binding = new Binding();
        dynamics = new DynamicsOpenMM();
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
    public void testDynamicsOpenMMStochastic(){
        
        if(!ffxOpenMM){
            return;
        }
        
        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-z", "1", "-t", "298.15", "-i", "Stochastic", "-b", "Adiabatic", "-r", "0.001", "src/main/java/" + filename};
        binding.setVariable("args", args);
        
        // Evaluate script
        dynamics.run();
        
        MolecularDynamicsOpenMM molDynOpenMM = dynamics.getMolecularDynamics();
        
        // Assert that the end energies are within the threshold for the dynamics trajectory.
        assertEquals(info + "End kinetic energy for OpenMM Langevin(Stochastic) integrator", endKineticEnergy, molDynOpenMM.getKineticEnergy(), kineticEnergyTolerance);
        
        assertEquals(info + "End potential energy for OpenMM Langevin(Stochastic) integrator", endPotentialEnergy, molDynOpenMM.getPotentialEnergy(), potentialEnergyTolerance);
        
        assertEquals(info + "End total energy for OpenMM Langevin(Stochastic) integrator", endTotalEnergy, molDynOpenMM.getTotalEnergy(), totalEnergyTolerance);
    }

}
