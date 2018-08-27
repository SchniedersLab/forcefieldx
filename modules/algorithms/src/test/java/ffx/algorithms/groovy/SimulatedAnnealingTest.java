/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms.groovy;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;

import ffx.algorithms.SimulatedAnnealing;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * @author hbernabe
 */
@RunWith(Parameterized.class)
public class SimulatedAnnealingTest extends ffx.algorithms.PJDependentTest {

    private String info;
    private String filename;
    private double endKineticEnergy;
    private double endPotentialEnergy;
    private double endTotalEnergy;
    private double endTemperature;
    private double tolerance = 1.0;

    private Binding binding;
    private Anneal anneal;

    private static final Logger logger = Logger.getLogger(SimulatedAnnealingTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "Acetamide with Stochastic integrator for Simulated Annealing Test", //info
                        "ffx/algorithms/structures/acetamide_annealing.xyz", //filename
                        3.2227, // endKineticEnergy
                        -29.5661, // endPotentialEnergy
                        -26.3434, // endTotalEnergy
                        120.13 // endTemperature
                }

        });
    }

    public SimulatedAnnealingTest(String info, String filename, double endKineticEnergy, double endPotentialEnergy, double endTotalEnergy,
                                  double endTemperature) {
        this.info = info;
        this.filename = filename;
        this.endKineticEnergy = endKineticEnergy;
        this.endPotentialEnergy = endPotentialEnergy;
        this.endTotalEnergy = endTotalEnergy;
        this.endTemperature = endTemperature;
    }

    @Before
    public void before() {
        binding = new Binding();
        anneal = new Anneal();
        anneal.setBinding(binding);
    }

    @After
    public void after() {
        anneal.destroyPotentials();
        System.gc();
    }

    @Test
    public void testSimulatedAnnealing() {

        // Set-up the input arguments for the script.
        String[] args = {"-n", "50", "-i", "Stochastic", "-r", "0.001", "--tl", "100", "--tu", "1000", "-W", "5", "src/main/java/" + filename};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("SaveAsPDB");
            anneal.setBaseDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        anneal.run();

        SimulatedAnnealing simulatedAnnealing = anneal.getAnnealing();

        // Assert that end energies are within the tolerance for the dynamics trajectory
        assertEquals(info + " Final kinetic energy", endKineticEnergy, simulatedAnnealing.getKineticEnergy(), tolerance);
        assertEquals(info + " Final potential energy", endPotentialEnergy, simulatedAnnealing.getPotentialEnergy(), tolerance);
        assertEquals(info + " Final total energy", endTotalEnergy, simulatedAnnealing.getTotalEnergy(), tolerance);
        assertEquals(info + " Final temperature", endTemperature, simulatedAnnealing.getTemperature(), tolerance);

        // Delate all created space grouop directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by SaveAsPDB.");
        }
    }

}
