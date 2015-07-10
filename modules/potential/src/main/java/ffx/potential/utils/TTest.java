/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.utils;

import org.junit.Test;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.nonbonded.ParticleMeshEwald;
import java.io.File;
import java.util.logging.Logger;
import static org.junit.Assert.fail;
import org.junit.Before;

/**
 *
 * @author slucore
 */
public class TTest {
    
    private static final Logger logger = Logger.getLogger(TTest.class.getName());
    
    private File structure;
    private MolecularAssembly molecularAssembly;
    private final double tolerance = 1.0e-3;
    private final double gradientTolerance = 1.0e-4;
    private ParticleMeshEwald.Polarization polarization;
    private ForceFieldEnergy forceFieldEnergy;
    private int nAtoms;
    
    public TTest() {}
    
    @Before
    public void setup() {
        String filename = System.getProperty("testFile");
        if (filename == null) {
            logger.warning("Couldn't find file: " + filename);
        }
        
        /**
         * Load the test system.
         */
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(filename);        
        PotentialsFileOpener opener = new PotentialsFileOpener(structure);
        opener.run();
        
        molecularAssembly = opener.getAssembly();
        nAtoms = molecularAssembly.getAtomArray().length;
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    }
    
    @Test
    public void testTimer() {
        logger.info(String.format("\n N-Body Test: "));
        int nEvals = 10;
        boolean gradient = true;
        boolean print = true;
        
        ForceFieldEnergy energy = molecularAssembly.getPotentialEnergy();

        long minTime = Long.MAX_VALUE;
        double sumTime2 = 0.0;
        int halfnEvals = (nEvals % 2 == 1) ? (nEvals/2) : (nEvals/2) - 1; // Halfway point
        for (int i=0; i<nEvals; i++) {
            long time = -System.nanoTime();
            energy.energy(gradient, print);
            time += System.nanoTime();
            minTime = time < minTime ? time : minTime;
            if (i >= (int) (nEvals/2)) {
                double time2 = time * 1.0E-9;
                sumTime2 += (time2*time2);
            }
        }
        ++halfnEvals;
        double rmsTime = Math.sqrt(sumTime2/halfnEvals);
        logger.info(String.format(" Minimum time: %14.5f (sec)", minTime * 1.0E-9));
        logger.info(String.format(" RMS time (latter half): %14.5f (sec)", rmsTime));
    }
}
