package ffx.algorithms.optimize;

import edu.rit.pj.Comm;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.PhReplicaExchange;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.misc.AlgorithmsTest;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.utils.PotentialsUtils;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.*;

public class PhRepexTest extends AlgorithmsTest {

    @Test
    public void testPhRepexRestarts() {
        System.setProperty("ffx.log", "severe");
        String structure = getResourcePath("PhRepexTestFiles/validRestart/LYS_penta.pdb");
        PotentialsUtils potentialUtils = new PotentialsUtils();
        MolecularAssembly molecularAssembly = potentialUtils.openQuietly(structure);
        ForceFieldEnergy potential = molecularAssembly.getPotentialEnergy();
        ExtendedSystem esvSystem = new ExtendedSystem(molecularAssembly, 7, null);
        potential.attachExtendedSystem(esvSystem);
        double[] pHLadder = new double[]{8.4, 11.4, 11.9};
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);
        File structureFile = new File(structure);
        int rank = 0;
        final String newMolAssemblyFile = structureFile.getParent() + File.separator + rank +
                File.separator + structureFile.getName();
        molecularAssembly.setFile(new File(newMolAssemblyFile));
        AlgorithmListener algorithmListener = new AlgorithmListener() {
            @Override
            public boolean algorithmUpdate(MolecularAssembly active) {
                return true;
            }
        };
        MolecularDynamics molecularDynamics = MolecularDynamics.dynamicsFactory(molecularAssembly, potential, algorithmListener,
                ThermostatEnum.ADIABATIC, IntegratorEnum.STOCHASTIC);
        try {
            // Tests restart feature
            PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 3);
            // Restarts based solely on the structure file, so thats all we need to change between tests
            structureFile = new File(getResourcePath("PhRepexTestFiles/validAfterPrimaryPhFail/LYS_penta.pdb"));
            pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 3);
            structureFile = new File(getResourcePath("PhRepexTestFiles/validAfterPrimaryCountFail/LYS_penta.pdb"));
            pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 3);
        } catch (Exception e) {
            // Fail if throws exception
            fail("Unexpected exception thrown from " + structureFile + ": " + e.getMessage());
        }
        try{
            structureFile = new File(getResourcePath("PhRepexTestFiles/invalidAfterPhFail/LYS_penta.pdb"));
            PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 3);
            fail("Expected exception not thrown from " + structureFile);
        } catch (Exception ignored) {}
        try{
            structureFile = new File(getResourcePath("PhRepexTestFiles/invalidAfterCountFail/LYS_penta.pdb"));
            PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 3);
            fail("Expected exception not thrown from " + structureFile);
        } catch (Exception ignored) {}
    }

    @Test
    public void testPhRepexExchangeValidity() {
        System.setProperty("ffx.log", "severe");
        String structure = getResourcePath("PhRepexTestFiles/validRestart/LYS_penta.pdb");
        PotentialsUtils potentialUtils = new PotentialsUtils();
        MolecularAssembly molecularAssembly = potentialUtils.openQuietly(structure);
        ForceFieldEnergy potential = molecularAssembly.getPotentialEnergy();
        ExtendedSystem esvSystem = new ExtendedSystem(molecularAssembly, 7, null);
        potential.attachExtendedSystem(esvSystem);
        double[] pHLadder = new double[]{8.4, 8.9, 9.4, 9.9, 10.4, 10.9, 11.4, 11.9};
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);
        File structureFile = new File(structure);
        int rank = 0;
        final String newMolAssemblyFile = structureFile.getParent() + File.separator + rank +
                File.separator + structureFile.getName();
        molecularAssembly.setFile(new File(newMolAssemblyFile));
        AlgorithmListener algorithmListener = new AlgorithmListener() {
            @Override
            public boolean algorithmUpdate(MolecularAssembly active) {
                return true;
            }
        };
        MolecularDynamics molecularDynamics = MolecularDynamics.dynamicsFactory(molecularAssembly, potential, algorithmListener,
                ThermostatEnum.ADIABATIC, IntegratorEnum.STOCHASTIC);
        try {
            PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, 7, pHLadder,
                    298, esvSystem, x, 8);
            int[] pHMap = pHReplicaExchange.setTestingParametersAndExchangeOnce();
            int[] expected = new int[]{3, 0, 1, 2, 7, 4, 5, 6,
                    1, 2, 3, 0, 5, 6, 7, 4};
            assertArrayEquals(expected, pHMap);
        } catch(Exception e) {
            fail("Unexpected exception thrown from " + structureFile + ": " + e.getMessage());
        }
    }
}
