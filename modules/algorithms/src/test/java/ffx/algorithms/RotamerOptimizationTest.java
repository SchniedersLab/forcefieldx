package ffx.algorithms;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsUtils;

import java.io.File;
import java.util.logging.Logger;

public class RotamerOptimizationTest {

    private static final Logger logger = Logger.getLogger(RotamerOptimizationTest.class.getName());

    String filename;
    File structure;
    MolecularAssembly molecularAssembly;
    ForceFieldEnergy forceFieldEnergy;

    public RotamerOptimizationTest(){
        filename = "ffx/algorithms/structures/trpcage.pdb";
    }

    public void load() {
        /**
         * Load the test system.
         */
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());
        PotentialsUtils potentialUtils = new PotentialsUtils();
        molecularAssembly = potentialUtils.openQuietly(structure.getAbsolutePath());
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    }

    @org.junit.Test
    public void testSelfEnergyElimination() {
        load();
        /**
        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(false);
        rotamerOptimization.setUseGoldstein(true);
        rotamerOptimization.setResidues(1, 20);
        rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
         */
    }

}
