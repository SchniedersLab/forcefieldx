package ffx.algorithms;

import edu.rit.pj.Comm;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.utils.PotentialsUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

public class RotamerOptimizationTest {

    private static final Logger logger = Logger.getLogger(RotamerOptimizationTest.class.getName());

    String filename;
    String restartName;
    File structure;
    File restartFile;
    MolecularAssembly molecularAssembly;
    ForceFieldEnergy forceFieldEnergy;

    public RotamerOptimizationTest() {
        filename = "ffx/algorithms/structures/trpcage.pdb";
        restartName = "ffx/algorithms/structures/trpcage.restart";
    }

    public void load() {
        /**
         * Load the test system.
         */
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());

        restartFile = new File(cl.getResource(restartName).getPath());

        PotentialsUtils potentialUtils = new PotentialsUtils();
        molecularAssembly = potentialUtils.openQuietly(structure.getAbsolutePath());
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    }

    @org.junit.Test
    public void testSelfEnergyElimination() {
        // Load the test system.
        load();

        // Initialize Parallel Java
        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
            message = String.format(" Skipping rotamer optimization test.");
            logger.log(Level.WARNING, message, e.toString());
            return;
        }

        // Run the optimization.
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(true);

        int counter = 1;
        int allStartResID = 1;
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = molecularAssembly.getChains();
        int nPolymers = polymers.length;
        for (int p = 0; p < nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            int nResidues = residues.size();
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers(rLib);
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    }
                    if (counter >= allStartResID) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }

        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(false);
        rotamerOptimization.setUseGoldstein(true);
        rotamerOptimization.setPruning(0);
        rotamerOptimization.setEnergyRestartFile(restartFile);
        rotamerOptimization.setResidues(residueList);
        rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
    }

}
