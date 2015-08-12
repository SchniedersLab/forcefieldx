/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFileOpener;
import java.io.File;
import java.util.ArrayList;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 *
 * @author slucore
 */
public class XeonProtonate {
    
    private static final Logger logger = Logger.getLogger(XeonProtonate.class.getName());
    private static File structure;
    private static MolecularAssembly molecularAssembly;
    private static ForceFieldEnergy forceFieldEnergy;
    private static int nAtoms;
    
    private static void main(String args[]) {
        System.setProperty("polarization", "none");
        System.setProperty("gkterm", "true");
        System.setProperty("cav_model", "born_solv");
        
        /**
         * Load the test system.
         */
        String testFile = System.getProperty("testFile");
        if (testFile != null) {
            structure = new File(testFile);
        } else {
            structure = new File("1omu.pdb");
        }
        
        PotentialsFileOpener opener = new PotentialsFileOpener(structure);
        opener.run();
        molecularAssembly = opener.getAssembly();
        nAtoms = molecularAssembly.getAtomArray().length;
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        
        logger.info(String.format("\n Running protonate: "));
        // default CLI options for 1omu.pdb
        int nSteps = 1000000;
        double timeStep = 1.0;
        double printInterval = 0.01;
        double saveInterval = 1.0;
        double temperature = 298.15;
        Thermostats thermostat = null;
        Integrators integrator = null;
        boolean initVelocities = true;
        double restartFrequency = 1000;
        String fileType = "PDB";
        int mcStepFrequency = 100;
        int rotamerStepFrequency = 995;
        double pH = 4.0;
        
        CompositeConfiguration properties = new CompositeConfiguration();

        // Single-residue titration option.
        ArrayList<String> resList = new ArrayList<>();
        resList.add("P7");
        resList.add("P10");
        resList.add("P13");
        resList.add("P19");
        resList.add("P27");
        resList.add("P29");
        resList.add("P34");
        resList.add("P55");

        // create the MD object
        MolecularDynamics molDyn = new MolecularDynamics(molecularAssembly, molecularAssembly.getPotentialEnergy(), properties, null, thermostat, integrator);
        molDyn.setFileType(fileType);
        molDyn.setRestartFrequency(restartFrequency);
        // create the Monte-Carlo listener and connect it to the MD
        Protonate mcProt = new Protonate(molecularAssembly, mcStepFrequency, rotamerStepFrequency, pH, molDyn.getThermostat());
        molDyn.addMCListener(mcProt);
        mcProt.addMolDyn(molDyn);
        // set residues to be titrated
        mcProt.chooseResID(resList);
        // finalize the Multi-Residue machinery
        mcProt.readyUp();
        // and away we go!
        molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, null);
    }
}
