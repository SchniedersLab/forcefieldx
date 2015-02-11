/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms;

import ffx.potential.ForceFieldEnergy;
import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.ForceField;
import java.util.List;
import java.util.Objects;

/**
 * @author S. LuCore
 */
public class Protonate implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(Protonate.class.getName());
    /**
     * The MD thermostat.
     */
    private final Thermostat thermostat;
    /**
     * Boltzmann's constant is kcal/mol/Kelvin.
     */
    private static final double boltzmann = 0.0019872041;
    /**
     * Energy of the system at initialization.
     */
    private final double systemReferenceEnergy;
    /**
     * Simulation pH.
     */
    private final double pH;
    /**
     * The current MD step.
     */
    private int stepCount = 0;
    /**
     * Number of simulation steps between MC move attempts.
     */
    private static int mcStepFrequency;
    /**
     * Number of accepted MD moves.
     */
    private int numMovesAccepted;
    /**
     * MultiResidue forms of entities from titratableResidues; ready to be (de-/)protonated.
     */
    private ArrayList<MultiResidue> titratingResidues = new ArrayList<>();
    private Random rng = new Random();
    private final ForceField forceField;
    private final ForceFieldEnergy forceFieldEnergy;

    /**
     * Construct a Monte-Carlo protonation state switching mechanism.
     *
     * @param molAss the molecular assembly
     * @param mcStepFrequency number of MD steps between switch attempts
     * @param pH the simulation pH
     * @param thermostat the MD thermostat
     */
    Protonate(MolecularAssembly molAss, int mcStepFrequency, double pH, Thermostat thermostat) {
        //initialize stepcount and the number of accepted moves
        numMovesAccepted = 0;

        this.forceField = molAss.getForceField();
        this.forceFieldEnergy = molAss.getPotentialEnergy();
        this.mcStepFrequency = mcStepFrequency;
        this.pH = pH;
        this.thermostat = thermostat;
        systemReferenceEnergy = molAss.getPotentialEnergy().getTotalEnergy();
        
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Running Protonate:\n"));
        sb.append(String.format("     mcStepFrequency: %4d\n", mcStepFrequency));
        sb.append(String.format("     referenceEnergy: %7.2f\n", systemReferenceEnergy));
        sb.append(String.format("     system pH:       %7.2f", pH));
        logger.info(sb.toString());

        // Identify titratable residues.
        List<Residue> titratableResidues = new ArrayList<>();
        Polymer polymers[] = molAss.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                if (isTitratable(residues.get(j).getName())) {
                    titratableResidues.add(residues.get(j));
                    logger.info(String.format(" Titratable: %s", residues.get(j)));
                }
            }
        }
        
        // Create MultiResidue objects to wrap titratables.
        for (Residue res : titratableResidues) {
            MultiResidue multiRes = new MultiResidue(res, forceField, forceFieldEnergy);
            Polymer polymer = findResiduePolymer(res, molAss);
            polymer.addMultiResidue(multiRes);
            String protFormName = Titratable.valueOf(res.getName()).protForm.toString();
            String deprotFormName = Titratable.valueOf(res.getName()).deprotForm.toString();
            int resNumber = res.getResidueNumber();
            ResidueType resType = res.getResidueType();
            if (!res.getName().equalsIgnoreCase(protFormName)) {
                multiRes.addResidue(new Residue(protFormName, resNumber, resType));
            } else {
                multiRes.addResidue(new Residue(deprotFormName, resNumber, resType));
            }
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(protFormName));
            forceFieldEnergy.reInit();
            titratingResidues.add(multiRes);
            logger.info(String.format(" Titrating: %s", multiRes));
        }
        
        forceFieldEnergy.reInit();
    }

    /**
     * True if passed residue name has multiple protonation states.
     *
     * @param residueName
     * @return if this residue has multiple protonation states
     */
    private boolean isTitratable(String residueName) {
        for (Titratable titrName : Titratable.values()) {
            if (residueName.equalsIgnoreCase(titrName.toString())) {
                return true;
            }
        }
        return false;
    }

    @Override
    public boolean mcUpdate(MolecularAssembly molAss) {
        stepCount++;
        if (stepCount % mcStepFrequency != 0) {
            return false;
        }        
        double referenceEnergy = molAss.getPotentialEnergy().getTotalEnergy();
        
        // Randomly choose a target titratable residue to attempt protonation switch.
        int random = rng.nextInt(titratingResidues.size());
        MultiResidue targetResidue = titratingResidues.get(random);

        // Switch titration state for chosen residue.
        switchProtonationState(targetResidue);

        String name = (titratingResidues.get(random)).getActive().getName();
        double pKaref = Titratable.valueOf(name).pKa;
        double dG_ref = Titratable.valueOf(name).refEnergy;
        double temperature = thermostat.getCurrentTemperature();
        double kT = boltzmann * temperature;
        double dG_elec = molAss.getPotentialEnergy().getTotalEnergy() - referenceEnergy;

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the
         * reference compound
         */
        double dG_MC = kT * (pH - pKaref) * Math.log(10) + dG_elec - dG_ref;
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC protonation step:\n"));
        sb.append(String.format("     pKaref:  %7.2f\n", pKaref));
        sb.append(String.format("     dG_ref:  %7.2f\n", dG_ref));
        sb.append(String.format("     dG_elec: %9.4f\n", dG_elec));
        sb.append(String.format("     dG_MC:   %9.4f\n", dG_MC));
        sb.append(String.format("     -----\n"));

        // Test Monte-Carlo criterion.
        if (dG_MC < 0) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        double boltzmann = exp(-dG_MC / kT);
        double metropolis = random();
        sb.append(String.format("     boltzmann:  %9.4f\n", boltzmann));
        sb.append(String.format("     metropolis: %9.4f\n", metropolis));
        if (metropolis < boltzmann) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        
        sb.append(String.format("     Denied."));
        logger.info(sb.toString());
        // Undo titration state change if criterion was not accepted.
        switchProtonationState(targetResidue);
        return false;
    }

    /**
     * Switch the protonation state of target residue and reinitialize FF.
     *
     * @param residue
     */
    private void switchProtonationState(MultiResidue multiRes) {
        String protFormName = Titratable.valueOf(multiRes.getActive().getName()).protForm.toString();
        String deprotFormName = Titratable.valueOf(multiRes.getActive().getName()).deprotForm.toString();
        if (multiRes.getName().equalsIgnoreCase(protFormName)) {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(deprotFormName));
        } else {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(protFormName));
        }
        forceFieldEnergy.reInit();            
        return;
    }

    /**
     * Get the current MC acceptance rate.
     *
     * @return the acceptance rate.
     */
    @Override
    public double getAcceptanceRate() {
        // Intentional integer division.
        int numTries = stepCount / mcStepFrequency;
        return (double) numMovesAccepted / numTries;
    }
    
    /**
     * Locate to which Polymer in the MolecularAssembly a given Residue belongs.
     * @param res
     * @param molAss
     * @return 
     */
    private Polymer findResiduePolymer(Residue res, MolecularAssembly molAss) {
        if (res.getChainID() == null) {
            logger.severe("No chain ID for residue " + res);
        }
        Polymer polymers[] = molAss.getChains();
        Polymer location = null;
        for (Polymer p : polymers) {
            if (p.getChainID() == res.getChainID()) {
                location = p;
            }
        }
        if (location == null) {
            logger.severe("Couldn't find polymer for residue " + res);
        }
        return location;
    }

    /**
     * Constant values for intrinsic pKa and reference energy of deprotonation.
     */
    public enum Titratable {

//        ARG(12.48, 1.00, AminoAcid3.ARD),
        // Standard Forms
        ASP(4.00, 1.00, AminoAcid3.ASH, AminoAcid3.ASP),
        GLU(4.25, 1.00, AminoAcid3.GLH, AminoAcid3.GLU),
        CYS(8.18, 1.00, AminoAcid3.CYS, AminoAcid3.CYD),
        HIS(6.00, 1.00, AminoAcid3.HIS, AminoAcid3.HID),
        LYS(10.53, 1.00, AminoAcid3.LYS, AminoAcid3.LYD),
        TYR(10.07, 1.00, AminoAcid3.TYR, AminoAcid3.TYD),
        // Protonated Forms
        ASH(4.00, 1.00, AminoAcid3.ASH, AminoAcid3.ASP),
        GLH(4.25, 1.00, AminoAcid3.GLH, AminoAcid3.GLU),
        // Deprotonated Forms
        CYD(8.18, 1.00, AminoAcid3.CYS, AminoAcid3.CYD),
        HID(6.00, 1.00, AminoAcid3.HIS, AminoAcid3.HID),
        LYD(10.53, 1.00, AminoAcid3.LYS, AminoAcid3.LYD),
        TYD(10.07, 1.00, AminoAcid3.TYR, AminoAcid3.TYD);

        public final double pKa;
        public final double refEnergy;
        public final AminoAcid3 protForm;
        public final AminoAcid3 deprotForm;

        Titratable(double pKa, double refEnergy, AminoAcid3 protForm, AminoAcid3 deprotForm) {
            this.pKa = pKa;
            this.refEnergy = refEnergy;
            this.protForm = protForm;
            this.deprotForm = deprotForm;
        }
    };
}
