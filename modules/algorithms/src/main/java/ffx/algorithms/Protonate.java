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
import static ffx.potential.bonded.AminoAcidUtils.cbType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedUtils;
import static ffx.potential.bonded.BondedUtils.buildHydrogen;
import static ffx.potential.bonded.BondedUtils.intxyz;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import java.io.File;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import org.apache.commons.io.FilenameUtils;

/**
 * @author S. LuCore
 */
public class Protonate implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(Protonate.class.getName());
    /**
     * The MoleularAssembly.
     */
    private final MolecularAssembly molAss;
    /**
     * The MolecularDynamics object controlling the simulation.
     */
    private MolecularDynamics molDyn;
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
    private int mcStepFrequency;
    /**
     * Number of simulation steps between rotamer move attempts.
     */
    private int rotamerStepFrequency = 0;
    /**
     * Number of accepted MD moves.
     */
    private int numMovesAccepted;
    /**
     * Residues selected by user.
     */
    private List<Residue> chosenResidues = new ArrayList<>();
    /**
     * MultiResidue forms of entities from chosenResidues; ready to be (de-/)protonated.
     */
    private List<MultiResidue> titratingResidues = new ArrayList<>();
    /**
     * Everyone's favorite.
     */
    private Random rng = new Random();
    /**
     * The forcefield being used.
     * Needed by MultiResidue constructor.
     */
    private final ForceField forceField;
    /**
     * The ForceFieldEnergy object being used by MD.
     * Needed by MultiResidue constructor and for reinitializing after a chemical change.
     */
    private final ForceFieldEnergy forceFieldEnergy;
    /**
     * Enum type to specify global override of MC acceptance criteria.
     */
    private MCOverride mcTitrationOverride = MCOverride.NONE;
    /**
     * Writes .s-[num] and .f-[num] files representing before/after MC move structures.
     * Note: The 'after' snapshot represents the change that was PROPOSED, regardless of accept/reject.
     */
    private SnapshotsType snapshotsType = SnapshotsType.NONE;
    /**
     * True once the titratingResidues list is ready.
     */
    private boolean finalized = false;
    /**
     * Debug mode: all reference energies set to zero.
     */
    private boolean zeroReferenceEnergies = false;
    
    /**
     * Construct a Monte-Carlo protonation state switching mechanism.
     *
     * @param molAss the molecular assembly
     * @param mcStepFrequency number of MD steps between switch attempts
     * @param pH the simulation pH
     * @param thermostat the MD thermostat
     */
    Protonate(MolecularAssembly molAss, int mcStepFrequency, int rotamerStepFrequency, double pH, Thermostat thermostat) {
        // process system flags
        String zeroReferenceEnergies = System.getProperty("zeroReferences");
        if (zeroReferenceEnergies != null) {
            this.zeroReferenceEnergies = Boolean.parseBoolean(zeroReferenceEnergies);
            if (this.zeroReferenceEnergies) {
                logger.info(" DEBUG: Zero-ing all reference energies.");
            }
        }
        String debugLogLevel = System.getProperty("debug");
        if (debugLogLevel != null) {
            this.debugLogLevel = Integer.parseInt(debugLogLevel);
        }
        String overrideFlag = System.getProperty("MCoverride");
        if (overrideFlag != null && overrideFlag.equalsIgnoreCase("accept")) {
            logger.info(" OVERRIDE: Accepting all MC moves.");
            mcTitrationOverride = MCOverride.ACCEPT;
        }
        if (overrideFlag != null && overrideFlag.equalsIgnoreCase("reject")) {
            logger.info(" OVERRIDE: Rejecting all MC moves.");
            mcTitrationOverride = MCOverride.REJECT;
        }
        String beforeAfter = System.getProperty("MCsnapshots");
        if (beforeAfter != null) {
            if (beforeAfter.equalsIgnoreCase("true") || beforeAfter.equalsIgnoreCase("separate")) {
                logger.info(" DEBUG: Writing before-and-after MC snapshots.");
                snapshotsType = SnapshotsType.SEPARATE;
            } else if (beforeAfter.equalsIgnoreCase("interleave") || beforeAfter.equalsIgnoreCase("interleaved")) {
                logger.info(" DEBUG: Writing interleaved MC snapshots.");
                snapshotsType = SnapshotsType.INTERLEAVED;
            }
        }

        numMovesAccepted = 0;
        // Set the rotamer library in case we do rotamer MC moves.
        RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        RotamerLibrary.setUseOrigCoordsRotamer(false);

        this.molAss = molAss;
        this.forceField = molAss.getForceField();
        this.forceFieldEnergy = molAss.getPotentialEnergy();
        this.mcStepFrequency = (mcStepFrequency == 0) ? Integer.MAX_VALUE : mcStepFrequency;
        this.rotamerStepFrequency = (rotamerStepFrequency == 0) ? Integer.MAX_VALUE : rotamerStepFrequency;
        this.pH = pH;
        this.thermostat = thermostat;
        systemReferenceEnergy = molAss.getPotentialEnergy().getTotalEnergy();
        
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Running Protonate:\n"));
        sb.append(String.format("     mcStepFrequency: %4d\n", mcStepFrequency));
        sb.append(String.format("     referenceEnergy: %7.2f\n", systemReferenceEnergy));
        sb.append(String.format("     system pH:       %7.2f", pH));
        logger.info(sb.toString());

        forceFieldEnergy.reInit();
    }
    
    /**
     * Identify titratable residues and choose them all.
     */
    private void chooseAllTitratables() {
        chosenResidues = new ArrayList<>();
        Polymer polymers[] = molAss.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                if (isTitratable(residues.get(j).getName())) {
                    chosenResidues.add(residues.get(j));
                    // logger.info(String.format(" Titratable: %s", residues.get(j)));
                }
            }
        }
    }
    
    /**
     * Choose titratables with intrinsic pKa inside (pH-window,pH+window).
     * @param pH
     * @param window 
     */
    private void chooseTitratablesInWindow(double pH, double window) {
        chosenResidues = new ArrayList<>();
        Polymer polymers[] = molAss.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                Residue res = residues.get(j);
                if (isTitratable(res.getName())) {
                    double pKa = Titratable.valueOf(res.getName()).pKa;
                    if (pKa >= pH - window && pKa <= pH + window) {
                        chosenResidues.add(residues.get(j));
                        // logger.info(String.format(" Titratable: %s", residues.get(j)));
                    }
                }
            }
        }
    }
        
    private void chooseResID(ArrayList<String> crIDs) {
        Polymer[] polymers = molAss.getChains();
        int n = 0;
        for (String s : crIDs) {
            Character chainID = s.charAt(0);
            int i = Integer.parseInt(s.substring(1));
            for (Polymer p : polymers) {
                if (p.getChainID() == chainID) {
                    List<Residue> rs = p.getResidues();
                    for (Residue r : rs) {
                        if (r.getResidueNumber() == i) {
                            chosenResidues.add(r);
                            // logger.info(String.format(" Chosen: %s", r));
                        }
                    }
                }
            }
        }
    }
    
    private void chooseResID(char chain, int resID) {
        Polymer polymers[] = molAss.getChains();
        for (Polymer polymer : polymers) {
            if (polymer.getChainID() == chain) {
                ArrayList<Residue> residues = polymer.getResidues();
                for (Residue residue : residues) {
                    if (residue.getResidueNumber() == resID) {
                        chosenResidues.add(residue);
                        logger.info(String.format(" Chosen: %s", residue));
                    }
                }
            }
        }
    }
    
    /**
     * Must be called after all titratable residues have been chosen, but before beginning MD.
     */
    private void readyUp() {
        // Create MultiResidue objects to wrap titratables.
        for (Residue res : chosenResidues) {
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
            multiRes.setActiveResidue(res);
            forceFieldEnergy.reInit();
            titratingResidues.add(multiRes);
            logger.info(String.format(" Titrating: %s", multiRes));
        }
        
        finalized = true;
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

    private int snapshotIndex = 0;
    
    /**
     * The primary driver. Called by the MD engine at each dynamics step.
     */
    @Override
    public boolean mcUpdate(MolecularAssembly molAss) {
        if (!finalized) {
            logger.severe("Monte-Carlo protonation engine was not finalized!");
        }
        
        propagateInactiveResidues(titratingResidues);
        
        stepCount++;
        // Decide on the type of step to be taken.
        StepType stepType;
        if (stepCount % mcStepFrequency == 0 && stepCount % rotamerStepFrequency == 0) {
            stepType = StepType.COMBO;
        } else if (stepCount % mcStepFrequency == 0) {
            stepType = StepType.TITRATE;
        } else if (stepCount % rotamerStepFrequency == 0) {
            stepType = StepType.ROTAMER;
        } else {
            // Not yet time for an MC step, return to MD.
            return false;
        }
        
        // Randomly choose a target titratable residue to attempt protonation switch.
        int random = rng.nextInt(titratingResidues.size());
        MultiResidue targetMulti = titratingResidues.get(random);
        
        // Check whether rotamer moves are possible for the selected residue.
        if (RotamerLibrary.getRotamers(targetMulti.getActive()) == null 
                || RotamerLibrary.getRotamers(targetMulti.getActive()).length <= 1) {
            if (stepType == StepType.ROTAMER) {
                return false;
            } else if (stepType == StepType.COMBO) {
                stepType = StepType.TITRATE;
            }
        }
        
        // Perform the MC move.
        boolean accepted = false;
        switch (stepType) {
            case TITRATE:
                accepted = tryTitrationStep(targetMulti);
                break;
            case ROTAMER:
                accepted = tryRotamerStep(targetMulti);
                break;
            case COMBO:
                accepted = tryComboStep(targetMulti);
                break;
        }
        
        // Increment the shared snapshot counter.
        snapshotIndex++;
        return accepted;
    }

    /**
     * Perform a titration MC move.
     * @param targetMulti
     * @return accept/reject
     */
    private boolean tryTitrationStep(MultiResidue targetMulti) {
        // Record the pre-change electrostatic energy.
        double previousElectrostaticEnergy = currentElectrostaticEnergy();
        
        // Write the pre-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);

        // Switch titration state for chosen residue.
        String startName = targetMulti.toString();
        switchProtonationState(targetMulti);

        // Write the post-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);
        
        String endName = targetMulti.getActive().getName();
        TitrationType type = determineTitrationType(startName, endName);
        
        // Test the MC criterion for a titration step.
        double pKaref = Titratable.valueOf(endName).pKa;
        double dG_ref = Titratable.valueOf(endName).refEnergy;
        double temperature = thermostat.getCurrentTemperature();
        double kT = boltzmann * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;
        
        if (zeroReferenceEnergies) {
            dG_ref = 0.0;
        }

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the reference compound
         */
        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEPROT) {
            prefix = -prefix;
        }
        double postfix = dG_elec - dG_ref;
        double dG_MC = prefix + postfix;
        
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC protonation step:\n"));
        sb.append(String.format("     %s --> %s\n", startName, targetMulti.toString()));
        sb.append(String.format("     pKaref:  %7.2f\n", pKaref));
        sb.append(String.format("     dG_ref:  %7.2f\n", dG_ref));
        sb.append(String.format("     dG_elec: %16.8f\n", dG_elec));
        sb.append(String.format("     dG_MC:   %16.8f\n", dG_MC));
        sb.append(String.format("     -----\n"));

        // Test Monte-Carlo criterion.
        if (dG_MC < 0 && mcTitrationOverride != MCOverride.REJECT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        double criterion = exp(-dG_MC / kT);
        double metropolis = random();
        sb.append(String.format("     criterion:  %9.4f\n", criterion));
        sb.append(String.format("     rng:        %9.4f\n", metropolis));
        if ((metropolis < criterion && mcTitrationOverride != MCOverride.REJECT) || mcTitrationOverride == MCOverride.ACCEPT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        sb.append(String.format("     Denied."));
        logger.info(sb.toString());
        
        // Move was rejected, undo the titration state change.
        switchProtonationState(targetMulti);
        return false;
    }
    
    private boolean tryRotamerStep(MultiResidue targetMulti) {
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();
        
        // Write the before-step snapshot.
        writeSnapshot(true, StepType.ROTAMER, snapshotsType);

        // Save coordinates so we can return to them if move is rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        double[][] origCoordinates = new double[atoms.size()][];
        for (int i = 0; i < atoms.size(); i++) {
            Atom atomi = atoms.get(i);
            origCoordinates[i] = new double[atomi.getXYZ().length];
            atomi.getXYZ(origCoordinates[i]);
        }
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origCoordinates, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);
        // Select a new rotamer and swap to it.
        Rotamer rotamers[] = RotamerLibrary.getRotamers(residue);
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);

        // Write the post-rotamer change snapshot.
        writeSnapshot(false, StepType.ROTAMER, snapshotsType);

        // Check the MC criterion.
        double temperature = thermostat.getCurrentTemperature();
        double kT = boltzmann * temperature;
        double postTotalEnergy = currentTotalEnergy();
        double dG_tot = postTotalEnergy - previousTotalEnergy;
        double criterion = exp(-dG_tot / kT);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC rotamer step:\n"));
        sb.append(String.format("     prev:   %16.8f\n", previousTotalEnergy));
        sb.append(String.format("     post:   %16.8f\n", postTotalEnergy));
        sb.append(String.format("     dG_tot: %16.8f\n", dG_tot));
        sb.append(String.format("     -----\n"));

        // Automatic acceptance if energy change is favorable.
        if (dG_tot < 0) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            propagateInactiveResidues(titratingResidues);
            return true;
        } else {
            // Conditional acceptance if energy change is positive.
            double metropolis = random();
            sb.append(String.format("     criterion:  %9.4f\n", criterion));
            sb.append(String.format("     rng:        %9.4f\n", metropolis));
            if (metropolis < criterion) {
                sb.append(String.format("     Accepted!"));
                logger.info(sb.toString());
                numMovesAccepted++;
                propagateInactiveResidues(titratingResidues);
                return true;
            } else {
                // Move was denied.
                sb.append(String.format("     Denied."));
                logger.info(sb.toString());

                // Undo the rejected move.
                RotamerLibrary.applyRotamer(residue, origCoordsRotamer);
                return false;
            }
        }
    }
    
    private boolean tryComboStep(MultiResidue targetMulti) {
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();
        double previousElectrostaticEnergy = currentElectrostaticEnergy();
        
        // Write the pre-combo snapshot.
        writeSnapshot(true, StepType.COMBO, snapshotsType);

        // Change titration state.
        String startName = targetMulti.toString();
        switchProtonationState(targetMulti);
        
        // Change rotamer state, but first save coordinates so we can return to them if rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        double[][] origCoordinates = new double[atoms.size()][];
        for (int i = 0; i < atoms.size(); i++) {
            Atom atomi = atoms.get(i);
            origCoordinates[i] = new double[atomi.getXYZ().length];
            atomi.getXYZ(origCoordinates[i]);
        }
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origCoordinates, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);
        
        // Swap to the new rotamer.
        Rotamer rotamers[] = RotamerLibrary.getRotamers(residue);
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);
        
        // Write the post-combo snapshot.
        writeSnapshot(false, StepType.COMBO, snapshotsType);
        
        // Evaluate both MC criteria.
        String endName = targetMulti.getActive().getName();
        TitrationType type = determineTitrationType(startName, endName);
        
        // Evaluate the titration probability of the step.
        double pKaref = Titratable.valueOf(endName).pKa;
        double dG_ref = Titratable.valueOf(endName).refEnergy;
        double temperature = thermostat.getCurrentTemperature();
        double kT = boltzmann * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;
        
        if (zeroReferenceEnergies) {
            dG_ref = 0.0;
        }
        
        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEPROT) {
            prefix = -prefix;
        }
        double postfix = dG_elec - dG_ref;
        double dG_titr = prefix + postfix;
        double titrCriterion = exp(-dG_titr / kT);
        
        // Evaluate the rotamer probability of the step.
        double dG_rota = currentTotalEnergy() - previousTotalEnergy;
        double rotaCriterion = exp(-dG_rota / kT);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC combo step:\n"));
        sb.append(String.format("     dG_elec: %16.8f\n", dG_elec));
        sb.append(String.format("     dG_titr: %16.8f\n", dG_titr));
        sb.append(String.format("     dG_rota: %16.8f\n", dG_rota));
        sb.append(String.format("     -----\n"));
        
        // Test the combined probability of this move.
        // Automatic acceptance if both energy changes are favorable.
        if (dG_titr < 0 && dG_rota < 0 && mcTitrationOverride != MCOverride.REJECT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            propagateInactiveResidues(titratingResidues);
            return true;
        } else {
            // Conditionally accept based on combined probabilities.
            if (dG_titr < 0 || mcTitrationOverride == MCOverride.ACCEPT) {
                titrCriterion = 1.0;
            }
            if (dG_rota < 0) {
                rotaCriterion = 1.0;
            }
            if (mcTitrationOverride == MCOverride.REJECT) {
                titrCriterion = 0.0;
            }
            double metropolis = random();
            double comboCriterion = titrCriterion * rotaCriterion;
            sb.append(String.format("     titrCrit:   %9.4f\n", titrCriterion));
            sb.append(String.format("     rotaCrit:   %9.4f\n", rotaCriterion));
            sb.append(String.format("     criterion:  %9.4f\n", comboCriterion));
            sb.append(String.format("     rng:        %9.4f\n", metropolis));
            if (metropolis < comboCriterion) {
                sb.append(String.format("     Accepted!"));
                logger.info(sb.toString());
                numMovesAccepted++;
                propagateInactiveResidues(titratingResidues);
                return true;
            } else {
                // Move was denied.
                sb.append(String.format("     Denied."));
                logger.info(sb.toString());

                // Undo both pieces of the rejected move IN THE RIGHT ORDER.
                RotamerLibrary.applyRotamer(residue, origCoordsRotamer);
                switchProtonationState(targetMulti);                
                return false;
            }
        }
    }
    
    /**
     * Switch the protonation state of target residue and reinitialize FF.
     *
     * @param residue
     */
    private void switchProtonationState(MultiResidue multiRes) {
        /* OLD WAY TO CALL MD_REINIT() == crazy slow
        Atom oldAtomArray[] = new Atom[molAss.getAtomArray().length];
        System.arraycopy(oldAtomArray, 0, molAss.getAtomArray(), 0, oldAtomArray.length);
        for (int i = 0; i < oldAtomArray.length; i++) {
            oldAtomArray[i] = molAss.getAtomArray()[i];
        } */
        
        List<Atom> oldAtoms = multiRes.getActive().getAtomList();        
        String protFormName = Titratable.valueOf(multiRes.getActive().getName()).protForm.toString();
        String deprotFormName = Titratable.valueOf(multiRes.getActive().getName()).deprotForm.toString();
        if (multiRes.getName().equalsIgnoreCase(protFormName)) {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(deprotFormName));
        } else {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(protFormName));
        }
        List<Atom> newAtoms = multiRes.getActive().getAtomList();
        
        // identify which atoms were actually inserted/removed
        List<Atom> removedAtoms = new ArrayList<>();
        List<Atom> insertedAtoms = new ArrayList<>();
        for (Atom oldAtom : oldAtoms) {
            boolean found = false;
            for (Atom newAtom : newAtoms) {
                if (newAtom == oldAtom || newAtom.toNameNumberString().equals(oldAtom.toNameNumberString())) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                removedAtoms.add(oldAtom);
            }
        }
        for (Atom newAtom : newAtoms) {
            boolean found = false;
            for (Atom oldAtom : oldAtoms) {
                if (newAtom == oldAtom || newAtom.toNameNumberString().equals(oldAtom.toNameNumberString())) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                insertedAtoms.add(newAtom);
            }
        }
        if (insertedAtoms.size() + removedAtoms.size() > 1) {
            logger.warning("Protonate: removed + inserted atom count > 1.");
        }
        
        forceFieldEnergy.reInit();
        molDyn.reInit(insertedAtoms, removedAtoms);
        
        StringBuilder sb = new StringBuilder();
        sb.append("Active:\n");
        for (Atom a : multiRes.getActive().getAtomList()) {
            sb.append(String.format("  %s\n", a));
        }
        sb.append("Inactive:\n");
        for (Atom a : multiRes.getInactive().get(0).getAtomList()) {
            sb.append(String.format("  %s\n", a));
        }
        debug(1, sb.toString());
        
        return;
    }
    
    /**
     * Copies atomic coordinates from each active residue to its inactive counterparts.
     * Assumes that these residues differ by only a hydrogen. If said hydrogen is 
     * in an inactive form, its coordinates are updated by geometry with the propagated heavies.
     * @param multiResidues 
     */
    private void propagateInactiveResidues(List<MultiResidue> multiResidues) {
        long startTime = System.nanoTime();
//        debug(3, " Begin multiResidue atomic coordinate propagation:");
//        debug(3, String.format(" multiResidues.size() = %d", multiResidues.size()));
        // Propagate all atom coordinates from active residues to their inactive counterparts.
        for (MultiResidue multiRes : multiResidues) {
            Residue active = multiRes.getActive();
//            debug(3, String.format(" active = %s", active.toString()));
            String activeResName = active.getName();
            List<Residue> inactives = multiRes.getInactive();
            for (Atom activeAtom : active.getAtomList()) {
//                debug(3, String.format(" activeAtom = %s", activeAtom.toString()));
                String activeName = activeAtom.getName();                
                for (Residue inactive : inactives) {
//                    debug(3, String.format(" inactive = %s", inactive.toString()));
//                    StringBuilder sb = new StringBuilder();
//                    sb.append("    inactiveAtomList: ");
//                    for (Atom test : inactive.getAtomList()) {
//                        sb.append(String.format("%s, ", test.getName()));
//                    }
//                    sb.append("\n");
//                    debug(3, sb.toString());
                    Atom inactiveAtom = (Atom) inactive.getAtomNode(activeName);
//                    debug(3, String.format(" inactiveAtom = %s", inactiveAtom));
                    if (inactiveAtom != null) {
                        debug(4, String.format(" Propagating %s\n          to %s.", activeAtom, inactiveAtom));
                        double activeXYZ[] = activeAtom.getXYZ();
                        double inactiveXYZ[] = inactiveAtom.getXYZ();
                        inactiveXYZ[0] = activeXYZ[0];
                        inactiveXYZ[1] = activeXYZ[1];
                        inactiveXYZ[2] = activeXYZ[2];
                        double grad[] = new double[3];
                        activeAtom.getXYZGradient(grad);
                        inactiveAtom.setXYZGradient(grad[0], grad[1], grad[2]);
                        debug(4, String.format("\n          to %s.", activeAtom, inactiveAtom));
                    } else {
                        if (activeName.equals("C") || activeName.equals("O") || activeName.equals("N") || activeName.equals("CA")
                                || activeName.equals("H") || activeName.equals("HA")) {
                            // Backbone atoms aren't supposed to exist in inactive multiResidue components; so no problem.
                        } else if ((activeResName.equals("LYS") && activeName.equals("HZ3"))
                                || (activeResName.equals("TYR") && activeName.equals("HH"))
                                || (activeResName.equals("CYS") && activeName.equals("HG"))
                                || (activeResName.equals("HIS") && (activeName.equals("HD1") || activeName.equals("HE2")))
                                || (activeResName.equals("ASH") && activeName.equals("HD2"))
                                || (activeResName.equals("GLH") && activeName.equals("HE2"))) {
                            // These titratable protons are handled below; so no problem.
                        } else {
                            // Now we have a problem.
                            logger.warning(String.format("Couldn't copy atom_xyz: %s: %s, %s", 
                                    multiRes, activeName, activeAtom.toString()));
                        }
                    }
                }
            }
        }
        
        // If inactive residue is a protonated form, move the stranded hydrogen to new coords (based on propagated heavies).
        for (MultiResidue multiRes : multiResidues) {
            Residue active = multiRes.getActive();
            List<Residue> inactives = multiRes.getInactive();
            for (Residue inactive : inactives) {
                switch (inactive.getName()) {
                    case "LYS": {
                        Atom HZ3 = (Atom) inactive.getAtomNode("HZ3");
                        Atom NZ = (Atom) inactive.getAtomNode("NZ");
                        Atom CE = (Atom) inactive.getAtomNode("CE");
                        Atom HZ1 = (Atom) inactive.getAtomNode("HZ1");
                        BondedUtils.intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HZ3));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HZ3 = buildHydrogen(inactive, "HZ3", NZ, 1.02, CE, 109.5, HZ1, 109.5, -1, k + 9, forceField, null);
                        break;
                    }
                    case "ASH": {
                        Atom HD2 = (Atom) inactive.getAtomNode("HD2");
                        Atom OD2 = (Atom) inactive.getAtomNode("OD2");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom OD1 = (Atom) inactive.getAtomNode("OD1");
                        BondedUtils.intxyz(HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HD2));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HD2 = buildHydrogen(residue, "HD2", OD2, 0.98, CG, 108.7, OD1, 0.0, 0, k + 5, forceField, bondList);
                        break;
                    }
                    case "GLH": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom OE2 = (Atom) inactive.getAtomNode("OE2");
                        Atom CD = (Atom) inactive.getAtomNode("CD");
                        Atom OE1 = (Atom) inactive.getAtomNode("OE1");
                        BondedUtils.intxyz(HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HE2));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HE2 = buildHydrogen(residue, "HE2", OE2, 0.98, CD, 108.7, OE1, 0.0, 0, k + 7, forceField, bondList);
                        break;
                    }
                    case "HIS": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom NE2 = (Atom) inactive.getAtomNode("NE2");
                        Atom CD2 = (Atom) inactive.getAtomNode("CD2");
                        Atom CE1 = (Atom) inactive.getAtomNode("CE1");
                        Atom HD1 = (Atom) inactive.getAtomNode("HD1");
                        Atom ND1 = (Atom) inactive.getAtomNode("ND1");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
                        BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HE2));
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HD1));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 10, forceField, bondList);
                        // Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4, forceField, bondList);
                        break;
                    }
                    case "CYS": {
                        Atom HG = (Atom) inactive.getAtomNode("HG");
                        Atom SG = (Atom) inactive.getAtomNode("SG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        Atom CA = (Atom) inactive.getAtomNode("CA");
                        BondedUtils.intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HG));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HG = buildHydrogen(residue, "HG", SG, 1.34, CB, 96.0, CA, 180.0, 0, k + 3, forceField, bondList);
                        break;
                    }
                    case "TYR": {
                        Atom HH = (Atom) inactive.getAtomNode("HH");
                        Atom OH = (Atom) inactive.getAtomNode("OH");
                        Atom CZ = (Atom) inactive.getAtomNode("CZ");
                        Atom CE2 = (Atom) inactive.getAtomNode("CE2");
                        BondedUtils.intxyz(HH, OH, 0.97, CZ, 108.0, CE2, 0.0, 0);
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HH));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HH = buildHydrogen(residue, "HH", OH, 0.97, CZ, 108.0, CE2, 0.0, 0, k + 9, forceField, bondList);
                        break;
                    }
                    default:
                }
            }
        }
        
        // Print out atomic comparisons.
        if (debugLogLevel >= 4) {
            for (MultiResidue multiRes : multiResidues) {
                Residue active = multiRes.getActive();
                List<Residue> inactives = multiRes.getInactive();
                for (Atom activeAtom : active.getAtomList()) {
                    for (Residue inactive : inactives) {
                        Atom inactiveAtom = (Atom) inactive.getAtomNode(activeAtom.getName());
                        StringBuilder sb = new StringBuilder();
                        sb.append(String.format(" %s\n %s\n", activeAtom, inactiveAtom));
                        debug(4, sb.toString());
                    }
                }
            }
        }
        
        long took = System.nanoTime() - startTime;
//        logger.info(String.format(" Propagating inactive residues took: %d ms", (long) (took * 1e-6)));
    }

    private TitrationType determineTitrationType(String startName, String endName) {
        AminoAcid3 start = AminoAcid3.valueOf(startName);
        AminoAcid3 end = AminoAcid3.valueOf(endName);
        AminoAcid3 protForm = Titratable.valueOf(startName).protForm;
        AminoAcid3 deprotForm = Titratable.valueOf(startName).deprotForm;
        if (protForm != Titratable.valueOf(endName).protForm || deprotForm != Titratable.valueOf(endName).deprotForm) {
            logger.severe(String.format("Titration reaction end states inconsistent: %s %s.", start, end));
            return null;
        }
        if (start == deprotForm && end == protForm) {
            return TitrationType.PROT;
        } else if (start == protForm && end == deprotForm) {
            return TitrationType.DEPROT;
        } else {
            logger.severe(String.format("Titration reaction end states inconsistent: %s %s.", start, end));
            return null;
        }
    }
    
    /**
     * Get the current MC acceptance rate.
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
     * Calculates the electrostatic energy at the current state.
     * @return Energy of the current state.
     */
    private double currentElectrostaticEnergy() {
        double x[] = new double[forceFieldEnergy.getNumberOfVariables() * 3];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x);
        return forceFieldEnergy.getTotalElectrostaticEnergy();
    }
    
    /**
     * Calculates the total energy at the current state.
     * @return Energy of the current state.
     */
    private double currentTotalEnergy() {
        double x[] = new double[forceFieldEnergy.getNumberOfVariables() * 3];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x);
        return forceFieldEnergy.getTotalEnergy();
    }
    
    private void writeSnapshot(boolean beforeChange, StepType stepType, SnapshotsType snapshotsType) {
        // Write the after-step snapshot.
        if (snapshotsType != SnapshotsType.NONE) {
            String postfixA = ".";
            switch (stepType) {
                case TITRATE:
                    postfixA = ".pro";
                    break;
                case ROTAMER:
                    postfixA = ".rot";
                    break;
                case COMBO:
                    postfixA = ".cbo";
                    break;
            }
            String postfixB = (beforeChange) ? "S-" : "F-";
            String filename = FilenameUtils.removeExtension(molAss.getFile().toString()) + postfixA + postfixB + snapshotIndex;
            if (snapshotsType == SnapshotsType.INTERLEAVED) {
                filename = molAss.getFile().getAbsolutePath();
                if (!filename.contains("dyn")) {
                    filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
                }
            }
            File afterFile = new File(filename);
            PDBFilter afterWriter = new PDBFilter(afterFile, molAss, null, null);
            afterWriter.writeFile(afterFile, false);
        }
    }
    
    private void addMolDyn(MolecularDynamics molDyn) {
        this.molDyn = molDyn;
    }
    
    private static int debugLogLevel = 0;
    private static void debug(int level, String message) {
        if (debugLogLevel >= level) {
            logger.info(message);
        }
    }
    
    /**
     * Constant values for intrinsic pKa and reference energy of a CHANGE IN protonation.
     */
    public enum Titratable {
        // NOTE: refEnergy is the energy you should subtract if this form is your PROPOSED TARGET.
        
//        ARG(12.48, 1.00, AminoAcid3.ARD),
        // Standard Forms
        ASP(3.90, -53.188, AminoAcid3.ASH, AminoAcid3.ASP),
        GLU(4.25, -59.390, AminoAcid3.GLH, AminoAcid3.GLU),
        CYS(8.18, 60.834, AminoAcid3.CYS, AminoAcid3.CYD),
        HIS(6.00, -42.920, AminoAcid3.HIS, AminoAcid3.HID),
        LYS(10.53, -50.440, AminoAcid3.LYS, AminoAcid3.LYD),
        TYR(10.07, 34.802, AminoAcid3.TYR, AminoAcid3.TYD),
        // Protonated Forms
        ASH(4.00, +53.188, AminoAcid3.ASH, AminoAcid3.ASP),
        GLH(4.25, +59.390, AminoAcid3.GLH, AminoAcid3.GLU),
        // Deprotonated Forms
        CYD(8.18, -60.834, AminoAcid3.CYS, AminoAcid3.CYD),
        HID(6.00, 42.920, AminoAcid3.HIS, AminoAcid3.HID),
        LYD(10.53, +50.440, AminoAcid3.LYS, AminoAcid3.LYD),
        TYD(10.07, -34.802, AminoAcid3.TYR, AminoAcid3.TYD);

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
    
    private enum MCOverride {
        ACCEPT, REJECT, NONE;
    }
    
    private enum StepType {
        TITRATE, ROTAMER, COMBO;
    }
    
    private enum SnapshotsType {
        NONE, SEPARATE, INTERLEAVED;
    }
    
    private enum TitrationType {
        PROT, DEPROT;
    }
}
