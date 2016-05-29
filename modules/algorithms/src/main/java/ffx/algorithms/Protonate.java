/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import ffx.algorithms.mc.RosenbluthCBMC;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedUtils;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.MultiTerminus;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.utils.PotentialsUtils;
import org.apache.commons.math3.util.FastMath;

/**
 * @author S. LuCore
 */
public class Protonate implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(Protonate.class.getName());
    private static final double NS_TO_SEC = 0.000000001;
    private final boolean logTimings = System.getProperty("cphmd-logTimings") != null ? true : false;
    private long startTime;
    /**
     * The MolecularAssembly.
     */
    private final MolecularAssembly mola;
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
    private static final double BOLTZMANN = 0.0019872041;
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
     * MultiResidue forms of entities from chosenResidues; ready to be
     * (de-/)protonated.
     */
    private List<MultiResidue> titratingResidues = new ArrayList<>();
    private List<MultiTerminus> titratingTermini = new ArrayList<>();
    /**
     * Maps Residue objects to their available Titration enumerations. Filled by
     * the readyUp() method during MultiResidue creation.
     */
    private HashMap<Residue, List<Titration>> titrationMap = new HashMap<>();
    /**
     * Whether to model histidine titration as three states or only two.
     */
    private HistidineMode histidineMode = HistidineMode.ALL;
    /**
     * Everyone's favorite.
     */
    private final Random rng = new Random();
    /**
     * The forcefield being used. Needed by MultiResidue constructor.
     */
    private final ForceField forceField;
    /**
     * The ForceFieldEnergy object being used by MD. Needed by MultiResidue
     * constructor and for reinitializing after a chemical change.
     */
    private final ForceFieldEnergy forceFieldEnergy;
    /**
     * Enum type to specify global override of MC acceptance criteria.
     */
    private MCOverride mcTitrationOverride = MCOverride.NONE;
    /**
     * Writes .s-[num] and .f-[num] files representing before/after MC move
     * structures. Note: The 'after' snapshot represents the change that was
     * PROPOSED, regardless of accept/reject.
     */
    private SnapshotsType snapshotsType = SnapshotsType.NONE;
    /**
     * Snapshot index for the [num] portion of filename above.
     */
    private int snapshotIndex = 0;
    /**
     * Target of the most recently accepted move.
     */
    private Residue previousTarget;
    /**
     * True once the titratingResidues list is ready.
     */
    private boolean finalized = false;
    /**
     * Debug mode: all reference energies set to zero.
     */
    private boolean zeroReferenceEnergies = false;
    private static final double DEFAULT_TEMPERATURE_CRISIS = 6000.0;
    
    private boolean refOverride = System.getProperty("cphmd-refOverride") != null ? true : false;
    private final double refOverrideValue = System.getProperty("cphmd-refOverride") != null ? 
            Double.parseDouble(System.getProperty("cphmd-refOverride")) : 0.0;
    private final double temperatureMonitor = System.getProperty("cphmd-tempMonitor") != null ?
            Double.parseDouble(System.getProperty("cphmd-tempMonitor")) : DEFAULT_TEMPERATURE_CRISIS;
    private final boolean titrateTermini = System.getProperty("cphmd-termini") != null ? true : false;
    private final int terminiOnly = System.getProperty("cphmd-terminiOnly") != null ?
            Integer.parseInt(System.getProperty("cphmd-terminiOnly")) : 0;

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
        if (refOverride) {
            logger.info(" (CpHMD) Reference_Override: " + refOverrideValue);
        }
        String zeroReferenceEnergies = System.getProperty("cphmd-zeroReferences");
        if (zeroReferenceEnergies != null) {
            if (refOverride) {
                logger.severe(" (CpHMD) Incompatible: refOverride,zeroReferenceEnergies");
            }
            this.zeroReferenceEnergies = true;
            logger.info(" OVERRIDE: Zero-ing all reference energies.");
        }
        String debugLogLevel = System.getProperty("debug");
        if (debugLogLevel != null) {
            this.debugLogLevel = Integer.parseInt(debugLogLevel);
        }
        String overrideFlag = System.getProperty("cphmd-override");
        if (overrideFlag != null && overrideFlag.equalsIgnoreCase("accept")) {
            logger.info(" OVERRIDE: Accepting all MC moves.");
            mcTitrationOverride = MCOverride.ACCEPT;
        }
        if (overrideFlag != null && overrideFlag.equalsIgnoreCase("reject")) {
            logger.info(" OVERRIDE: Rejecting all MC moves.");
            mcTitrationOverride = MCOverride.REJECT;
        }
        String beforeAfter = System.getProperty("cphmd-snapshots");
        if (beforeAfter != null) {
            if (beforeAfter.equalsIgnoreCase("true") || beforeAfter.equalsIgnoreCase("separate")) {
                logger.info(" DEBUG: Writing before-and-after MC snapshots.");
                snapshotsType = SnapshotsType.SEPARATE;
            } else if (beforeAfter.equalsIgnoreCase("interleave") || beforeAfter.equalsIgnoreCase("interleaved")) {
                logger.info(" DEBUG: Writing interleaved MC snapshots.");
                snapshotsType = SnapshotsType.INTERLEAVED;
            }
        }
        String histidineMode = System.getProperty("cphmd-histidineMode");
        if (histidineMode != null) {
            if (histidineMode.equalsIgnoreCase("HIE-only")) {
                logger.info(" MC: Histidine mode set to HIE-only.");
                this.histidineMode = HistidineMode.HIE_ONLY;
            } else if (histidineMode.equalsIgnoreCase("HID-only")) {
                logger.info(" MC: Histidine mode set to HID-only.");
                this.histidineMode = HistidineMode.HID_ONLY;
            }
        }

        numMovesAccepted = 0;
        // Set the rotamer library in case we do rotamer MC moves.
        RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        RotamerLibrary.setUseOrigCoordsRotamer(false);

        this.mola = molAss;
        this.forceField = molAss.getForceField();
        this.forceFieldEnergy = molAss.getPotentialEnergy();
        this.mcStepFrequency = (mcStepFrequency == 0) ? Integer.MAX_VALUE : mcStepFrequency;
        this.rotamerStepFrequency = (rotamerStepFrequency == 0) ? Integer.MAX_VALUE : rotamerStepFrequency;
        this.pH = pH;
        this.thermostat = thermostat;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Running Constant-pH MD:\n"));
        sb.append(String.format("     Protonation Step Freq:  %4d\n", mcStepFrequency));
        sb.append(String.format("     Conformation Step Freq: %4d\n", rotamerStepFrequency));
        sb.append(String.format("     system pH:       %7.2f", this.pH));
        logger.info(sb.toString());

        forceFieldEnergy.reInit();
    }

    /**
     * Identify titratable residues and choose them all.
     */
    private void chooseAllTitratables() {
        chosenResidues = new ArrayList<>();
        Polymer polymers[] = mola.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                Residue res = residues.get(j);
                List<Titration> avail = mapTitrations(res, false);
                if (avail.size() > 0) {
                    chosenResidues.add(residues.get(j));
                    // logger.info(String.format(" Titratable: %s", residues.get(j)));
                }
            }
        }
    }

    /**
     * Choose titratables with intrinsic pKa inside (pH-window,pH+window).
     *
     * @param pH
     * @param window
     */
    private void chooseTitratablesInWindow(double pH, double window) {
        chosenResidues = new ArrayList<>();
        Polymer polymers[] = mola.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                Residue res = residues.get(j);
                List<Titration> avail = mapTitrations(res, false);
                for (Titration titration : avail) {
                    double pKa = titration.pKa;
                    if (pKa >= pH - window && pKa <= pH + window) {
                        chosenResidues.add(residues.get(j));
                        // logger.info(String.format(" Titratable: %s", residues.get(j)));
                    }
                }
            }
        }
    }
    
    /**
     * Selecting titrating residues by a list of names, i.e. "LYS,TYR,HIS" will get all K/k/Y/y/H/U/D.
     * @param names 
     */
    public void chooseByName(String names) {
        chosenResidues = new ArrayList<>();
        String tok[] = names.split(",");
        for (int k = 0; k < tok.length; k++) {
            String name = tok[k];
            AminoAcid3 aa3 = AminoAcid3.valueOf(name);
            Polymer polymers[] = mola.getChains();
            for (int i = 0; i < polymers.length; i++) {
                ArrayList<Residue> residues = polymers[i].getResidues();
                for (int j = 0; j < residues.size(); j++) {
                    Residue res = residues.get(j);
                    List<Titration> avail = mapTitrations(res, false);
                    for (Titration titration : avail) {
                        AminoAcid3 from = titration.source;
                        AminoAcid3 to = titration.target;
                        if (aa3 == from || aa3 == to) {
                            chosenResidues.add(res);
                        }                        
                    }
                }
            }
        }
    }

    public void chooseResID(ArrayList<String> crIDs) {
        Polymer[] polymers = mola.getChains();
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

    public void chooseResID(char chain, int resID) {
        Polymer polymers[] = mola.getChains();
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
     * DEPRECATED.
     */
    private void old_readyUp() {
        // Create MultiResidue objects to wrap titratables.
        for (Residue res : chosenResidues) {
            MultiResidue multiRes = new MultiResidue(res, forceField, forceFieldEnergy);
            Polymer polymer = findResiduePolymer(res, mola);
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
     * Must be called after all titratable residues have been chosen, but before
     * beginning MD.
     */
    public void readyUp() {
        // Create MultiTerminus objects to wrap termini.
        if (titrateTermini) {
            for (Residue res : mola.getResidueList()) {
                if (res.getPreviousResidue() == null || res.getNextResidue() == null) {
                    MultiTerminus multiTerminus = new MultiTerminus(res, forceField, forceFieldEnergy, mola);
                    Polymer polymer = findResiduePolymer(res, mola);
                    polymer.addMultiTerminus(res, multiTerminus);
                    forceFieldEnergy.reInit();
                    titratingTermini.add(multiTerminus);
                    logger.info(String.format(" Titrating: %s", multiTerminus));
                }
            }
        }
        // Create MultiResidue objects to wrap titratables.
        for (Residue res : chosenResidues) {
            MultiResidue multiRes = new MultiResidue(res, forceField, forceFieldEnergy);
            Polymer polymer = findResiduePolymer(res, mola);
            polymer.addMultiResidue(multiRes);

            /* OLD WAY
             int resNumber = res.getResidueNumber();
             ResidueType resType = res.getResidueType();
             AminoAcid3 resAA = AminoAcid3.valueOf(res.getName());
             List<Titration> titrs = mapTitrations(res, true);

             // For each available titration, add a target to this MultiResidue.
             for (Titration titr : titrs) {
             String targetName = titr.target.toString();
             Residue targetRes = new Residue(targetName, resNumber, resType);
             multiRes.addResidue(targetRes);
             // Also map the reverse titration.
             List<Titration> reversals = mapTitrations(targetRes, true);
             //  For Histidine, if the original residue was U then mapping titrations
             //  for the reversal will yield H->Z as well as H->U.
             //  We have to make sure that Z gets added as a MultiResidue option in this case.
             for (Titration reversal : reversals) {
             if (reversal.target != resAA) {
             String altReversalName = reversal.target.toString();
             Residue altReversalRes = new Residue(altReversalName, resNumber, resType);
             multiRes.addResidue(altReversalRes);
             }
             }
             }
             */
            // NEW WAY USING RECURSIVE CALL
            recursiveMap(res, multiRes);

            // Switch back to the original form and ready the ForceFieldEnergy.
            multiRes.setActiveResidue(res);
            forceFieldEnergy.reInit();
            titratingResidues.add(multiRes);
            logger.info(String.format(" Titrating: %s", multiRes));
        }

        finalized = true;
    }

    /**
     * Recursively maps Titration events and adds target Residues to a
     * MultiResidue object.
     *
     * @param member
     * @param multiRes
     */
    private void recursiveMap(Residue member, MultiResidue multiRes) {
        if (finalized) {
            logger.severe("Programming error: improper function call.");
        }
        // Map titrations for this member.
        List<Titration> titrs = mapTitrations(member, true);

        // For each titration, check whether it needs added as a MultiResidue option.
        for (Titration titr : titrs) {
            // Allow manual override of Histidine treatment.
            if ((titr.target == AminoAcid3.HID && histidineMode == HistidineMode.HIE_ONLY)
                    || (titr.target == AminoAcid3.HIE && histidineMode == HistidineMode.HID_ONLY)) {
                continue;
            }
            // Find all the choices currently available to this MultiResidue.
            List<String> choices = new ArrayList<>();
            for (Residue choice : multiRes.getConsideredResidues()) {
                choices.add(choice.getName());
            }
            // If this Titration target is not a choice for the MultiResidue, then add it.
            String targetName = titr.target.toString();
            if (!choices.contains(targetName)) {
                int resNumber = member.getResidueNumber();
                ResidueType resType = member.getResidueType();
                Residue newChoice = new Residue(targetName, resNumber, resType);
                multiRes.addResidue(newChoice);
                // Recursively call this method on each added choice.
                recursiveMap(newChoice, multiRes);
            }
        }
    }

    /**
     * Maps available Titration enums to a given Residue; used to fill the
     * titrationMap field.
     *
     * @param res
     * @param store add identified Titrations to the HashMap
     * @return list of Titrations available for the given residue
     */
    private List<Titration> mapTitrations(Residue res, boolean store) {
        if (finalized) {
            logger.severe("Programming error: improper function call.");
        }
        if (store && titrationMap.containsKey(res)) {
            logger.warning(String.format("Titration map already contained key for residue: %s", res));
            return new ArrayList<>();
        }
        AminoAcid3 source = AminoAcid3.valueOf(res.getName());
        List<Titration> avail = new ArrayList<>();
        for (Titration titr : Titration.values()) {
            // Allow manual override of Histidine treatment.
            if ((titr.target == AminoAcid3.HID && histidineMode == HistidineMode.HIE_ONLY)
                    || (titr.target == AminoAcid3.HIE && histidineMode == HistidineMode.HID_ONLY)) {
                continue;
            }
            if (titr.source == source) {
                avail.add(titr);
            }
        }
        if (store) {
            if (avail.size() == 0) {
                logger.severe(String.format("Chosen residue couldn't map to any Titration object: %s", res));
            }
            titrationMap.put(res, avail);
        }
        return avail;
    }
    
    private void meltdown() {
        writeSnapshot(".meltdown-");
        forceFieldEnergy.energy(false, true);
        if (forceFieldEnergy.getBondEnergy() > 1000) {
            for (ROLS rols : previousTarget.getBondList()) {
                ((Bond) rols).log();
            }
        }
        if (forceFieldEnergy.getAngleEnergy() > 1000) {
            for (ROLS rols : previousTarget.getAngleList()) {
                ((Angle) rols).log();
            }
        }
        if (forceFieldEnergy.getVanDerWaalsEnergy() > 1000) {
            for (Atom a1 : previousTarget.getAtomList()) {
                for (Atom a2 : mola.getAtomArray()) {
                    if (a1 == a2 || a1.getBond(a2) != null) {
                        continue;
                    }
                    double dist = FastMath.sqrt(
                            FastMath.pow((a1.getX()-a2.getX()),2) +
                            FastMath.pow((a1.getY()-a2.getY()),2) +
                            FastMath.pow((a1.getZ()-a2.getZ()),2));
                    if (dist < 0.8*(a1.getVDWR() + a2.getVDWR())) {
                        logger.warning(String.format("Close vdW contact for atoms: \n   %s\n   %s", a1, a2));
                    }
                }
            }
        }
        logger.severe(String.format("Temperature above critical threshold: %f", thermostat.getCurrentTemperature()));
    }
    
    /**
     * The primary driver. Called by the MD engine at each dynamics step.
     */
    @Override
    public boolean mcUpdate(MolecularAssembly molAss) {
        startTime = System.nanoTime();
        if (!finalized) {
            logger.severe("Monte-Carlo protonation engine was not finalized!");
        }
        if (thermostat.getCurrentTemperature() > temperatureMonitor) {
            meltdown();
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
//            if (logTimings) {
//                long took = System.nanoTime() - startTime;
//                logger.info(String.format(" CpHMD propagation time: %.6f", took * NS_TO_SEC));
//            }
            return false;
        }

        // Randomly choose a target titratable residue to attempt protonation switch.
        int random = rng.nextInt(titratingResidues.size() + titratingTermini.size());
        
        if (terminiOnly == 0) {
            random = rng.nextInt(titratingTermini.size()) + titratingResidues.size();
        } else if (terminiOnly == 1) {
            random = titratingResidues.size();
        } else if (terminiOnly == 2) {
            random = titratingResidues.size() + 1;
        }
        
        if (random >= titratingResidues.size()) {
            Residue target = titratingTermini.get(random - titratingResidues.size());
            boolean accepted = tryTerminusTitration((MultiTerminus) target);
            snapshotIndex++;
            if (accepted) {
                molDyn.reInit();
                previousTarget = target;
            }
            return accepted;
        }
        MultiResidue targetMulti = titratingResidues.get(random);

        // Check whether rotamer moves are possible for the selected residue.
        Residue targetMultiActive = targetMulti.getActive();
        Rotamer[] targetMultiRotamers = targetMultiActive.getRotamers();
        if (targetMultiRotamers == null || targetMultiRotamers.length <= 1) {
            if (stepType == StepType.ROTAMER) {
                return false;
            } else if (stepType == StepType.COMBO) {
                stepType = StepType.TITRATE;
            }
        }

        // forceFieldEnergy.checkAtoms();
        // Perform the MC move.
        boolean accepted = false;
        switch (stepType) {
            case TITRATE:
                accepted = tryTitrationStep(targetMulti);
                break;
            case ROTAMER:
//                accepted = tryRotamerStep(targetMulti);
                accepted = tryCBMCStep(targetMulti);
                break;
            case COMBO:
//                accepted = tryComboStep(targetMulti);
                accepted = tryCBMCStep(targetMulti);
                accepted = accepted || tryTitrationStep(targetMulti);
                break;
        }

        // forceFieldEnergy.checkAtoms();
        // Increment the shared snapshot counter.
//        if (logTimings) {
//            long took = System.nanoTime() - startTime;
//            logger.info(String.format(" CpHMD step time:        %.6f", took * NS_TO_SEC));
//        }
        snapshotIndex++;
        if (accepted) {
            previousTarget = targetMulti;
        }
        return accepted;
    }

    /**
     * Perform a titration MC move.
     *
     * @param targetMulti
     * @return accept/reject
     */
    private boolean tryTitrationStep(Residue target) {
        boolean terminus = false;
        MultiResidue targetMulti = null;
        MultiTerminus targetTerm = null;
        if (target instanceof MultiResidue) {
            targetMulti = (MultiResidue) target;
            terminus = false;
        } else if (target instanceof MultiTerminus) {
            targetTerm = (MultiTerminus) target;
            terminus = true;
        } else {
            logger.warning("Improper method call.");
        }
        // Record the pre-change electrostatic energy.
        double previousElectrostaticEnergy = currentElectrostaticEnergy();

        // Write the pre-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);
        String startString = target.toString();
        String startName = target.getName();

        double pKaref = 0;
        double dG_ref = 0;
        Titration titration = null;
        TitrationType type = null;
        
        if (terminus) {
            if (targetTerm.end == MultiTerminus.END.NTERM) {
                pKaref = 10.0;
                dG_ref = 0.0;
            } else {
                pKaref = 3.0;
                dG_ref = 0.0;
            }
            targetTerm.titrateTerminus_v1(thermostat.getCurrentTemperature());
        } else {
            // Choose from the list of available titrations for the active residue.
            List<Titration> avail = titrationMap.get(targetMulti.getActive());
            titration = avail.get(rng.nextInt(avail.size()));
            type = titration.type;

            // Perform the chosen titration.
            performTitration(targetMulti, titration);

            // Test the MC criterion for a titration step.
            pKaref = titration.pKa;
            dG_ref = titration.refEnergy;
        }
        // Write the post-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);
        
        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;

        if (zeroReferenceEnergies) {
            dG_ref = 0.0;
        }
        if (refOverride) {
            dG_ref = refOverrideValue;
        }

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the
         * reference compound
         */
        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEP) {
            prefix = -prefix;
        }
        // Change this to use a single value for reference and then switch based on reaction.
        // Either positive ref == deprotonation or == standard -> nonstandard transition.
        if (type == TitrationType.PROT) {
            dG_ref = -dG_ref;
        }
        double postfix = dG_elec - dG_ref;
        double dG_MC = prefix + postfix;

//        StringBuilder sb = new StringBuilder();
//        sb.append(String.format(" Assessing possible MC protonation step:\n"));
//        sb.append(String.format("     %s --> %s\n", startString, targetMulti.toString()));
//        sb.append(String.format("     pKaref:  %7.2f\n", pKaref));
//        sb.append(String.format("     dG_ref:  %7.2f\n", dG_ref));
//        sb.append(String.format("     dG_elec: %16.8f\n", dG_elec));
//        sb.append(String.format("     dG_MC:   %16.8f\n", dG_MC));
//        sb.append(String.format("     -----\n"));
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC protonation step:\n"));
        sb.append(String.format("     %s --> %s\n", startString, target.toString()));
        sb.append(String.format("     dG_ref:  %7.2f                pKaref:  %7.2f\n", dG_ref, pKaref));
        sb.append(String.format("     pH_term: %9.4f              elec_term: %10.4f\n", prefix, postfix));
        sb.append(String.format("     dG_elec: %9.4f              dG_MC:     %10.4f\n", dG_elec, dG_MC));
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
        sb.append(String.format("     crit:    %9.4f              rng:       %10.4f\n", criterion, metropolis));
        if ((metropolis < criterion && mcTitrationOverride != MCOverride.REJECT) || mcTitrationOverride == MCOverride.ACCEPT) {
            numMovesAccepted++;
            long took = System.nanoTime() - startTime;
            sb.append(String.format("     Accepted!                                                %1.3f", took * NS_TO_SEC));
            logger.info(sb.toString());
            return true;
        }

        // Move was rejected, undo the titration state change.
        if (terminus) {
            
        } else {
            Titration reverse = inverseReactions.get(titration);
            performTitration(targetMulti, reverse);
        }
        long took = System.nanoTime() - startTime;
        sb.append(String.format("     Denied.                                                  %1.3f", took * NS_TO_SEC));
        logger.info(sb.toString());
        return false;
    }
    
    /**
     * Perform a titration MC move.
     *
     * @param targetMulti
     * @return accept/reject
     */
    private boolean tryTerminusTitration(MultiTerminus target) {
        // Record the pre-change electrostatic energy.
        double previousElectrostaticEnergy = currentElectrostaticEnergy();

        // Write the pre-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);
        String startString = target.toString();
        String startName = target.getName();

        double pKaref = 0;
        double dG_ref = 0;
        Titration titration = null;
        TitrationType type = null;
        
        if (target.end == MultiTerminus.END.NTERM) {
            pKaref = 10.0;
            dG_ref = 0.0;
        } else if (target.end == MultiTerminus.END.CTERM) {
            pKaref = 3.0;
            dG_ref = 0.0;
        }
        boolean beganCharged = target.isCharged;
        target.titrateTerminus_v1(thermostat.getCurrentTemperature());
        // Write the post-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, snapshotsType);
        
        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;

        if (zeroReferenceEnergies) {
            dG_ref = 0.0;
        }
        if (refOverride) {
            dG_ref = refOverrideValue;
        }

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the
         * reference compound
         */
        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEP) {
            prefix = -prefix;
        }
        // Change this to use a single value for reference and then switch based on reaction.
        // Either positive ref == deprotonation or == standard -> nonstandard transition.
        if (type == TitrationType.PROT) {
            dG_ref = -dG_ref;
        }
        double postfix = dG_elec - dG_ref;
        double dG_MC = prefix + postfix;

//        StringBuilder sb = new StringBuilder();
//        sb.append(String.format(" Assessing possible MC protonation step:\n"));
//        sb.append(String.format("     %s --> %s\n", startString, targetMulti.toString()));
//        sb.append(String.format("     pKaref:  %7.2f\n", pKaref));
//        sb.append(String.format("     dG_ref:  %7.2f\n", dG_ref));
//        sb.append(String.format("     dG_elec: %16.8f\n", dG_elec));
//        sb.append(String.format("     dG_MC:   %16.8f\n", dG_MC));
//        sb.append(String.format("     -----\n"));
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC protonation step:\n"));
        if (beganCharged) {
            sb.append(String.format("     %sc --> %sn\n", startString, target.toString()));
        } else {
            sb.append(String.format("     %sn --> %sc\n", startString, target.toString()));
        }
        sb.append(String.format("     dG_ref:  %7.2f                pKaref:  %7.2f\n", dG_ref, pKaref));
        sb.append(String.format("     pH_term: %9.4f              elec_term: %10.4f\n", prefix, postfix));
        sb.append(String.format("     dG_elec: %9.4f              dG_MC:     %10.4f\n", dG_elec, dG_MC));
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
        sb.append(String.format("     crit:    %9.4f              rng:       %10.4f\n", criterion, metropolis));
        if ((metropolis < criterion && mcTitrationOverride != MCOverride.REJECT) || mcTitrationOverride == MCOverride.ACCEPT) {
            numMovesAccepted++;
            long took = System.nanoTime() - startTime;
            sb.append(String.format("     Accepted!                                                %1.3f", took * NS_TO_SEC));
            logger.info(sb.toString());
            return true;
        }

        // Move was rejected, undo the titration state change.
        target.titrateTerminus_v1(thermostat.getCurrentTemperature());
        long took = System.nanoTime() - startTime;
        sb.append(String.format("     Denied.                                                  %1.3f", took * NS_TO_SEC));
        logger.info(sb.toString());
        return false;
    }
    
    private boolean tryCBMCStep(MultiResidue targetMulti) {
        List<Residue> targets = new ArrayList<>();
        targets.add(targetMulti.getActive());
        int trialSetSize = 5;   // cost still scales with this, unfortunately
        int mcFrequency = 1;    // irrelevant for manual step call
        boolean writeSnapshots = false;
        System.setProperty("cbmc-type", "CHEAP");
        RosenbluthCBMC cbmc = new RosenbluthCBMC(mola, mola.getPotentialEnergy(), null,
            targets, mcFrequency, trialSetSize, writeSnapshots);
        boolean accepted = cbmc.cbmcStep();
        if (logTimings) {
            long took = System.nanoTime() - startTime;
            logger.info(String.format(" CBMC time: %1.3f", took * NS_TO_SEC));
        }
        return accepted;
    }

    /**
     * Attempt a rotamer MC move.
     *
     * @param targetMulti
     * @return accept/reject
     */
    private boolean tryRotamerStep(MultiResidue targetMulti) {
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();

        // Write the before-step snapshot.
        writeSnapshot(true, StepType.ROTAMER, snapshotsType);

        // Save coordinates so we can return to them if move is rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        ResidueState origState = residue.storeState();
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origState, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);
        // Select a new rotamer and swap to it.
        //Rotamer rotamers[] = residue.getRotamers();
        Rotamer[] rotamers = residue.getRotamers();
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);

        // Write the post-rotamer change snapshot.
        writeSnapshot(false, StepType.ROTAMER, snapshotsType);

        // Check the MC criterion.
        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
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

    /**
     * Attempt a combination titration/rotamer MC move.
     *
     * @param targetMulti
     * @return accept/reject
     */
    private boolean tryComboStep(MultiResidue targetMulti) {
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();
        double previousElectrostaticEnergy = currentElectrostaticEnergy();

        // Write the pre-combo snapshot.
        writeSnapshot(true, StepType.COMBO, snapshotsType);
        String startString = targetMulti.toString();
        String startName = targetMulti.getActive().getName();

        // Choose from the list of available titrations for the active residue.
        List<Titration> avail = titrationMap.get(targetMulti.getActive());
        Titration titration = avail.get(rng.nextInt(avail.size()));
        TitrationType type = titration.type;

        // Perform the chosen titration.
        performTitration(targetMulti, titration);

        // Change rotamer state, but first save coordinates so we can return to them if rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        ResidueState origState = residue.storeState();
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origState, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);

        // Swap to the new rotamer.
        //Rotamer rotamers[] = residue.getRotamers();
        Rotamer[] rotamers = residue.getRotamers();
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);

        // Write the post-combo snapshot.
        writeSnapshot(false, StepType.COMBO, snapshotsType);

        // Evaluate both MC criteria.
        String endName = targetMulti.getActive().getName();

        // Evaluate the titration probability of the step.
        double pKaref = titration.pKa;
        double dG_ref = titration.refEnergy;
        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;

        if (zeroReferenceEnergies) {
            dG_ref = 0.0;
        }

        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEP) {
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
                Titration reverse = inverseReactions.get(titration);
                performTitration(targetMulti, reverse);
                return false;
            }
        }
    }

    /**
     * Perform the requested titration on the given MultiResidue and
     * reinitialize the FF.
     *
     * @param multiRes
     * @param titration
     */
    private void performTitration(MultiResidue multiRes, Titration titration) {
        if (titration.source != AminoAcid3.valueOf(multiRes.getActive().getName())) {
            logger.severe(String.format("Requested titration source didn't match target MultiResidue: %s", multiRes.toString()));
        }

        List<Atom> oldAtoms = multiRes.getActive().getAtomList();
        boolean success = multiRes.requestSetActiveResidue(titration.target);
        if (!success) {
            logger.severe(String.format("Couldn't perform requested titration for MultiRes: %s", multiRes.toString()));
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
        molDyn.reInit();

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
     * Copies atomic coordinates from each active residue to its inactive
     * counterparts. Assumes that these residues differ by only a hydrogen. If
     * said hydrogen is in an inactive form, its coordinates are updated by
     * geometry with the propagated heavies.
     *
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
                        // Propagate position and gradient.
                        double activeXYZ[] = activeAtom.getXYZ(null);
                        inactiveAtom.setXYZ(activeXYZ);
                        double grad[] = new double[3];
                        activeAtom.getXYZGradient(grad);
                        inactiveAtom.setXYZGradient(grad[0], grad[1], grad[2]);
                        // Propagate velocity, acceleration, and previous acceleration.
                        double activeVelocity[] = new double[3];
                        activeAtom.getVelocity(activeVelocity);
                        inactiveAtom.setVelocity(activeVelocity);
                        double activeAccel[] = new double[3];
                        activeAtom.getAcceleration(activeAccel);
                        inactiveAtom.setAcceleration(activeAccel);
                        double activePrevAcc[] = new double[3];
                        activeAtom.getPreviousAcceleration(activePrevAcc);
                        inactiveAtom.setPreviousAcceleration(activePrevAcc);
                        debug(4, String.format("\n          to %s.", activeAtom, inactiveAtom));
                    } else {
                        if (activeName.equals("C") || activeName.equals("O") || activeName.equals("N") || activeName.equals("CA")
                                || activeName.equals("H") || activeName.equals("HA")) {
                            // Backbone atoms aren't supposed to exist in inactive multiResidue components; so no problem.
                        } else if ((activeResName.equals("LYS") && activeName.equals("HZ3"))
                                || (activeResName.equals("TYR") && activeName.equals("HH"))
                                || (activeResName.equals("CYS") && activeName.equals("HG"))
                                || (activeResName.equals("HIS") && (activeName.equals("HD1") || activeName.equals("HE2")))
                                || (activeResName.equals("HID") && activeName.equals("HD1"))
                                || (activeResName.equals("HIE") && activeName.equals("HE2"))
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
        // Also give the stranded hydrogen a maxwell velocity and remove its accelerations.
        for (MultiResidue multiRes : multiResidues) {
            Residue active = multiRes.getActive();
            List<Residue> inactives = multiRes.getInactive();
            for (Residue inactive : inactives) {
                Atom resetMe = null;
                switch (inactive.getName()) {
                    case "LYS": {
                        Atom HZ3 = (Atom) inactive.getAtomNode("HZ3");
                        Atom NZ = (Atom) inactive.getAtomNode("NZ");
                        Atom CE = (Atom) inactive.getAtomNode("CE");
                        Atom HZ1 = (Atom) inactive.getAtomNode("HZ1");
                        BondedUtils.intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
                        resetMe = HZ3;
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
                        resetMe = HD2;
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
                        resetMe = HE2;
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
                        // Manual reset since we gotta reset two of 'em.
                        HE2.setXYZGradient(0, 0, 0);
                        HE2.setVelocity(thermostat.maxwellIndividual(HE2.getMass()));
                        HE2.setAcceleration(new double[]{0, 0, 0});
                        HE2.setPreviousAcceleration(new double[]{0, 0, 0});
                        HD1.setXYZGradient(0, 0, 0);
                        HD1.setVelocity(thermostat.maxwellIndividual(HD1.getMass()));
                        HD1.setAcceleration(new double[]{0, 0, 0});
                        HD1.setPreviousAcceleration(new double[]{0, 0, 0});
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HE2));
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HD1));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 10, forceField, bondList);
                        // Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4, forceField, bondList);
                        break;
                    }
                    case "HID": {
                        Atom HD1 = (Atom) inactive.getAtomNode("HD1");
                        Atom ND1 = (Atom) inactive.getAtomNode("ND1");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                        resetMe = HD1;
                        // Parameters from AminoAcidUtils, line:
                        // Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4, forceField, bondList);
                        break;
                    }
                    case "HIE": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom NE2 = (Atom) inactive.getAtomNode("NE2");
                        Atom CD2 = (Atom) inactive.getAtomNode("CD2");
                        Atom CE1 = (Atom) inactive.getAtomNode("CE1");
                        BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
                        resetMe = HE2;
                        // Parameters from AminoAcidUtils, line:
                        // Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 9, forceField, bondList);
                        break;
                    }
                    case "CYS": {
                        Atom HG = (Atom) inactive.getAtomNode("HG");
                        Atom SG = (Atom) inactive.getAtomNode("SG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        Atom CA = (Atom) inactive.getAtomNode("CA");
                        BondedUtils.intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
                        resetMe = HG;
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
                        resetMe = HH;
                        debug(4, String.format(" Moved 'stranded' hydrogen %s.", HH));
                        // Parameters from AminoAcidUtils, line:
                        // Atom HH = buildHydrogen(residue, "HH", OH, 0.97, CZ, 108.0, CE2, 0.0, 0, k + 9, forceField, bondList);
                        break;
                    }
                    default:
                }
                if (resetMe != null) {
                    resetMe.setXYZGradient(0, 0, 0);
                    resetMe.setVelocity(thermostat.maxwellIndividual(resetMe.getMass()));
                    resetMe.setAcceleration(new double[]{0, 0, 0});
                    resetMe.setPreviousAcceleration(new double[]{0, 0, 0});
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

    /**
     * Get the current MC acceptance rate.
     *
     * @return the acceptance rate.
     */
    public double getAcceptanceRate() {
        // Intentional integer division.
        int numTries = stepCount / mcStepFrequency;
        return (double) numMovesAccepted / numTries;
    }

    /**
     * Locate to which Polymer in the MolecularAssembly a given Residue belongs.
     *
     * @param residue
     *
     * @param molecularAssembly
     *
     * @return the Polymer where the passed Residue is located.
     */
    private Polymer findResiduePolymer(Residue residue,
            MolecularAssembly molecularAssembly) {
        if (residue.getChainID() == null) {
            logger.severe("No chain ID for residue " + residue);
        }
        Polymer polymers[] = molecularAssembly.getChains();
        Polymer location = null;
        for (Polymer p : polymers) {
            if (p.getChainID() == residue.getChainID()) {
                location = p;
            }
        }
        if (location == null) {
            logger.severe("Couldn't find polymer for residue " + residue);
        }
        return location;
    }

    /**
     * Calculates the electrostatic energy at the current state.
     *
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
     *
     * @return Energy of the current state.
     */
    private double currentTotalEnergy() {
        double x[] = new double[forceFieldEnergy.getNumberOfVariables() * 3];
        forceFieldEnergy.getCoordinates(x);
        forceFieldEnergy.energy(x);
        return forceFieldEnergy.getTotalEnergy();
    }

    private void writeSnapshot(String extension) {
        String filename = FilenameUtils.removeExtension(mola.getFile().toString()) + extension + snapshotIndex;
        if (snapshotsType == SnapshotsType.INTERLEAVED) {
            filename = mola.getFile().getAbsolutePath();
            if (!filename.contains("dyn")) {
                filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
            }
        }
        File file = new File(filename);
        PDBFilter writer = new PDBFilter(file, mola, null, null);
        writer.writeFile(file, false);
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
            String filename = FilenameUtils.removeExtension(mola.getFile().toString()) + postfixA + postfixB + snapshotIndex;
            if (snapshotsType == SnapshotsType.INTERLEAVED) {
                filename = mola.getFile().getAbsolutePath();
                if (!filename.contains("dyn")) {
                    filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
                }
            }
            File afterFile = new File(filename);
            PDBFilter afterWriter = new PDBFilter(afterFile, mola, null, null);
            afterWriter.writeFile(afterFile, false);
        }
    }

    public void addMolDyn(MolecularDynamics molDyn) {
        this.molDyn = molDyn;
    }

    /**
     * DEPRECATED. Constant values for intrinsic pKa and reference energy of a
     * CHANGE IN protonation. NOTE: refEnergy is the energy you should subtract
     * if this form is your PROPOSED TARGET.
     */
    private enum Titratable {

        // Standard Forms
        ASP(3.90, -53.188, AminoAcid3.ASH, AminoAcid3.ASP),
        GLU(4.25, -59.390, AminoAcid3.GLH, AminoAcid3.GLU),
        CYS(8.18, 60.834, AminoAcid3.CYS, AminoAcid3.CYD),
        HIS(6.00, -42.920, AminoAcid3.HIS, AminoAcid3.HID),
        LYS(10.53, -50.440, AminoAcid3.LYS, AminoAcid3.LYD),    // new: 48.6928
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
    }

    /**
     * Enumerated titration reactions for source/target amino acid pairs.
     */
    public enum Titration {

        Ctoc(8.18, -60.168, TitrationType.DEP, AminoAcid3.CYS, AminoAcid3.CYD),
        ctoC(8.18, +60.168, TitrationType.PROT, AminoAcid3.CYD, AminoAcid3.CYS),
        Dtod(3.90, +53.188, TitrationType.PROT, AminoAcid3.ASP, AminoAcid3.ASH),
        dtoD(3.90, -53.188, TitrationType.DEP, AminoAcid3.ASH, AminoAcid3.ASP),
        Etoe(4.25, +59.390, TitrationType.PROT, AminoAcid3.GLU, AminoAcid3.GLH),
        etoE(4.25, -59.390, TitrationType.DEP, AminoAcid3.GLH, AminoAcid3.GLU),
        Ktok(10.53, +50.440, TitrationType.DEP, AminoAcid3.LYS, AminoAcid3.LYD),    // new dG_elec: 48.6928
        ktoK(10.53, -50.440, TitrationType.PROT, AminoAcid3.LYD, AminoAcid3.LYS),
        Ytoy(10.07, -34.961, TitrationType.DEP, AminoAcid3.TYR, AminoAcid3.TYD),
        ytoY(10.07, +34.961, TitrationType.PROT, AminoAcid3.TYD, AminoAcid3.TYR),
        HtoU(6.00, +42.923, TitrationType.DEP, AminoAcid3.HIS, AminoAcid3.HID),
        UtoH(6.00, -42.923, TitrationType.PROT, AminoAcid3.HID, AminoAcid3.HIS),
        HtoZ(6.00, +00.000, TitrationType.DEP, AminoAcid3.HIS, AminoAcid3.HIE),
        ZtoH(6.00, +00.000, TitrationType.PROT, AminoAcid3.HIE, AminoAcid3.HIS);

        public final double pKa, refEnergy;
        public final TitrationType type;
        public final AminoAcid3 source, target;

        Titration(double pKa, double refEnergy, TitrationType type,
                AminoAcid3 source, AminoAcid3 target) {
            this.pKa = pKa;
            this.refEnergy = refEnergy;
            this.type = type;
            this.source = source;
            this.target = target;
        }
    }

    /**
     * Maps each Titration reaction to its inverse for the purpose of reverting
     * failed MC steps.
     */
    private static final HashMap<Titration, Titration> inverseReactions = new HashMap<>();

    static {
        inverseReactions.put(Titration.Ctoc, Titration.ctoC);
        inverseReactions.put(Titration.ctoC, Titration.Ctoc);
        inverseReactions.put(Titration.Dtod, Titration.dtoD);
        inverseReactions.put(Titration.dtoD, Titration.Dtod);
        inverseReactions.put(Titration.Etoe, Titration.etoE);
        inverseReactions.put(Titration.etoE, Titration.Etoe);
        inverseReactions.put(Titration.Ktok, Titration.ktoK);
        inverseReactions.put(Titration.ktoK, Titration.Ktok);
        inverseReactions.put(Titration.Ytoy, Titration.ytoY);
        inverseReactions.put(Titration.ytoY, Titration.Ytoy);
        inverseReactions.put(Titration.HtoU, Titration.UtoH);
        inverseReactions.put(Titration.UtoH, Titration.HtoU);
        inverseReactions.put(Titration.HtoZ, Titration.ZtoH);
        inverseReactions.put(Titration.ZtoH, Titration.HtoZ);
    }

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

        PROT, DEP;
    }

    public enum HistidineMode {

        ALL, HID_ONLY, HIE_ONLY;
    }

    private static int debugLogLevel = 0;

    private static void debug(int level, String message) {
        if (debugLogLevel >= level) {
            logger.info(message);
        }
    }
}
