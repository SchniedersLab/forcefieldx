/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.algorithms.MolecularDynamics.MonteCarloNotification;
import ffx.algorithms.mc.RosenbluthCBMC;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
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
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.extended.TitrationESV;
import ffx.potential.extended.TitrationUtils;
import ffx.potential.extended.TitrationUtils.HistidineMode;
import ffx.potential.extended.TitrationUtils.MCOverride;
import ffx.potential.extended.TitrationUtils.Snapshots;
import ffx.potential.extended.TitrationUtils.Titration;
import ffx.potential.extended.TitrationUtils.TitrationConfig;
import ffx.potential.extended.TitrationUtils.TitrationType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;

import static ffx.potential.extended.TitrationUtils.findResiduePolymer;
import static ffx.potential.extended.TitrationUtils.inactivateResidue;
import static ffx.potential.extended.TitrationUtils.performTitration;
import static ffx.potential.extended.TitrationUtils.propagateInactiveResidues;

/**
 * @author S. LuCore
 */
public class PhMD implements MonteCarloListener {
    private static final Logger logger = Logger.getLogger(PhMD.class.getName());
    
    /**
     * Disallow code paths other than that which is being validated.
     * DEBUG: reset to false.
     */
    private static final boolean CAUTIOUS = true;
    
    private static final double NS_TO_SEC = 0.000000001;
    private long startTime;
    private StringBuilder discountLogger;
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
    private int rotamerStepFrequency = Integer.MAX_VALUE;
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
    private List<MultiResidue> titratingMultis = new ArrayList<>();
    private List<MultiTerminus> titratingTermini = new ArrayList<>();
    private List<ExtendedVariable> titratingESVs = new ArrayList<>();
    /**
     * Maps Residue objects to their available Titration enumerations. Filled by
     * the readyup() method during MultiResidue creation.
     */
    private HashMap<Residue, List<Titration>> titrationMap = new HashMap<>();
    /**
     * Everyone's favorite.
     */
    private final Random rng = new Random();
    /**
     * The forcefield being used. Needed by MultiResidue constructor.
     */
    private final ForceField ff;
    /**
     * The ForceFieldEnergy object being used by MD. Needed by MultiResidue
     * constructor and for reinitializing after a chemical change.
     */
    private final ForceFieldEnergy ffe;
    /**
     * Snapshot index for the [num] portion of filename above.
     */
    private int snapshotIndex = 0;
    /**
     * Target of the most recently accepted move.
     */
    private Residue previousTarget;

    private ExtendedSystem esvSystem;
    private Object[] mdOptions;
    private DynamicsLauncher mdLauncher;
    private RotamerLibrary library;
    private final TitrationConfig config;
    
    private final Distribution distribution;
    public enum Distribution {
        DISCRETE, CONTINUOUS;
    }

    /**
     * Construct a Monte-Carlo protonation state switching mechanism.
     *
     * @param mola the molecular assembly
     * @param mcStepFrequency number of MD steps between switch attempts
     * @param pH the simulation pH
     * @param thermostat the MD thermostat
     */
    public PhMD(Distribution distribution, MolecularAssembly mola, MolecularDynamics molDyn,
            List<Residue> titrating, double pH, int mcStepFrequency) {
        this.config = new TitrationConfig();
        this.distribution = distribution;
        this.mola = mola;
        this.molDyn = molDyn;
        this.thermostat = molDyn.getThermostat();
        this.ff = mola.getForceField();
        this.ffe = mola.getPotentialEnergy();
        this.mcStepFrequency = (mcStepFrequency == 0) ? Integer.MAX_VALUE : mcStepFrequency;
        this.pH = pH;
        this.numMovesAccepted = 0;
        this.chosenResidues = titrating;
        
        // Set the rotamer library in case we do rotamer MC moves.
        library = RotamerLibrary.getDefaultLibrary();
        library.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        library.setUseOrigCoordsRotamer(false);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Running Constant-pH MD:\n"));
        sb.append(String.format("     Protonation Step Freq:  %4d\n", mcStepFrequency));
        sb.append(String.format("     Conformation Step Freq: %4d\n", rotamerStepFrequency));
        sb.append(String.format("     system pH:       %7.2f", this.pH));
        logger.info(sb.toString());

        reInitialize(true, false);
        readyup();
    }
    
    public void setRotamerStepFrequency(int frequency) {
        this.rotamerStepFrequency = (frequency <= 0) ? Integer.MAX_VALUE : frequency;
    }
    
    private void readyup() {
        // Create MultiTerminus objects to wrap termini.
        if (config.titrateTermini) {
            for (Residue res : mola.getResidueList()) {
                if (res.getPreviousResidue() == null || res.getNextResidue() == null) {
                    MultiTerminus multiTerminus = new MultiTerminus(res, ff, ffe, mola);
                    Polymer polymer = findResiduePolymer(res, mola);
                    polymer.addMultiTerminus(res, multiTerminus);
                    reInitialize(true, false);
                    titratingTermini.add(multiTerminus);
                    logger.info(String.format(" Titrating: %s", multiTerminus));
                }
            }
        }
        
        /* Create containers for titratables: MultiResidues for discrete, ExtendedVariables for continuous. */
        if (distribution == Distribution.CONTINUOUS) {
            esvSystem = new ExtendedSystem(mola);
			esvSystem.setConstantPh(pH);
            for (Residue res : chosenResidues) {
                MultiResidue multi = TitrationUtils.titratingMultiresidueFactory(mola, res);
                TitrationESV esv = new TitrationESV(esvSystem, multi);
                titratingESVs.add(esv);
                for (Residue background : multi.getInactive()) {
                    inactivateResidue(background);
                }
                esvSystem.addVariable(esv);
            }
            ffe.attachExtendedSystem(esvSystem);
            logger.info(format(" Continuous pHMD readied with %d residues.", titratingESVs.size()));
        } else {
            for (Residue res : chosenResidues) {
                // Create MultiResidue objects to wrap titratables.
               MultiResidue multiRes = new MultiResidue(res, ff, ffe);
               Polymer polymer = findResiduePolymer(res, mola);
               polymer.addMultiResidue(multiRes);
               recursiveMap(res, multiRes);

               // Switch back to the original form and ready the ForceFieldEnergy.
               multiRes.setActiveResidue(res);
               reInitialize(true, false);
               titratingMultis.add(multiRes);
               logger.info(String.format(" Titrating: %s", multiRes));
            }
            logger.info(format(" Discrete MCMD readied with %d residues.", titratingMultis.size()));
        }
        
        switch (distribution) {
            default:
            case DISCRETE:
                molDyn.setMonteCarloListener(this, MonteCarloNotification.EACH_STEP);
                break;
            case CONTINUOUS:
                ffe.attachExtendedSystem(esvSystem);
                molDyn.attachExtendedSystem(esvSystem, 100);
                break;
        }
    }

    /**
     * Finds Titration definitions for the given Residue and adds them to the given MultiResidue.
     * For three-state transitions, simply populate the enumeration with multiple titrations
     * from a shared state and this will include them in MultiResidue construction.
     */
    private void recursiveMap(Residue member, MultiResidue multiRes) {
        // Map titrations for this member.
        Titration[] titrations = Titration.multiLookup(member);
        titrationMap.put(member, Arrays.asList(titrations));

        // For each titration, check whether it needs added as a MultiResidue option.
        for (Titration titration : titrations) {
            // Allow manual override of Histidine treatment.
            if ((titration.deprotForm == AminoAcid3.HID && config.histidineMode == HistidineMode.HIE_ONLY)
                    || (titration.deprotForm == AminoAcid3.HIE && config.histidineMode == HistidineMode.HID_ONLY)) {
                continue;
            }
            // Find all the choices currently available to this MultiResidue.
            List<AminoAcid3> choices = new ArrayList<>();
            for (Residue choice : multiRes.getConsideredResidues()) {
                choices.add(choice.getAminoAcid3());
            }
            // If this Titration target is not a choice for the MultiResidue, then add it.
            if (!choices.contains(titration.protForm) || !(choices.contains(titration.deprotForm))) {
                String targetName = (member.getAminoAcid3() == titration.protForm)
                        ? titration.deprotForm.toString() : titration.protForm.toString();
                int resNumber = member.getResidueNumber();
                ResidueType resType = member.getResidueType();
                Residue newChoice = new Residue(targetName, resNumber, resType);
                multiRes.addResidue(newChoice);
                titrationMap.put(newChoice, Arrays.asList(Titration.multiLookup(newChoice)));
            }
        }
    }

    private void meltdown() {
        writeSnapshot(".meltdown-");
        ffe.energy(false, true);
        if (ffe.getBondEnergy() > 1000) {
            for (ROLS rols : previousTarget.getBondList()) {
                ((Bond) rols).log();
            }
        }
        if (ffe.getAngleEnergy() > 1000) {
            for (ROLS rols : previousTarget.getAngleList()) {
                ((Angle) rols).log();
            }
        }
        if (ffe.getVanDerWaalsEnergy() > 1000) {
            for (Atom a1 : previousTarget.getAtomList()) {
                for (Atom a2 : mola.getAtomArray()) {
                    if (a1 == a2 || a1.getBond(a2) != null) {
                        continue;
                    }
                    double dist = sqrt(
                            pow((a1.getX()-a2.getX()),2) +
                            pow((a1.getY()-a2.getY()),2) +
                            pow((a1.getZ()-a2.getZ()),2));
                    if (dist < 0.8*(a1.getVDWR() + a2.getVDWR())) {
                        logger.warning(String.format("Close vdW contact for atoms: \n   %s\n   %s", a1, a2));
                    }
                }
            }
        }
        logger.severe(String.format("Temperature above critical threshold: %f", thermostat.getCurrentTemperature()));
    }
    
    /**
     * Wraps reinitialization calls so as to provide un-fucked atom numbering.
     */
    private void reInitialize(boolean initFFE, boolean initMolDyn) {
//        renumberAtoms(mola);	// TODO Determine if+why this is necessary.
        if (initFFE) {
            ffe.reInit();
        }
        if (initMolDyn) {
            molDyn.reInit();
        }
    }

    public void setDynamicsLauncher(DynamicsLauncher launcher) {
        this.mdLauncher = launcher;
    }
    public void setDynamicsLauncher(Object[] opt) {
        this.mdOptions = opt;
    }

    public void launch(MolecularDynamics md, int mcStepFrequency, double timeStep,
            double printInterval, double saveInterval,
            double temperature, boolean initVelocities, File dyn) {
        this.molDyn = md;
        this.mcStepFrequency = mcStepFrequency;
    }

    /**
     * The primary driver. Called by the MD engine at each dynamics step.
     */
    @Override
    public boolean mcUpdate(double temperature) {
        startTime = System.nanoTime();
        if (thermostat.getCurrentTemperature() > config.meltdownTemperature) {
            meltdown();
		}
        
        if (thermostat.getCurrentTemperature() > config.warningTemperature) {
            Atom[] atoms = mola.getAtomArray();
            logger.info(format(" System heating! Dumping atomic velocities for %d D.o.F.:", ffe.getNumberOfVariables()));
            double[] velocity = new double[3];
            for (Atom atom : atoms) {
                atom.getVelocity(velocity);
                logger.info(format(" %s: %s", atom.describe(Atom.Descriptions.Trim), Arrays.toString(velocity)));
            }
        }
		esvSystem.setTemperature(temperature);
        
        propagateInactiveResidues(titratingMultis);
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
            if (config.logTimings) {
                long took = System.nanoTime() - startTime;
                logger.info(String.format(" CpHMD propagation time: %.6f", took * NS_TO_SEC));
            }
            return false;
        }

        logger.info(format("TitratingMultis: %d", titratingMultis.size()));

        // Randomly choose a target titratable residue to attempt protonation switch.
        int random = (config.titrateTermini)
                ? rng.nextInt(titratingMultis.size() + titratingTermini.size())
                : rng.nextInt(titratingMultis.size());

        if (random >= titratingMultis.size()) {
            Residue target = titratingTermini.get(random - titratingMultis.size());
            boolean accepted = tryTerminusTitration((MultiTerminus) target);
            snapshotIndex++;
            if (accepted) {
                molDyn.reInit();
                previousTarget = target;
            }
            return accepted;
        }
        MultiResidue targetMulti = titratingMultis.get(random);

        // Check whether rotamer moves are possible for the selected residue.
        Residue targetMultiActive = targetMulti.getActive();
        Rotamer[] targetMultiRotamers = targetMultiActive.getRotamers(library);
        if (targetMultiRotamers == null || targetMultiRotamers.length <= 1) {
            if (stepType == StepType.ROTAMER) {
                return false;
            } else if (stepType == StepType.COMBO) {
                stepType = StepType.TITRATE;
            }
        }

        // Perform the MC move.
        boolean accepted;
        switch (stepType) {
            case TITRATE:
                accepted = tryTitrationStep(targetMulti);
                break;
            case ROTAMER:
                accepted = (config.useConformationalBias) 
                        ? tryCBMCStep(targetMulti)
                        : tryRotamerStep(targetMulti);
                break;
            case COMBO:
                accepted = (config.useConformationalBias)
                        ? tryCBMCStep(targetMulti) || tryTitrationStep(targetMulti)
                        : tryComboStep(targetMulti);
                break;
            default:
                accepted = false;
                throw new IllegalStateException();
        }

        snapshotIndex++;
        if (accepted) {
            previousTarget = targetMulti;
        }
        if (config.logTimings) {
            long took = System.nanoTime() - startTime;
            logger.info(String.format(" CpHMD step time:        %.6f", took * NS_TO_SEC));
        }
        return accepted;
    }

    private void log() {
        logger.info(discountLogger.toString());
        discountLogger = new StringBuilder();
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
        ffe.energy(false, false);
        final double previousElectrostaticEnergy = ffe.getElectrostaticEnergy();

        // Write the pre-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, config.snapshots);
        String startString = target.toString();
        String startName = target.getName();

        double pKaref = 0;
        double dG_ref = 0;
        Titration titration = null;
        final TitrationType titrationType;

        if (terminus) {
            if (targetTerm.end == MultiTerminus.END.NTERM) {
                pKaref = 10.0;
                dG_ref = 0.0;
            } else {
                pKaref = 3.0;
                dG_ref = 0.0;
            }
            titrationType = targetTerm.titrateTerminus_v1(thermostat.getCurrentTemperature());
        } else {
            logger.info(format("targetMulti:  %s", targetMulti.toString()));
            logger.info(format("getActive:    %s", targetMulti.getActive().toString()));
            logger.info(format("titrationMap: %s", Arrays.toString(titrationMap.get(targetMulti.getActive()).toArray())));
            // Choose from the list of available titrations for the active residue.
            List<Titration> avail = titrationMap.get(targetMulti.getActive());
            titration = avail.get(rng.nextInt(avail.size()));

            // Perform the chosen titration.
            titrationType = performTitration(targetMulti, titration, config.inactivateBackground);
            reInitialize(true, true);

            // Test the MC criterion for a titration step.
            pKaref = titration.pKa;
            dG_ref = titration.refEnergy;
        }
        // Write the post-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, config.snapshots);

        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        ffe.energy(false, false);
        final double currentElectrostaticEnergy = ffe.getElectrostaticEnergy();
        final double dG_elec = currentElectrostaticEnergy - previousElectrostaticEnergy;

        if (config.zeroReferenceEnergies) {
            dG_ref = 0.0;
        }
        if (config.refOverride.isPresent()) {
            dG_ref = config.refOverride.getAsDouble();
        }

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the
         * reference compound
         */
        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (titrationType == TitrationType.DEPROT) {
            prefix = -prefix;
        }
        // Change this to use a single value for reference and then switch based on reaction.
        // Either positive ref == deprotonation or == standard -> nonstandard transition.
        if (titrationType == TitrationType.PROT) {
            dG_ref = -dG_ref;
        }
        double postfix = dG_elec - dG_ref;
        double dG_MC = prefix + postfix;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC protonation step:\n"));
        sb.append(String.format("     %s --> %s\n", startString, target.toString()));
        sb.append(String.format("     dG_ref:  %7.2f                pKaref:  %7.2f\n", dG_ref, pKaref));
        sb.append(String.format("     pH_term: %9.4f              elec_term: %10.4f\n", prefix, postfix));
        sb.append(String.format("     dG_elec: %9.4f              dG_MC:     %10.4f\n", dG_elec, dG_MC));
        sb.append(String.format("     -----\n"));

        // Test Monte-Carlo criterion.
        if (dG_MC < 0 && config.mcOverride != MCOverride.REJECT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        double criterion = exp(-dG_MC / kT);
        double metropolis = random();
        sb.append(String.format("     crit:    %9.4f              rng:       %10.4f\n", criterion, metropolis));
        if ((metropolis < criterion && config.mcOverride != MCOverride.REJECT) || config.mcOverride == MCOverride.ACCEPT) {
            numMovesAccepted++;
            molDyn.reInit();
            long took = System.nanoTime() - startTime;
            sb.append(String.format("     Accepted!                                                %1.3f", took * NS_TO_SEC));
            logger.info(sb.toString());
            return true;
        }

        // Move was rejected, undo the titration state change.
        performTitration(targetMulti, titration, config.inactivateBackground);
        reInitialize(true, true);
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
        if (CAUTIOUS) {
            throw new UnsupportedOperationException();
        }
        
        // Record the pre-change electrostatic energy.
        double previousElectrostaticEnergy = currentElectrostaticEnergy();

        // Write the pre-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, config.snapshots);
        String startString = target.toString();
        String startName = target.getName();

        double pKaref = 0;
        double dG_ref = 0;
        Titration titration = null;
        TitrationType type = null;

        if (target.end == MultiTerminus.END.NTERM) {
            pKaref = 8.23;
            dG_ref = 85.4929;
            type = target.isCharged ? TitrationType.DEPROT : TitrationType.PROT;
        } else if (target.end == MultiTerminus.END.CTERM) {
            pKaref = 3.55;
            dG_ref = -61.3825;
            type = target.isCharged ? TitrationType.PROT : TitrationType.DEPROT;
        }
        boolean beganCharged = target.isCharged;
        target.titrateTerminus_v1(thermostat.getCurrentTemperature());
        // Write the post-titration change snapshot.
        writeSnapshot(true, StepType.TITRATE, config.snapshots);

        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;

        if (config.zeroReferenceEnergies) {
            dG_ref = 0.0;
        }
        if (config.refOverride.isPresent()) {
            dG_ref = config.refOverride.getAsDouble();
        }

        /**
         * dG_elec = electrostatic energy component of the titratable residue
         * dG_ref = electrostatic component of the transition energy for the
         * reference compound
         */
        double pHterm = Math.log(10) * kT * (pH - pKaref);
        if (type == TitrationType.DEPROT) {
            pHterm = -pHterm;
        }
        // Change this to use a single value for reference and then switch based on reaction.
        // Either positive ref == deprotonation or == standard -> nonstandard transition.
        if (type == TitrationType.PROT) {
            dG_ref = -dG_ref;
        }
        double ddGterm = dG_elec - dG_ref;
        double dG_MC = pHterm + ddGterm;

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
        sb.append(String.format("     pH_term: %9.4f              elec_term: %10.4f\n", pHterm, ddGterm));
        sb.append(String.format("     dG_elec: %9.4f              dG_MC:     %10.4f\n", dG_elec, dG_MC));
        sb.append(String.format("     -----\n"));

        // Test Monte-Carlo criterion.
        if (dG_MC < 0 && config.mcOverride != MCOverride.REJECT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        double criterion = exp(-dG_MC / kT);
        double metropolis = random();
        sb.append(String.format("     crit:    %9.4f              rng:       %10.4f\n", criterion, metropolis));
        if ((metropolis < criterion && config.mcOverride != MCOverride.REJECT)
                || config.mcOverride == MCOverride.ACCEPT
                || config.mcOverride == MCOverride.ONCE) {
            numMovesAccepted++;
            long took = System.nanoTime() - startTime;
            sb.append(String.format("     Accepted!                                                %1.3f", took * NS_TO_SEC));
            logger.info(sb.toString());
            if (config.mcOverride == MCOverride.ONCE) {
                config.mcOverride = MCOverride.REJECT;
            }
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
        if (CAUTIOUS) {
            throw new UnsupportedOperationException();
        }
        
        List<Residue> targets = new ArrayList<>();
        targets.add(targetMulti.getActive());
        int trialSetSize = 5;   // cost still scales with this, unfortunately
        int mcFrequency = 1;    // irrelevant for manual step call
        boolean writeSnapshots = false;
        System.setProperty("cbmc-type", "CHEAP");
        RosenbluthCBMC cbmc = new RosenbluthCBMC(mola, mola.getPotentialEnergy(), null,
            targets, mcFrequency, trialSetSize, writeSnapshots);
        boolean accepted = cbmc.cbmcStep();
        if (config.logTimings) {
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
        if (CAUTIOUS) {
            throw new UnsupportedOperationException();
        }
        
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();

        // Write the before-step snapshot.
        writeSnapshot(true, StepType.ROTAMER, config.snapshots);

        // Save coordinates so we can return to them if move is rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        ResidueState origState = residue.storeState();
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origState, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);
        // Select a new rotamer and swap to it.
        Rotamer[] rotamers = residue.getRotamers(library);
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);

        // Write the post-rotamer change snapshot.
        writeSnapshot(false, StepType.ROTAMER, config.snapshots);

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
            propagateInactiveResidues(titratingMultis, true);
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
                propagateInactiveResidues(titratingMultis, true);
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
        if (CAUTIOUS) {
            throw new UnsupportedOperationException();
        }
        
        // Record the pre-change total energy.
        double previousTotalEnergy = currentTotalEnergy();
        double previousElectrostaticEnergy = currentElectrostaticEnergy();

        // Write the pre-combo snapshot.
        writeSnapshot(true, StepType.COMBO, config.snapshots);
        String startString = targetMulti.toString();
        String startName = targetMulti.getActive().getName();

        // Choose from the list of available titrations for the active residue.
        List<Titration> avail = titrationMap.get(targetMulti.getActive());
        Titration titration = avail.get(rng.nextInt(avail.size()));

        // Perform the chosen titration.
        TitrationType titrationType = performTitration(targetMulti, titration, config.inactivateBackground);

        // Change rotamer state, but first save coordinates so we can return to them if rejected.
        Residue residue = targetMulti.getActive();
        ArrayList<Atom> atoms = residue.getAtomList();
        ResidueState origState = residue.storeState();
        double chi[] = new double[4];
        RotamerLibrary.measureAARotamer(residue, chi, false);
        AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
        Rotamer origCoordsRotamer = new Rotamer(aa, origState, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);

        // Swap to the new rotamer.
        Rotamer[] rotamers = residue.getRotamers(library);
        int rotaRand = rng.nextInt(rotamers.length);
        RotamerLibrary.applyRotamer(residue, rotamers[rotaRand]);

        // Write the post-combo snapshot.
        writeSnapshot(false, StepType.COMBO, config.snapshots);

        // Evaluate both MC criteria.
        String endName = targetMulti.getActive().getName();

        // Evaluate the titration probability of the step.
        double pKaref = titration.pKa;
        double dG_ref = titration.refEnergy;
        double temperature = thermostat.getCurrentTemperature();
        double kT = BOLTZMANN * temperature;
        double dG_elec = currentElectrostaticEnergy() - previousElectrostaticEnergy;

        if (config.zeroReferenceEnergies) {
            dG_ref = 0.0;
        }

        double prefix = Math.log(10) * kT * (pH - pKaref);
        if (titrationType == TitrationType.DEPROT) {
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
        if (dG_titr < 0 && dG_rota < 0 && config.mcOverride != MCOverride.REJECT) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            propagateInactiveResidues(titratingMultis, false);
            return true;
        } else {
            // Conditionally accept based on combined probabilities.
            if (dG_titr < 0 || config.mcOverride == MCOverride.ACCEPT) {
                titrCriterion = 1.0;
            }
            if (dG_rota < 0) {
                rotaCriterion = 1.0;
            }
            if (config.mcOverride == MCOverride.REJECT) {
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
                propagateInactiveResidues(titratingMultis, false);
                return true;
            } else {
                // Move was denied.
                sb.append(String.format("     Denied."));
                logger.info(sb.toString());

                // Undo both pieces of the rejected move IN THE RIGHT ORDER.
                RotamerLibrary.applyRotamer(residue, origCoordsRotamer);
                performTitration(targetMulti, titration, config.inactivateBackground);
                ffe.reInit();
                molDyn.reInit();
                return false;
            }
        }
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
     * Calculates the electrostatic energy at the current state.
     *
     * @return Energy of the current state.
     */
    private double currentElectrostaticEnergy() {
        double x[] = new double[ffe.getNumberOfVariables() * 3];
        ffe.getCoordinates(x);
        ffe.energy(x);
        return ffe.getTotalElectrostaticEnergy();
    }
    
    /**
     * Calculates the total energy at the current state.
     *
     * @return Energy of the current state.
     */
    private double currentTotalEnergy() {
        double x[] = new double[ffe.getNumberOfVariables() * 3];
        ffe.getCoordinates(x);
        ffe.energy(x);
        return ffe.getTotalEnergy();
    }

    private void writeSnapshot(String extension) {
        String filename = FilenameUtils.removeExtension(mola.getFile().toString()) + extension + snapshotIndex;
        if (config.snapshots == Snapshots.INTERLEAVED) {
            filename = mola.getFile().getAbsolutePath();
            if (!filename.contains("dyn")) {
                filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
            }
        }
        File file = new File(filename);
        PDBFilter writer = new PDBFilter(file, mola, null, null);
        writer.writeFile(file, false);
    }

    private void writeSnapshot(boolean beforeChange, StepType stepType, Snapshots snapshotsType) {
        // Write the after-step snapshot.
        if (snapshotsType != Snapshots.NONE) {
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
            if (snapshotsType == Snapshots.INTERLEAVED) {
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

    public class DynamicsLauncher {
        private final int nSteps;
        private final boolean initVelocities;
        private final double timeStep, print, save, restart, temperature;
        private final String fileType;
        private final File dynFile;
        public DynamicsLauncher(MolecularDynamics md, Object[] opt) {
            molDyn = md;
            nSteps = (int) opt[0];
            timeStep = (double) opt[1]; print = (double) opt[2];
            save = (double) opt[3]; restart = (double) opt[4]; initVelocities = (boolean) opt[5];
            fileType = (String) opt[6]; temperature = (double) opt[7];  dynFile = (File) opt[8];
        }
        public DynamicsLauncher(MolecularDynamics md,
                int nSteps, double timeStep,
                double print, double save,
                double temperature, boolean initVelocities,
                String fileType, double restart,
                File dynFile) {
            molDyn = md;
            this.nSteps = nSteps; this.initVelocities = initVelocities;
            this.timeStep = timeStep; this.print = print; this.save = save; this.restart = restart;
            this.temperature = temperature; this.fileType = fileType; this.dynFile = dynFile;
        }
        public DynamicsLauncher(MolecularDynamics md,
                int nSteps, double timeStep,
                double print, double save,
                double temperature, boolean initVelocities,
                File dynFile) {
            molDyn = md;
            this.nSteps = nSteps; this.initVelocities = initVelocities;
            this.timeStep = timeStep; this.print = print; this.save = save; this.restart = 0.1;
            this.temperature = temperature; this.fileType = "PDB"; this.dynFile = dynFile;
        }
        public void launch() {
            launch(nSteps);
        }
        public void launch(int nSteps) {
            /* For reference:
            molDyn.init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);
            molDyn.dynamic(nSteps, timeStep,
                  printInterval, saveInterval,
                  temperature, initVelocities,
                  fileType, restartFrequency, dyn);   */
            discountLogger.append(format("    dynamic launcher parameters: %d %g %g %g %g %s %g\n",
                    mdOptions[0], mdOptions[1], mdOptions[2], mdOptions[3],
                    mdOptions[4], mdOptions[5], mdOptions[7]));
//            discountLogger.append(format("    terminating current md...\n"));
//            log();
//            molDyn.terminate();
            discountLogger.append(format("    launching new md process...\n"));
            log();
            molDyn.dynamic(nSteps, timeStep, print, save, temperature, initVelocities, fileType, restart, dynFile);
        }
    }

    private enum StepType {
        TITRATE, ROTAMER, COMBO, CONTINUOUS_DYNAMICS;
    }
    
    private void debug(int level, String message) {
        if (config.debugLogLevel >= level) {
            logger.info(message);
        }
    }
}
