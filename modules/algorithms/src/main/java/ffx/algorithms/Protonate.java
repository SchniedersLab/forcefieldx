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
import ffx.potential.parameters.ForceField;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;

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
        String debugLogLevel = System.getProperty("debug");
        if (debugLogLevel != null) {
            this.debugLogLevel = Integer.parseInt(debugLogLevel);
        }

        //initialize stepcount and the number of accepted moves
        numMovesAccepted = 0;

        this.molAss = molAss;
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
            multiRes.setActiveResidue(res);
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
        propagateInactiveResidues(titratingResidues);
        
        stepCount++;
        if (stepCount % mcStepFrequency != 0) {
            return false;
        }        
        double referenceEnergy = molAss.getPotentialEnergy().getTotalEnergy();
        
        // Randomly choose a target titratable residue to attempt protonation switch.
        int random = rng.nextInt(titratingResidues.size());
//        MultiResidue targetMulti = titratingResidues.get(random);
//        String startingName = targetMulti.toString();
        MultiResidue targetMulti = null;
        
        for (MultiResidue mR : titratingResidues) {
            if (mR.getActive().getResidueNumber() == 27) {
                targetMulti = mR;
            }
        }
        String startingName = targetMulti.toString();
        
        // Switch titration state for chosen residue.
        switchProtonationState(targetMulti);

        String name = targetMulti.getActive().getName();
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
        sb.append(String.format("     %s --> %s\n", startingName, targetMulti.toString()));
        sb.append(String.format("     pKaref:  %7.2f\n", pKaref));
        sb.append(String.format("     dG_ref:  %7.2f\n", dG_ref));
        sb.append(String.format("     dG_elec: %9.4f\n", dG_elec));
        sb.append(String.format("     dG_MC:   %9.4f\n", dG_MC));
        sb.append(String.format("     -----\n"));
        
        String acceptMC = System.getProperty("acceptMC");
        if (acceptMC != null && acceptMC.equalsIgnoreCase("true")) {
            sb.append("     Accept override: ");
            dG_MC = -1;
        }
        
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
        switchProtonationState(targetMulti);
        return false;
    }

    // TESTING ONLY
    private MolecularDynamics molDyn;
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
     * Switch the protonation state of target residue and reinitialize FF.
     *
     * @param residue
     */
    private void switchProtonationState(MultiResidue multiRes) {
        Atom oldAtomArray[] = new Atom[molAss.getAtomArray().length];
        System.arraycopy(oldAtomArray, 0, molAss.getAtomArray(), 0, oldAtomArray.length);
        for (int i = 0; i < oldAtomArray.length; i++) {
            oldAtomArray[i] = molAss.getAtomArray()[i];
        }
        
        String protFormName = Titratable.valueOf(multiRes.getActive().getName()).protForm.toString();
        String deprotFormName = Titratable.valueOf(multiRes.getActive().getName()).deprotForm.toString();
        if (multiRes.getName().equalsIgnoreCase(protFormName)) {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(deprotFormName));
        } else {
            multiRes.requestSetActiveResidue(AminoAcid3.valueOf(protFormName));
        }
        forceFieldEnergy.reInit();
        molDyn.reInit(oldAtomArray);
        
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
