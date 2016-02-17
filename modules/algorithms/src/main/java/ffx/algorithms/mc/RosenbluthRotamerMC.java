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
package ffx.algorithms.mc;

import ffx.algorithms.MonteCarloListener;
import ffx.algorithms.Thermostat;
import ffx.algorithms.mc.MCMove;
import ffx.algorithms.mc.RosenbluthRotamerMove;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Torsion;
import ffx.potential.parsers.PDBFilter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;

/**
 * As described by Frenkel/Smit in "Understanding Molecular Simulation" Chapter 13.1.2.
 * This uses the "orientational biasing" method to select chi[0] moves that are frequently accepted.
 * @author S. LuCore
 */
public class RosenbluthRotamerMC implements MonteCarloListener {
    private static final Logger logger = Logger.getLogger(RosenbluthRotamerMC.class.getName());
    
    private final double BOLTZMANN = 0.0019872041; // In kcal/(mol*K)
    private final MolecularAssembly mola;
    private final ForceFieldEnergy forceFieldEnergy;
    private final Thermostat thermostat;
    /**
     * At each move, one of these residues will be chosen as the target.
     */
    private final List<Residue> targets;
    /**
     * Number of mcUpdate() calls (e.g. MD steps) between move proposals.
     */
    private final int mcFrequency;
    /**
     * Keeps track of calls to mcUpdate (e.g. MD steps).
     */
    private int steps = 0;
    /**
     * Rosenbluth factor for the forward move.
     */
    private double Wn;
    /**
     * Rosenbluth factor for the backward move.
     */
    private double Wo;
    /**
     * Trial moveset for new configuration, bn.
     */
    private List<MCMove> newTrialSet;
    /**
     * Trial moveset for old configuration, bo.
     */
    private List<MCMove> oldTrialSet;
    /**
     * Size of the trial sets, k.
     */
    private final int trialSetSize;
    /**
     * Counters for proposed and accepted moves.
     */
    private int numMovesProposed = 0;
    private int numMovesAccepted = 0;
    /**
     * Creates verbose output.
     */
    private StringBuilder report = new StringBuilder();
    
    public RosenbluthRotamerMC(MolecularAssembly mola, ForceFieldEnergy ffe, Thermostat thermostat, 
            List<Residue> targets, int mcFrequency, int trialSetSize) {
        this.targets = targets;
        this.mcFrequency = mcFrequency;
        this.trialSetSize = trialSetSize;
        this.mola = mola;
        this.forceFieldEnergy = ffe;
        this.thermostat = thermostat;
    }
    
    @Override
    public boolean mcUpdate(MolecularAssembly mola) {
        steps++;
        if (steps % mcFrequency == 0) {
            return mcStep();
        }
        return false;
    }
    
    /**
     * Does all the work for a move.
     * Moveset is a continuous 360 degree spin of the chi[0] torsion.
     * U_or in Frenkel's notation (uDep here) is the associated torsion energy.
     * Evaluation criterion: P_accept = Min( 1, (Wn/Wo)*exp(-beta(U[n]-U[o]) )
     */
    private boolean mcStep() {
        numMovesProposed++;
        boolean accepted;
        
        // Select a target residue.
        int which = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(which);
        ResidueState state = ResidueState.storeAllCoordinates(target);
        Torsion chi0 = getChiZeroTorsion(target);
        
        /* Create old and new trial sets, calculate Wn and Wo, and choose a move bn.
            When doing strictly chi[0] moves, Frenkel/Smit's 'old' and 'new' configurations
            are the same state. The distinction is made here only to aid in future generalization.
        */
        createTrialSets(target, state, state);
        report = new StringBuilder();
        report.append(String.format(" Rosenbluth Rotamer MC Move: \n"));
        report.append(String.format("    residue:   %s\n", target.toString()));
        report.append(String.format("    chi0:      %s\n", chi0.toString()));
        MCMove proposal = calculateRosenbluthFactors(target, state, state, chi0);
        
        /* Calculate the independent portion of the total old-conf energy.
            Then apply the move and calculate the independent total new-conf energy.
        */
        double uIndO = getTotalEnergy() - getTorsionEnergy(chi0);
        proposal.move();
        double uIndN = getTotalEnergy() - getTorsionEnergy(chi0);
        
        // Apply acceptance criterion.
        double temperature = thermostat.getCurrentTemperature();
        double beta = 1.0 / (BOLTZMANN * temperature);
        double dInd = uIndN - uIndO;
        double dIndE = FastMath.exp(-beta*dInd);
        double criterion = (Wn / Wo) * FastMath.exp(-beta*(uIndN - uIndO));
        double metropolis = Math.min(1, criterion);
        double rng = ThreadLocalRandom.current().nextDouble();
        
        report.append(String.format("    theta:     %3.2f\n", ((RosenbluthRotamerMove) proposal).theta));
        report.append(String.format("    criterion: %1.4f\n", criterion));
        report.append(String.format("       Wn/Wo:     %.2f\n", Wn/Wo));
        report.append(String.format("       uIndN,O:  %7.2f\t%7.2f\n", uIndN, uIndO));
        report.append(String.format("       dInd(E):  %7.2f\t%7.2f\n", dInd, dIndE));
        report.append(String.format("    rng:       %1.4f\n", rng));
        if (rng < metropolis) {
            numMovesAccepted++;
            report.append(String.format(" Accepted.\n"));
            accepted = true;
        } else {
            proposal.revertMove();
            report.append(String.format(" Denied.\n"));
            accepted = false;
        }
        logger.info(report.toString());
        
        // Cleanup.
        oldTrialSet = null;
        newTrialSet = null;
        Wn = 0.0;
        Wo = 0.0;
        return accepted;
    }
    
    /**
     * Generates trial movesets around new and old configurations.
     * This involves loading the move-dependent energy component
     * into each member of the trial sets.
     */
    private void createTrialSets(Residue target, ResidueState oldConf, ResidueState newConf) {
        List<MCMove> oldMoves = new ArrayList<>();
        List<MCMove> newMoves = new ArrayList<>();
        // Trial set around old configuration is of size (k-1).
        ResidueState.revertAllCoordinates(target, oldConf);
        for (int i = 0; i < trialSetSize - 1; i++) {
            oldMoves.add(new RosenbluthRotamerMove(target));
        }
        ResidueState.revertAllCoordinates(target, newConf);
        for (int i = 0; i < trialSetSize; i++) {
            newMoves.add(new RosenbluthRotamerMove(target));
        }
        oldTrialSet = oldMoves;
        newTrialSet = newMoves;
    }
    
    /**
     * Chooses a move, bn, from amongst the new trial set, {b}k, based on 
     * the Boltzmann-weighted orientational energy, U_or[n].
     * Also calculates Rosenbluth factors for both sets, Wn and Wo.
     * Wn = sum(i=1:k, exp(-beta * uDep[i]))
     * Wo = exp(-beta * uDep[current]) + sum(i=2:k, exp(-beta * uDep[i]))
     */
    private MCMove calculateRosenbluthFactors(Residue target,
            ResidueState oldConf, ResidueState newConf, Torsion chi0) {
        if (oldTrialSet == null || newTrialSet == null) {
            logger.severe("Improper function call.");
        }
        double temperature = thermostat.getCurrentTemperature();
        double beta = 1.0 / (BOLTZMANN * temperature);
        
        // Initialize and add up Wo.
        Wo = FastMath.exp(-beta * getTorsionEnergy(chi0));
        report.append(String.format("    TestSet (Old): %5s\t%7s\t\t%7s\n", "uDepO", "uDepOe", "Sum(Wo)"));
        for (int i = 0; i < oldTrialSet.size(); i++) {
            ResidueState.revertAllCoordinates(target, oldConf);
            MCMove move = oldTrialSet.get(i);
            move.move();
            double uDepO = getTorsionEnergy(chi0);
            double uDepOe = FastMath.exp(-beta*uDepO);
            Wo += uDepOe;
            report.append(String.format("       Prop %d:   %7.4f\t%7.4f\t\t%7.4f\n", i, uDepO, uDepOe, Wo));
        }
        
        // Initialize and add up Wn.  Record dependent energy of each trial set member.
        Wn = 0.0;
        double uDepN[] = new double[newTrialSet.size()];
        double uDepNe[] = new double[newTrialSet.size()];
        report.append(String.format("    TestSet (New): %5s\t%7s\t\t%7s\n", "uDepN", "uDepNe", "Sum(Wn)"));
        for (int i = 0; i < newTrialSet.size(); i++) {
            ResidueState.revertAllCoordinates(target, newConf);
            MCMove move = newTrialSet.get(i);
            move.move();
            uDepN[i] = getTorsionEnergy(chi0);
            uDepNe[i] = FastMath.exp(-beta*uDepN[i]);
            Wn += uDepNe[i];
            report.append(String.format("       Prop %d:   %7.4f\t%7.4f\t\t%7.4f\n", i, uDepN[i], uDepNe[i], Wn));
        }
        ResidueState.revertAllCoordinates(target, oldConf);
        
        // Choose a proposal move from the new trial set.
        MCMove proposal = null;
        double rng = ThreadLocalRandom.current().nextDouble(Wn);
        double running = 0.0;
        for (int i = 0; i < newTrialSet.size(); i++) {
            running += uDepNe[i];
            if (rng < running) {
                proposal = newTrialSet.get(i);
                double prob = uDepNe[i] / Wn * 100;
                report.append(String.format("       Chose %d   %7.4f\t%7.4f\t  %4.1f%%\n", i, uDepN[i], uDepNe[i], prob));
                break;
            }
        }
        if (proposal == null) {
            logger.severe("Programming error.");
        }
        return proposal;
    }
    
    private double getTotalEnergy() {
        double x[] = new double[forceFieldEnergy.getNumberOfVariables() * 3];
        forceFieldEnergy.getCoordinates(x);
        return forceFieldEnergy.energy(x);
    }
    
    private double getTorsionEnergy(Torsion torsion) {
        return torsion.energy(false);
    }
    
    private Torsion getChiZeroTorsion(Residue residue) {
        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
        ArrayList<ROLS> torsions = residue.getTorsionList();
        switch (name) {
            case VAL: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG1)) {
                        return torsion;
                    }
                }
                break;
            }
            case ILE: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG1)) {
                        return torsion;
                    }
                }
                break;
            }
            case SER: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom OG = (Atom) residue.getAtomNode("OG");
                Atom HG = (Atom) residue.getAtomNode("HG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, OG)) {
                        return torsion;
                    }
                }
                break;
            }
            case THR: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom OG1 = (Atom) residue.getAtomNode("OG1");
                Atom HG1 = (Atom) residue.getAtomNode("HG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, OG1)) {
                        return torsion;
                    }
                }
                break;
            }
            case CYX: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom SG = (Atom) residue.getAtomNode("SG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, SG)) {
                        return torsion;
                    }
                }
                break;
            }
            case CYD: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom SG = (Atom) residue.getAtomNode("SG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, SG)) {
                        return torsion;
                    }
                }
                break;
            }
            default: {  // All other residues' chi[0] are defined by N,CA,CB,CG.
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        return torsion;
                    }
                }
                logger.info("Couldn't find chi[0] for residue " + residue.toString());
                return null;
            }
        }
        logger.info("Couldn't find chi[0] for residue " + residue.toString());
        return null;
    }
    
    private final boolean snapshotInterleaving = false;
    private void writeSnapshot(String suffix) {
        String filename = FilenameUtils.removeExtension(mola.getFile().toString()) + "." + suffix + "-" + numMovesProposed;
        if (snapshotInterleaving) {
            filename = mola.getFile().getAbsolutePath();
            if (!filename.contains("dyn")) {
                filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
            }
        }
        File file = new File(filename);
        PDBFilter writer = new PDBFilter(file, mola, null, null);
        writer.writeFile(file, false);
    }
}
