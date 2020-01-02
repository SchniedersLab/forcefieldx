//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.mc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Torsion;
import ffx.potential.parsers.PDBFilter;
import static ffx.utilities.Constants.R;

/**
 * Orientational Biased Monte Carlo (as applied to chi0 torsion of peptide side-chains.
 * <p>
 * As described by Frenkel/Smit in "Understanding Molecular Simulation" Chapter 13.1.2.
 * This uses the "orientational biasing" method to select chi[0] moves that are frequently accepted.
 *
 * @author Stephen D. LuCore
 */
public class RosenbluthOBMC implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(RosenbluthOBMC.class.getName());

    /**
     * MolecularAssembly to operate on.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * Force field energy to use.
     */
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
     * Size of the trial sets, k.
     */
    private final int trialSetSize;
    /**
     * Counters for proposed and accepted moves.
     */
    private int numMovesProposed = 0;
    /**
     * Creates verbose output.
     */
    private StringBuilder report = new StringBuilder();
    /**
     * Writes PDBs of each trial set and original/proposed configurations.
     */
    private boolean writeSnapshots = false;

    /**
     * RRMC constructor.
     *
     * @param targets           Residues to undergo RRMC.
     * @param mcFrequency       Number of MD steps between RRMC proposals.
     * @param trialSetSize      Larger values cost more but increase acceptance.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forceFieldEnergy  a {@link ffx.potential.ForceFieldEnergy} object.
     * @param thermostat        a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
     */
    public RosenbluthOBMC(MolecularAssembly molecularAssembly, ForceFieldEnergy forceFieldEnergy, Thermostat thermostat,
                          List<Residue> targets, int mcFrequency, int trialSetSize) {
        this.targets = targets;
        this.mcFrequency = mcFrequency;
        this.trialSetSize = trialSetSize;
        this.molecularAssembly = molecularAssembly;
        this.forceFieldEnergy = forceFieldEnergy;
        this.thermostat = thermostat;
    }

    /**
     * <p>Constructor for RosenbluthOBMC.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param ffe               a {@link ffx.potential.ForceFieldEnergy} object.
     * @param thermostat        a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
     * @param targets           a {@link java.util.List} object.
     * @param mcFrequency       a int.
     * @param trialSetSize      a int.
     * @param writeSnapshots    a boolean.
     */
    public RosenbluthOBMC(MolecularAssembly molecularAssembly, ForceFieldEnergy ffe, Thermostat thermostat,
                          List<Residue> targets, int mcFrequency, int trialSetSize, boolean writeSnapshots) {
        this(molecularAssembly, ffe, thermostat, targets, mcFrequency, trialSetSize);
        this.writeSnapshots = writeSnapshots;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean mcUpdate(double temperature) {
        steps++;
        if (steps % mcFrequency == 0) {
            return mcStep();
        }
        return false;
    }

    /**
     * Does all the work for a move. Moveset is a continuous 360 degree spin of
     * the chi[0] torsion. U_or in Frenkel's notation (uDep here) is the
     * associated torsion energy. Evaluation criterion: P_accept = Min( 1,
     * (Wn/Wo)*exp(-beta(U[n]-U[o]) )
     */
    private boolean mcStep() {
        numMovesProposed++;
        boolean accepted;

        // Select a target residue.
        int index = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(index);
        ResidueState origState = target.storeState();
        Torsion chi0 = getChiZeroTorsion(target);
        writeSnapshot("orig");

        /*
            Create old and new trial sets, calculate Wn and Wo, and choose a move bn.
            When doing strictly chi[0] moves, Frenkel/Smit's 'old' and 'new' configurations
            are the same state. The distinction is made here only to aid in future generalization.
         */
        List<MCMove> oldTrialSet = createTrialSet(target, origState, trialSetSize - 1);
        List<MCMove> newTrialSet = createTrialSet(target, origState, trialSetSize);
        report = new StringBuilder();
        report.append(format(" Rosenbluth Rotamer MC Move: %4d\n", numMovesProposed));
        report.append(format("    residue:   %s\n", target.toString()));
        report.append(format("    chi0:      %s\n", chi0.toString()));
        MCMove proposal = calculateRosenbluthFactors(target, chi0,
                origState, oldTrialSet, origState, newTrialSet);

        /*
            Calculate the independent portion of the total old-conf energy.
            Then apply the move and calculate the independent total new-conf energy.
         */
        setState(target, origState);
        writeSnapshot("uIndO");
        double uIndO = getTotalEnergy() - getTorsionEnergy(chi0);
        proposal.move();
        writeSnapshot("uIndN");
        double uIndN = getTotalEnergy() - getTorsionEnergy(chi0);

        // Apply acceptance criterion.
        double temperature = thermostat.getCurrentTemperature();
        double beta = 1.0 / (R * temperature);
        double dInd = uIndN - uIndO;
        double dIndE = exp(-beta * dInd);
        double criterion = (Wn / Wo) * exp(-beta * (uIndN - uIndO));
        double metropolis = min(1, criterion);
        double rng = ThreadLocalRandom.current().nextDouble();

        report.append(format("    theta:     %3.2f\n", ((RosenbluthChi0Move) proposal).theta));
        report.append(format("    criterion: %1.4f\n", criterion));
        report.append(format("       Wn/Wo:     %.2f\n", Wn / Wo));
        report.append(format("       uIndN,O:  %7.2f\t%7.2f\n", uIndN, uIndO));
        report.append(format("       dInd(E):  %7.2f\t%7.2f\n", dInd, dIndE));
        report.append(format("    rng:       %1.4f\n", rng));
        if (rng < metropolis) {
            report.append(" Accepted.\n");
            accepted = true;
        } else {
            proposal.revertMove();
            report.append(" Denied.\n");
            accepted = false;
        }
        logger.info(report.toString());

        // Cleanup.
        Wn = 0.0;
        Wo = 0.0;
        return accepted;
    }

    /**
     * Generates trial movesets around new and old configurations. This involves
     * loading the move-dependent energy component into each member of the trial
     * sets.
     */
    private List<MCMove> createTrialSet(Residue target, ResidueState state, int setSize) {
        List<MCMove> moves = new ArrayList<>();
        // Trial set around old configuration is of size (k-1).
        setState(target, state);
        for (int i = 0; i < setSize; i++) {
            moves.add(new RosenbluthChi0Move(target));
        }
        return moves;
    }

    /**
     * Chooses a move, bn, from amongst the new trial set, {b}k, based on the
     * Boltzmann-weighted orientational energy, U_or[n]. Also calculates
     * Rosenbluth factors for both sets, Wn and Wo. Wn = sum(i=1:k, exp(-beta *
     * uDep[i])) Wo = exp(-beta * uDep[current]) + sum(i=2:k, exp(-beta *
     * uDep[i]))
     */
    private MCMove calculateRosenbluthFactors(Residue target, Torsion chi0,
                                              ResidueState oldConf, List<MCMove> oldTrialSet,
                                              ResidueState newConf, List<MCMove> newTrialSet) {
        double temperature = thermostat.getCurrentTemperature();
        double beta = 1.0 / (R * temperature);

        // Initialize and add up Wo.
        Wo = exp(-beta * getTorsionEnergy(chi0));
        report.append(format("    TestSet (Old): %5s\t%7s\t\t%7s\n", "uDepO", "uDepOe", "Sum(Wo)"));
        report.append(format("       Orig %d:   %7.4f\t%7.4f\t\t%7.4f\n",
                0, getTorsionEnergy(chi0), exp(-beta * getTorsionEnergy(chi0)), Wo));
        for (int i = 0; i < oldTrialSet.size(); i++) {
            setState(target, oldConf);
            MCMove move = oldTrialSet.get(i);
            move.move();
            double uDepO = getTorsionEnergy(chi0);
            double uDepOe = exp(-beta * uDepO);
            Wo += uDepOe;
            if (i < 5 || i >= oldTrialSet.size() - 5) {
                report.append(format("       Prop %d:   %7.4f\t%7.4f\t\t%7.4f\n", i + 1, uDepO, uDepOe, Wo));
                writeSnapshot("ots");
            } else if (i == 5) {
                report.append("        ... \n");
            }
        }

        // Initialize and add up Wn.  Record dependent energy of each trial set member.
        Wn = 0.0;
        double[] uDepN = new double[newTrialSet.size()];
        double[] uDepNe = new double[newTrialSet.size()];
        report.append(format("    TestSet (New): %5s\t%7s\t\t%7s\n", "uDepN", "uDepNe", "Sum(Wn)"));
        for (int i = 0; i < newTrialSet.size(); i++) {
            setState(target, newConf);
            MCMove move = newTrialSet.get(i);
            move.move();
            uDepN[i] = getTorsionEnergy(chi0);
            uDepNe[i] = exp(-beta * uDepN[i]);
            Wn += uDepNe[i];
            if (i < 5 || i >= newTrialSet.size() - 5) {
                report.append(format("       Prop %d:   %7.4f\t%7.4f\t\t%7.4f\n", i, uDepN[i], uDepNe[i], Wn));
                writeSnapshot("nts");
            } else if (i == 5) {
                report.append("        ... \n");
            }
        }
        setState(target, oldConf);

        // Choose a proposal move from the new trial set.
        MCMove proposal = null;
        double rng = ThreadLocalRandom.current().nextDouble(Wn);
        double running = 0.0;
        for (int i = 0; i < newTrialSet.size(); i++) {
            running += uDepNe[i];
            if (rng < running) {
                proposal = newTrialSet.get(i);
                double prob = uDepNe[i] / Wn * 100;
                report.append(format("       Chose %d   %7.4f\t%7.4f\t  %4.1f%%\n", i, uDepN[i], uDepNe[i], prob));
                break;
            }
        }
        if (proposal == null) {
            logger.severe("Programming error.");
        }
        return proposal;
    }

    private double getTotalEnergy() {
        double[] x = new double[forceFieldEnergy.getNumberOfVariables() * 3];
        forceFieldEnergy.getCoordinates(x);
        return forceFieldEnergy.energy(x);
    }

    private double getTorsionEnergy(Torsion torsion) {
        return torsion.energy(false);
    }

    private Torsion getChiZeroTorsion(Residue residue) {
        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
        ArrayList<Torsion> torsions = residue.getTorsionList();
        switch (name) {
            case VAL: {
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (Torsion torsion : torsions) {
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
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (Torsion torsion : torsions) {
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
                for (Torsion torsion : torsions) {
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
                for (Torsion torsion : torsions) {
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
                for (Torsion torsion : torsions) {
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
                for (Torsion torsion : torsions) {
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
                for (Torsion torsion : torsions) {
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

    /**
     * Calls through to residue.revertState() but also updates the Torsion
     * objects associated with that residue (so they contain appropriate chi
     * values).
     */
    private void setState(Residue target, ResidueState state) {
        target.revertState(state);
        for (Torsion torsion : target.getTorsionList()) {
            torsion.update();
        }
    }

    private void writeSnapshot(String suffix) {
        if (!writeSnapshots) {
            return;
        }
        String filename = FilenameUtils.removeExtension(molecularAssembly.getFile().toString()) + "." + suffix + "-" + numMovesProposed;
        File file = new File(filename);
        PDBFilter writer = new PDBFilter(file, molecularAssembly, null, null);
        writer.writeFile(file, false);
    }
}
