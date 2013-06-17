/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */
package ffx.algorithms;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.Rotamer;
import ffx.potential.RotamerLibrary;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;

/**
 * Optimize protein side-chain conformations using rotamers.
 *
 * @author Michael J. Schnieders
 */
public class RotamerOptimization implements Terminatable {

    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());

    public enum Algorithm {

        INDEPENDENT, GLOBAL, SLIDING_WINDOW
    }

    public enum Direction {

        FORWARD, BACKWARD
    }
    MolecularAssembly molecularAssembly;
    AlgorithmListener algorithmListener;
    ForceFieldEnergy potential;
    Polymer[] polymers;
    int startResID = -1;
    int finalResID = -1;
    int windowSize = 3;
    Algorithm algorithm = null;
    Direction direction = null;
    boolean terminate = false;
    boolean done = true;
    /**
     * Self-energy of each residue for each rotamer. [residue][rotamer]
     */
    double selfEnergy[][];
    /**
     * Pair-energies for each pair of residue and pair of rotamers.
     * [residue1][rotamer1][residue2][rotamer2]
     */
    double pairEnergy[][][][];
    /**
     * Eliminated rotamers. [residue][rotamer]
     */
    boolean eliminatedRotamers[][];
    /**
     * Eliminated rotamer pairs. [residue1][rotamer1][residue2][rotamer2]
     */
    boolean eliminatedRotamerPairs[][][][];
    /**
     * Private internal memory for applying eliminating rotamers and rotamter
     * pairs.
     */
    private double minEnergyPairs[][];
    private double maxEnergyPairs[][];
    private double minEnergyTriples[][][][];
    private double maxEnergyTriples[][][][];

    /**
     * RotamerOptimization constructor.
     *
     * @param molecularAssembly The MolecularAssembly to search rotamers for.
     * @param algorithmListener AlgorithmListener to update the GUI.
     */
    public RotamerOptimization(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener) {
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
    }

    public RotamerOptimization(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener,
            int startResID, int finalResID, Algorithm algorithm) {
        this(molecularAssembly, algorithmListener);
        this.startResID = startResID;
        this.finalResID = finalResID;
        this.algorithm = algorithm;
    }

    public double optimize() {
        double e = Double.MAX_VALUE;
        done = false;
        terminate = false;
        if (finalResID >= startResID) {
            potential = (ForceFieldEnergy) molecularAssembly.getPotentialEnergy();
            polymers = molecularAssembly.getChains();

            Residue residues[] = new Residue[finalResID - startResID + 1];
            for (int i = startResID, j = 0; i <= finalResID; i++, j++) {
                residues[j] = polymers[0].getResidue(i);
            }

            initDEE(residues);

            switch (algorithm) {
                case INDEPENDENT:
                    e = independent();
                    break;
                case GLOBAL:
                    e = global();
                    break;
                case SLIDING_WINDOW:
                    e = slidingWindow(windowSize, direction);
                    break;
            }
        }
        terminate = false;
        done = true;
        return e;
    }

    private void logRotamers() {
        int rotamerCount = 0;
        int pairCount = 0;
        int eliminated = 0;
        int eliminatedPairs = 0;
        int nres = eliminatedRotamers.length;

        for (int i = 0; i < nres; i++) {
            if (eliminatedRotamers[i] != null) {
                rotamerCount += eliminatedRotamers[i].length;
                for (int ri = 0; ri < eliminatedRotamers[i].length; ri++) {
                    if (eliminatedRotamers[i][ri]) {
                        eliminated++;
                    }
                    for (int j = i + 1; j < nres; j++) {
                        if (eliminatedRotamerPairs[i][ri][j] != null) {
                            pairCount += eliminatedRotamerPairs[i][ri][j].length;
                            for (int rj = 0; rj < eliminatedRotamerPairs[i][ri][j].length; rj++) {
                                if (eliminatedRotamerPairs[i][ri][j][rj]) {
                                    eliminatedPairs++;
                                }
                            }
                        }
                    }
                }
            }
        }
        logger.info(String.format(" %d out of %d rotamers eliminated.", eliminated, rotamerCount));
        logger.info(String.format(" %d out of %d rotamer pairs eliminated.", eliminatedPairs, pairCount));
    }

    public void initDEE(Residue residues[]) {
        selfAndPairEnergies(residues);
        allocateDEEMemory(residues);
        logRotamers();
        while (applyDEEConditions(residues)) {
            logRotamers();
        }
    }

    /**
     * Turn off non-bonded contributions from all residues except for one.
     * Compute the self-energy for each residue relative to the backbone
     * contribution.
     *
     * @param residues A list of residues that we undergo rotamer optimization.
     * @return template energy
     */
    private double selfAndPairEnergies(Residue residues[]) {
        if (residues == null || residues.length < 2) {
            logger.warning(" Attempt to compute self energy for an empty array of residues.");
            return 0;
        }

        int nResidues = residues.length;
        int nAtoms = molecularAssembly.getAtomArray().length;

        selfEnergy = new double[nResidues][];
        pairEnergy = new double[nResidues][][][];

        /**
         * Initialize all atoms to be used.
         */
        boolean use[] = new boolean[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            use[i] = true;
        }

        logger.info(" Turning off side-chain atoms.");
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer rotamers[] = RotamerLibrary.getRotamers(name);
            if (rotamers == null) {
                continue;
            }
            // Place each residue into its zeroth rotamer.
            RotamerLibrary.applyRotamer(name, residue, rotamers[0]);

            // Turn off all side-chain atoms that will be optimized.
            ArrayList<Atom> atoms = residue.getSideChainAtoms();
            for (Atom atom : atoms) {
                logger.info(atom.toString());
                use[atom.getXYZIndex() - 1] = false;
            }
        }
        ForceFieldEnergy energy = molecularAssembly.getPotentialEnergy();
        energy.setUse(use);

        boolean print = false;

        // Compute the backbone energy.
        double backboneEnergy = energy.energy(false, print);
        logger.info(String.format(" Backbone energy: %16.8f\n", backboneEnergy));

        // Compute the self-energy for each rotamer of each residue
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer rotamers[] = RotamerLibrary.getRotamers(name);
            if (rotamers == null) {
                continue;
            }

            // Turn on this residue's side-chain atoms
            ArrayList<Atom> atoms = residue.getSideChainAtoms();
            for (Atom atom : atoms) {
                use[atom.getXYZIndex()] = true;
            }
            energy.setUse(use);

            // Create space for this residue's rotamer self-energies.
            selfEnergy[i] = new double[rotamers.length];

            // Loop over rotamers computing self-energies.
            for (int ri = 0; ri < rotamers.length; ri++) {
                Rotamer rotamer = rotamers[ri];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
                selfEnergy[i][ri] = energy.energy(false, print) - backboneEnergy;
                logger.info(String.format(" Self-energy %s %d %16.8f", name, ri, selfEnergy[i][ri]));
            }

            // Reset the residue to its zeroth rotamer.
            RotamerLibrary.applyRotamer(name, residue, rotamers[0]);

            // Turn off the side-chain
            for (Atom atom : atoms) {
                use[atom.getXYZIndex() - 1] = false;
            }
            energy.setUse(use);
        }

        // Compute the pair-energy for each pair of rotamers
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residues[i];
            AminoAcid3 namei = AminoAcid3.valueOf(residuei.getName());
            Rotamer rotamersi[] = RotamerLibrary.getRotamers(namei);
            if (rotamersi == null) {
                continue;
            }
            ArrayList<Atom> atomsi = residuei.getSideChainAtoms();

            // Turn on residue i
            for (Atom atom : atomsi) {
                use[atom.xyzIndex - 1] = true;
            }
            int ni = rotamersi.length;
            pairEnergy[i] = new double[ni][][];
            for (int ri = 0; ri < ni; ri++) {
                Rotamer rotameri = rotamersi[ri];
                RotamerLibrary.applyRotamer(namei, residuei, rotameri);
                //int npairs = residues.length - (i + 1);
                // TODO: reduce memory use.
                pairEnergy[i][ri] = new double[nResidues][];
                for (int j = i + 1; j < nResidues; j++) {
                    Residue residuej = residues[j];
                    AminoAcid3 namej = AminoAcid3.valueOf(residuej.getName());
                    Rotamer rotamersj[] = RotamerLibrary.getRotamers(namej);
                    if (rotamersj == null) {
                        continue;
                    }
                    ArrayList<Atom> atomsj = residuej.getSideChainAtoms();
                    // Turn on residue j
                    for (Atom atom : atomsj) {
                        use[atom.getXYZIndex() - 1] = true;
                    }
                    energy.setUse(use);

                    // Loop over residue j's rotamers and compute pair-wise energies.
                    int nj = rotamersj.length;
                    pairEnergy[i][ri][j] = new double[nj];
                    for (int rj = 0; rj < nj; rj++) {
                        Rotamer rotamerj = rotamersj[rj];
                        RotamerLibrary.applyRotamer(namej, residuej, rotamerj);
                        pairEnergy[i][ri][j][rj] = energy.energy(false, print)
                                - selfEnergy[i][ri] - selfEnergy[j][rj] - backboneEnergy;
                        logger.info(String.format(" Pair energy %s-%d %d %s-%d %d %16.8f",
                                namei, residuei.getResidueNumber(), ri,
                                namej, residuej.getResidueNumber(), rj, pairEnergy[i][ri][j][rj]));
                    }
                    // Reset the residue to rotamer zero.
                    RotamerLibrary.applyRotamer(namej, residuej, rotamersj[0]);
                    // Turn off residue j
                    for (Atom atom : atomsj) {
                        use[atom.getXYZIndex() - 1] = false;
                    }
                }
            }
            // Reset the residue to rotamer zero.
            RotamerLibrary.applyRotamer(namei, residuei, rotamersi[0]);
            // Turn off residue i
            for (Atom atom : atomsi) {
                use[atom.getXYZIndex() - 1] = false;
            }
            energy.setUse(use);
        }

        // Turn on all atoms.
        for (int i = 0; i < nAtoms; i++) {
            use[i] = true;
        }
        energy.setUse(use);

        return backboneEnergy;
    }

    private void allocateDEEMemory(Residue[] residues) {

        int nres = residues.length;
        eliminatedRotamers = new boolean[nres][];
        eliminatedRotamerPairs = new boolean[nres][][][];
        /**
         * Find the min/max energy for each { residue (i), rotamer (ri)} set by
         * summing the min/max pairwise energy for each residue (j).
         */
        minEnergyPairs = new double[nres][];
        maxEnergyPairs = new double[nres][];
        /**
         * Find the min/max energy for each {residue (i), rotamer (ri), residue
         * (j), rotamer (rj)} set by summing the min/max pairwise energy for
         * each residue (k).
         */
        minEnergyTriples = new double[nres][][][];
        maxEnergyTriples = new double[nres][][][];
        // Loop over residues.
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            AminoAcid3 namei = AminoAcid3.valueOf(residuei.getName());
            Rotamer rotamersi[] = RotamerLibrary.getRotamers(namei);
            if (rotamersi == null) {
                continue;
            }
            int lenri = rotamersi.length;
            eliminatedRotamers[i] = new boolean[lenri];
            eliminatedRotamerPairs[i] = new boolean[lenri][][];
            minEnergyPairs[i] = new double[lenri];
            maxEnergyPairs[i] = new double[lenri];
            minEnergyTriples[i] = new double[lenri][][];
            maxEnergyTriples[i] = new double[lenri][][];
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // int npairs = nres - (i + 1);
                // TODO - reduce memory by half.
                minEnergyTriples[i][ri] = new double[nres][];
                maxEnergyTriples[i][ri] = new double[nres][];
                eliminatedRotamerPairs[i][ri] = new boolean[nres][];
                eliminatedRotamers[i][ri] = false;
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    AminoAcid3 namej = AminoAcid3.valueOf(residuej.getName());
                    Rotamer rotamersj[] = RotamerLibrary.getRotamers(namej);
                    // Some residues do not have 2 or more rotamers.
                    if (rotamersj == null || rotamersj.length == 1) {
                        continue;
                    }
                    int lenrj = rotamersj.length;
                    minEnergyTriples[i][ri][j] = new double[lenrj];
                    maxEnergyTriples[i][ri][j] = new double[lenrj];
                    eliminatedRotamerPairs[i][ri][j] = new boolean[lenrj];
                    for (int rj = 0; rj < lenrj; rj++) {
                        eliminatedRotamerPairs[i][ri][j][rj] = false;
                    }
                }
            }
        }
    }

    /**
     * Elimination of rotamers.
     */
    private boolean applyDEEConditions(Residue[] residues) {
        int nres = residues.length;
        boolean eliminated = false;
        // Loop over residues.
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            AminoAcid3 namei = AminoAcid3.valueOf(residuei.getName());
            Rotamer rotamersi[] = RotamerLibrary.getRotamers(namei);
            if (rotamersi == null) {
                continue;
            }
            int lenri = rotamersi.length;
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // Check for an eliminated single.
                if (eliminatedRotamers[i][ri]) {
                    continue;
                }
                // Start the min/max summation with the self-energy.
                minEnergyPairs[i][ri] = selfEnergy[i][ri];
                maxEnergyPairs[i][ri] = minEnergyPairs[i][ri];
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    AminoAcid3 namej = AminoAcid3.valueOf(residuej.getName());
                    Rotamer rotamersj[] = RotamerLibrary.getRotamers(namej);
                    // Some residues do not have 2 or more rotamers.
                    if (rotamersj == null || rotamersj.length == 1) {
                        continue;
                    }
                    int lenrj = rotamersj.length;
                    double minPairE = Double.MAX_VALUE;
                    double maxPairE = Double.MIN_VALUE;
                    // Loop over residue j's rotamers.
                    for (int rj = 0; rj < lenrj; rj++) {
                        // Check for an eliminated single.
                        if (eliminatedRotamers[j][rj]) {
                            continue;
                        }
                        // Start the min/max summation with the "pair" self-energy.
                        minEnergyTriples[i][ri][j][rj] = selfEnergy[i][ri] + selfEnergy[j][rj] + pairEnergy[i][ri][j][rj];
                        maxEnergyTriples[i][ri][j][rj] = minEnergyTriples[i][ri][j][rj];
                        double current = pairEnergy[i][ri][j][rj];
                        if (current < minPairE) {
                            minPairE = current;
                        }
                        if (current > maxPairE) {
                            maxPairE = current;
                        }

                        // Loop over the third residue.
                        for (int k = j + 1; k < nres; k++) {
                            Residue residuek = residues[k];
                            AminoAcid3 namek = AminoAcid3.valueOf(residuek.getName());
                            Rotamer rotamersk[] = RotamerLibrary.getRotamers(namek);
                            // Some residues do not have 2 or more rotamers.
                            if (rotamersk == null || rotamersk.length == 1) {
                                continue;
                            }
                            int lenrk = rotamersk.length;
                            double minTripleE = Double.MAX_VALUE;
                            double maxTripleE = Double.MIN_VALUE;
                            // Loop over the third residues' rotamers.
                            for (int rk = 0; rk < lenrk; rk++) {
                                if (eliminatedRotamerPairs[i][ri][k][rk] || eliminatedRotamerPairs[j][rj][k][rk]) {
                                    continue;
                                }
                                double currentTripleE = pairEnergy[i][ri][k][rk] + pairEnergy[j][rj][k][rk];
                                if (currentTripleE < minTripleE) {
                                    minTripleE = currentTripleE;
                                }
                                if (currentTripleE > maxTripleE) {
                                    maxTripleE = currentTripleE;
                                }
                            }
                            minEnergyTriples[i][ri][j][rj] += minTripleE;
                            maxEnergyTriples[i][ri][j][rj] += maxTripleE;
                        }
                    }
                    // Sum the min/max contributions.
                    minEnergyPairs[i][ri] += minPairE;
                    maxEnergyPairs[i][ri] += maxPairE;
                }
            }

            /**
             * Apply the singles elimination criteria to rotamers of residue i
             * by determining the most favorable maximum energy.
             */
            double eliminationEnergy = Double.MAX_VALUE;
            for (int ri = 0; ri < lenri; ri++) {
                if (eliminatedRotamers[i][ri]) {
                    continue;
                }
                if (maxEnergyPairs[i][ri] < eliminationEnergy) {
                    eliminationEnergy = maxEnergyPairs[i][ri];
                }
            }

            /**
             * Eliminate rotamers whose minimum energy is greater than the worst
             * case for another rotamer.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (eliminatedRotamers[i][ri]) {
                    continue;
                }
                if (minEnergyPairs[i][ri] > eliminationEnergy) {
                    eliminatedRotamers[i][ri] = true;
                    logger.info(String.format(" Eliminating single %d %d", i, ri));
                    eliminated = true;
                    if (eliminatedRotamerPairs[i][ri] == null) {
                        continue;
                    }
                    int lenj = eliminatedRotamerPairs[i][ri].length;
                    for (int j = i + 1; j < lenj; j++) {
                        if (eliminatedRotamerPairs[i][ri][j] == null) {
                            continue;
                        }
                        int lenrj = eliminatedRotamerPairs[i][ri][j].length;
                        for (int rj = 0; rj < lenrj; rj++) {
                            eliminatedRotamerPairs[i][ri][j][rj] = true;
                            logger.info(String.format(" Eliminating pair %d %d %d %d", i, ri, j, rj));
                        }
                    }
                }
            }

            /**
             * Apply the double elimination criteria to rotamer pairs for
             * residue i, rotamer ri, residue j and rotamer rj by determining
             * the most favorable maximum energy.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (eliminatedRotamers[i][ri]) {
                    continue;
                }
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    AminoAcid3 namej = AminoAcid3.valueOf(residuej.getName());
                    Rotamer rotamersj[] = RotamerLibrary.getRotamers(namej);
                    // Some residues do not have 2 or more rotamers.
                    if (rotamersj == null || rotamersj.length == 1) {
                        continue;
                    }
                    int lenrj = rotamersj.length;
                    double pairEliminationEnergy = Double.MAX_VALUE;
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (eliminatedRotamerPairs[i][ri][j][rj]) {
                            continue;
                        }
                        if (maxEnergyTriples[i][ri][j][rj] < pairEliminationEnergy) {
                            pairEliminationEnergy = maxEnergyTriples[i][ri][j][rj];
                        }
                    }
                    /**
                     * Eliminate rotamer pairs whose minimum energy is higher
                     * than the worst case for an alternative pair.
                     */
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (eliminatedRotamerPairs[i][ri][j][rj]) {
                            continue;
                        }
                        if (minEnergyTriples[i][ri][j][rj] > pairEliminationEnergy) {
                            eliminatedRotamerPairs[i][ri][j][rj] = true;
                            logger.info(String.format(" Eliminating pair %d %d %d %d", i, ri, j, rj));
                            eliminated = true;
                        }
                    }
                }
            }
        }
        return eliminated;
    }

    public double optimize(int startResID, int finalResID, Algorithm algorithm) {
        this.startResID = startResID;
        this.finalResID = finalResID;
        this.algorithm = algorithm;
        return optimize();
    }

    public double optimize(int startResID, int finalResID, int windowSize, Direction direction, Algorithm algorithm) {
        this.startResID = startResID;
        this.finalResID = finalResID;
        this.windowSize = windowSize;
        this.direction = direction;
        this.algorithm = algorithm;
        return optimize();
    }

    private double independent() {
        double e = Double.MAX_VALUE;
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            logger.info(String.format(" Optimizing %s side-chain.", residue));
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            if (rotamers == null) {
                continue;
            }
            e = potential.energy(false, true);
            int bestRotamer = -1;
            for (int j = 0; j < rotamers.length; j++) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                double newE = potential.energy(false, true);
                if (newE < e) {
                    bestRotamer = j;
                }
            }
            if (bestRotamer > -1) {
                Rotamer rotamer = rotamers[bestRotamer];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
            }
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            if (terminate) {
                logger.info(String.format("\n Terminating after residue %s.\n", residue));
                break;
            }
        }
        return e;
    }

    private double global() {
        ArrayList<Residue> residues = new ArrayList<Residue>();
        int permutations = 1;
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            residues.add(residue);
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            if (rotamers != null) {
                permutations *= rotamers.length;
            }
        }
        logger.info(String.format(" Number of permutations: %d.", permutations));
        ArrayList<Integer> optimum = new ArrayList<Integer>();
        double e = RotamerLibrary.rotamerOptimization(molecularAssembly, residues, Double.MAX_VALUE, optimum);
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            int j = optimum.remove(0);
            if (rotamers != null) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
            }
        }
        return e;
    }

    private double slidingWindow(int windowSize, Direction direction) {
        double e = Double.MAX_VALUE;

        if ((finalResID - startResID) < windowSize - 1) {
            logger.warning("StartResID and FinalResID too close for sliding window size.");
        }

        switch (direction) {
            case FORWARD:
                for (int startWindow = startResID; startWindow + (windowSize - 1) <= finalResID; startWindow++) {

                    if (polymers[0].getResidue(startWindow + (windowSize - 1)) == null) {
                        logger.warning("FinalResID references non-existent residue; terminating at end of chain.");
                        break;
                    }
                    ArrayList<Residue> residues = new ArrayList<Residue>();
                    int permutations = 1;
                    for (int i = startWindow; i <= startWindow + (windowSize - 1); i++) {
                        Residue residue = polymers[0].getResidue(i);
                        residues.add(residue);
                        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
                        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
                        if (rotamers != null) {
                            permutations *= rotamers.length;
                        }
                    }
                    logger.info(String.format(" Number of permutations: %d.", permutations));
                    ArrayList<Integer> optimum = new ArrayList<Integer>();
                    e = RotamerLibrary.rotamerOptimization(molecularAssembly, residues, Double.MAX_VALUE, optimum);
                    for (int i = startWindow; i <= startWindow + (windowSize - 1); i++) {
                        Residue residue = polymers[0].getResidue(i);
                        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
                        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
                        int j = optimum.remove(0);
                        if (rotamers != null) {
                            Rotamer rotamer = rotamers[j];
                            RotamerLibrary.applyRotamer(name, residue, rotamer);
                        }
                    }
                }
                break;
            case BACKWARD:
                for (int endWindow = finalResID; endWindow - (windowSize - 1) >= startResID; endWindow--) {

                    if (polymers[0].getResidue(endWindow - windowSize) == null) {
                        logger.warning("StartResID references non-existent residue; terminating at beginning of chain.");
                        break;
                    }
                    ArrayList<Residue> residues = new ArrayList<Residue>();
                    int permutations = 1;
                    for (int i = endWindow; i >= endWindow - (windowSize - 1); i--) {
                        Residue residue = polymers[0].getResidue(i);
                        residues.add(residue);
                        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
                        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
                        if (rotamers != null) {
                            permutations *= rotamers.length;
                        }
                    }
                    logger.info(String.format(" Number of permutations: %d.", permutations));
                    ArrayList<Integer> optimum = new ArrayList<Integer>();
                    e = RotamerLibrary.rotamerOptimization(molecularAssembly, residues, Double.MAX_VALUE, optimum);
                    for (int i = endWindow; i >= endWindow - (windowSize - 1); i--) {
                        Residue residue = polymers[0].getResidue(i);
                        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
                        Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
                        int j = optimum.remove(0);
                        if (rotamers != null) {
                            Rotamer rotamer = rotamers[j];
                            RotamerLibrary.applyRotamer(name, residue, rotamer);
                        }
                    }
                }
                break;
        }
        return e;
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating rotamer optimization.\n", e);
                }
            }
        }
    }
}
