//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.optimize.manybody;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerRegion;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;

/**
 * Compute 2-Body energy values in parallel across nodes.
 */
public class TwoBodyEnergyRegion extends WorkerRegion {

    private static final Logger logger = Logger.getLogger(TwoBodyEnergyRegion.class.getName());

    private RotamerOptimization rO;
    private DistanceMatrix dM;
    private RotamerOptimization.EnergyExpansion eE;
    private final Residue[] residues;
    private Set<Integer> keySet;
    /**
     * A list of all residues being optimized. Note that Box and Window
     * optimizations operate on subsets of this list.
     */
    private ArrayList<Residue> allResiduesList;
    /**
     * RotamerLibrary instance.
     */
    private RotamerLibrary library;
    /**
     * Map of self-energy values to compute.
     */
    private HashMap<Integer, Integer[]> twoBodyEnergyMap;
    /**
     * Writes energies to restart file.
     */
    private BufferedWriter energyWriter;
    /**
     * World Parallel Java communicator.
     */
    private Comm world;
    /**
     * Number of Parallel Java processes.
     */
    private int numProc;
    /**
     * Flag to prune clashes.
     */
    private boolean prunePairClashes;
    /**
     * If a pair of residues have two atoms closer together than the
     * superposition threshold, the energy is set to NaN.
     */
    private double superpositionThreshold;
    /**
     * Flag to indicate if this is the master process.
     */
    private boolean master;
    /**
     * Rank of this process.
     */
    private int rank;
    /**
     * Flag to indicate verbose logging.
     */
    private boolean verbose;
    /**
     * If true, write out an energy restart file.
     */
    private boolean writeEnergyRestart;
    /**
     * Sets whether files should be printed; true for standalone applications,
     * false for some applications which use rotamer optimization as part of a
     * larger process.
     */
    private boolean printFiles;

    public TwoBodyEnergyRegion(RotamerOptimization rotamerOptimization, DistanceMatrix dM, RotamerOptimization.EnergyExpansion eE,
                               Residue[] residues,
                               ArrayList<Residue> allResiduesList, RotamerLibrary library,
                               HashMap<Integer, Integer[]> twoBodyEnergyMap, BufferedWriter energyWriter,
                               Comm world, int numProc, boolean prunePairClashes, double superpositionThreshold,
                               boolean master, int rank, boolean verbose, boolean writeEnergyRestart, boolean printFiles) {
        this.rO = rotamerOptimization;
        this.dM = dM;
        this.eE = eE;
        this.residues = residues;
        this.allResiduesList = allResiduesList;
        this.library = library;
        this.twoBodyEnergyMap = twoBodyEnergyMap;
        this.energyWriter = energyWriter;
        this.world = world;
        this.numProc = numProc;
        this.prunePairClashes = prunePairClashes;
        this.superpositionThreshold = superpositionThreshold;
        this.master = master;
        this.rank = rank;
        this.verbose = verbose;
        this.writeEnergyRestart = writeEnergyRestart;
        this.printFiles = printFiles;
    }

    @Override
    public void start() {
        int numPair = twoBodyEnergyMap.size();
        int remainder = numPair % numProc;

        // Set padded residue and rotamer to less than zero.
        Integer[] padding = {-1, -1, -1, -1};

        int padKey = numPair;
        while (remainder != 0) {
            twoBodyEnergyMap.put(padKey++, padding);
            remainder = twoBodyEnergyMap.size() % numProc;
        }

        numPair = twoBodyEnergyMap.size();
        if (numPair % numProc != 0) {
            logger.severe(" Logic error padding pair energies.");
        }

        // Load the keySet of pair energies.
        keySet = twoBodyEnergyMap.keySet();
    }

    @Override
    public void run() throws Exception {
        if (!keySet.isEmpty()) {
            execute(0, keySet.size() - 1, new TwoBodyEnergyLoop());
        }
    }

    @Override
    public void finish() {
        // Pre-Prune if pair-energy is Double.NaN.
        rO.prePrunePairs(residues);

        // Prune each rotamer that clashes with all rotamers from a 2nd residue.
        if (prunePairClashes) {
            rO.prunePairClashes(residues);
        }

        // Print what we've got so far.
        if (master && verbose) {
            for (int i = 0; i < residues.length; i++) {
                Residue resi = residues[i];
                Rotamer[] roti = resi.getRotamers(library);
                for (int ri = 0; ri < roti.length; ri++) {
                    if (rO.check(i, ri)) {
                        continue;
                    }
                    for (int j = i + 1; j < residues.length; j++) {
                        Residue resj = residues[j];
                        Rotamer[] rotj = resj.getRotamers(library);
                        for (int rj = 0; rj < rotj.length; rj++) {
                            if (rO.check(j, rj) ||
                                    rO.check(i, ri, j, rj)) {
                                continue;
                            }
                            logger.info(format(" Pair energy %8s %-2d, %8s %-2d: %s",
                                    residues[i].toFormattedString(false, true),
                                    ri, residues[j].toFormattedString(false, true),
                                    rj, rO.formatEnergy(
                                            eE.get2Body(i, ri, j, rj))));
                        }
                    }
                }
            }
        }
    }

    private class TwoBodyEnergyLoop extends WorkerIntegerForLoop {

        DoubleBuf[] resultBuffer;
        DoubleBuf myBuffer;

        TwoBodyEnergyLoop() {
            resultBuffer = new DoubleBuf[numProc];
            for (int i = 0; i < numProc; i++) {
                resultBuffer[i] = DoubleBuf.buffer(new double[5]);
            }
            myBuffer = resultBuffer[rank];
        }

        @Override
        public IntegerSchedule schedule() {
            // The schedule must be fixed.
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) {
            for (int key = lb; key <= ub; key++) {
                long time = -System.nanoTime();
                Integer[] job = twoBodyEnergyMap.get(key);
                int i = job[0];
                int ri = job[1];
                int j = job[2];
                int rj = job[3];

                myBuffer.put(0, i);
                myBuffer.put(1, ri);
                myBuffer.put(2, j);
                myBuffer.put(3, rj);
                myBuffer.put(4, 0.0);

                // Initialize result.
                if (i >= 0 && ri >= 0 && j >= 0 && rj >= 0) {
                    if (!rO.check(i, ri) ||
                            !rO.check(j, rj) ||
                            !rO.check(i, ri, j, rj)) {
                        Residue residueI = residues[i];
                        Residue residueJ = residues[j];
                        int indexI = allResiduesList.indexOf(residueI);
                        int indexJ = allResiduesList.indexOf(residueJ);
                        double resDist = dM.getResidueDistance(indexI, ri, indexJ, rj);
                        String resDistString = "large";
                        if (resDist < Double.MAX_VALUE) {
                            resDistString = format("%5.3f", resDist);
                        }

                        double dist = dM.checkDistMatrix(indexI, ri, indexJ, rj);
                        String distString = "     large";
                        if (dist < Double.MAX_VALUE) {
                            distString = format("%10.3f", dist);
                        }

                        double twoBodyEnergy;
                        if (dist < superpositionThreshold) {
                            // Set the energy to NaN for superposed atoms.
                            twoBodyEnergy = Double.NaN;
                            logger.info(format(" Pair %8s %-2d, %8s %-2d:\t    NaN at %13.6f Ang (%s Ang by residue) < %5.3f Ang",
                                    residueI.toFormattedString(false, true),
                                    ri, residueJ.toFormattedString(false, true),
                                    rj, dist, resDist, superpositionThreshold));
                        } else if (dM.checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                            // Set the two-body energy to 0.0 for separation distances larger than the two-body cutoff.
                            twoBodyEnergy = 0.0;
                            time += System.nanoTime();
                            logger.info(format(" Pair %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true),
                                    ri, residueJ.toFormattedString(false, true),
                                    rj, rO.formatEnergy(twoBodyEnergy),
                                    distString, resDistString, time * 1.0e-9));
                        } else {
                            try {
                                twoBodyEnergy = eE.compute2BodyEnergy(residues, i, ri, j, rj);
                                time += System.nanoTime();
                                logger.info(format(" Pair %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                        residueI.toFormattedString(false, true),
                                        ri, residueJ.toFormattedString(false, true),
                                        rj, rO.formatEnergy(twoBodyEnergy),
                                        distString, resDistString, time * 1.0e-9));
                            } catch (ArithmeticException ex) {
                                twoBodyEnergy = Double.NaN;
                                time += System.nanoTime();
                                logger.info(format(" Pair %8s %-2d, %8s %-2d:              NaN at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                        residueI.toFormattedString(false, true),
                                        ri, residueJ.toFormattedString(false, true),
                                        rj, distString, resDistString, time * 1.0e-9));
                            }
                        }
                        myBuffer.put(4, twoBodyEnergy);
                    }
                }

                // All to All communication
                if (numProc > 1) {
                    try {
                        world.allGather(myBuffer, resultBuffer);
                    } catch (Exception e) {
                        logger.log(Level.SEVERE, " Exception communicating pair energies.", e);
                    }
                }

                // Process the two-body energy received from each process.
                for (DoubleBuf doubleBuf : resultBuffer) {
                    int resi = (int) doubleBuf.get(0);
                    int roti = (int) doubleBuf.get(1);
                    int resj = (int) doubleBuf.get(2);
                    int rotj = (int) doubleBuf.get(3);
                    double energy = doubleBuf.get(4);
                    // Skip for padded result.
                    if (resi >= 0 && roti >= 0 && resj >= 0 && rotj >= 0) {
                        if (!Double.isFinite(energy)) {
                            logger.info(" Rotamer pair eliminated: " + resi + ", " + roti + ", " + resj + ", " + rotj);
                            rO.eliminateRotamerPair(residues, resi, roti, resj, rotj, false);
                        }
                        eE.set2Body(resi, roti, resj, rotj, energy);
                        if (rank == 0 && writeEnergyRestart && printFiles) {
                            try {
                                energyWriter.append(format("Pair %d %d, %d %d: %16.8f", resi, roti, resj, rotj, energy));
                                energyWriter.newLine();
                                energyWriter.flush();
                            } catch (IOException ex) {
                                logger.log(Level.SEVERE, " Exception writing energy restart file.", ex);
                            }
                        }
                    }
                }
            }
        }
    }
}