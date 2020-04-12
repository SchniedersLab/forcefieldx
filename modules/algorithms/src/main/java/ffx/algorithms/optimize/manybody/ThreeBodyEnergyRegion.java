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
package ffx.algorithms.optimize.manybody;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.min;

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
 * Compute 3-Body energy values in parallel across nodes.
 */
public class ThreeBodyEnergyRegion extends WorkerRegion {

    private static final Logger logger = Logger.getLogger(ThreeBodyEnergyRegion.class.getName());
    private final Residue[] residues;
    private final RotamerOptimization rO;
    private final DistanceMatrix dM;
    private final EnergyExpansion eE;
    private final EliminatedRotamers eR;
    /**
     * A list of all residues being optimized. Note that Box and Window
     * optimizations operate on subsets of this list.
     */
    private final List<Residue> allResiduesList;
    /**
     * RotamerLibrary instance.
     */
    private final RotamerLibrary library;
    /**
     * Map of 3-body energy values to compute.
     */
    private final Map<Integer, Integer[]> threeBodyEnergyMap;
    /**
     * Writes energies to restart file.
     */
    private final BufferedWriter energyWriter;
    /**
     * World Parallel Java communicator.
     */
    private final Comm world;
    /**
     * Number of Parallel Java processes.
     */
    private final int numProc;
    /**
     * If a pair of residues have two atoms closer together than the
     * superposition threshold, the energy is set to NaN.
     */
    private final double superpositionThreshold;
    /**
     * Flag to indicate if this is the master process.
     */
    private final boolean master;
    /**
     * Rank of this process.
     */
    private final int rank;
    /**
     * Flag to indicate verbose logging.
     */
    private final boolean verbose;
    /**
     * If true, write out an energy restart file.
     */
    private final boolean writeEnergyRestart;
    /**
     * Sets whether files should be printed; true for standalone applications,
     * false for some applications which use rotamer optimization as part of a
     * larger process.
     */
    private final boolean printFiles;
    private Set<Integer> keySet;

    public ThreeBodyEnergyRegion(RotamerOptimization rotamerOptimization, DistanceMatrix dM,
                                 EnergyExpansion eE, EliminatedRotamers eR,
                                 Residue[] residues, List<Residue> allResiduesList, RotamerLibrary library,
                                 BufferedWriter energyWriter, Comm world, int numProc, double superpositionThreshold,
                                 boolean master, int rank, boolean verbose, boolean writeEnergyRestart, boolean printFiles) {
        this.rO = rotamerOptimization;
        this.dM = dM;
        this.eE = eE;
        this.eR = eR;
        this.residues = residues;
        this.allResiduesList = allResiduesList;
        this.library = library;
        this.energyWriter = energyWriter;
        this.world = world;
        this.numProc = numProc;
        this.superpositionThreshold = superpositionThreshold;
        this.master = master;
        this.rank = rank;
        this.verbose = verbose;
        this.writeEnergyRestart = writeEnergyRestart;
        this.printFiles = printFiles;

        this.threeBodyEnergyMap = eE.getThreeBodyEnergyMap();
        logger.info(format(" Number of 3-Body energies to calculate: %d", threeBodyEnergyMap.size()));
    }

    @Override
    public void finish() {
        // Print what we've got so far.
        if (master && verbose) {
            for (int i = 0; i < residues.length; i++) {
                Residue resi = residues[i];
                Rotamer[] roti = resi.getRotamers(library);
                for (int ri = 0; ri < roti.length; ri++) {
                    if (eR.check(i, ri)) {
                        continue;
                    }
                    for (int j = i + 1; j < residues.length; j++) {
                        Residue resj = residues[j];
                        Rotamer[] rotj = resj.getRotamers(library);
                        for (int rj = 0; rj < rotj.length; rj++) {
                            if (eR.check(j, rj) || eR.check(i, ri, j, rj)) {
                                continue;
                            }
                            for (int k = j + 1; k < residues.length; k++) {
                                Residue resk = residues[k];
                                Rotamer[] rotk = resk.getRotamers(library);
                                for (int rk = 0; rk < rotk.length; rk++) {
                                    if (eR.check(k, rk) || eR.check(i, ri, k, rk) || eR.check(j, rj, k, rk)) {
                                        continue;
                                    }
                                    logger.info(format(" 3-Body energy %8s %-2d, %8s %-2d, %8s %-2d: %s",
                                            resi.toFormattedString(false, true),
                                            ri, resj.toFormattedString(false, true),
                                            rj, resk.toFormattedString(false, true),
                                            rk, rO.formatEnergy(eE.get3Body(residues, i, ri, j, rj, k, rk))));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    @Override
    public void run() throws Exception {
        if (!keySet.isEmpty()) {
            execute(0, keySet.size() - 1, new ThreeBodyEnergyLoop());
        }
    }

    @Override
    public void start() {
        int numTriple = threeBodyEnergyMap.size();
        int remainder = numTriple % numProc;
        // Set padded residue and rotamer to less than zero.
        Integer[] padding = {-1, -1, -1, -1, -1, -1};

        int padKey = numTriple;
        while (remainder != 0) {
            threeBodyEnergyMap.put(padKey++, padding);
            remainder = threeBodyEnergyMap.size() % numProc;
        }

        numTriple = threeBodyEnergyMap.size();
        if (numTriple % numProc != 0) {
            logger.severe(" Logic error padding triple energies.");
        }

        // Load the keySet of triple energies.
        keySet = threeBodyEnergyMap.keySet();
    }

    private class ThreeBodyEnergyLoop extends WorkerIntegerForLoop {

        DoubleBuf[] resultBuffer;
        DoubleBuf myBuffer;

        ThreeBodyEnergyLoop() {
            resultBuffer = new DoubleBuf[numProc];
            for (int i = 0; i < numProc; i++) {
                resultBuffer[i] = DoubleBuf.buffer(new double[7]);
            }
            myBuffer = resultBuffer[rank];
        }

        @Override
        public void run(int lb, int ub) {
            for (int key = lb; key <= ub; key++) {
                long time = -System.nanoTime();
                Integer[] job = threeBodyEnergyMap.get(key);
                int i = job[0];
                int ri = job[1];
                int j = job[2];
                int rj = job[3];
                int k = job[4];
                int rk = job[5];

                myBuffer.put(0, i);
                myBuffer.put(1, ri);
                myBuffer.put(2, j);
                myBuffer.put(3, rj);
                myBuffer.put(4, k);
                myBuffer.put(5, rk);
                myBuffer.put(6, 0.0);

                // Initialize result.
                if (i >= 0 && ri >= 0 && j >= 0 && rj >= 0 && k >= 0 && rk >= 0) {
                    if ((!eR.check(i, ri) || !eR.check(j, rj) || !eR.check(k, rk) ||
                            !eR.check(i, ri, j, rj) || !eR.check(i, ri, k, rk) || !eR.check(j, rj, k, rk))) {

                        Residue residueI = residues[i];
                        Residue residueJ = residues[j];
                        Residue residueK = residues[k];

                        int indexI = allResiduesList.indexOf(residueI);
                        int indexJ = allResiduesList.indexOf(residueJ);
                        int indexK = allResiduesList.indexOf(residueK);

                        double rawDist = dM.getRawNBodyDistance(indexI, ri, indexJ, rj, indexK, rk);
                        double dIJ = dM.checkDistMatrix(indexI, ri, indexJ, rj);
                        double dIK = dM.checkDistMatrix(indexI, ri, indexK, rk);
                        double dJK = dM.checkDistMatrix(indexJ, rj, indexK, rk);
                        double minDist = min(min(dIJ, dIK), dJK);

                        double resDist = dM.get3BodyResidueDistance(indexI, ri, indexJ, rj, indexK, rk);
                        String resDistString = "     large";
                        if (resDist < Double.MAX_VALUE) {
                            resDistString = format("%5.3f", resDist);
                        }

                        String distString = "     large";
                        if (rawDist < Double.MAX_VALUE) {
                            distString = format("%10.3f", rawDist);
                        }

                        double threeBodyEnergy;
                        if (minDist < superpositionThreshold) {
                            threeBodyEnergy = Double.NaN;
                            logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d:\t    NaN      at %13.6f Ang (%s Ang by residue) < %5.3f Ang.",
                                    residueI.toFormattedString(false, true), ri,
                                    residueJ.toFormattedString(false, true), rj,
                                    residueK.toFormattedString(false, true), rk,
                                    minDist, resDistString, superpositionThreshold));
                        } else if (dM.checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
                            // Set the two-body energy to 0.0 for separation distances larger than the two-body cutoff.
                            threeBodyEnergy = 0.0;
                            time += System.nanoTime();
                            logger.fine(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true), ri,
                                    residueJ.toFormattedString(false, true), rj,
                                    residueK.toFormattedString(false, true), rk,
                                    rO.formatEnergy(threeBodyEnergy), distString, resDistString, time * 1.0e-9));
                        } else {
                            try {
                                threeBodyEnergy = eE.compute3BodyEnergy(residues, i, ri, j, rj, k, rk);
                                time += System.nanoTime();
                                logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                        residueI.toFormattedString(false, true), ri,
                                        residueJ.toFormattedString(false, true), rj,
                                        residueK.toFormattedString(false, true), rk,
                                        rO.formatEnergy(threeBodyEnergy), distString, resDistString, time * 1.0e-9));
                            } catch (ArithmeticException ex) {
                                threeBodyEnergy = Double.NaN;
                                time += System.nanoTime();
                                logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d:\t    NaN      at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                        residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, residueK.toFormattedString(false, true), rk, distString, resDistString, time * 1.0e-9));
                            }
                        }
                        myBuffer.put(6, threeBodyEnergy);
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

                // Process the three-body energy received from each process.
                for (DoubleBuf doubleBuf : resultBuffer) {
                    int resi = (int) doubleBuf.get(0);
                    int roti = (int) doubleBuf.get(1);
                    int resj = (int) doubleBuf.get(2);
                    int rotj = (int) doubleBuf.get(3);
                    int resk = (int) doubleBuf.get(4);
                    int rotk = (int) doubleBuf.get(5);
                    double energy = doubleBuf.get(6);
                    // Skip for padded result.
                    if (resi >= 0 && roti >= 0 && resj >= 0 && rotj >= 0 && resk >= 0 && rotk >= 0) {
                        if (!Double.isFinite(energy)) {
                            logger.info(" Rotamer pair eliminated: " + resi + ", " + roti + ", " + resj + ", " + rotj);
                            eR.eliminateRotamerPair(residues, resi, roti, resj, rotj, false);
                        }
                        eE.set3Body(residues, resi, roti, resj, rotj, resk, rotk, energy);
                        if (rank == 0 && writeEnergyRestart && printFiles) {
                            try {
                                energyWriter.append(format("Triple %d %d, %d %d, %d %d: %16.8f", resi, roti, resj, rotj, resk, rotk, energy));
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

        @Override
        public IntegerSchedule schedule() {
            // The schedule must be fixed.
            return IntegerSchedule.fixed();
        }
    }
}
