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
import java.util.Collection;
import java.util.HashMap;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.MultipleParallelException;
import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerRegion;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.Utilities;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;

/**
 * Compute residue self-energy values in parallel across nodes.
 */
public class SelfEnergyRegion extends WorkerRegion {

    private static final Logger logger = Logger.getLogger(SelfEnergyRegion.class.getName());

    private RotamerOptimization rotamerOptimization;
    private final Residue[] residues;
    private Set<Integer> keySet;
    /**
     * RotamerLibrary instance.
     */
    private RotamerLibrary library;
    /**
     * Map of self-energy values to compute.
     */
    private HashMap<Integer, Integer[]> selfEnergyMap;
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
    private boolean pruneClashes;
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

    /**
     * This needs to be returned to the RotOpt instance.
     */
    private double backboneEnergy;

    public SelfEnergyRegion(RotamerOptimization rotamerOptimization, Residue[] residues, RotamerLibrary library,
                            HashMap<Integer, Integer[]> selfEnergyMap, BufferedWriter energyWriter,
                            Comm world, int numProc, boolean pruneClashes, boolean master,
                            int rank, boolean verbose, boolean writeEnergyRestart, boolean printFiles) {
        this.rotamerOptimization = rotamerOptimization;
        this.residues = residues;
        this.library = library;
        this.selfEnergyMap = selfEnergyMap;
        this.energyWriter = energyWriter;
        this.world = world;
        this.numProc = numProc;
        this.pruneClashes = pruneClashes;
        this.master = master;
        this.rank = rank;
        this.verbose = verbose;
        this.writeEnergyRestart = writeEnergyRestart;
        this.printFiles = printFiles;
    }

    public double getBackboneEnergy() {
        return backboneEnergy;
    }

    @Override
    public void start() {

        int numSelf = selfEnergyMap.size();
        int remainder = numSelf % numProc;
        // Set padded residue and rotamer to less than zero.
        Integer[] padding = {-1, -1};

        int padKey = numSelf;
        while (remainder != 0) {
            selfEnergyMap.put(padKey++, padding);
            remainder = selfEnergyMap.size() % numProc;
        }

        numSelf = selfEnergyMap.size();
        if (numSelf % numProc != 0) {
            logger.severe(" Logic error padding self energies.");
        }

        // Load the keySet of self energies.
        keySet = selfEnergyMap.keySet();

        // Compute backbone energy.
        try {
            backboneEnergy = rotamerOptimization.computeBackboneEnergy(residues);
        } catch (ArithmeticException ex) {
            logger.severe(format(" Error in calculation of backbone energy %s", ex.getMessage()));
        }
        rotamerOptimization.logIfMaster(format(" Backbone energy:  %s\n",
                rotamerOptimization.formatEnergy(backboneEnergy)));
    }

    @Override
    public void run() throws Exception {
        if (!keySet.isEmpty()) {
            try {
                execute(0, keySet.size() - 1, new SelfEnergyLoop());
            } catch (MultipleParallelException mpx) {
                Collection<Throwable> subErrors = mpx.getExceptionMap().values();
                logger.info(format(" MultipleParallelException caught: %s\n Stack trace:\n%s",
                        mpx, Utilities.stackTraceToString(mpx)));
                for (Throwable subError : subErrors) {
                    logger.info(format(" Exception %s\n Stack trace:\n%s",
                            subError, Utilities.stackTraceToString(subError)));
                }
                throw mpx; // Or logger.severe.
            } catch (Throwable t) {
                Throwable cause = t.getCause();
                logger.info(format(" Throwable caught: %s\n Stack trace:\n%s",
                        t, Utilities.stackTraceToString(t)));
                if (cause != null) {
                    logger.info(format(" Cause: %s\n Stack trace:\n%s",
                            cause, Utilities.stackTraceToString(cause)));
                }
                throw t;
            }
        }
    }

    @Override
    public void finish() {
        // Pre-Prune if self-energy is Double.NaN.
        rotamerOptimization.prePruneSelves(residues);

        // Prune clashes for all singles (not just the ones this node did).
        if (pruneClashes) {
            rotamerOptimization.pruneSingleClashes(residues);
        }

        // Print what we've got so far.
        if (master && verbose) {
            for (int i = 0; i < residues.length; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                for (int ri = 0; ri < rotamers.length; ri++) {
                    logger.info(format(" Self energy %8s %-2d: %s",
                            residues[i].toFormattedString(false, true),
                            ri, rotamerOptimization.formatEnergy(
                                    rotamerOptimization.getSelf(i, ri))));
                }
            }
        }
    }

    private class SelfEnergyLoop extends WorkerIntegerForLoop {

        DoubleBuf[] resultBuffer;
        DoubleBuf myBuffer;

        SelfEnergyLoop() {
            resultBuffer = new DoubleBuf[numProc];
            for (int i = 0; i < numProc; i++) {
                resultBuffer[i] = DoubleBuf.buffer(new double[3]);
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
                Integer[] job = selfEnergyMap.get(key);
                int i = job[0];
                int ri = job[1];
                // Initialize result.
                myBuffer.put(0, i);
                myBuffer.put(1, ri);
                myBuffer.put(2, 0.0);

                if (i >= 0 && ri >= 0) {
                    if (!rotamerOptimization.check(i, ri)) {
                        long time = -System.nanoTime();
                        double selfEnergy;
                        try {
                            selfEnergy = rotamerOptimization.computeSelfEnergy(residues, i, ri);
                            time += System.nanoTime();
                            logger.info(format(" Self %8s %-2d: %s in %6.4f (sec).", residues[i].toFormattedString(false, true), ri,
                                    rotamerOptimization.formatEnergy(selfEnergy), time * 1.0e-9));
                        } catch (ArithmeticException ex) {
                            selfEnergy = Double.NaN;
                            time += System.nanoTime();
                            logger.info(format(" Self %8s %-2d:\t    pruned in %6.4f (sec).", residues[i].toFormattedString(false, true), ri, time * 1.0e-9));
                        }
                        myBuffer.put(2, selfEnergy);
                    }
                } else {
                    // allGather parallel command below requires resultBuffer
                    // to have no null elements. Therefore, the padded energies that
                    // enter this else statement must be given an energy of 0.
                    myBuffer.put(2, 0.0);
                }

                // All to All communication
                if (numProc > 1) {
                    try {
                        world.allGather(myBuffer, resultBuffer);
                    } catch (Exception e) {
                        logger.log(Level.SEVERE, " Exception communicating self energies.", e);
                    }
                }

                // Process the self energy received from each process.
                for (DoubleBuf doubleBuf : resultBuffer) {
                    int resi = (int) doubleBuf.get(0);
                    int roti = (int) doubleBuf.get(1);
                    double energy = doubleBuf.get(2);
                    // Skip for padded result.
                    if (resi >= 0 && roti >= 0) {
                        if (Double.isNaN(energy)) {
                            logger.info(" Rotamer  eliminated: " + resi + ", " + roti);
                            rotamerOptimization.eliminateRotamer(residues, resi, roti, false);
                        }
                        rotamerOptimization.setSelf(resi, roti, energy);
                        if (rank == 0 && writeEnergyRestart && printFiles) {
                            try {
                                energyWriter.append(format("Self %d %d: %16.8f", resi, roti, energy));
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
