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
package ffx.algorithms.thermodynamics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Optional;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.LongConsumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import edu.rit.mp.LongBuf;
import edu.rit.pj.Comm;

import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.cli.OSTOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.Constants;

/**
 * An implementation of RepEx between Orthogonal Space Tempering potentials.
 *
 * @author Michael J. Schnieders
 * @author Jacob Litman
 * @since 1.0
 */
public class RepExOST {
    private static final Logger logger = Logger.getLogger(RepExOST.class.getName());

    private final OrthogonalSpaceTempering ost;
    private final OrthogonalSpaceTempering.Histogram[] allHistograms;
    private final SynchronousSend[] sends;
    private final LongConsumer algoRun;
    private final MolecularDynamics molDyn;
    private final DynamicsOptions dynamics;
    private final String fileType;
    private final MonteCarloOST mcOST;
    private final long stepsBetweenExchanges;
    private final Comm world;
    private final int rank;
    private final int numPairs;
    private final int[] rankToHisto;
    private final int[] histoToRank;
    private final boolean isMC;
    private final Random random;
    private final double invKT;
    private final String basePath;
    private final String[] allFilenames;
    private final File dynFile;
    private final String extension;

    private final long[] totalSwaps;
    private final long[] acceptedSwaps;

    // Message tags to use.
    private static final int lamTag = 42;
    private static final int mainLoopTag = 2020;

    private boolean reinitVelocities = true;
    private int currentHistoIndex;

    private double currentLambda;

    private boolean automaticWriteouts = true; // False if the repex OST is not responsible for writing files out.

    /**
     * Private constructor used here to centralize shared logic.
     *
     * @param ost           An OrthogonalSpaceTempering for each repex rung.
     * @param mcOST         A MonteCarloOST for each repex rung, or null (for MD).
     * @param dyn           A MolecularDynamics for each repex rung (never null).
     * @param oType         Type of OST to run (MD, MC 1-step, MC 2-step).
     * @param dynamics      DynamicsOptions to apply universally.
     * @param ostOpts       OST options to apply.
     * @param fileType      File type to save to.
     * @param repexInterval Interval in psec between repex attempts.
     */
    private RepExOST(OrthogonalSpaceTempering ost, MonteCarloOST mcOST, MolecularDynamics dyn, OstType oType,
                     DynamicsOptions dynamics, OSTOptions ostOpts, CompositeConfiguration properties,
                     String fileType, double repexInterval) throws IOException {
        this.ost = ost;
        switch (oType) {
            case MD:
                algoRun = this::runMD;
                isMC = false;
                break;
            case MC_ONESTEP:
                algoRun = this::runMCOneStep;
                isMC = true;
                break;
            case MC_TWOSTEP:
                algoRun = this::runMCTwoStep;
                isMC = true;
                break;
            default:
                throw new IllegalArgumentException(" Could not recognize whether this is supposed to be MD, MC 1-step, or MC 2-step!");
        }
        this.molDyn = dyn;
        molDyn.setAutomaticWriteouts(false);
        this.dynamics = dynamics;
        this.fileType = fileType;
        this.mcOST = mcOST;
        if (mcOST != null) {
            mcOST.setAutomaticWriteouts(false);
        }
        this.extension = WriteoutOptions.toArchiveExtension(fileType);

        this.world = Comm.world();
        this.rank = world.rank();
        int size = world.size();

        MolecularAssembly[] allAssemblies = molDyn.getAssemblies();
        allFilenames = Arrays.stream(allAssemblies).
                map(MolecularAssembly::getFile).
                map(File::getName).
                map(FilenameUtils::getBaseName).
                toArray(String[]::new);

        File firstFile = allAssemblies[0].getFile();
        basePath = FilenameUtils.getFullPath(firstFile.getAbsolutePath()) + File.separator;
        String baseFileName = FilenameUtils.getBaseName(firstFile.getAbsolutePath());
        dynFile = new File(String.format("%s%d%s%s.dyn", basePath, rank, File.separator, baseFileName));
        molDyn.setFallbackDynFile(dynFile);

        File lambdaFile = new File(String.format("%s%d%s%s.lam", basePath, rank, File.separator, baseFileName));
        currentHistoIndex = rank;
        if (lambdaFile.exists()) {
            try (LambdaReader lr = new LambdaReader(new BufferedReader(new FileReader(lambdaFile)))) {
                lr.readLambdaFile(false);
                currentHistoIndex = lr.getHistogramIndex();
            }
        }
        
        allHistograms = ost.getAllHistograms();
        // TODO: Possibly de-comment this... though HistogramSettings makes it obsolete.
        //Arrays.stream(allHistograms).forEach((OrthogonalSpaceTempering.Histogram h) -> h.setIndependentWrites(true));

        this.numPairs = size - 1;
        this.invKT = -1.0 / (Constants.R * dynamics.getTemp());

        long seed;
        // TODO: Set this per-process individually if we move back to sending accept/reject messages up/down the chain.
        LongBuf seedBuf = LongBuf.buffer(0L);
        if (rank == 0) {
            seed = properties.getLong("randomseed", ThreadLocalRandom.current().nextLong());
            seedBuf.put(0, seed);
            world.broadcast(0, seedBuf);
        } else {
            world.broadcast(0, seedBuf);
            seed = seedBuf.get(0);
        }
        this.random = new Random(seed);

        double timestep = dynamics.getDt() * Constants.FSEC_TO_PSEC;
        stepsBetweenExchanges = Math.max(1, (int) (repexInterval / timestep));

        sends = Arrays.stream(allHistograms).
                map(OrthogonalSpaceTempering.Histogram::getSynchronousSend).
                map(Optional::get).
                toArray(SynchronousSend[]::new);
        if (sends.length < 1) {
            throw new IllegalArgumentException(" No SynchronousSend objects were found!");
        }

        rankToHisto = IntStream.range(0, size).toArray();
        // TODO: Properly back-copy instead of assuming everything is in order at the start.
        histoToRank = Arrays.copyOf(rankToHisto, size);

        Arrays.stream(sends).forEach((SynchronousSend ss) -> ss.setHistograms(allHistograms, rankToHisto));

        totalSwaps = new long[numPairs];
        acceptedSwaps = new long[numPairs];
        Arrays.fill(totalSwaps, 0);
        Arrays.fill(acceptedSwaps, 0);

        setFiles();
        setHistogram(rank);
    }

    private void setHistogram(int index) {
        currentHistoIndex = index;
        ost.switchHistogram(index);
    }

    /**
     * Construct a RepExOST for Monte Carlo orthogonal space tempering.
     *
     * @param ost           An OrthogonalSpaceTempering for each repex rung.
     * @param mcOST         A MonteCarloOST for each repex rung
     * @param dynamics      DynamicsOptions to apply universally.
     * @param ostOpts       OST options to apply.
     * @param fileType      File type to save to.
     * @param twoStep       Whether to use the 2-step MC algorithm (instead of the 1-step).
     * @param repexInterval Interval in psec between repex attempts.
     * @return A RepExOST.
     */
    public static RepExOST repexMC(OrthogonalSpaceTempering ost, MonteCarloOST mcOST,
                                   DynamicsOptions dynamics, OSTOptions ostOpts, CompositeConfiguration properties,
                                   String fileType, boolean twoStep, double repexInterval) throws IOException {
        MolecularDynamics md = mcOST.getMD();
        OstType type = twoStep ? OstType.MC_TWOSTEP : OstType.MC_ONESTEP;
        return new RepExOST(ost, mcOST, md, type, dynamics, ostOpts, properties, fileType, repexInterval);
    }

    /**
     * Construct a RepExOST for Molecular Dynamics orthogonal space tempering.
     *
     * @param ost           An OrthogonalSpaceTempering for each repex rung.
     * @param dyn           A MolecularDynamics for each repex rung.
     * @param dynamics      DynamicsOptions to apply universally.
     * @param ostOpts       OST options to apply.
     * @param fileType      File type to save to.
     * @param repexInterval Interval in psec between repex attempts.
     * @return A RepExOST.
     */
    public static RepExOST repexMD(OrthogonalSpaceTempering ost, MolecularDynamics dyn,
                                   DynamicsOptions dynamics, OSTOptions ostOpts, CompositeConfiguration properties,
                                   String fileType, double repexInterval) throws IOException {
        return new RepExOST(ost, null, dyn, OstType.MD, dynamics, ostOpts, properties, fileType, repexInterval);
    }

    /**
     * Executes the main loop of RepExOST.
     *
     * @throws IOException Possible from Parallel Java.
     */
    public void mainLoop(long numTimesteps, boolean equilibrate) throws IOException {
        if (isMC) {
            mcOST.setEquilibration(equilibrate);
        }
        currentLambda = ost.getLambda();

        Arrays.fill(totalSwaps, 0);
        Arrays.fill(acceptedSwaps, 0);

        if (equilibrate) {
            logger.info(String.format(" Equilibrating repex OST without exchanges on histogram %d.", currentHistoIndex));
            algoRun.accept(numTimesteps);
            reinitVelocities = false;
        } else {
            long numExchanges = numTimesteps / stepsBetweenExchanges;
            for (int i = 0; i < numExchanges; i++) {
                logger.info(String.format(" Beginning of repex loop %d of %d, operating on histogram %d", (i + 1), numExchanges, currentHistoIndex));
                world.barrier(mainLoopTag);
                algoRun.accept(stepsBetweenExchanges);
                ost.logOutputFiles(currentHistoIndex);
                world.barrier(mainLoopTag);
                proposeSwaps((i % 2), 2);
                setFiles();

                long mdMoveNum = i * stepsBetweenExchanges;
                currentLambda = ost.getLambda();
                boolean forceSnapshot = currentLambda >= ost.lambdaWriteOut;
                if (automaticWriteouts) {
                    EnumSet<MolecularDynamics.WriteActions> written = molDyn.writeFilesForStep(mdMoveNum, forceSnapshot, true);
                    if (written.contains(MolecularDynamics.WriteActions.RESTART)) {
                        ost.writeAdditionalRestartInfo(false);
                    }
                }

                reinitVelocities = false;
            }
        }

        logger.info(" Final rank-to-histogram mapping: " + Arrays.toString(rankToHisto));
    }

    // Below logIf methods are intended to reduce redundant logging

    /**
     * Logs a message if this is the master (rank 0) process.
     *
     * @param message Message to log (at INFO)
     */
    private void logIfMaster(String message) {
        logIfMaster(Level.INFO, message);
    }

    /**
     * Logs a message if this is the master (rank 0) process.
     *
     * @param message Message to log
     * @param level   Logging level to log at
     */
    private void logIfMaster(Level level, String message) {
        if (rank == 0) {
            logger.log(level, message);
        }
    }

    /**
     * Logs a message if it is the lowest-ranked process computing a given swap.
     * Currently, this just wraps logIfMaster, as all processes compute all swaps.
     *
     * @param message Message to log (at INFO)
     */
    private void logIfSwapping(String message) {
        logIfMaster(message);
    }
    /**
     *
     * Logs a message if it is the lowest-ranked process computing a given swap.
     * Currently, this just wraps logIfMaster, as all processes compute all swaps.
     *
     * @param message Message to log
     * @param level   Logging level to log at
     */
    private void logIfSwapping(Level level, String message) {
        logIfMaster(level, message);
    }

    private void setFiles() {
        File[] trajFiles = Arrays.stream(allFilenames).
                map((String fn) -> String.format("%s%d%s%s.%s", basePath, currentHistoIndex, File.separator, fn, extension)).
                map(File::new).
                toArray(File[]::new);
        molDyn.setTrajectoryFiles(trajFiles);
    }

    /**
     * Main loop for consistent PRNG-based repex (i.e. every process tests every swap independently).
     * Typically, to create an odd-even staggered schedule (i.e. each pair is tested every other cycle),
     * offset is either 0 or 1, and stride is 2.
     *
     * @param offset Index of the first pair to test swaps for.
     * @param stride Test every nth pair.
     */
    private void proposeSwaps(final int offset, final int stride) {
        for (int i = offset; i < numPairs; i += stride) {
            int rankLow = histoToRank[i];
            int rankHigh = histoToRank[i + 1];
            OrthogonalSpaceTempering.Histogram histoLow = allHistograms[i];
            OrthogonalSpaceTempering.Histogram histoHigh = allHistograms[i+1];
            
            double lamLow = histoLow.getLastReceivedLambda();
            double dUdLLow = histoLow.getLastReceivedDUDL();
            double lamHigh = histoHigh.getLastReceivedLambda();
            double dUdLHigh = histoHigh.getLastReceivedDUDL();

            double eii = histoLow.computeBiasEnergy(lamLow, dUdLLow);
            double eij = histoLow.computeBiasEnergy(lamHigh, dUdLHigh);
            double eji = histoHigh.computeBiasEnergy(lamLow, dUdLLow);
            double ejj = histoHigh.computeBiasEnergy(lamHigh, dUdLHigh);

            logIfSwapping(String.format("\n Proposing exchange between histograms %d (rank %d) and %d (rank %d).\n" +
                            " Li: %.6f dU/dLi: %.6f Lj: %.6f dU/dLj: %.6f",
                    i, rankLow, i + 1, rankHigh,
                    lamLow, dUdLLow, lamHigh, dUdLHigh));

            double e1 = eii + ejj;
            double e2 = eji + eij;
            boolean accept = BoltzmannMC.evaluateMove(random, invKT, e1, e2);
            double acceptChance = BoltzmannMC.acceptChance(invKT, e1, e2);

            String desc = accept ? "Accepted" : "Rejected";
            logIfSwapping(String.format(" %s exchange with probability %.5f based on Eii %.6f, Ejj %.6f, Eij %.6f, Eji %.6f kcal/mol",
                    desc, acceptChance, eii, ejj, eij, eji));

            ++totalSwaps[i];
            if (accept) {
                ++acceptedSwaps[i];
                switchHistos(rankLow, rankHigh, i);
            }

            double acceptRate = ((double) acceptedSwaps[i]) / ((double) totalSwaps[i]);
            logIfSwapping(String.format(" Replica exchange acceptance rate for pair %d-%d is %.3f%%", i, (i+1), acceptRate * 100));
        }
    }

    private void switchHistos(int rankLow, int rankHigh, int histoLow) {
        int histoHigh = histoLow + 1;
        rankToHisto[rankLow] = histoHigh;
        rankToHisto[rankHigh] = histoLow;
        histoToRank[histoLow] = rankHigh;
        histoToRank[histoHigh] = rankLow;
        setHistogram(rankToHisto[rank]);

        ost.setLambda(currentLambda);
        /* TODO: If there is ever a case where an algorithm will not update coordinates itself at the start, we have to
         * update coordinates here (from the OST we used to be running on to the new OST). */

        for (SynchronousSend send : sends) {
            send.updateRanks(rankToHisto);
        }
    }

    /**
     * Run 1-step MC-OST for the specified number of MD steps.
     *
     * @param numSteps Number of MD steps (not MC cycles) to run.
     */
    private void runMCOneStep(long numSteps) {
        mcOST.setRunLength(numSteps);
        mcOST.sampleOneStep();
    }

    /**
     * Run 2-step MC-OST for the specified number of MD steps.
     *
     * @param numSteps Number of MD steps (not MC cycles) to run.
     */
    private void runMCTwoStep(long numSteps) {
        mcOST.setRunLength(numSteps);
        mcOST.sampleTwoStep();
    }

    /**
     * Run MD for the specified number of steps.
     *
     * @param numSteps MD steps to run.
     */
    private void runMD(long numSteps) {
        molDyn.dynamic(numSteps, dynamics.getDt(), dynamics.getReport(), dynamics.getSnapshotInterval(), dynamics.getTemp(),
                reinitVelocities, fileType, dynamics.getCheckpoint(), dynFile);
    }

    public OrthogonalSpaceTempering getOST() {
        return ost;
    }

    private enum OstType {
        MD, MC_ONESTEP, MC_TWOSTEP
    }
}
