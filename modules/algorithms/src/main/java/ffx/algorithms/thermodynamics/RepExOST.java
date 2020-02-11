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

import edu.rit.mp.BooleanBuf;
import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;

import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.utilities.Constants;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Optional;
import java.util.Random;
import java.util.function.IntConsumer;
import java.util.function.LongConsumer;

/**
 * An implementation of RepEx between Orthogonal Space Tempering potentials.
 *
 * @author Michael J. Schnieders
 * @author Jacob Litman
 * @since 1.0
 */
public class RepExOST {
    private final OrthogonalSpaceTempering[] osts;
    private final SynchronousSend[] sends;
    private final LongConsumer algoRun;
    private final MolecularDynamics[] molDyns;
    private final File[] dynFiles;
    private final DynamicsOptions dynamics;
    private final String fileType;
    private final MonteCarloOST[] mcOSTs;
    private final long stepsBetweenExchanges;
    private final Comm world;
    private final int rank;
    private final int numPairs;
    private int currentOST;
    private int[] rankToHisto;
    private final boolean isMC;
    private boolean reinitVelocities = true;

    private final Random random;
    private final double invKT;
    
    private DoubleBuf lamBuf = DoubleBuf.buffer(0);
    private BooleanBuf swapBuf = BooleanBuf.buffer(false);
    // Message tags to use.
    private static final int lamTag = 42;
    private static final int mainLoopTag = 2020;

    /**
     * Private constructor used here to centralize shared logic.
     *
     * @param osts          An OrthogonalSpaceTempering for each repex rung.
     * @param mcOSTs        A MonteCarloOST for each repex rung, or null (for MD).
     * @param dyns          A MolecularDynamics for each repex rung (never null).
     * @param oType         Type of OST to run (MD, MC 1-step, MC 2-step).
     * @param dynamics      DynamicsOptions to apply universally.
     * @param fileType      File type to save to.
     * @param repexInterval Interval in psec between repex attempts.
     */
    private RepExOST(OrthogonalSpaceTempering[] osts, MonteCarloOST[] mcOSTs, MolecularDynamics[] dyns, OstType oType,
                     DynamicsOptions dynamics, String fileType, double repexInterval) {
        this.osts = osts;
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
        this.molDyns = dyns;
        this.dynFiles = Arrays.stream(molDyns).map(MolecularDynamics::getDynFile).toArray(File[]::new);
        this.dynamics = dynamics;
        this.fileType = fileType;
        this.mcOSTs = mcOSTs;

        this.world = Comm.world();
        this.rank = world.rank();
        int size = world.size();

        this.numPairs = size - 1;
        this.random = new Random();
        this.invKT = -1.0 / (Constants.R * dynamics.getTemp());

        double timestep = dynamics.getDt() * Constants.FSEC_TO_PSEC;
        stepsBetweenExchanges = Math.max(1, (int) (repexInterval / timestep));

        sends = Arrays.stream(osts).
                map(OrthogonalSpaceTempering::getHistogram).
                map(OrthogonalSpaceTempering.Histogram::getSynchronousSend).
                map(Optional::get).
                toArray(SynchronousSend[]::new);
        rankToHisto = sends[0].getRankToHistogramMap();
    }

    /**
     * Construct a RepExOST for Monte Carlo orthogonal space tempering.
     *
     * @param osts          An OrthogonalSpaceTempering for each repex rung.
     * @param mcOSTs        A MonteCarloOST for each repex rung
     * @param dynamics      DynamicsOptions to apply universally.
     * @param fileType      File type to save to.
     * @param twoStep       Whether to use the 2-step MC algorithm (instead of the 1-step).
     * @param repexInterval Interval in psec between repex attempts.
     * @return              A RepExOST.
     */
    public static RepExOST repexMC(OrthogonalSpaceTempering[] osts, MonteCarloOST[] mcOSTs,
                                   DynamicsOptions dynamics, String fileType, boolean twoStep, double repexInterval) {
        MolecularDynamics[] dyns = Arrays.stream(mcOSTs).map(MonteCarloOST::getMD).toArray(MolecularDynamics[]::new);
        OstType type = twoStep ? OstType.MC_TWOSTEP : OstType.MC_ONESTEP;
        return new RepExOST(osts, mcOSTs, dyns, type, dynamics, fileType, repexInterval);
    }

    /**
     * Construct a RepExOST for Molecular Dynamics orthogonal space tempering.
     *
     * @param osts          An OrthogonalSpaceTempering for each repex rung.
     * @param dyns          A MolecularDynamics for each repex rung.
     * @param dynamics      DynamicsOptions to apply universally.
     * @param fileType      File type to save to.
     * @param repexInterval Interval in psec between repex attempts.
     * @return              A RepExOST.
     */
    public static RepExOST repexMD(OrthogonalSpaceTempering[] osts, MolecularDynamics[] dyns,
                                   DynamicsOptions dynamics, String fileType, double repexInterval) {
        return new RepExOST(osts, null, dyns, OstType.MD, dynamics, fileType, repexInterval);
    }

    /**
     * Executes the main loop of RepExOST.
     * @throws IOException Possible from Parallel Java.
     */
    public void mainLoop(long numTimesteps, boolean equilibrate) throws IOException {
        long numExchanges = numTimesteps / stepsBetweenExchanges;
        if (isMC) {
            Arrays.stream(mcOSTs).forEach((MonteCarloOST mco) -> mco.setEquilibration(equilibrate));
        }
        for (int i = 0; i < numExchanges; i++) {
            world.barrier(mainLoopTag);
            algoRun.accept(stepsBetweenExchanges);
            for (int j = 0; j < numPairs; j++) {
                if (j == rank) {
                    proposeSwapUp();
                } else if (j == (rank - 1)) {
                    proposeSwapDown();
                } else {
                    listenSwap(j);
                }
            }
            reinitVelocities = false;
        }
    }

    /**
     * Attempts to exchange histograms via a Boltzmann trial. Implemented such
     * that every replica has a copy of every histogram, so just a pointer must
     * be updated.
     * @throws IOException If Parallel Java has an issue
     */
    private void proposeSwapUp() throws IOException {
        int rankUp = rank + 1;

        world.receive(rankUp, lamTag, lamBuf);
        double otherLam = lamBuf.get(0);
        world.receive(rankUp, lamTag, lamBuf);
        double otherDUDL = lamBuf.get(0);

        int ostIndex = rankToHisto[rank];
        OrthogonalSpaceTempering currOST = osts[ostIndex];
        double currentLambda = currOST.getLambda();
        double currentDUDL = currOST.getdEdL();
        OrthogonalSpaceTempering.Histogram currentHistogram = osts[ostIndex].getHistogram();

        double eii = currentHistogram.computeBiasEnergy(currentLambda, currentDUDL);
        double eij = currentHistogram.computeBiasEnergy(otherLam, otherDUDL);

        ostIndex = rankToHisto[rankUp];
        OrthogonalSpaceTempering.Histogram otherHisto = osts[ostIndex].getHistogram();
        double eji = otherHisto.computeBiasEnergy(currentLambda, currentDUDL);
        double ejj = otherHisto.computeBiasEnergy(otherLam, otherDUDL);

        testSwap(eii, eij, eji, ejj, rankUp);
    }

    /**
     * Sends information to the rank below for a proposed replica exchange.
     * @throws IOException If Parallel Java has an issue
     */
    private void proposeSwapDown() throws IOException {
        int rankDown = rank - 1;

        int ostIndex = rankToHisto[rank];
        OrthogonalSpaceTempering currOST = osts[ostIndex];
        double currentLambda = currOST.getLambda();
        double currentDUDL = currOST.getdEdL();

        lamBuf.put(0, currentLambda);
        world.send(rankDown, lamTag, lamBuf);
        lamBuf.put(0, currentDUDL);
        world.send(rankDown, lamTag, lamBuf);

        listenSwap(rankDown);
    }

    /**
     * Listen to rootRank to check if the swap succeeded or is rejected.
     *
     * @param rootRank     Rank checking if the swap succeeded (lower half)
     * @throws IOException For Parallel Java
     */
    private void listenSwap(int rootRank) throws IOException {
        world.broadcast(rootRank, swapBuf);
        if (swapBuf.get(0)) {
            acceptSwap(rootRank);
        }
    }

    /**
     * Apply the Metropolis criterion to a proposed replica exchange.
     *
     * @param eii          Bias energy i of L/FL i
     * @param eij          Bias energy i of L/FL j
     * @param eji          Bias energy j of L/FL i
     * @param ejj          Bias energy j of L/FL j
     * @param rankUp       Rank above this process
     * @throws IOException For Parallel Java
     */
    private void testSwap(double eii, double eij, double eji, double ejj, int rankUp) throws IOException {
        double e1 = eii + ejj;
        double e2 = eji + eij;
        // I'm assuming this class extends BoltzmannMC here.
        boolean accept = BoltzmannMC.evaluateMove(random, invKT, e1, e2);
        swapBuf.put(0, accept);

        // All other processes should be calling comm.broadcast via the listenSwap method.
        world.broadcast(rank, swapBuf);
        if (accept) {
            acceptSwap(rank);
        }
    }

    /**
     * Accept a swap: called by every process if a swap is accepted.
     *
     * @param rootRank Lower end of the swap.
     */
    private void acceptSwap(int rootRank) {
        int histoDown = rankToHisto[rootRank];
        int histoUp = rankToHisto[rootRank + 1];
        rankToHisto[rootRank] = histoUp;
        rankToHisto[rootRank + 1] = histoDown;
        currentOST = rankToHisto[rank];

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
        mcOSTs[currentOST].setRunLength(numSteps);
        mcOSTs[currentOST].sampleOneStep();
    }

    /**
     * Run 2-step MC-OST for the specified number of MD steps.
     *
     * @param numSteps Number of MD steps (not MC cycles) to run.
     */
    private void runMCTwoStep(long numSteps) {
        MonteCarloOST currMC = mcOSTs[currentOST];
        currMC.setRunLength(numSteps);
        currMC.sampleTwoStep();
    }

    /**
     * Run MD for the specified number of steps.
     * @param numSteps MD steps to run.
     */
    private void runMD(long numSteps) {
        MolecularDynamics molDyn = molDyns[currentOST];
        File dyn = dynFiles[currentOST];
        molDyn.dynamic(numSteps, dynamics.getDt(), dynamics.getReport(), dynamics.getSnapshotInterval(), dynamics.getTemp(),
                reinitVelocities, fileType, dynamics.getCheckpoint(), dyn);
    }

    public OrthogonalSpaceTempering getCurrentOST() {
        return osts[currentOST];
    }

    int getCurrentIndex() {
        return currentOST;
    }

    public OrthogonalSpaceTempering[] getAllOST() {
        return Arrays.copyOf(osts, osts.length);
    }

    private enum OstType {
        MD, MC_ONESTEP, MC_TWOSTEP
    }
}
