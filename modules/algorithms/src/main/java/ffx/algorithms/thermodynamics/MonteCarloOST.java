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

import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.lang.System.nanoTime;

import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.abs;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.LambdaMove;
import ffx.algorithms.mc.MDMove;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.EnergyException;
import static ffx.utilities.Constants.NS2SEC;

/**
 * Sample a thermodynamic path using the OST method, with the time-dependent
 * bias built up using Metropolis Monte Carlo steps.
 * <p>
 * The algorithm generates coordinate (X) MC moves using molecular dynamics at a
 * fixed lambda value (i.e. using OpenMM), followed by MC lambda moves.
 * <p>
 * 1.) At a fixed Lambda, run a defined length MD trajectory to "move"
 * coordinates and dU/dL on an approximate potential U* (i.e. no OST Bias).
 * <p>
 * 2.) Accept / Reject the MD move with probability exp[-Beta(dU - dU*)] where
 * dU is the change in AMOEBA + Bias energy and dU* is the change in AMOEBA +
 * Kinetic energy from the MD.
 * <p>
 * 3.) Randomly change the value of Lambda.
 * <p>
 * 4.) Accept / Reject the Lambda move using the AMOEBA + OST Bias energy.
 * <p>
 * 5.) Add to the time dependent 2D bias using the current values of Lambda and
 * dU/dL.
 *
 * @author Michael J. Schnieders
 * @author Hernan Beranbe
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class MonteCarloOST extends BoltzmannMC {

    /**
     * Logger object to print out information for this class.
     */
    private static final Logger logger = Logger.getLogger(MonteCarloOST.class.getName());
    /**
     * Energy conservation during MD moves should generally be within ~0.1
     * kcal/mol. A change in total energy of 1.0 kcal/mol or more is of
     * significant concern that the time step is too large, or lambda moves are
     * too aggressive.
     */
    private static final double ENERGY_CONSERVATION_TOLERANCE = 10.0;
    /**
     * Potential object used to retrieve the coordinates for the system.
     */
    private final Potential potential;
    /**
     * OST object used to retrieve OST energy throughout the simulation.
     */
    private final OrthogonalSpaceTempering orthogonalSpaceTempering;
    /**
     * MDMove object for completing MC-OST molecular dynamics moves.
     */
    private final MDMove mdMove;
    /**
     * Total number of steps to take for MC-OST sampling.
     */
    private long totalSteps;
    /**
     * Number of steps to take per MC-OST round.
     */
    private long stepsPerMove;
    /**
     * Lambda move object for completing MC-OST lambda moves.
     */
    private LambdaMove lambdaMove;
    /**
     * Double that keeps track of our lambda value.
     */
    private double lambda = 1.0;
    /**
     * Boolean that tells algorithm that we are in the equilibration phase of MC-OST.
     */
    private boolean equilibration = false;

    /**
     * Controls the effect of verbose by logging at FINE vs. INFO.
     */
    private final Level verboseLoggingLevel;
    /**
     * How verbose MD should be.
     */
    private final MolecularDynamics.VerbosityLevel mdVerbosityLevel;
    /**
     * Deposit a bias once every N MC cycles. Defaults to 1.
     */
    private final int biasDepositionFrequency;
    /**
     * Only print out logging information every N MC cycles. Defaults to 1.
     */
    private final int loggingFrequency;

    /**
     * <p>
     * Constructor for MonteCarloOST.</p>
     *
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param orthogonalSpaceTempering                a {@link OrthogonalSpaceTempering} object.
     * @param molecularAssembly   a {@link ffx.potential.MolecularAssembly} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param dynamics            CLI object containing key information.
     * @param verbose             Whether to be verbose.
     * @param cycleLength         Length of an MC cycle in MD steps.
     */
    public MonteCarloOST(Potential potentialEnergy, OrthogonalSpaceTempering orthogonalSpaceTempering,
                         MolecularAssembly molecularAssembly, CompositeConfiguration properties,
                         AlgorithmListener listener, DynamicsOptions dynamics, boolean verbose,
                         int cycleLength) {
        this.potential = potentialEnergy;
        this.orthogonalSpaceTempering = orthogonalSpaceTempering;
        verboseLoggingLevel = verbose ? Level.INFO : Level.FINE;
        mdVerbosityLevel = verbose ? MolecularDynamics.VerbosityLevel.QUIET : MolecularDynamics.VerbosityLevel.SILENT;
        stepsPerMove = cycleLength;
        totalSteps = dynamics.getNumSteps();

        ThermostatEnum tstat = dynamics.thermostat;
        if (!tstat.equals(ThermostatEnum.ADIABATIC)) {
            logger.warning(format("MC-OST requires the ADIABATIC thermostat, found %s.", tstat));
            dynamics.setThermostat(ThermostatEnum.ADIABATIC);
        }

        IntegratorEnum integ = dynamics.integrator;
        if (!integ.knownReversible || !integ.knownDeterministic) {
            throw new IllegalArgumentException(format("MC-OST requires " +
                    "a reversible deterministic integrator (e.g. VERLET, RESPA), found %s!", integ));
        }

        // Create the MC MD and Lambda moves.
        mdMove = new MDMove(molecularAssembly, potential, properties, listener, dynamics, stepsPerMove);
        if (properties.containsKey("randomseed")) {
            int randomSeed = properties.getInt("randomseed", 0);
            logger.info(format(" Setting random seed for lambdaMove to %d ", randomSeed));
            lambdaMove = new LambdaMove(randomSeed, orthogonalSpaceTempering);
            setRandomSeed(randomSeed);
        } else {
            lambdaMove = new LambdaMove(orthogonalSpaceTempering);
        }

        // Changing the value of lambda will be handled by this class, as well as adding the time dependent bias.
        orthogonalSpaceTempering.setPropagateLambda(false);
        biasDepositionFrequency = properties.getInt("mc-ost-biasf", 1);
        if (biasDepositionFrequency < 1) {
            throw new IllegalArgumentException("The property mc-ost-biasf must be a positive integer, found " + biasDepositionFrequency + " !");
        } else if (biasDepositionFrequency > 1) {
            logger.info(format(" MC-OST will deposit a bias only once per %d MC cycles (mc-ost-biasf).", biasDepositionFrequency));
        }

        loggingFrequency = properties.getInt("mc-ost-logf", 1);
    }

    /**
     * Takes in parameters and calls the MDMove method setMDParameters to update
     * the stepsPerMove and timeStep parameters to the current value in this
     * class
     *
     * @param totalSteps   a int.
     * @param stepsPerMove a int.
     * @param mcMDE        a boolean
     */
    public void setMDMoveParameters(long totalSteps, int stepsPerMove, boolean mcMDE) {

        if (mcMDE) {
            if (equilibration) {
                this.stepsPerMove = (int) Math.round(stepsPerMove * 0.1);
            } else {
                mdMove.setMDIntervalSteps(stepsPerMove);
            }
        }
        this.totalSteps = totalSteps;
    }

    /**
     * Calls on LambdaMove class method setLambdaStdDev to update the lambda
     * standard deviation to the current value in this class
     *
     * @param stdDev a double.
     */
    public void setLambdaStdDev(double stdDev) {
        lambdaMove.setStdDev(stdDev);
    }

    /**
     * Sets the value of the boolean equilibration variables to true or false to
     * either allow an equilibration step or skip it.
     *
     * @param equilibration a boolean.
     */
    public void setEquilibration(boolean equilibration) {
        this.equilibration = equilibration;
    }

    /**
     * Calls on the OST method set lambda to update lambda to the current value
     * in this class
     *
     * @param lambda a double.
     */
    public void setLambda(double lambda) {
        this.lambda = lambda;
        orthogonalSpaceTempering.setLambda(lambda);
    }

    /**
     * Returns the current value of lambda
     *
     * @return lambda
     */
    public double getLambda() {
        return lambda;
    }

    public void setLambdaWriteOut(double lambdaWriteOut) {
        orthogonalSpaceTempering.setLambdaWriteOut(lambdaWriteOut);
    }

    /**
     * Only log stuff every N steps.
     * @param i       MC step
     * @param level   Logging level
     * @param message Message to log
     */
    private void logForStep(int i, Level level, String message) {
        if (i % loggingFrequency == 0) {
            logger.log(level, message);
        }
    }

    /**
     * The goal is to sample lambda and coordinates (X) separately to converge
     * the ensemble average dU/dL for every state (lambda) along the
     * thermodynamic path.
     * <p>
     * 1.) At a fixed lambda, run a defined length MD trajectory to "move"
     * coordinates and dU/dL.
     * <p>
     * 2.) Accept / Reject the MD move using the total Hamiltonian (Kinetic
     * energy + OST energy).
     * <p>
     * 3.) Randomly change the value of Lambda.
     * <p>
     * 4.) Accept / Reject the Lambda move using the OST energy.
     * <p>
     * 5.) Add to the bias.
     */
    public void sampleTwoStep() {

        int n = potential.getNumberOfVariables();
        double[] gradient = new double[n];
        double[] currentCoordinates = new double[n];
        double[] proposedCoordinates = new double[n];
        long numMoves = totalSteps / stepsPerMove;
        int acceptLambda = 0;
        int acceptMD = 0;

        // Initialize the current coordinates.
        potential.getCoordinates(currentCoordinates);

        // Update time dependent bias.
        Histogram histogram = orthogonalSpaceTempering.getHistogram();

        // Compute the current OST potential energy.
        double currentOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);

        // Collect the current dU/dL, Force Field Energy and Bias Energy.
        double currentdUdL = orthogonalSpaceTempering.getForceFielddEdL();
        double currentForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
        double currentBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

        // Initialize MC move instances.
        for (int imove = 0; imove < numMoves; imove++) {
            long totalMoveTime = -nanoTime();
            long mdMoveAndEvalTime = -nanoTime();

            if (equilibration) {
                logForStep(imove, Level.INFO, format("\n MD Equilibration Round %d", imove + 1));
            } else {
                logForStep(imove, Level.INFO, format("\n MC Orthogonal Space Sampling Round %d: Independent Steps", imove + 1));
            }

            // Run MD in an approximate potential U* (U star) that does not include the OST bias.
            long mdMoveTime = -nanoTime();
            mdMove.move(mdVerbosityLevel);
            mdMoveTime += nanoTime();
            logForStep(imove, verboseLoggingLevel, format("  Total time for MD move: %6.3f", mdMoveTime * NS2SEC));

            // Get the starting and final kinetic energy for the MD move.
            double currentKineticEnergy = mdMove.getInitialKinetic();
            double proposedKineticEnergy = mdMove.getKineticEnergy();

            // Get the new coordinates.
            potential.getCoordinates(proposedCoordinates);

            // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
            long proposedOSTEnergyTime = -nanoTime();

            double proposedOSTEnergy;
            try {
                proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(proposedCoordinates, gradient);
            } catch (EnergyException e) {
                mdMove.revertMove();
                logForStep(imove, Level.INFO, " Unstable MD Move skipped.");
                continue;
            }
            proposedOSTEnergyTime += nanoTime();

            logForStep(imove, Level.FINE, format("  Time to complete MD OST energy method call %6.3f", proposedOSTEnergyTime * NS2SEC));

            // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
            double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
            double proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
            double proposedBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

            // The Metropolis criteria is based on the sum of the OST Energy and Kinetic Energy.
            double currentTotalEnergy = currentOSTEnergy + currentKineticEnergy;
            double proposedTotalEnergy = proposedOSTEnergy + proposedKineticEnergy;

            logForStep(imove, verboseLoggingLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logForStep(imove, verboseLoggingLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logForStep(imove, verboseLoggingLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logForStep(imove, verboseLoggingLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            double energyChange = mdMove.getEnergyChange();
            if (abs(energyChange) > ENERGY_CONSERVATION_TOLERANCE) {
                mdMove.revertMove();
                logger.warning(" MC Move skipped due to lack of MD energy conservation");
                continue;
            }

            if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(), proposeddUdL) &&
                    evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                // Accept MD move.
                acceptMD++;
                double percent = (acceptMD * 100.0) / (imove + 1);
                logForStep(imove, Level.INFO, format(" Accept [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                currentOSTEnergy = proposedOSTEnergy;
                currentdUdL = proposeddUdL;
                currentForceFieldEnergy = proposedForceFieldEnergy;
                currentBiasEnergy = proposedBiasEnergy;
                arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
            } else {
                double percent = (acceptMD * 100.0) / (imove + 1);
                logForStep(imove, Level.INFO, format(" Reject [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                mdMove.revertMove();
            }
            mdMoveAndEvalTime += nanoTime();

            logForStep(imove, Level.FINE, format("\n  Total time to run and evaluate MD move: %6.3f", mdMoveAndEvalTime * NS2SEC));

            // During equilibration, do not change Lambda or contribute to the OST bias.
            if (!equilibration) {
                // Update Lambda.
                logForStep(imove, Level.INFO, " MC Lambda Step");

                long lambdaMoveTime = -nanoTime();
                double currentLambda = orthogonalSpaceTempering.getLambda();
                lambdaMove.move();
                double proposedLambda = orthogonalSpaceTempering.getLambda();

                // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
                long proposedOSTEnergyTime2 = -nanoTime();
                proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);
                proposedOSTEnergyTime2 += nanoTime();

                logForStep(imove, verboseLoggingLevel, format("  Time to complete Lambda OST energy method call %6.3f ", proposedOSTEnergyTime2 * NS2SEC));

                // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
                proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
                proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();

                logForStep(imove, verboseLoggingLevel, format("\n  Current  OST     %12.3f at L=%5.3f.", currentOSTEnergy, currentLambda));
                logForStep(imove, verboseLoggingLevel, format("  Proposed OST     %12.3f at L=%5.3f.", proposedOSTEnergy, proposedLambda));
                logForStep(imove, verboseLoggingLevel, format("  MC Energy change: %12.3f (kcal/mol).", proposedOSTEnergy - currentOSTEnergy));

                if (orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL) &&
                        evaluateMove(currentOSTEnergy, proposedOSTEnergy)) {
                    acceptLambda++;
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format("  Accept [ L=%8.3f,E=%12.4f]   -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSTEnergy, proposedLambda, proposedOSTEnergy, percent));
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    currentdUdL = proposeddUdL;
                    lambda = proposedLambda;
                } else {
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format("  Reject [ L=%8.3f,E=%12.4f]   -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSTEnergy, proposedLambda, proposedOSTEnergy, percent));
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                }

                lambdaMoveTime += nanoTime();
                logForStep(imove, verboseLoggingLevel, format("  Lambda move completed in %6.3f", lambdaMoveTime * NS2SEC));

                if (imove % biasDepositionFrequency == 0) {
                    histogram.addBias(currentdUdL, currentCoordinates, null);
                } else {
                    // TODO: Step down to FINE when we know this works.
                    logForStep(imove, Level.INFO, format(" Cycle %d: skipping bias deposition.", imove));
                }

                logForStep(imove, verboseLoggingLevel, format("  Added Bias at [L=%5.3f, FL=%9.3f]", lambda, currentdUdL));

                // Compute the updated OST bias.
                currentBiasEnergy = histogram.computeBiasEnergy(lambda, currentdUdL);

                // Update the current OST Energy to be the sum of the current Force Field Energy and updated OST Bias.
                currentOSTEnergy = currentForceFieldEnergy + currentBiasEnergy;

                if (lambda >= orthogonalSpaceTempering.lambdaWriteOut) {
                    long mdMoveNum = imove * stepsPerMove;
                    mdMove.writeFilesForStep(mdMoveNum);
                }
            }

            totalMoveTime += nanoTime();
            logForStep(imove, Level.INFO, format(" Round complete in %6.3f sec.", totalMoveTime * NS2SEC));
        }
    }

    private double singleStepLambda() {
        lambdaMove.move();
        double proposedLambda = orthogonalSpaceTempering.getLambda();
        logger.log(verboseLoggingLevel, format(" Proposed lambda: %5.3f.", proposedLambda));
        return proposedLambda;
    }

    private void singleStepMD() {
        // Run MD in an approximate potential U* (U star) that does not include the OST bias.
        long mdMoveTime = -nanoTime();
        mdMove.move(mdVerbosityLevel);
        mdMoveTime += nanoTime();
        logger.log(verboseLoggingLevel, format(" Total time for MD move: %6.3f", mdMoveTime * NS2SEC));
    }

    /**
     * The goal is to sample lambda and coordinates (X) simultaneously to
     * converge the ensemble average dU/dL for every state (lambda) along the
     * thermodynamic path.
     * <p>
     * 1.) Randomly change the value of Lambda.
     * <p>
     * 2.) At the proposed lambda, run a defined length MD trajectory to "move"
     * coordinates and dU/dL.
     * <p>
     * 3.) Accept / Reject the Lambda + MD move using the total Hamiltonian
     * (Kinetic energy + OST energy).
     * <p>
     * 4.) Add to the bias.
     */
    public void sampleOneStep() {

        int n = potential.getNumberOfVariables();
        double[] gradient = new double[n];
        double[] currentCoordinates = new double[n];
        double[] proposedCoordinates = new double[n];
        long numMoves = totalSteps / stepsPerMove;
        int acceptMD = 0;
        int acceptMCOST = 0;

        // Initialize the current coordinates.
        potential.getCoordinates(currentCoordinates);

        // Update time-dependent bias.
        Histogram histogram = orthogonalSpaceTempering.getHistogram();

        // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
        double currentOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);

        // Retrieve the computed dU/dL, Force Field Energy and Bias Energy.
        double currentdUdL = orthogonalSpaceTempering.getForceFielddEdL();
        double currentForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
        double currentBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

        // Initialize MC move instances.
        for (int imove = 0; imove < numMoves; imove++) {
            long totalMoveTime = -nanoTime();

            double currentLambda = orthogonalSpaceTempering.getLambda();
            double proposedLambda = currentLambda;

            boolean lambdaFirst = random.nextBoolean();

            if (equilibration) {
                logForStep(imove, Level.INFO, format("\n Equilibration Round %d", imove + 1));
            } else {
                String moveType = lambdaFirst ? "Single-step lambda plus MD move." : " Single-step MD plus lambda move.";
                logForStep(imove, Level.INFO, format("\n MC-OST Round %d: %s", imove + 1, moveType));
            }

            logForStep(imove, Level.FINE, format(" Starting force field energy for move %16.8f", currentForceFieldEnergy));

            if (equilibration) {
                singleStepMD();
            } else {
                if (lambdaFirst) {
                    proposedLambda = singleStepLambda();
                    singleStepMD();
                } else {
                    singleStepMD();
                    proposedLambda = singleStepLambda();
                }
            }

            // Get the starting and final kinetic energy for the MD move.
            double currentKineticEnergy = mdMove.getInitialKinetic();
            double proposedKineticEnergy = mdMove.getKineticEnergy();

            // Get the new coordinates.
            potential.getCoordinates(proposedCoordinates);

            long proposedOSTEnergyTime = -nanoTime();

            // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
            double proposedOSTEnergy;
            try {
                proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(proposedCoordinates, gradient);
            } catch (EnergyException e) {
                mdMove.revertMove();
                mdMove.writeErrorFiles();
                if (!equilibration) {
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                }
                logger.log(Level.WARNING, " Unstable MD Move skipped.");
                continue;
            }

            proposedOSTEnergyTime += nanoTime();
            logForStep(imove, Level.FINE, format(" Time to complete MD OST energy method call %6.3f", proposedOSTEnergyTime * NS2SEC));

            // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
            double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
            double proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
            double proposedBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

            // The Metropolis criteria is based on the sum of the OST Energy and Kinetic Energy.
            double currentTotalEnergy = currentOSTEnergy + currentKineticEnergy;
            double proposedTotalEnergy = proposedOSTEnergy + proposedKineticEnergy;

            logForStep(imove, verboseLoggingLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logForStep(imove, verboseLoggingLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logForStep(imove, verboseLoggingLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logForStep(imove, verboseLoggingLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            double energyChange = mdMove.getEnergyChange();
            if (abs(energyChange) > ENERGY_CONSERVATION_TOLERANCE) {
                mdMove.revertMove();
                if (!equilibration) {
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                }
                logger.warning(" MC Move skipped due to lack of MD energy conservation.");
                continue;
            }

            if (equilibration) {
                if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(), proposeddUdL) &&
                        evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                    // Accept MD.
                    acceptMD++;
                    double percent = (acceptMD * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format(" Accept [ FL=%8.3f, E=%12.4f]  -> [ FL=%8.3f, E=%12.4f] (%5.1f%%)",
                            currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                    currentOSTEnergy = proposedOSTEnergy;
                    currentdUdL = proposeddUdL;
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    currentBiasEnergy = proposedBiasEnergy;
                    arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
                } else {
                    double percent = (acceptMD * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format(" Reject [ FL=%8.3f, E=%12.4f]  -> [ FL=%8.3f, E=%12.4f] (%5.1f%%)",
                            currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                    mdMove.revertMove();
                }
            } else {
                if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(), proposeddUdL) &&
                        evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                    acceptMCOST++;
                    double percent = (acceptMCOST * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format(" Accept [ L=%5.3f, FL=%8.3f, E=%12.4f]  -> [ L=%5.3f, FL=%8.3f, E=%12.4f] (%5.1f%%)",
                            currentLambda, currentdUdL, currentOSTEnergy, proposedLambda, proposeddUdL, proposedOSTEnergy, percent));
                    lambda = proposedLambda;
                    currentdUdL = proposeddUdL;
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
                } else {
                    double percent = (acceptMCOST * 100.0) / (imove + 1);
                    logForStep(imove, Level.INFO, format(" Reject [ L=%5.3f, FL=%8.3f, E=%12.4f]  -> [ L=%5.3f, FL=%8.3f, E=%12.4f] (%5.1f%%)",
                            currentLambda, currentdUdL, currentOSTEnergy, proposedLambda, proposeddUdL, proposedOSTEnergy, percent));
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                    mdMove.revertMove();
                }

                if (imove % biasDepositionFrequency == 0) {
                    histogram.addBias(currentdUdL, currentCoordinates, null);
                } else {
                    logForStep(imove, Level.FINE, format(" Cycle %d: skipping bias deposition.", imove));
                }

                logForStep(imove, Level.FINE, format(" Added Bias at [ L=%5.3f, FL=%9.3f]", lambda, currentdUdL));

                // Compute the updated OST bias.
                currentBiasEnergy = histogram.computeBiasEnergy(lambda, currentdUdL);

                // Update the current OST Energy to be the sum of the current Force Field Energy and updated OST Bias.
                currentOSTEnergy = currentForceFieldEnergy + currentBiasEnergy;

                if (lambda >= orthogonalSpaceTempering.lambdaWriteOut) {
                    long mdMoveNum = imove * stepsPerMove;
                    EnumSet<MolecularDynamics.WriteActions> written = mdMove.writeFilesForStep(mdMoveNum);
                    if (written.contains(MolecularDynamics.WriteActions.RESTART)) {
                        orthogonalSpaceTempering.writeAdditionalRestartInfo(false);
                    }
                }
            }

            totalMoveTime += nanoTime();
            logForStep(imove, Level.INFO, format(" Round complete in %6.3f sec.", totalMoveTime * NS2SEC));
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double currentEnergy() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void storeState() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertStep() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
