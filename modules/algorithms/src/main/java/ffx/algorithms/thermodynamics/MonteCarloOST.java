// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.algorithms.thermodynamics;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MDVerbosity;
import ffx.algorithms.dynamics.MDWriteAction;
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
import org.apache.commons.configuration2.CompositeConfiguration;

import java.io.File;
import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.lang.System.nanoTime;
import static org.apache.commons.math3.util.FastMath.abs;

/**
 * Sample a thermodynamic path using the OST method, with the time-dependent bias built up using
 * Metropolis Monte Carlo steps.
 *
 * <p>The algorithm generates coordinate (X) MC moves using molecular dynamics at a fixed lambda
 * value (i.e. using OpenMM), followed by MC lambda moves.
 *
 * <p>1.) At a fixed Lambda, run a defined length MD trajectory to "move" coordinates and dU/dL on
 * an approximate potential U* (i.e. no OST Bias).
 *
 * <p>2.) Accept / Reject the MD move with probability exp[-Beta(dU - dU*)] where dU is the change
 * in AMOEBA + Bias energy and dU* is the change in AMOEBA + Kinetic energy from the MD.
 *
 * <p>3.) Randomly change the value of Lambda.
 *
 * <p>4.) Accept / Reject the Lambda move using the AMOEBA + OST Bias energy.
 *
 * <p>5.) Add to the time dependent 2D bias using the current values of Lambda and dU/dL.
 *
 * @author Michael J. Schnieders
 * @author Hernan Beranbe
 * @author Mallory R. Tollefson
 * @author Jacob Litman
 * @since 1.0
 */
public class MonteCarloOST extends BoltzmannMC {

  /**
   * Logger object to print out information for this class.
   */
  private static final Logger logger = Logger.getLogger(MonteCarloOST.class.getName());

  /**
   * Controls the effect of verbose by logging at FINE vs. INFO.
   */
  private final Level mcLogLevel;
  /**
   * How verbose MD should be.
   */
  private final MDVerbosity mdVerbosity;

  /**
   * The MD moves are only valid if energy is conserved. For this reason, energy drift is monitored.
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
   * Deposit a bias once every N MC cycles. Defaults to 1.
   */
  private final int biasDepositionFrequency;
  /**
   * Total number of steps to take for MC-OST sampling.
   */
  private long totalSteps;
  /**
   * Number of MD steps to take per MC-OST round.
   */
  private final long stepsPerMove;
  /**
   * Lambda move object for completing MC-OST lambda moves.
   */
  private final LambdaMove lambdaMove;
  /**
   * Double that keeps track of our lambda value.
   */
  private double lambda = 1.0;
  /**
   * Boolean that tells algorithm that we are in the equilibration phase of MC-OST.
   */
  private boolean equilibration = false;
  /**
   * True if MC-OST should handle writing out files.
   */
  private boolean automaticWriteouts = true;

  /**
   * Constructor for MonteCarloOST.
   *
   * @param potentialEnergy          a {@link ffx.numerics.Potential} object.
   * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
   * @param molecularAssembly        a {@link ffx.potential.MolecularAssembly} object.
   * @param properties               a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *                                 object.
   * @param listener                 a {@link ffx.algorithms.AlgorithmListener} object.
   * @param dynamics                 CLI object containing key information.
   * @param verbose                  Whether to be verbose.
   * @param cycleLength              Length of an MC cycle in MD steps.
   */
  public MonteCarloOST(Potential potentialEnergy, OrthogonalSpaceTempering orthogonalSpaceTempering,
                       MolecularAssembly molecularAssembly, CompositeConfiguration properties,
                       AlgorithmListener listener, DynamicsOptions dynamics,
                       boolean verbose, int cycleLength, File dynRestartFile) {
    this.potential = potentialEnergy;
    this.orthogonalSpaceTempering = orthogonalSpaceTempering;
    mcLogLevel = verbose ? Level.INFO : Level.FINE;
    mdVerbosity = verbose ? MDVerbosity.QUIET : MDVerbosity.SILENT;
    stepsPerMove = cycleLength;
    totalSteps = dynamics.getNumSteps();

    ThermostatEnum thermostat = dynamics.thermostat;
    if (!thermostat.equals(ThermostatEnum.ADIABATIC)) {
      logger.info(format(" MC-OST MD moves will use the Adiabatic thermostat (%s ignored).", thermostat));
      dynamics.setThermostat(ThermostatEnum.ADIABATIC);
    }

    IntegratorEnum integrator = dynamics.integrator;
    if (!integrator.reversible || !integrator.deterministic) {
      logger.info(format(" MC-OST MD moves require a reversible, deterministic integrator (Verlet replaced %s).", integrator));
      dynamics.setIntegrator(IntegratorEnum.VERLET);
    }

    // Create the MC MD and Lambda moves.
    boolean useOST = properties.getBoolean("mdmove-full", false);
    Potential mdPotential = useOST ? orthogonalSpaceTempering : potential;
    mdMove = new MDMove(molecularAssembly, mdPotential, listener, dynamics, stepsPerMove, dynRestartFile);
    if (properties.containsKey("randomseed")) {
      int randomSeed = properties.getInt("randomseed", 0);
      logger.info(format(" Setting random seed for lambdaMove to %d ", randomSeed));
      lambdaMove = new LambdaMove(randomSeed, orthogonalSpaceTempering);
      setRandomSeed(randomSeed);
    } else {
      lambdaMove = new LambdaMove(orthogonalSpaceTempering);
    }

    // Configure discrete Lambda.
    boolean discreteLambda = orthogonalSpaceTempering.getHistogram().hd.discreteLambda;
    if (discreteLambda) {
      lambdaMove.setContinuous(false);
      // Set the move size to the Lambda bin width.
      lambdaMove.setMoveSize(orthogonalSpaceTempering.getHistogram().hd.lambdaBinWidth);
    }

    // The OST class will handle adding the time dependent bias.
    orthogonalSpaceTempering.setPropagateLambda(false);
    int frequency = properties.getInt("mc-ost-bias-frequency", 1);
    if (frequency <= 1) {
      biasDepositionFrequency = 1;
    } else {
      biasDepositionFrequency = frequency;
      logger.info(format(" MC-OST will deposit a bias only once per %d MC cycles.", biasDepositionFrequency));
    }
  }

  /**
   * Returns the current value of lambda
   *
   * @return lambda
   */
  public double getLambda() {
    return lambda;
  }

  /**
   * Calls on the OST method set lambda to update lambda to the current value in this class
   *
   * @param lambda a double.
   */
  public void setLambda(double lambda) {
    this.lambda = lambda;
    orthogonalSpaceTempering.setLambda(lambda);
  }

  public MolecularDynamics getMD() {
    return mdMove.getMD();
  }

  /**
   * The goal is to sample lambda and coordinates (X) simultaneously to converge the ensemble average
   * dU/dL for every state (lambda) along the thermodynamic path. Note that the order of 1 & 2 below
   * can be swapped (i.e. run MD and then change lambda). Here the order is random for each trial.
   * <p>
   * 1.) Randomly change the value of Lambda.
   * <p>
   * 2.) At the proposed lambda, run a am MD trajectory to "move" coordinates and dU/dL.
   * <p>
   * 3.) Accept / Reject the Lambda + MD move using the total Hamiltonian (Kinetic energy + OST energy).
   * <p>
   * 4.) Add to the bias.
   */
  public void sampleOneStep() {
    // Validate the starting value of lambda.
    lambda = orthogonalSpaceTempering.getLambda();
    lambda = lambdaMove.validateLambda(lambda);
    orthogonalSpaceTempering.setLambda(lambda);

    int n = potential.getNumberOfVariables();
    double[] gradient = new double[n];
    double[] currentCoordinates = new double[n];
    double[] proposedCoordinates = new double[n];
    long numMoves = totalSteps / stepsPerMove;
    int acceptMCOST = 0;

    // Initialize the current coordinates.
    potential.getCoordinates(currentCoordinates);

    // The Histogram maintains the time-dependent bias.
    Histogram histogram = orthogonalSpaceTempering.getHistogram();

    // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
    double currentOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);

    // Retrieve the computed dU/dL, Force Field Energy and Bias Energy.
    double currentdUdL = orthogonalSpaceTempering.getForceFielddEdL();
    double currentForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
    double currentBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

    // Loop over the number of requested MC Moves.
    for (int imove = 0; imove < numMoves; imove++) {
      long totalMoveTime = -nanoTime();

      // Change lambda before or after the MD trajectory.
      boolean lambdaBeforeMD = random.nextBoolean();

      if (logger.isLoggable(Level.FINE)) {
        if (equilibration) {
          logger.fine(format("\n Equilibration Round %d", imove + 1));
        } else {
          String moveType = lambdaBeforeMD ? "Single-step lambda plus MD move." : "Single-step MD plus lambda move.";
          logger.fine(format("\n MC-OST Round %d: %s", imove + 1, moveType));
        }
        logger.fine(format(" Starting force field energy for move %16.8f", currentForceFieldEnergy));
      }

      // Set the current and proposed lambda.
      double currentLambda = orthogonalSpaceTempering.getLambda();
      double proposedLambda;
      if (equilibration) {
        // The proposed lambda equals the current lambda because its frozen for equilibration.
        proposedLambda = currentLambda;
        singleStepMD();
      } else {
        if (lambdaBeforeMD) {
          // Change lambda before MD.
          proposedLambda = singleStepLambda();
          singleStepMD();
        } else {
          // Change lambda after MD.
          singleStepMD();
          proposedLambda = singleStepLambda();
        }
      }

      // Get the starting and final kinetic energy for the MD move.
      double currentKineticEnergy = mdMove.getInitialKinetic();
      double proposedKineticEnergy = mdMove.getKineticEnergy();

      // Get the new coordinates.
      potential.getCoordinates(proposedCoordinates);

      // Compute the Total OST energy and gradient as the sum of the force field energy and bias energy.
      long proposedOSTEnergyTime = -nanoTime();
      double proposedOSTEnergy;
      try {
        proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(proposedCoordinates, gradient);
      } catch (EnergyException e) {
        mdMove.revertMove();
        if (!equilibration) {
          lambdaMove.revertMove();
          lambda = currentLambda;
        }
        logger.log(Level.INFO, " Unstable MD Move skipped.");
        continue;
      }
      proposedOSTEnergyTime += nanoTime();

      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Time to complete MD OST energy method call %6.3f", proposedOSTEnergyTime * NS2SEC));
      }

      // Retrieve the proposed dU/dL, force field energy and bias energy.
      double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
      double proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
      double proposedBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

      // The Metropolis criteria is based on the sum of the OST Energy and Kinetic Energy.
      double currentTotalEnergy = currentOSTEnergy + currentKineticEnergy;
      double proposedTotalEnergy = proposedOSTEnergy + proposedKineticEnergy;

      // Log energy differences for the proposed move.
      if (logger.isLoggable(mcLogLevel)) {
        logger.log(mcLogLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
        logger.log(mcLogLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f", currentKineticEnergy,
                currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
        logger.log(mcLogLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f", proposedKineticEnergy,
                proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
        logger.log(mcLogLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
            proposedKineticEnergy - currentKineticEnergy,
            proposedForceFieldEnergy - currentForceFieldEnergy, proposedBiasEnergy - currentBiasEnergy,
            proposedTotalEnergy - currentTotalEnergy));
      }

      // Monitor energy conservation.
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

      // Initialize the result string to reflect rejection.
      String result = "    ";
      // Energy change.
      double dE = proposedTotalEnergy - currentTotalEnergy;
      if (equilibration) {
        if (orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL)
            && evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
          // Accept Equilibration MD move.
          result = " -> ";
          acceptMCOST++;
          currentOSTEnergy = proposedOSTEnergy;
          currentdUdL = proposeddUdL;
          currentForceFieldEnergy = proposedForceFieldEnergy;
          currentBiasEnergy = proposedBiasEnergy;
          arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
        } else {
          mdMove.revertMove();
        }
      } else {
        if (orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL)
            && evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
          // Accept MD = Lambda move.
          result = " -> ";
          acceptMCOST++;
          lambda = proposedLambda;
          currentdUdL = proposeddUdL;
          currentForceFieldEnergy = proposedForceFieldEnergy;
          arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
        } else {
          lambdaMove.revertMove();
          lambda = currentLambda;
          mdMove.revertMove();
        }
        if (imove % biasDepositionFrequency == 0) {
          histogram.addBias(currentdUdL);
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Added Bias at move %d: [ L=%5.3f, FL=%9.3f] ", imove + 1, lambda, currentdUdL));
          }
        } else {
          logger.info(format(" Cycle %d: skipping bias deposition.", imove + 1));
        }

        // Compute the updated OST bias.
        currentBiasEnergy = histogram.computeBiasEnergy(lambda, currentdUdL);

        // Update the current OST Energy to be the sum of the current force field energy and updated OST Bias.
        currentOSTEnergy = currentForceFieldEnergy + currentBiasEnergy;

        boolean snapShot = lambda >= orthogonalSpaceTempering.getLambdaWriteOut();
        long mdMoveNum = (imove + 1) * stepsPerMove;
        orthogonalSpaceTempering.getHistogram().ld.stepsTaken += stepsPerMove;
        EnumSet<MDWriteAction> written = mdMove.writeFilesForStep(mdMoveNum, snapShot, true);
        if (written.contains(MDWriteAction.RESTART)) {
          orthogonalSpaceTempering.writeAdditionalRestartInfo(false);
        }
      }

      // Log the move.
      totalMoveTime += nanoTime();
      double percent = (acceptMCOST * 100.0) / (imove + 1);
      logger.info(format(" %4d [ L=%5.3f, dU/dL=%8.3f]%s[ L=%5.3f, dU/dL=%8.3f] ΔE=%12.4f (%5.1f%%) in %6.3f sec.",
          imove + 1, currentLambda, currentdUdL, result,
          proposedLambda, proposeddUdL, dE, percent, totalMoveTime * NS2SEC));
    }
  }

  /**
   * The goal is to sample lambda and coordinates (X) separately to converge the ensemble average
   * dU/dL for every state (lambda) along the thermodynamic path.
   * <p>
   * 1.) At a fixed lambda, run a defined length MD trajectory to "move" coordinates and dU/dL.
   * <p>
   * 2.) Accept / Reject the MD move using the total Hamiltonian (Kinetic energy + OST energy).
   * <p>
   * 3.) Randomly change the value of Lambda.
   * <p>
   * 4.) Accept / Reject the Lambda move using the OST energy.
   * <p>
   * 5.) Add a hill to the histogram.
   */
  public void sampleTwoStep() {
    // Validate the starting value of lambda.
    lambda = orthogonalSpaceTempering.getLambda();
    lambda = lambdaMove.validateLambda(lambda);
    orthogonalSpaceTempering.setLambda(lambda);

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

      if (logger.isLoggable(Level.FINE)) {
        if (equilibration) {
          logger.fine(format("\n MD Equilibration Round %d", imove + 1));
        } else {
          logger.fine(format("\n MC Orthogonal Space Sampling Round %d: Independent Steps", imove + 1));
        }
      }

      // Run MD in an approximate potential U* (U star) that does not include the OST bias.
      long mdMoveTime = -nanoTime();
      mdMove.move(mdVerbosity);
      mdMoveTime += nanoTime();
      logger.log(mcLogLevel, format("  Total time for MD move: %6.3f", mdMoveTime * NS2SEC));

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
        logger.log(Level.INFO, " Unstable MD Move skipped.");
        continue;
      }
      proposedOSTEnergyTime += nanoTime();

      logger.fine(format("  Time to complete MD OST energy method call %6.3f", proposedOSTEnergyTime * NS2SEC));

      // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
      double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
      double proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
      double proposedBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

      // The Metropolis criteria is based on the sum of the OST Energy and Kinetic Energy.
      double currentTotalEnergy = currentOSTEnergy + currentKineticEnergy;
      double proposedTotalEnergy = proposedOSTEnergy + proposedKineticEnergy;

      if (logger.isLoggable(mcLogLevel)) {
        logger.log(mcLogLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
        logger.log(mcLogLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f", currentKineticEnergy,
            currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
        logger.log(mcLogLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f", proposedKineticEnergy,
            proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
        logger.log(mcLogLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
            proposedKineticEnergy - currentKineticEnergy,
            proposedForceFieldEnergy - currentForceFieldEnergy, proposedBiasEnergy - currentBiasEnergy,
            proposedTotalEnergy - currentTotalEnergy));
      }

      double energyChange = mdMove.getEnergyChange();
      if (abs(energyChange) > ENERGY_CONSERVATION_TOLERANCE) {
        mdMove.revertMove();
        logger.warning(" MC Move skipped due to lack of MD energy conservation");
        continue;
      }

      String result = "    ";
      double dE = proposedTotalEnergy - currentTotalEnergy;
      if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(),
          proposeddUdL) && evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
        // Accept MD move.
        acceptMD++;
        result = " -> ";
        currentOSTEnergy = proposedOSTEnergy;
        currentdUdL = proposeddUdL;
        currentForceFieldEnergy = proposedForceFieldEnergy;
        currentBiasEnergy = proposedBiasEnergy;
        arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
      } else {
        mdMove.revertMove();
      }

      // Log the move.
      mdMoveAndEvalTime += nanoTime();
      double percent = (acceptMD * 100.0) / (imove + 1);
      logger.info(format(" %4d MD [ L=%5.3f, dU/dL=%8.3f]%s[ L=%5.3f, dU/dL=%8.3f] ΔE=%12.4f (%5.1f%%) in %6.3f sec.", imove + 1,
          lambda, currentdUdL, result, lambda, proposeddUdL, dE, percent, mdMoveAndEvalTime * NS2SEC));

      // During equilibration, do not change Lambda or contribute to the OST bias.
      if (!equilibration) {
        // Update Lambda.
        long lambdaMoveTime = -nanoTime();
        double currentLambda = orthogonalSpaceTempering.getLambda();
        lambdaMove.move();
        double proposedLambda = orthogonalSpaceTempering.getLambda();

        // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
        long proposedOSTEnergyTime2 = -nanoTime();
        proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);
        proposedOSTEnergyTime2 += nanoTime();

        // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
        proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
        proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();

        if (logger.isLoggable(mcLogLevel)) {
          logger.log(mcLogLevel, format("\n  Time to complete Lambda OST energy method call %6.3f ", proposedOSTEnergyTime2 * NS2SEC));
          logger.log(mcLogLevel, format("  Current  OST     %12.3f at L=%5.3f.", currentOSTEnergy, currentLambda));
          logger.log(mcLogLevel, format("  Proposed OST     %12.3f at L=%5.3f.", proposedOSTEnergy, proposedLambda));
          logger.log(mcLogLevel, format("  MC Energy change: %12.3f (kcal/mol).", proposedOSTEnergy - currentOSTEnergy));
        }

        dE = proposedOSTEnergy - currentOSTEnergy;
        if (orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL)
            && evaluateMove(currentOSTEnergy, proposedOSTEnergy)) {
          acceptLambda++;
          result = " -> ";
          currentForceFieldEnergy = proposedForceFieldEnergy;
          currentdUdL = proposeddUdL;
          lambda = proposedLambda;
        } else {
          result = "    ";
          lambdaMove.revertMove();
          lambda = currentLambda;
        }

        lambdaMoveTime += nanoTime();
        percent = (acceptLambda * 100.0) / (imove + 1);
        logger.info(format("      L  [ L=%5.3f, dU/dL=%8.3f]%s[ L=%5.3f, dU/dL=%8.3f] ΔE=%12.4f (%5.1f%%) in %6.3f sec.",
            currentLambda, currentdUdL, result, proposedLambda, proposeddUdL, dE, percent, lambdaMoveTime * NS2SEC));

        if (imove % biasDepositionFrequency == 0) {
          histogram.addBias(currentdUdL);
          logger.fine(format(" Added Bias at move %d: [ L=%5.3f, FL=%9.3f] ", imove + 1, lambda, currentdUdL));
        } else {
          logger.info(format(" Cycle %d: skipping bias deposition.", imove + 1));
        }

        // Compute the updated OST bias.
        currentBiasEnergy = histogram.computeBiasEnergy(lambda, currentdUdL);

        // Update the current OST Energy to be the sum of the current Force Field Energy and updated
        // OST Bias.
        currentOSTEnergy = currentForceFieldEnergy + currentBiasEnergy;

        if (automaticWriteouts) {
          long mdMoveNum = (imove + 1) * stepsPerMove;
          boolean trySnapshot = lambda >= orthogonalSpaceTempering.getLambdaWriteOut();
          EnumSet<MDWriteAction> written = mdMove.writeFilesForStep(mdMoveNum, trySnapshot, true);
          if (written.contains(MDWriteAction.RESTART)) {
            orthogonalSpaceTempering.writeAdditionalRestartInfo(false);
          }
        }
      }

      if (logger.isLoggable(Level.FINE)) {
        totalMoveTime += nanoTime();
        logger.fine(format(" Round complete in %6.3f sec.", totalMoveTime * NS2SEC));
      }
    }
  }

  public void setAutomaticWriteouts(boolean autoWrite) {
    automaticWriteouts = autoWrite;
  }

  /**
   * Sets the value of the boolean equilibration variables to true or false to either allow an
   * equilibration step or skip it.
   *
   * @param equilibration a boolean.
   */
  public void setEquilibration(boolean equilibration) {
    this.equilibration = equilibration;
  }

  /**
   * Calls on LambdaMove class method setLambdaStdDev to update the lambda standard deviation to the
   * current value in this class
   *
   * @param stdDev a double.
   */
  public void setLambdaStdDev(double stdDev) {
    if (!lambdaMove.isContinuous()) {
      if (stdDev != orthogonalSpaceTempering.getHistogram().hd.lambdaBinWidth) {
        logger.info(
            format(" Requested Lambda step size change %6.4f is not equal to OST bin width %6.4f.",
                stdDev, orthogonalSpaceTempering.getHistogram().hd.lambdaBinWidth));
        return;
      }
    }
    lambdaMove.setMoveSize(stdDev);
  }

  /**
   * Sets the next number of steps MC-OST should run.
   *
   * @param totalSteps The length of the next MC-OST run.
   */
  public void setTotalSteps(long totalSteps) {
    this.totalSteps = totalSteps;
  }

  /**
   * Propose a lambda move.
   *
   * @return The proposed lambda.
   */
  private double singleStepLambda() {
    lambdaMove.move();
    double proposedLambda = orthogonalSpaceTempering.getLambda();
    logger.log(mcLogLevel, format(" Proposed lambda: %5.3f.", proposedLambda));
    return proposedLambda;
  }

  /**
   * Run MD in an approximate potential U* (U star) that does not include the OST bias.
   */
  private void singleStepMD() {
    long mdMoveTime = -nanoTime();
    mdMove.move(mdVerbosity);
    mdMoveTime += nanoTime();
    logger.log(mcLogLevel, format(" Total time for MD move: %6.3f", mdMoveTime * NS2SEC));
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void revertStep() {
    throw new UnsupportedOperationException("Not supported yet.");
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
}
