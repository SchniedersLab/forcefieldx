/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms;

import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.LambdaMove;
import ffx.algorithms.mc.MDMove;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import static ffx.algorithms.MolecularDynamics.NS2SEC;

/**
 * Sample a thermodynamic path using the OSRW method, with the time-dependent
 * bias built up using Metropolis Monte Carlo steps.
 * <p>
 * The algorithm generates coordinate (X) MC moves using molecular dynamics at a
 * fixed lambda value (i.e. using OpenMM), followed by MC lambda moves.
 * <p>
 * 1.) At a fixed Lambda, run a defined length MD trajectory to "move"
 * coordinates and dU/dL on an approximate potential U* (i.e. no OSRW Bias).
 * <p>
 * 2.) Accept / Reject the MD move with probability exp[-Beta(dU - dU*)] where
 * dU is the change in AMOEBA + Bias energy and dU* is the change in AMOEBA +
 * Kinetic energy from the MD.
 * <p>
 * 3.) Randomly change the value of Lambda.
 * <p>
 * 4.) Accept / Reject the Lambda move using the AMOEBA + OSRW Bias energy.
 * <p>
 * 5.) Add to the time dependent 2D bias using the current values of Lambda and
 * dU/dL.
 *
 * @author Michael J. Schnieders
 * @author Hernan Beranbe
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class MonteCarloOSRW extends BoltzmannMC {

    /**
     * Logger object to print out information for this class.
     */
    private static final Logger logger = Logger.getLogger(MonteCarloOSRW.class.getName());
    /**
     * Potential object used to retrieve the coordinates for the system.
     */
    private final Potential potential;
    /**
     * OSRW object used to retrieve OSRW energy throughout the simulation.
     */
    private final AbstractOSRW osrw;
    /**
     * Double that keeps track of our lambda value.
     */
    private double lambda = 0.0;
    /**
     * Keeps track of the percentage of lambda moves accepted.
     */
    private int acceptLambda = 0;
    /**
     * Keeps track of the percentage of MD moves accepted.
     */
    private int acceptMD = 0;
    /**
     * Keeps track of the percentage of MCOSRW moves accepted.
     */
    private int acceptMCOSRW = 0;
    /**
     * MDMove object for completing MC-OSRW molecular dynamics moves.
     */
    private MDMove mdMove;
    /**
     * Total number of steps to take for MC-OSRW sampling.
     */
    private int totalSteps = 10000000;
    /**
     * Number of steps to take per MC-OSRW round.
     */
    private int stepsPerMove = 50;
    /**
     * Total time required to complete one MC-OSRW round.
     */
    private long totalMoveTime = 0;
    /**
     * Total time required to complete a lambda move.
     */
    private long lambdaMoveTime = 0;
    /**
     * Total time required to complete and MD move.
     */
    private long mdMoveTime = 0;
    /**
     * Total time required to run and evaluate MD move.
     */
    private long mdMoveAndEvalTime = 0;
    /**
     * Lambda move object for completing MC-OSRW lambda moves.
     */
    private LambdaMove lambdaMove;
    /**
     * Current lambda at the start of the MCOSRW move.
     */
    private double currentLambda;
    /**
     * Lambda proposed during MCOSRW move.
     */
    private double proposedLambda;
    /**
     * Boolean that tells algorithm that we are in the equilibration phase of
     * MC-OSRW.
     */
    private boolean equilibration = false;

    /**
     * <p>
     * Constructor for MonteCarloOSRW.</p>
     *
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param osrw a {@link ffx.algorithms.AbstractOSRW} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param properties a
     * {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     * {@link ffx.algorithms.thermostats.ThermostatEnum} object.
     * @param requestedIntegrator a
     * {@link ffx.algorithms.integrators.IntegratorEnum} object.
     */
    public MonteCarloOSRW(Potential potentialEnergy, AbstractOSRW osrw,
            MolecularAssembly molecularAssembly, CompositeConfiguration properties,
            AlgorithmListener listener, ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator) {
        this.potential = potentialEnergy;
        this.osrw = osrw;

        /**
         * Create the MC MD and Lambda moves.
         */
        mdMove = new MDMove(molecularAssembly, potential, properties, listener, requestedThermostat, requestedIntegrator);
        if (properties.containsKey("randomseed")) {
            int randomSeed = properties.getInt("randomseed", 0);
            logger.info(String.format(" Setting random seed for lambdaMove to %d ", randomSeed));
            lambdaMove = new LambdaMove(randomSeed, lambda, osrw);
            setRandomSeed(randomSeed);
        } else {
            lambdaMove = new LambdaMove(lambda, osrw);
        }

        /**
         * Changing the value of lambda will be handled by this class, as well
         * as adding the time dependent bias.
         */
        osrw.setPropagateLambda(false);

    }

    /**
     * Takes in parameters and calls the MDMove method setMDParameters to update
     * the stepsPerMove and timeStep parameters to the current value in this
     * class
     *
     * @param totalSteps a int.
     * @param stepsPerMove a int.
     * @param timeStep a double.
     */
    public void setMDMoveParameters(int totalSteps, int stepsPerMove, double timeStep) {
        this.totalSteps = totalSteps;
        this.stepsPerMove = stepsPerMove;
        mdMove.setMDParameters(stepsPerMove, timeStep);
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
     * Calls on the OSRW method set lambda to update lambda to the current value
     * in this class
     *
     * @param lambda a double.
     */
    public void setLambda(double lambda) {
        this.lambda = lambda;
        osrw.setLambda(lambda);
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
     * The goal is to sample coordinates (X) and converge "dU/dL" for every
     * state (lambda) along the thermodynamic path.
     * <p>
     * 1.) At a fixed lambda, run a defined length MD trajectory to "move"
     * coordinates and dU/dL.
     * <p>
     * 2.) Accept / Reject the MD move using the OSRW energy.
     * <p>
     * 3.) Randomly change the value of Lambda.
     * <p>
     * 4.) Accept / Reject the Lambda move using the OSRW energy.
     * <p>
     * 5.) Add to the bias.
     */
    public void sample() {

        int n = potential.getNumberOfVariables();
        double[] gradient = new double[n];
        double[] currentCoordinates = new double[n];
        double[] proposedCoordinates = new double[n];
        int numMoves = totalSteps / stepsPerMove;
        acceptLambda = 0;
        acceptMD = 0;

        // Initialize the current Force Field Energy, Bias Energy, Total OSRW Energy and dU/dL.
        potential.getCoordinates(currentCoordinates);
        double currentOSRWEnergy = osrw.energyAndGradient(currentCoordinates, gradient);
        double currentdUdL = osrw.getForceFielddEdL();
        double currentForceFieldEnergy = osrw.getForceFieldEnergy();
        double currentBiasEnergy = osrw.getBiasEnergy();

        /**
         * Initialize MC move instances.
         */
        for (int imove = 0; imove < numMoves; imove++) {
            totalMoveTime = -System.nanoTime();
            mdMoveAndEvalTime = -System.nanoTime();

            if (equilibration) {
                logger.info(String.format("\n MD Equilibration Round %d", imove + 1));
            } else {
                logger.info(String.format("\n MCOSRW Round %d", imove + 1));
            }

            logger.info(String.format(" Starting force field energy for move %f", currentForceFieldEnergy));

            /**
             * Run MD in an approximate potential U* (U star) that does not
             * include the OSRW bias.
             */
            mdMoveTime = -System.nanoTime();
            mdMove.move();
            mdMoveTime += System.nanoTime();
            logger.info(String.format(" Total time for MD move: %6.3f", mdMoveTime * NS2SEC));

            // Get the starting and final kinetic energy for the MD move.
            double currentKineticEnergy = mdMove.getStartingKineticEnergy();
            double proposedKineticEnergy = mdMove.getKineticEnergy();

            // Get the new coordinates.
            potential.getCoordinates(proposedCoordinates);

            // ToDo: Get rid of redundant call to OpenMM in "osrw.energy" to compute the force field energy we already have from MD.
            // ToDo: Instead, refactor OSRW to compute only the Bias energy (i.e. not re-compute the force field energy).
            long proposedOSRWEnergyTime = 0;
            proposedOSRWEnergyTime = -System.nanoTime();
            double proposedOSRWEnergy = osrw.energyAndGradient(proposedCoordinates, gradient);
            proposedOSRWEnergyTime += System.nanoTime();
            logger.info(String.format(" Time to complete MD OSRW energy method call %6.3f", proposedOSRWEnergyTime * NS2SEC));
            double proposeddUdL = osrw.getForceFielddEdL();
            double proposedForceFieldEnergy = osrw.getForceFieldEnergy();
            double proposedBiasEnergy = osrw.getBiasEnergy();

            double currentTotalEnergy = currentOSRWEnergy + currentKineticEnergy;
            double proposedTotalEnergy = proposedOSRWEnergy + proposedKineticEnergy;

            logger.info(format("  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logger.info(format("  Current  %12.3f %12.3f %12.3f %12.3f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logger.info(format("  Proposed %12.3f %12.3f %12.3f %12.3f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logger.info(format("  Delta    %12.3f %12.3f %12.3f %12.3f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            if (evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                /**
                 * Accept MD move.
                 */
                acceptMD++;
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(String.format(" MCMD step   :      Accepted [FL=%8.3f,E=%12.4f] -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSRWEnergy, proposeddUdL, proposedOSRWEnergy, percent));
                currentTotalEnergy = proposedTotalEnergy;
                currentOSRWEnergy = proposedOSRWEnergy;
                currentdUdL = proposeddUdL;
                currentForceFieldEnergy = proposedForceFieldEnergy;
                currentBiasEnergy = proposedBiasEnergy;
                System.arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
            } else {
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(String.format(" MCMD step   :      Rejected [FL=%8.3f,E=%12.4f] -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSRWEnergy, proposeddUdL, proposedOSRWEnergy, percent));
                long revertTime = 0;
                revertTime = -System.nanoTime();
                mdMove.revertMove();
                revertTime += System.nanoTime();
                logger.info(String.format(" Reverted move in %6.3f", revertTime * NS2SEC));
            }
            mdMoveAndEvalTime += System.nanoTime();
            logger.info(String.format("\n Total time to run and evaluate MD move: %6.3f", mdMoveAndEvalTime * NS2SEC));

            /**
             * During equilibration, do not change Lambda or contribute to the
             * OSRW bias.
             */
            if (!equilibration) {
                /**
                 * Update Lambda.
                 */
                lambdaMoveTime = -System.nanoTime();
                long proposeLambdaTime = 0;
                proposeLambdaTime = -System.nanoTime();
                currentLambda = osrw.getLambda();
                lambdaMove.move();
                proposedLambda = osrw.getLambda();
                proposeLambdaTime += System.nanoTime();
                logger.info(String.format(" Time spent generating and assigning lambda %6.3f", proposeLambdaTime * NS2SEC));

                long proposedOSRWEnergyTime2 = 0;
                proposedOSRWEnergyTime2 = -System.nanoTime();
                proposedOSRWEnergy = osrw.energyAndGradient(currentCoordinates, gradient);
                proposedOSRWEnergyTime2 += System.nanoTime();
                logger.info(String.format(" Time to complete Lambda OSRW energy method call %6.3f ", proposedOSRWEnergyTime2 * NS2SEC));
                proposedForceFieldEnergy = osrw.getForceFieldEnergy();
                proposedBiasEnergy = osrw.getBiasEnergy();
                proposeddUdL = osrw.getForceFielddEdL();

                logger.info(format("\n Current  OSRW     %12.3f at L=%5.3f.", currentOSRWEnergy, currentLambda));
                logger.info(format(" Proposed OSRW     %12.3f at L=%5.3f.", proposedOSRWEnergy, proposedLambda));
                logger.info(format(" MC Energy change: %12.3f (kcal/mol).", proposedOSRWEnergy - currentOSRWEnergy));

                if (evaluateMove(currentOSRWEnergy, proposedOSRWEnergy)) {
                    acceptLambda++;
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logger.info(String.format(" MC Lambda step   : Accepted [ L=%8.3f,E=%12.4f] -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSRWEnergy, proposedLambda, proposedOSRWEnergy, percent));
                    currentOSRWEnergy = proposedOSRWEnergy;
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    currentBiasEnergy = proposedBiasEnergy;
                    currentdUdL = proposeddUdL;
                    lambda = proposedLambda;
                } else {
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logger.info(String.format(" MC Lambda step   : Rejected [ L=%8.3f,E=%12.4f] -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSRWEnergy, proposedLambda, proposedOSRWEnergy, percent));
                    long revertLambdaTime = 0;
                    revertLambdaTime = -System.nanoTime();
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                    revertLambdaTime += System.nanoTime();
                    logger.info(String.format(" Time to revert lambda move: %6.3f", revertLambdaTime * NS2SEC));
                }
                lambdaMoveTime += System.nanoTime();
                logger.info(String.format(" Lambda move completed in %6.3f", lambdaMoveTime * NS2SEC));

                /**
                 * Update time dependent bias.
                 */
                osrw.addBias(currentdUdL, currentCoordinates, null);
                logger.info(format(" Added Bias at [L=%5.3f, FL=%9.3f]", lambda, currentdUdL));
            }

            totalMoveTime += System.nanoTime();
            logger.info(format(" Total MC-OSRW Round Time: %6.3f sec.", totalMoveTime * NS2SEC));
        }
    }

    public void simplifiedSample() {

        int n = potential.getNumberOfVariables();
        double[] gradient = new double[n];
        double[] currentCoordinates = new double[n];
        double[] proposedCoordinates = new double[n];
        int numMoves = totalSteps / stepsPerMove;
        acceptMD = 0;

        // Initialize the current Force Field Energy, Bias Energy, Total OSRW Energy and dU/dL.
        potential.getCoordinates(currentCoordinates);
        double currentOSRWEnergy = osrw.energyAndGradient(currentCoordinates, gradient);
        double currentdUdL = osrw.getForceFielddEdL();
        double currentForceFieldEnergy = osrw.getForceFieldEnergy();
        double currentBiasEnergy = osrw.getBiasEnergy();

        /**
         * Initialize MC move instances.
         */
        for (int imove = 0; imove < numMoves; imove++) {
            totalMoveTime = -System.nanoTime();

            if (equilibration) {
                logger.info(String.format("\n MD Equilibration Round %d", imove + 1));
            } else {
                logger.info(String.format("\n MCOSRW Round %d", imove + 1));
                lambdaMoveTime = -System.nanoTime();
                long proposeLambdaTime = 0;
                proposeLambdaTime = -System.nanoTime();
                currentLambda = osrw.getLambda();
                lambdaMove.move();
                proposedLambda = osrw.getLambda();
                proposeLambdaTime += System.nanoTime();
                logger.info(String.format(" Time spent generating and assigning lambda %6.3f", proposeLambdaTime * NS2SEC));
            }

            logger.info(String.format(" Starting force field energy for move %f", currentForceFieldEnergy));

            /**
             * Run MD in an approximate potential U* (U star) that does not
             * include the OSRW bias.
             */
            mdMoveTime = -System.nanoTime();
            mdMove.move();
            mdMoveTime += System.nanoTime();
            logger.info(String.format(" Total time for MD move: %6.3f", mdMoveTime * NS2SEC));

            // Get the starting and final kinetic energy for the MD move.
            double currentKineticEnergy = mdMove.getStartingKineticEnergy();
            double proposedKineticEnergy = mdMove.getKineticEnergy();

            // Get the new coordinates.
            potential.getCoordinates(proposedCoordinates);

            // ToDo: Get rid of redundant call to OpenMM in "osrw.energy" to compute the force field energy we already have from MD.
            // ToDo: Instead, refactor OSRW to compute only the Bias energy (i.e. not re-compute the force field energy).
            long proposedOSRWEnergyTime = 0;
            proposedOSRWEnergyTime = -System.nanoTime();
            double proposedOSRWEnergy = osrw.energyAndGradient(proposedCoordinates, gradient);
            proposedOSRWEnergyTime += System.nanoTime();
            logger.info(String.format(" Time to complete MD OSRW energy method call %6.3f", proposedOSRWEnergyTime * NS2SEC));
            double proposeddUdL = osrw.getForceFielddEdL();
            double proposedForceFieldEnergy = osrw.getForceFieldEnergy();
            double proposedBiasEnergy = osrw.getBiasEnergy();
            
            logger.info(String.format(" Proposed dUdL: %f", proposeddUdL));

            double currentTotalEnergy = currentOSRWEnergy + currentKineticEnergy;
            double proposedTotalEnergy = proposedOSRWEnergy + proposedKineticEnergy;

            logger.info(format("  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logger.info(format("  Current  %12.3f %12.3f %12.3f %12.3f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logger.info(format("  Proposed %12.3f %12.3f %12.3f %12.3f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logger.info(format("  Delta    %12.3f %12.3f %12.3f %12.3f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            if (equilibration) {
                if (evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                    /**
                     * Accept MD.
                     */
                    acceptMD++;
                    double percent = (acceptMD * 100.0) / (imove + 1);
                    logger.info(String.format(" MCMD step   :      Accepted [FL=%8.3f,E=%12.4f] -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentdUdL, currentOSRWEnergy, proposeddUdL, proposedOSRWEnergy, percent));
                    currentOSRWEnergy = proposedOSRWEnergy;
                    currentdUdL = proposeddUdL;
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    currentBiasEnergy = proposedBiasEnergy;
                    System.arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
                } else {
                    double percent = (acceptMD * 100.0) / (imove + 1);
                    logger.info(String.format(" MCMD step   :      Rejected [FL=%8.3f,E=%12.4f] -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentdUdL, currentOSRWEnergy, proposeddUdL, proposedOSRWEnergy, percent));
                    //long revertTime = 0;
                    //revertTime = -System.nanoTime();
                    mdMove.revertMove();
                    //revertTime += System.nanoTime();
                    //logger.info(String.format(" Reverted move in %6.3f", revertTime * NS2SEC));
                }
            } else {
                if (evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                    acceptMCOSRW++;
                    double percent = (acceptMCOSRW * 100.0) / (imove + 1);
                    logger.info(String.format(" MCOSRW step :      Accepted [ L=%8.3f,FL=%8.3f,E=%12.4f] -> [ L=%8.3f,FL=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentdUdL, currentOSRWEnergy, proposedLambda, proposeddUdL, proposedOSRWEnergy, percent));
                    lambda = proposedLambda;
                    currentdUdL = proposeddUdL;
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    System.arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
                } else {
                    double percent = (acceptMCOSRW * 100.0) / (imove + 1);
                    logger.info(String.format(" MCOSRW step :      Rejected [ L=%8.3f,FL=%8.3f,E=%12.4f] -> [ L=%8.3f,FL=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentdUdL, currentOSRWEnergy, proposedLambda, proposeddUdL, proposedOSRWEnergy, percent));
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                    //long revertTime = 0;
                    //revertTime = -System.nanoTime();
                    mdMove.revertMove();
                    //revertTime += System.nanoTime();
                    //logger.info(String.format(" Reverted move in %6.3f", revertTime * NS2SEC));                
                }
                /**
                 * Update time dependent bias.
                 */
                osrw.addBias(currentdUdL, currentCoordinates, null);
                logger.info(format(" Added Bias at [L=%5.3f, FL=%9.3f]", lambda, currentdUdL));

                currentBiasEnergy = osrw.computeBiasEnergy(currentLambda, currentdUdL);
                currentOSRWEnergy = currentForceFieldEnergy + currentBiasEnergy;

            }

            totalMoveTime += System.nanoTime();
            logger.info(format(" Total MC-OSRW Round Time: %6.3f sec.", totalMoveTime * NS2SEC));
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
