/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.LambdaMove;
import ffx.algorithms.mc.MDMove;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;

/**
 * Sample a thermodynamic path using the OSRW method, with the time-dependent
 * bias built up using Metropolis Monte Carlo steps.
 *
 * The goal is generate coordinate (X) MC moves using molecular dynamics at a
 * fixed lambda value, following by MC lambda moves.
 *
 * 1.) At a fixed Lambda, run a defined length MD trajectory to "move"
 * coordinates and dU/dL.
 *
 * 2.) Accept / Reject the MD move using the OSRW energy.
 *
 * 3.) Randomly change the value of Lambda.
 *
 * 4.) Accept / Reject the Lambda move using the OSRW energy.
 *
 * 5.) Add a biasing potential to the current value of Lambda.
 *
 * @author Michael J. Schnieders
 *
 * @author Hernan Beranbe
 *
 * @author Mallory R. Tollefson
 *
 * @since 1.0
 */
public class MonteCarloOSRW extends BoltzmannMC {

    private static final Logger logger = Logger.getLogger(MonteCarloOSRW.class.getName());

    private final Potential potential;
    private final AbstractOSRW osrw;
    private double lambda = 0.0;

    private int acceptLambda = 0;
    private int acceptMD = 0;

    private MDMove mdMove;
    private int totalSteps = 10000000;
    private int stepsPerMove = 50;

    private LambdaMove lambdaMove;

    private boolean equilibration = false;

    /**
     * @param potentialEnergy
     * @param osrw
     * @param molecularAssembly
     * @param properties
     * @param listener
     * @param requestedThermostat
     * @param requestedIntegrator
     */
    public MonteCarloOSRW(Potential potentialEnergy, AbstractOSRW osrw,
            MolecularAssembly molecularAssembly, CompositeConfiguration properties,
            AlgorithmListener listener, Thermostats requestedThermostat, Integrators requestedIntegrator) {
        this.potential = potentialEnergy;
        this.osrw = osrw;

        /**
         * Create the MC MD and Lambda moves.
         */
        mdMove = new MDMove(molecularAssembly, potential, properties, listener, requestedThermostat, requestedIntegrator);
        lambdaMove = new LambdaMove(lambda, osrw);

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
     * @param totalSteps
     * @param stepsPerMove
     * @param timeStep
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
     * @param stdDev
     */
    public void setLambdaStdDev(double stdDev) {
        lambdaMove.setStdDev(stdDev);
    }

    /**
     * Sets the value of the boolean equilibration variables to true or false to
     * either allow an equilibration step or skip it.
     *
     * @param equilibration
     */
    public void setEquilibration(boolean equilibration) {
        this.equilibration = equilibration;
    }

    /**
     * Calls on the OSRW method set lambda to update lambda to the current value
     * in this class
     *
     * @param lambda
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
     *
     * 1.) At a fixed lambda, run a defined length MD trajectory to "move"
     * coordinates and dU/dL.
     *
     * 2.) Accept / Reject the MD move using the OSRW energy.
     *
     * 3.) Randomly change the value of Lambda.
     *
     * 4.) Accept / Reject the Lambda move using the OSRW energy.
     *
     * 5.) Add to the bias.
     */
    public void sample() {
        int n = potential.getNumberOfVariables();
        double[] coordinates = new double[n];
        double[] proposedCoordinates = new double[n];
        double[] newCoordinates;
        double[] gradient = new double[n];
        int numMoves = totalSteps / stepsPerMove;
        acceptLambda = 0;
        acceptMD = 0;

        /**
         * Initialize MC move instances.
         */
        for (int imove = 0; imove < numMoves; imove++) {

            potential.getCoordinates(coordinates);
            double currentEnergy = osrw.energyAndGradient(coordinates, gradient);
            double currentdUdL = osrw.getForceFielddEdL();
            currentEnergy += mdMove.getKineticEnergy();

            /**
             * Run MD.
             */
            mdMove.move();
            potential.getCoordinates(proposedCoordinates);
            double proposedEnergy = osrw.energyAndGradient(proposedCoordinates, gradient);
            double proposeddUdL = osrw.getForceFielddEdL();
            proposedEnergy += mdMove.getKineticEnergy();

            if (evaluateMove(currentEnergy, proposedEnergy)) {
                /**
                 * Accept MD move.
                 */
                acceptMD++;
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(String.format(" MC MD step %d:     Accepted E(%8.3f)=%12.6f -> E(%8.3f)=%12.6f (%5.1f)",
                        imove + 1, currentdUdL, currentEnergy, proposeddUdL, proposedEnergy, percent));
                currentEnergy = proposedEnergy;
                newCoordinates = proposedCoordinates;

            } else {
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(String.format(" MC MD step %d:     Rejected E(%8.3f)=%12.6f -> E(%8.3f)=%12.6f (%5.1f)",
                        imove + 1, currentdUdL, currentEnergy, proposeddUdL, proposedEnergy, percent));
                mdMove.revertMove();
                newCoordinates = coordinates;
            }

            /**
             * During equilibration, do not change Lambda or contribute to the
             * OSRW bias.
             */
            if (equilibration) {
                continue;
            }

            /**
             * Update Lambda.
             */
            currentEnergy = osrw.energyAndGradient(newCoordinates, gradient);
            currentdUdL = osrw.getForceFielddEdL();
            double currentLambda = osrw.getLambda();
            lambdaMove.move();
            proposedEnergy = osrw.energyAndGradient(newCoordinates, gradient);
            proposeddUdL = osrw.getForceFielddEdL();
            double proposedLambda = osrw.getLambda();

            if (evaluateMove(currentEnergy, proposedEnergy)) {
                acceptLambda++;
                double percent = (acceptLambda * 100.0) / (imove + 1);
                logger.info(String.format(" MC Lambda step %d: Accepted E(%8.3f)=%12.6f -> E(%8.3f)=%12.6f (%5.1f)",
                        imove + 1, currentLambda, currentEnergy, proposedLambda, proposedEnergy, percent));
                currentdUdL = proposeddUdL;
            } else {
                double percent = (acceptLambda * 100.0) / (imove + 1);
                logger.info(String.format(" MC Lambda step %d: Rejected E(%8.3f)=%12.6f -> E(%8.3f)=%12.6f (%5.1f)",
                        imove + 1, currentLambda, currentEnergy, proposedLambda, proposedEnergy, percent));
                lambdaMove.revertMove();
                potential.getCoordinates(coordinates);
            }
            lambda = osrw.getLambda();

            /**
             * Update time dependent bias.
             */
            logger.info(format(" Adding bias at L=%6.4f and dU/dL=%16.8f.", lambda, currentdUdL));
            osrw.addBias(currentdUdL, coordinates, gradient);
        }
    }

    @Override
    protected double currentEnergy() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected void storeState() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void revertStep() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
