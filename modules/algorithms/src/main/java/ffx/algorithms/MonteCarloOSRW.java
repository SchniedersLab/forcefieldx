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
 * 1.) At a fixed lambda, run a defined length MD trajectory to "move"
 * coordinates and dU/dL.
 *
 * 2.) Accept / Reject the MD move using the OSRW energy.
 *
 * 3.) Randomly change the value of Lambda.
 *
 * 4.) Accept / Reject the Lambda move using the OSRW energy.
 *
 * @author Mallory R. Tollefson
 */
public class MonteCarloOSRW extends BoltzmannMC {

    private static final Logger logger = Logger.getLogger(MonteCarloOSRW.class.getName());

    private final Potential potential;
    private final AbstractOSRW osrw;
    private final MolecularAssembly molecularAssembly;
    private final CompositeConfiguration properties;
    private final AlgorithmListener listener;
    private final Thermostats requestedThermostat;
    private final Integrators requestedIntegrator;
    //private LambdaInterface linter;
    private double lambda = 0.0;
    private int totalMDSteps = 10000000;
    private final int mdSteps = 100;
    private int lambdaAccept = 0;
    private int dUdLAccept = 0;

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
        this.molecularAssembly = molecularAssembly;
        this.properties = properties;
        this.listener = listener;
        this.requestedThermostat = requestedThermostat;
        this.requestedIntegrator = requestedIntegrator;

        /**
         * Changing the value of lambda will be handled by this class, as well
         * as adding the time dependent bias.
         */
        osrw.setPropagateLambda(false);

    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
        osrw.setLambda(lambda);
    }

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
        double[] gradient = new double[n];
        int numMoves = totalMDSteps / mdSteps;

        /**
         * Initialize MC move instances.
         */
        MDMove mdMove = new MDMove(molecularAssembly, potential, properties, listener, requestedThermostat, requestedIntegrator);
        LambdaMove lambdaMove = new LambdaMove(lambda, osrw);

        for (int imove = 0; imove < numMoves; imove++) {

            potential.getCoordinates(coordinates);
            double currentEnergy = osrw.energyAndGradient(coordinates, gradient);
            double currentdUdL = osrw.getForceFielddEdL();
            double currentFFEnergy = osrw.getForceFieldEnergy();
            currentEnergy += mdMove.getKineticEnergy();

            /**
             * Run MD.
             */
            mdMove.move();
            potential.getCoordinates(coordinates);
            double proposedEnergy = osrw.energyAndGradient(coordinates, gradient);
            double proposeddUdL = osrw.getForceFielddEdL();
            double proposedFFEnergy = osrw.getForceFieldEnergy();
            proposedEnergy += mdMove.getKineticEnergy();

            logger.info(String.format(" Current %16.8f and Proposed %16.8f force field energies.", currentFFEnergy, proposedFFEnergy));

            if (evaluateMove(currentEnergy, proposedEnergy)) {
                dUdLAccept++;
                double percent = (dUdLAccept * 100.0) / (imove + 1);
                logger.info(String.format(" X,dU/dL MC step accepted: e(%8.3f)=%12.6f -> e(%8.3f)=%12.6f (%5.1f)",
                        currentdUdL, currentEnergy, proposeddUdL, proposedEnergy, percent));
                currentEnergy = proposedEnergy;
                //Accept MD move using OSRW energy
            } else {
                double percent = (dUdLAccept * 100.0) / (imove + 1);
                logger.info(String.format(" X,dU/dL MC step rejected: e(%8.3f)=%12.6f -> e(%8.3f)=%12.6f (%5.1f)",
                        currentdUdL, currentEnergy, proposeddUdL, proposedEnergy, percent));
                mdMove.revertMove();
            }

            /**
             * Update Lambda.
             */
            potential.getCoordinates(coordinates);
            currentEnergy = osrw.energyAndGradient(coordinates, gradient);
            currentdUdL = osrw.getForceFielddEdL();
            double currentLambda = osrw.getLambda();
            lambdaMove.move();
            proposedEnergy = osrw.energyAndGradient(coordinates, gradient);
            proposeddUdL = osrw.getForceFielddEdL();
            double proposedLambda = osrw.getLambda();

            if (evaluateMove(currentEnergy, proposedEnergy)) {
                lambdaAccept++;
                double percent = (lambdaAccept * 100.0) / (imove + 1);
                logger.info(String.format(" Lambda MC step accepted: e(%8.3f)=%12.6f -> e(%8.3f)=%12.6f (%5.1f)",
                        currentLambda, currentEnergy, proposedLambda, proposedEnergy, percent));
                currentdUdL = proposeddUdL;
            } else {
                double percent = (lambdaAccept * 100.0) / (imove + 1);
                logger.info(String.format(" Lambda MC step rejected: e(%8.3f)=%12.6f -> e(%8.3f)=%12.6f (%5.1f)",
                        currentLambda, currentEnergy, proposedLambda, proposedEnergy, percent));
                lambdaMove.revertMove();
            }
            lambda = osrw.getLambda();

            /**
             * Update time dependent bias.
             */
            double freeEnergy = osrw.updateFLambda(true);
            osrw.addBias(currentdUdL, freeEnergy);

            logger.info(String.format(" MC free energy %10.6f", freeEnergy));
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
