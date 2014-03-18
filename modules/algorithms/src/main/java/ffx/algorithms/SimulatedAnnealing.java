/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.algorithms;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;

/**
 * Run NVT molecular dynamics at a series of temperatures.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class SimulatedAnnealing implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(SimulatedAnnealing.class.getName());
    private final MolecularDynamics molecularDynamics;
    private double highTemperature;
    private double lowTemperature;
    private int annealingSteps;
    private int mdSteps;
    private double timeStep;
    private boolean done, terminate;

    /**
     * <p>Constructor for SimulatedAnnealing.</p>
     *
     * @param assembly a {@link ffx.potential.bonded.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     */
    public SimulatedAnnealing(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener) {

        this(assembly, potentialEnergy, properties, listener,
                Thermostats.BERENDSEN, Integrators.BEEMAN);
    }

    /**
     * <p>Constructor for SimulatedAnnealing.</p>
     *
     * @param assembly a {@link ffx.potential.bonded.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     * {@link ffx.algorithms.Thermostat.Thermostats}
     * @param requestedIntegrator a
     * {@link ffx.algorithms.Integrator.Integrators}
     */
    public SimulatedAnnealing(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener,
            Thermostats requestedThermostat,
            Integrators requestedIntegrator) {

        molecularDynamics = new MolecularDynamics(assembly,
                potentialEnergy, properties,
                listener,
                requestedThermostat,
                requestedIntegrator);
        done = true;
    }

    /**
     * <p>anneal</p>
     *
     * @param highTemperature a double.
     * @param lowTemperature a double.
     * @param annealingSteps a int.
     * @param mdSteps a int.
     */
    public void anneal(double highTemperature,
            double lowTemperature,
            int annealingSteps,
            int mdSteps) {
        anneal(highTemperature, lowTemperature, annealingSteps, mdSteps, 1.0);
    }

    /**
     * <p>anneal</p>
     *
     * @param highTemperature a double.
     * @param lowTemperature a double.
     * @param annealingSteps a int.
     * @param mdSteps a int.
     * @param timeStep a double
     */
    public void anneal(double highTemperature,
            double lowTemperature,
            int annealingSteps,
            int mdSteps,
            double timeStep) {

        /**
         * Return if already running; Could happen if two threads call dynamic
         * on the same SimulatedAnnealing instance.
         */
        if (!done) {
            logger.warning(" Programming error - a thread invoked anneal when it was already running.");
            return;
        }
        done = false;
        logger.info(" Simulated annealing starting up");

        if (annealingSteps <= 0) {
            annealingSteps = 1;
        }
        this.annealingSteps = annealingSteps;
        if (highTemperature < 0) {
            highTemperature = 400.0;
        }
        this.highTemperature = highTemperature;
        if (annealingSteps == 1) {
            lowTemperature = highTemperature;
        } else if (lowTemperature < 0.0 || lowTemperature > highTemperature) {
            lowTemperature = 100;
        }
        this.lowTemperature = lowTemperature;

        if (mdSteps <= 0) {
            mdSteps = 100;
        }
        this.mdSteps = mdSteps;

        if (timeStep <= 0) {
            timeStep = 1.0;
        }
        this.timeStep = timeStep;

        logger.info(String.format(" Initial temperature:    %8.3f (Kelvin)", highTemperature));
        logger.info(String.format(" Final temperature:      %8.3f (Kelvin)", lowTemperature));
        logger.info(String.format(" Annealing steps:        %8d", annealingSteps));
        logger.info(String.format(" MD steps/temperature:   %8d", mdSteps));
        logger.info(String.format(" MD time step:           %8.3f (fs)", timeStep));

        Thread annealingThread = new Thread(this);
        annealingThread.start();
        synchronized (this) {
            try {
                while (annealingThread.isAlive()) {
                    wait(100);
                }
            } catch (Exception e) {
                String message = "Simualted annealing interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }

    }

    /**
     * {@inheritDoc}
     *
     * This method should only be invoked within the SimulatedAnnealing
     * instance.
     */
    @Override
    public void run() {
        done = false;
        terminate = false;

        double dt = (highTemperature - lowTemperature) / (annealingSteps - 1);
        for (int i = 0; i < annealingSteps; i++) {
            double temperature = highTemperature - dt * i;
            molecularDynamics.dynamic(mdSteps, timeStep, 0.001, 0.002, temperature, true, null);
            if (terminate) {
                logger.info(String.format("\n Terminating at temperature %8.3f.\n", temperature));
                break;
            }
        }

        if (!terminate) {
            logger.info(String.format(" Completed %8d annealing steps\n", annealingSteps));
        }

        done = true;
        terminate = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating annealing.\n", e);
                }
            }
        }
    }
}
