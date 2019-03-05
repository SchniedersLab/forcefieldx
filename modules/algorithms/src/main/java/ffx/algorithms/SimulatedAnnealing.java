/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;

/**
 * Run NVT molecular dynamics at a series of temperatures.
 *
 * @author Michael J. Schnieders
 * @since 1.0
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
    private boolean targetTemperaturesPresent = false;
    private double[] targetTemperatures;
    private double printInterval = 0.01;

    /**
     * <p>
     * Constructor for SimulatedAnnealing.</p>
     *
     * @param assembly        a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties      a
     *                        {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener        a {@link ffx.algorithms.AlgorithmListener} object.
     */
    public SimulatedAnnealing(MolecularAssembly assembly,
                              Potential potentialEnergy,
                              CompositeConfiguration properties,
                              AlgorithmListener listener) {

        this(assembly, potentialEnergy, properties, listener,
                ThermostatEnum.BERENDSEN, IntegratorEnum.VERLET);
    }

    /**
     * <p>
     * Constructor for SimulatedAnnealing.</p>
     *
     * @param assembly            a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param properties          a
     *                            {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     *                            {@link ffx.algorithms.thermostats.ThermostatEnum}
     * @param requestedIntegrator a
     *                            {@link ffx.algorithms.integrators.IntegratorEnum}
     */
    public SimulatedAnnealing(MolecularAssembly assembly,
                              Potential potentialEnergy,
                              CompositeConfiguration properties,
                              AlgorithmListener listener,
                              ThermostatEnum requestedThermostat,
                              IntegratorEnum requestedIntegrator) {

        molecularDynamics = new MolecularDynamics(assembly,
                potentialEnergy, properties,
                listener,
                requestedThermostat,
                requestedIntegrator);
        done = true;
    }

    /**
     * <p>
     * anneal</p>
     *
     * @param highTemperature a double.
     * @param lowTemperature  a double.
     * @param annealingSteps  a int.
     * @param mdSteps         a int.
     */
    public void anneal(double highTemperature,
                       double lowTemperature,
                       int annealingSteps,
                       int mdSteps) {
        anneal(highTemperature, lowTemperature, annealingSteps, mdSteps, 1.0);
    }

    /**
     * <p>
     * anneal</p>
     *
     * @param highTemperature a double.
     * @param lowTemperature  a double.
     * @param annealingSteps  a int.
     * @param mdSteps         a int.
     * @param timeStep        a double
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
     * <p>annealToTargetValues.</p>
     *
     * @param targetTemperatures an array of {@link double} objects.
     * @param mdSteps            a int.
     * @param timeStep           a double.
     */
    public void annealToTargetValues(double[] targetTemperatures,
                                     int mdSteps,
                                     double timeStep) {

        targetTemperaturesPresent = true;
        this.targetTemperatures = targetTemperatures;

        if (mdSteps <= 0) {
            mdSteps = 100;
        }
        this.mdSteps = mdSteps;

        if (timeStep <= 0) {
            timeStep = 1.0;
        }
        this.timeStep = timeStep;

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
     * <p>
     * This method should only be invoked within the SimulatedAnnealing
     * instance.
     */
    @Override
    public void run() {
        done = false;
        terminate = false;

        if (!targetTemperaturesPresent) {
            double dt = (highTemperature - lowTemperature) / (annealingSteps - 1);
            for (int i = 0; i < annealingSteps; i++) {
                double temperature = highTemperature - dt * i;
                molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 10.0, temperature, true, null);
                if (terminate) {
                    logger.info(String.format("\n Terminating at temperature %8.3f.\n", temperature));
                    break;
                }
            }
        } else {
            for (int i = 0; i < targetTemperatures.length; i++) {
                double temperature = targetTemperatures[i];
                molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 10.0, temperature, true, null);
                if (terminate) {
                    logger.info(String.format("\n Terminating at temperature %8.3f.\n", temperature));
                    break;
                }
            }
        }
        if (!terminate) {
            logger.info(String.format(" Completed %8d annealing steps\n", annealingSteps));
        }

        done = true;
        terminate = false;
    }

    //Set print interval to report thermodyanamics (psec)

    /**
     * <p>Setter for the field <code>printInterval</code>.</p>
     *
     * @param printInterval a double.
     */
    public void setPrintInterval(double printInterval) {
        this.printInterval = printInterval;
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

    /**
     * <p>getKineticEnergy.</p>
     *
     * @return a double.
     */
    public double getKineticEnergy() {
        return molecularDynamics.getKineticEnergy();
    }

    /**
     * <p>getPotentialEnergy.</p>
     *
     * @return a double.
     */
    public double getPotentialEnergy() {
        return molecularDynamics.getPotentialEnergy();
    }

    /**
     * <p>getTotalEnergy.</p>
     *
     * @return a double.
     */
    public double getTotalEnergy() {
        return molecularDynamics.getTotalEnergy();
    }

    /**
     * <p>getTemperature.</p>
     *
     * @return a double.
     */
    public double getTemperature() {
        return molecularDynamics.getTemperature();
    }
}
