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
package ffx.algorithms.optimize;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;

/**
 * Run NVT molecular dynamics at a series of temperatures to optimize a structure.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class SimulatedAnnealing implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(SimulatedAnnealing.class.getName());
    /**
     * The MolecularDynamics instance used for Simulated Annealing.
     */
    private final MolecularDynamics molecularDynamics;
    /**
     * High temperature annealing limit.
     */
    private double highTemperature;
    /**
     * Low temperature annealing limit.
     */
    private double lowTemperature;
    /**
     * Number of annealing windows.
     */
    private int annealingSteps;
    /**
     * Flag to indicate a temperature ladder was provided.
     */
    private boolean targetTemperaturesPresent = false;
    /**
     * Array of target simulated annealing target temperatures.
     */
    private double[] targetTemperatures;
    /**
     * Number of MD steps per annealing window.
     */
    private int mdSteps;
    /**
     * Integration time step.
     */
    private double timeStep;
    /**
     * Interval to print updates to the screen.
     */
    private double printInterval = 0.01;
    /**
     * Flag to indicate the algorithm is done.
     */
    private boolean done;
    /**
     * Flag to indicate the UI has requested th algorithm terminate.
     */
    private boolean terminate;

    /**
     * <p>
     * Constructor for SimulatedAnnealing.</p>
     *
     * @param molecularAssembly The Molecular Assembly to operate on.
     * @param potentialEnergy   The potential to anneal against.
     * @param properties        The system properties to use.
     * @param listener          The algorithm listener is a callback to UI.
     */
    public SimulatedAnnealing(MolecularAssembly molecularAssembly,
                              Potential potentialEnergy,
                              CompositeConfiguration properties,
                              AlgorithmListener listener) {

        this(molecularAssembly, potentialEnergy, properties, listener,
                ThermostatEnum.BERENDSEN, IntegratorEnum.VERLET);
    }

    /**
     * <p>
     * Constructor for SimulatedAnnealing.</p>
     *
     * @param assembly            The Molecular Assembly to operate on.
     * @param potentialEnergy     The potential to anneal against.
     * @param properties          The system properties to use.
     * @param listener            The algorithm listener is a callback to UI.
     * @param requestedThermostat The requested thermostat.
     * @param requestedIntegrator The requested integrator.
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
     * @param highTemperature High temperature annealing limit.
     * @param lowTemperature  Low temperature annealing limit.
     * @param annealingSteps  The number of annealing steps.
     * @param mdSteps         The number of MD steps.
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
     * @param highTemperature High temperature annealing limit.
     * @param lowTemperature  Low temperature annealing limit.
     * @param annealingSteps  The number of annealing steps.
     * @param mdSteps         The number of MD steps.
     * @param timeStep        The MD time step (fsec).
     */
    public void anneal(double highTemperature,
                       double lowTemperature,
                       int annealingSteps,
                       int mdSteps,
                       double timeStep) {

        // Return if already running; Could happen if two threads call anneal
        // on the same SimulatedAnnealing instance.
        if (!done) {
            logger.warning(" Programming error - a thread invoked anneal when it was already running.");
            return;
        }
        done = false;
        logger.info(" Simulated annealing starting up");

        setAnnealingSteps(annealingSteps);
        setHighTemperature(highTemperature);
        setLowTemperature(lowTemperature);
        setMdSteps(mdSteps);
        setTimeStep(timeStep);
        begin();
    }

    /**
     * <p>annealToTargetValues.</p>
     *
     * @param targetTemperatures The array of annealing temperatures.
     * @param mdSteps            The number of MD steps.
     * @param timeStep           The MD time step (fsec).
     */
    public void annealToTargetValues(double[] targetTemperatures,
                                     int mdSteps,
                                     double timeStep) {
        if (targetTemperatures == null) {
            anneal(400, 100, 10, mdSteps, timeStep);
            return;
        }
        targetTemperaturesPresent = true;
        this.targetTemperatures = targetTemperatures;
        this.annealingSteps = targetTemperatures.length;
        this.highTemperature = targetTemperatures[0];
        this.lowTemperature = targetTemperatures[annealingSteps - 1];
        setMdSteps(mdSteps);
        setTimeStep(timeStep);
        begin();
    }

    private void begin() {
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
                molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 10.0,
                        temperature, true, null);
                if (terminate) {
                    logger.info(String.format("\n Terminating at temperature %8.3f.\n", temperature));
                    break;
                }
            }
        } else {
            for (double temperature : targetTemperatures) {
                molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 10.0,
                        temperature, true, null);
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

    /**
     * Set the low temperature annealing limit.
     *
     * @param highTemperature The low temperature limit (K).
     */
    public void setHighTemperature(double highTemperature) {
        if (highTemperature < 0) {
            highTemperature = 400.0;
        }
        this.highTemperature = highTemperature;
    }

    /**
     * Set the low temperature annealing limit.
     *
     * @param lowTemperature The low temperature limit (K).
     */
    public void setLowTemperature(double lowTemperature) {
        if (annealingSteps == 1) {
            lowTemperature = highTemperature;
        } else if (lowTemperature < 0.0 || lowTemperature > highTemperature) {
            lowTemperature = 100;
        }
        this.lowTemperature = lowTemperature;
    }

    /**
     * Set the number of annealing steps.
     *
     * @param annealingSteps The number of annealing steps.
     */
    public void setAnnealingSteps(int annealingSteps) {
        if (annealingSteps <= 0) {
            annealingSteps = 1;
        }
        this.annealingSteps = annealingSteps;
    }

    /**
     * Set the MD time step.
     *
     * @param timeStep The time step (fsec).
     */
    public void setTimeStep(double timeStep) {
        if (timeStep <= 0) {
            timeStep = 1.0;
        }
        this.timeStep = timeStep;
    }

    /**
     * Set the number of MD steps.
     *
     * @param mdSteps THe number of MD steps.
     */
    public void setMdSteps(int mdSteps) {
        if (mdSteps <= 0) {
            mdSteps = 100;
        }
        this.mdSteps = mdSteps;
    }

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
