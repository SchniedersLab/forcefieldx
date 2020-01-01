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
package ffx.algorithms.optimize.anneal;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
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
 * @author Jacob M. Litman
 * @since 1.0
 */
public class SimulatedAnnealing implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(SimulatedAnnealing.class.getName());
    /**
     * The MolecularDynamics instance used for Simulated Annealing.
     */
    private final MolecularDynamics molecularDynamics;
    /**
     * Schedule for annealing.
     */
    private final AnnealingSchedule schedule;
    /**
     * Number of MD steps per annealing window.
     */
    private final int mdSteps;
    /**
     * Integration time step.
     */
    private final double timeStep;
    /**
     * Whether to reinitialize velocities at the start of each timestep.
     */
    private final boolean reinitV;
    /**
     * Minimum length of a window in psec.
     */
    private final double minSimLength;
    /**
     * Interval to print updates to the screen.
     */
    private double printInterval = 0.01;
    /**
     * Flag to indicate the algorithm is done.
     */
    private boolean done = true;
    /**
     * Flag to indicate the UI has requested the algorithm terminate.
     */
    private boolean terminate;
    /**
     * Number of MD steps per OpenMM cycle (assuming OpenMM is used!).
     */
    private int trajSteps = 1;
    private double saveFrequency = 0.1;
    /**
     * Restart file.
     */
    private File dynFile = null;

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
     * @param schedule            Schedule of temperatures to simulate at.
     * @param mdSteps             Steps per SA window.
     * @param timeStep            Timestep for MD in psec.
     */
    public SimulatedAnnealing(MolecularAssembly assembly,
                              Potential potentialEnergy,
                              CompositeConfiguration properties,
                              AlgorithmListener listener,
                              ThermostatEnum requestedThermostat,
                              IntegratorEnum requestedIntegrator,
                              AnnealingSchedule schedule,
                              int mdSteps,
                              double timeStep,
                              boolean reinitVelocities,
                              File dynFile) {

        molecularDynamics = MolecularDynamics.dynamicsFactory(assembly, potentialEnergy,
                properties, listener, requestedThermostat, requestedIntegrator);
        this.schedule = schedule;
        this.mdSteps = mdSteps;
        this.timeStep = timeStep;
        this.reinitV = reinitVelocities;
        minSimLength = mdSteps * schedule.minWindowLength() * timeStep;
        this.dynFile = dynFile;
    }

    /**
     * <p>
     * anneal</p>
     */
    public void anneal() {
        // Return if already running; Could happen if two threads call anneal
        // on the same SimulatedAnnealing instance.
        if (!done) {
            logger.warning(" Programming error - a thread invoked anneal when it was already running.");
            return;
        }
        done = false;
        logger.info(" Beginning simulated annealing");
        begin();
    }

    private void begin() {
        logger.info(String.format(" Initial temperature:    %8.3f (Kelvin)", schedule.getHighTemp()));
        logger.info(String.format(" Final temperature:      %8.3f (Kelvin)", schedule.getLowTemp()));
        logger.info(String.format(" Annealing steps:        %8d", schedule.getNumWindows()));
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
                String message = "Simulated annealing interrupted.";
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

        int minMdSteps = (int) (mdSteps * schedule.minWindowLength());
        if (minMdSteps < trajSteps) {
            logger.warning(String.format(" Minimum number of MD steps per annealing cycle %d was less than steps per OpenMM MD cycle %d! Setting steps per MD cycle to %d", minMdSteps, trajSteps, minMdSteps));
            setTrajectorySteps(minMdSteps);
        }

        int nWindows = schedule.getNumWindows();
        boolean forceFirstReinit = (dynFile == null);

        for (int i = 0; i < nWindows; i++) {
            double temperature = schedule.getTemperature(i);
            int nSteps = (int) (schedule.windowLength(i) * mdSteps);
            logger.info(String.format(" Annealing window %d: %d steps at %9.4g K", (i+1), nSteps, temperature));
            molecularDynamics.dynamic(nSteps, timeStep, printInterval, saveFrequency, temperature, (reinitV || forceFirstReinit), dynFile);
            if (dynFile == null) {
                dynFile = molecularDynamics.getDynFile();
            }
            forceFirstReinit = false;
            if (terminate) {
                logger.info(String.format("\n Terminating at temperature %8.3f.\n", temperature));
                break;
            }
        }
        if (!terminate) {
            logger.info(String.format(" Completed %8d annealing steps\n", nWindows));
        }

        done = true;
        terminate = false;
    }

    /**
     * Sets the number of steps to use per OpenMM cycle.
     * @param trajectorySteps Steps per OpenMM cycle.
     */
    public void setTrajectorySteps(int trajectorySteps) {
        molecularDynamics.setIntervalSteps(trajectorySteps);
    }

    /**
     * Sets the frequency of writing to the trajectory file.
     *
     * @param save Frequency (psec^-1) to write out the trajectory.
     */
    public void setSaveFrequency(double save) {
        this.saveFrequency = save;
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
     * Method to set the Restart Frequency.
     *
     * @param restart the time between writing restart files.
     * @throws java.lang.IllegalArgumentException If restart frequency is not a
     *                                            positive number
     */
    public void setRestartFrequency(double restart) throws IllegalArgumentException {
        if (Double.isFinite(restart) && restart > 0) {
            molecularDynamics.setRestartFrequency(restart);
        } else {
            throw new IllegalArgumentException(String.format(" Restart frequency must be positive finite, was %10.4g", restart));
        }
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

    /**
     * Functional interface corresponding to constructors of non-composite AnnealingSchedules.
     */
    private interface ScheduleConstructor {
        AnnealingSchedule asConstruct(int nWindows, double tLow, double tHigh);
    }

    /**
     * Represents non-composite AnnealingSchedules known (i.e. not FlatEndAnnealSchedule).
     */
    public enum Schedules {
        EXP(ExpAnnealSchedule::new, "EXP", "EXPONENTIAL"), LINEAR(LinearAnnealSchedule::new, "LINEAR");

        private final ScheduleConstructor sc;
        private final Set<String> aliases;

        Schedules(ScheduleConstructor sc, String... names) {
            this.sc = sc;
            aliases = Collections.unmodifiableSet(new HashSet<>(Arrays.asList(names)));
        }

        /**
         * Creates an AnnealingSchedule corresponding to this enum and provided values.
         * @param nWindows Number of annealing windows.
         * @param tLow     Final temperature.
         * @param tHigh    Starting temperature.
         * @return         An AnnealingSchedule.
         */
        public AnnealingSchedule generate(int nWindows, double tLow, double tHigh) {
            return sc.asConstruct(nWindows, tLow, tHigh);
        }

        /**
         * Attempt to parse a String to a Schedules in a case-insensitive, alias-recognizing fashion.
         *
         * @param name Name of a schedule.
         * @return     A Schedules enum.
         */
        public static Schedules parse(String name) {
            name = name.toUpperCase();
            for (Schedules s : values()) {
                if (s.aliases.contains(name)) {
                    return s;
                }
            }
            return valueOf(name);
        }
    }
}
