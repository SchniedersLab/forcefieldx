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
package ffx.algorithms.dynamics;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import static java.lang.String.format;

import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;

import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getVelocities;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.extended.ExtendedSystem;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.NS2SEC;
import static ffx.utilities.Constants.kB;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Hernan V. Bernabe
 */
public class MolecularDynamicsOpenMM extends MolecularDynamics {

    private static final Logger logger = Logger.getLogger(MolecularDynamicsOpenMM.class.getName());

    /**
     * OpenMM ForceFieldEnergy.
     */
    private ForceFieldEnergyOpenMM forceFieldEnergyOpenMM;
    /**
     * OpenMM Context.
     */
    private PointerByReference context;
    /**
     * OpenMM Integrator.
     */
    private PointerByReference integrator;
    /**
     * Number of OpenMM Particles (typically the number of FFX atoms).
     */
    private int numParticles;
    /**
     * Number of OpenMM Degrees of Freedom.
     */
    private int dof;
    /**
     * Integrator Type.
     */
    private final IntegratorEnum integratorType;
    /**
     * Thermostat Type.
     */
    private final ThermostatEnum thermostatType;
    /**
     * Integrator String.
     */
    private String integratorString;
    /**
     * Number of OpenMM MD steps per iteration.
     */
    private int intervalSteps;
    /**
     * Flag to indicate OpenMM MD interactions are running.
     */
    private boolean running;
    /**
     * Run time.
     */
    private long time;
    /**
     * Obtain all variables with each update (i.e. include velocities, gradients).
     */
    private boolean getAllVars = true;
    /**
     * Method to run on update for obtaining variables. Will either grab
     * everything (default) or energies + positions (MC-OST).
     */
    private Runnable obtainVariables = this::getAllOpenMMVariables;

    /**
     * Constructs an MolecularDynamicsOpenMM object, to perform molecular
     * dynamics using native OpenMM routines, avoiding the cost of communicating
     * coordinates, gradients, and energies back and forth across the PCI bus.
     *
     * @param assembly     MolecularAssembly to operate on
     * @param potential    Either a ForceFieldEnergyOpenMM, or a Barostat.
     * @param properties   Associated properties
     * @param listener     a {@link ffx.algorithms.AlgorithmListener} object.
     * @param thermostat   May have to be slightly modified for native OpenMM routines
     * @param integratorMD May have to be slightly modified for native OpenMM routines
     */
    public MolecularDynamicsOpenMM(MolecularAssembly assembly, Potential potential,
                                   CompositeConfiguration properties, AlgorithmListener listener,
                                   ThermostatEnum thermostat, IntegratorEnum integratorMD) {
        super(assembly, potential, properties, listener, thermostat, integratorMD);

        // Initialization specific to MolecularDynamicsOpenMM
        running = false;
        List<Potential> potentialStack = new ArrayList<>(potential.getUnderlyingPotentials());
        potentialStack.add(potential);

        List<ForceFieldEnergyOpenMM> feOMM = potentialStack.stream().
                filter((Potential p) -> p instanceof ForceFieldEnergyOpenMM).
                map((Potential p) -> (ForceFieldEnergyOpenMM) p).
                collect(Collectors.toList());
        if (feOMM.size() != 1) {
            logger.severe(String.format(" Attempting to create a MolecularDynamicsOpenMM with %d OpenMM force field energies: this presently only allows one!", feOMM.size()));
        }
        forceFieldEnergyOpenMM = feOMM.get(0);

        List<Barostat> barostats = potentialStack.stream().
                filter((Potential p) -> p instanceof Barostat).
                map((Potential p) -> (Barostat) p).
                collect(Collectors.toList());
        if (barostats.isEmpty()) {
            constantPressure = false;
        } else if (barostats.size() > 1) {
            logger.severe(String.format(" Attempting to create a MolecularDynamicsOpenMM with %d barostats: this presently only allows 0-1!", barostats.size()));
        } else {
            barostat = barostats.get(0);
            barostat.setActive(false);
        }

        numParticles = forceFieldEnergyOpenMM.getNumParticles();
        forceFieldEnergyOpenMM.addCOMMRemover(false);
        thermostatType = thermostat;
        integratorType = integratorMD;
        integratorToString(integratorType);

        // Pseudo-random number generator used to seed the OpenMM velocity generator method.
        Random random = new Random();
        if (properties.containsKey("velRandomSeed")) {
            random.setSeed(properties.getInt("velRandomSeed", 0));
        } else {
            random.setSeed(0);
        }

        logger.info(" Molecular Dynamics OpenMM instance created.");
    }

    /**
     * Sets whether or not to obtain all variables (velocities, gradients) from OpenMM,
     * or just positions and energies.
     *
     * @param obtainVA If true, obtain all variables from OpenMM each update.
     */
    @Override
    public void setObtainVelAcc(boolean obtainVA) {
        // TODO: Make this more generic by letting it obtain any weird combination of variables.
        getAllVars = obtainVA;
        obtainVariables = obtainVA ? this::getAllOpenMMVariables : this::getOpenMMEnergiesAndPositions;
    }

    @Override
    protected void appendSnapshot(String[] extraLines) {
        if (!getAllVars) {
            // If !getAllVars, need to ensure coordinates are synced before writing a snapshot.
            getOpenMMEnergiesAndPositions();
        }
        super.appendSnapshot(extraLines);
    }

    @Override
    public void writeRestart() {
        if (!getAllVars) {
            // If !getAllVars, need to ensure all variables are synced before writing the restart.
            getAllOpenMMVariables();
        }
        super.writeRestart();
    }

    /**
     * Get the ForceFieldEnergyOpenMM instance used to run MD.
     *
     * @return a {@link ffx.potential.ForceFieldEnergyOpenMM} object.
     */
    public ForceFieldEnergyOpenMM getForceFieldEnergyOpenMM() {
        return forceFieldEnergyOpenMM;
    }

    /**
     * <p>
     * Setter for the field <code>intervalSteps</code>.</p>
     *
     * @param intervalSteps a int.
     */
    @Override
    public void setIntervalSteps(int intervalSteps) {
        this.intervalSteps = intervalSteps;
        logger.info(String.format(" Interval Steps set at %d", intervalSteps));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void init(long numSteps, double timeStep, double loggingInterval, double trajectoryInterval,
                     String fileType, double restartInterval, double temperature, boolean initVelocities, File dyn) {
        super.init(numSteps, timeStep, loggingInterval, trajectoryInterval,
                fileType, restartInterval, temperature, initVelocities, dyn);

        boolean isLangevin = integratorType.equals(IntegratorEnum.STOCHASTIC);
        switch (thermostatType) {
            case BUSSI:
            case BERENDSEN: {
                if (!isLangevin) {
                    logger.log(basicLogging, String.format(" Replacing FFX thermostat %s with OpenMM Andersen thermostat", thermostatType));
                    forceFieldEnergyOpenMM.addAndersenThermostat(targetTemperature);
                } else {
                    logger.log(basicLogging, " Langevin/Stochastic dynamics already has temperature control, will not be adding thermostat!");
                }
            }
            break;
            // Stochastic integrator gets updated when the context is updated.
            // Adiabatic constant-energy simulations don't need updated temperature.
        }

        if (constantPressure) {
            setMonteCarloBarostat(targetTemperature);
        }

        updateContext();

        dof = forceFieldEnergyOpenMM.calculateDegreesOfFreedom();
        // I'm curious why this was done...
        forceFieldEnergyOpenMM.setLambda(forceFieldEnergyOpenMM.getLambda());
    }

    void preRunOps() {
        // Update the time step in Picoseconds.
        OpenMM_Integrator_setStepSize(integrator, dt);
        super.preRunOps();
        setOpenMMState();
    }

    private void initializeEnergies() {
        // Compute the current potential energy.
        currentPotentialEnergy = forceFieldEnergyOpenMM.energyAndGradient(x, gradient);

        // Initialize current and previous accelerations.
        boolean respa = integratorType.equals(IntegratorEnum.RESPA);
        if (!loadRestart || initialized || respa) {
            // For the Respa integrator, initial accelerations are from the slowly varying forces.
            if (respa) {
                forceFieldEnergyOpenMM.setEnergyTermState(Potential.STATE.SLOW);
                forceFieldEnergyOpenMM.energyAndGradient(x, gradient);
            }

            for (int i = 0; i < numberOfVariables; i++) {
                a[i] = -KCAL_TO_GRAM_ANG2_PER_PS2 * gradient[i] / mass[i];
            }

            if (aPrevious != null) {
                System.arraycopy(a, 0, aPrevious, 0, numberOfVariables);
            }
        }

        updateContext();

        getOpenMMEnergies();
        initialKinetic = currentKineticEnergy;
        initialPotential = currentPotentialEnergy;
        initialTotal = currentTotalEnergy;
        initialTemp = currentTemperature;
    }

    void postInitEnergies() {
        super.postInitEnergies();
        running = true;
    }

    private void mainLoop(long numSteps) {
        long i = 0;
        time = System.nanoTime();

        while (i < numSteps) {

            // Take MD steps in OpenMM.
            long takeStepsTime = -System.nanoTime();
            takeOpenMMSteps(intervalSteps);
            takeStepsTime += System.nanoTime();
            logger.fine(String.format("\n Took steps in %6.3f", takeStepsTime * NS2SEC));
            totalSimTime += intervalSteps * dt;

            // Update the total step count.
            i += intervalSteps;

            long secondUpdateTime = -System.nanoTime();
            updateFromOpenMM(i, running);
            secondUpdateTime += System.nanoTime();

            logger.fine(String.format("\n Update finished in %6.3f", secondUpdateTime * NS2SEC));
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUpdate and takeOpenMMSteps to run the molecular dynamics
     * simulation.
     */
    @Override
    public void dynamic(long numSteps, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {
        // Return if already running;
        // Could happen if two threads call dynamic on the same MolecularDynamics instance.
        if (!done) {
            logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        init(numSteps, timeStep, printInterval, saveInterval, fileType, restartInterval, temperature, initVelocities, dyn);

        if (intervalSteps == 0 || intervalSteps > numSteps) {
            // Safe cast: if intervalSteps > numSteps, then numSteps must be less than Integer.MAX_VALUE.
            intervalSteps = (int) numSteps;
        }

        try {
            preRunOps();
        } catch (IllegalStateException ise) {
            return;
        }
        initializeEnergies();
        postInitEnergies(); // Not over-ridden.
        mainLoop(numSteps);
        postRun(); // Not over-ridden.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTimeStep() {
        return dt;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getIntervalSteps() {
        return intervalSteps;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * <p>
     * updateContext.</p>
     */
    private void updateContext() {
        if (context == null) {
            forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
        } else {
            String currentIntegrator = forceFieldEnergyOpenMM.getIntegratorString();
            double currentTimeStp = forceFieldEnergyOpenMM.getTimeStep();
            double currentTemperature = forceFieldEnergyOpenMM.getTemperature();
            if (currentTemperature != targetTemperature || currentTimeStp != dt
                    || !currentIntegrator.equalsIgnoreCase(integratorString)) {
                forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
            }
        }
        context = forceFieldEnergyOpenMM.getContext();
        integrator = forceFieldEnergyOpenMM.getIntegrator();
    }

    /**
     * takeOpenMMSteps moves the simulation forward in time a user defined
     * number of steps and integrates the equations of motion for each step.
     * This method ensures that the algorithm reports back only when the time
     * interval (steps) specified by the user is completed.
     *
     * @param intervalSteps Number of MD steps to take.
     */
    private void takeOpenMMSteps(int intervalSteps) {
        OpenMM_Integrator_step(integrator, intervalSteps);
    }

    @Override
    public void revertState() throws Exception {
        super.revertState();
        setOpenMMState();
    }

    /**
     * Get OpenMM Energies.
     */
    private void getOpenMMEnergies() {
        context = forceFieldEnergyOpenMM.getContext();

        PointerByReference state = OpenMM_Context_getState(context, OpenMM_State_Energy, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);

        OpenMM_State_destroy(state);
    }

    /**
     * Get OpenMM Energies and Positions.
     */
    private void getOpenMMEnergiesAndPositions() {
        context = forceFieldEnergyOpenMM.getContext();

        int infoMask = OpenMM_State_Positions + OpenMM_State_Energy;
        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);

        PointerByReference positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles * 3, x);
        forceFieldEnergyOpenMM.getPeriodicBoxVectors(state);

        OpenMM_State_destroy(state);
    }

    /**
     * Get OpenMM energies, positions, velocities, and accelerations.
     */
    private void getAllOpenMMVariables() {
        context = forceFieldEnergyOpenMM.getContext();
        int infoMask = OpenMM_State_Positions + OpenMM_State_Energy + OpenMM_State_Velocities + OpenMM_State_Forces;
        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);

        int nVars = numParticles * 3;
        PointerByReference positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, nVars, x);
        forceFieldEnergyOpenMM.getPeriodicBoxVectors(state);

        PointerByReference vels = OpenMM_State_getVelocities(state);
        forceFieldEnergyOpenMM.getOpenMMVelocities(vels, nVars, v);

        PointerByReference forces = OpenMM_State_getForces(state);
        forceFieldEnergyOpenMM.getOpenMMAccelerations(forces, nVars, mass, a);

        OpenMM_State_destroy(state);
    }

    private void setOpenMMState() {
        forceFieldEnergyOpenMM.setOpenMMPositions(x, numberOfVariables);
        forceFieldEnergyOpenMM.setOpenMMPeriodicBoxVectors();
        forceFieldEnergyOpenMM.setOpenMMVelocities(v, numberOfVariables);
        forceFieldEnergyOpenMM.setAcceleration(a);
        forceFieldEnergyOpenMM.setPreviousAcceleration(aPrevious);
    }

    /**
     * updateFromOpenMM obtains the state of the simulation from OpenMM,
     * completes some logging, and saves restart files.
     *
     * @param i       Number of OpenMM MD rounds.
     * @param running True if OpenMM MD rounds have begun running.
     */
    private void updateFromOpenMM(long i, boolean running) {

        double priorPE = currentPotentialEnergy;

        obtainVariables.run();

        double defaultDeltaPEThresh = 1.0E6;
        detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

        if (running) {
            if (i == 0) {
                logger.log(basicLogging, format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
                logger.log(basicLogging, format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
                logger.log(basicLogging, format("  %8s %12.4f %12.4f %12.4f %8.2f",
                        "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));
            }
            time = logThermoForTime(i, time);

            if (automaticWriteouts) {
                writeFilesForStep(i, true, true);
            }
        }
    }

    /**
     * <p>
     * integratorToString.</p>
     *
     * @param integrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
     *                   object.
     */
    private void integratorToString(IntegratorEnum integrator) {
        if (integrator == null) {
            integratorString = "VERLET";
            logger.info(" No specified integrator, will use Verlet");
        } else {
            switch (integratorType) {
                case STOCHASTIC:
                    integratorString = "LANGEVIN";
                    break;
                case VERLET:
                case VELOCITYVERLET:
                    integratorString = "VERLET";
                    break;
                case RESPA:
                    integratorString = "RESPA";
                    break;
                default:
                    integratorString = "VERLET";
                    logger.warning(String.format(" Integrator %s incompatible with "
                            + "OpenMM MD integration; defaulting to %s", integratorType, integratorString));
                    break;
            }
        }
    }

    /**
     * <p>
     * setMonteCarloBarostat.</p>
     *
     * @param temperature Temperature in Kelvin.
     */
    private void setMonteCarloBarostat(double temperature) {
        // Pressure is converted atm -> bar a few calls inwards.
        double pressure = barostat.getPressure();
        int frequency = barostat.getMeanBarostatInterval();
        forceFieldEnergyOpenMM.addMonteCarloBarostat(pressure, temperature, frequency);
    }
}
