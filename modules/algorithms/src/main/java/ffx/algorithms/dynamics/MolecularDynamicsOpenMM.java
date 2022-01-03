// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.algorithms.dynamics;

import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.ForceFieldEnergyOpenMM.Context;
import ffx.potential.ForceFieldEnergyOpenMM.State;
import ffx.potential.MolecularAssembly;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import ffx.potential.extended.ExtendedSystem;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Michael J. Schnieders
 */
public class MolecularDynamicsOpenMM extends MolecularDynamics {

  private static final Logger logger = Logger.getLogger(MolecularDynamicsOpenMM.class.getName());
  /** Integrator Type. */
  private final IntegratorEnum integratorType;
  /** Thermostat Type. */
  private final ThermostatEnum thermostatType;
  /** OpenMM ForceFieldEnergy. */
  private final ForceFieldEnergyOpenMM forceFieldEnergyOpenMM;
  /** Integrator String. */
  private String integratorString;
  /** Number of OpenMM MD steps per iteration. */
  private int intervalSteps;
  /** Flag to indicate OpenMM MD interactions are running. */
  private boolean running;
  /** Run time. */
  private long time;
  /** Obtain all variables with each update (i.e. include velocities, gradients). */
  private boolean getAllVars = true;
  /**
   * Method to run on update for obtaining variables. Will either grab everything (default) or
   * energies + positions (MC-OST).
   */
  private Runnable obtainVariables = this::getAllOpenMMVariables;

  /**
   * Constructs an MolecularDynamicsOpenMM object, to perform molecular dynamics using native OpenMM
   * routines, avoiding the cost of communicating coordinates, gradients, and energies back and
   * forth across the PCI bus.
   *
   * @param assembly MolecularAssembly to operate on
   * @param potential Either a ForceFieldEnergyOpenMM, or a Barostat.
   * @param properties Associated properties
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param thermostat May have to be slightly modified for native OpenMM routines
   * @param integrator May have to be slightly modified for native OpenMM routines
   */
  public MolecularDynamicsOpenMM(
      MolecularAssembly assembly,
      Potential potential,
      CompositeConfiguration properties,
      AlgorithmListener listener,
      ThermostatEnum thermostat,
      IntegratorEnum integrator) {
    super(assembly, potential, properties, listener, thermostat, integrator);

    logger.info("\n Initializing OpenMM molecular dynamics.");

    // Initialization specific to MolecularDynamicsOpenMM
    running = false;
    List<Potential> potentialStack = new ArrayList<>(potential.getUnderlyingPotentials());
    potentialStack.add(potential);

    List<ForceFieldEnergyOpenMM> feOMM =
        potentialStack.stream()
            .filter((Potential p) -> p instanceof ForceFieldEnergyOpenMM)
            .map((Potential p) -> (ForceFieldEnergyOpenMM) p)
            .collect(Collectors.toList());
    if (feOMM.size() != 1) {
      logger.severe(
          format(
              " Attempting to create a MolecularDynamicsOpenMM with %d OpenMM force field energies: this presently only allows one!",
              feOMM.size()));
    }
    forceFieldEnergyOpenMM = feOMM.get(0);

    List<Barostat> barostats =
        potentialStack.stream()
            .filter((Potential p) -> p instanceof Barostat)
            .map((Potential p) -> (Barostat) p)
            .collect(Collectors.toList());
    if (barostats.isEmpty()) {
      constantPressure = false;
    } else if (barostats.size() > 1) {
      logger.severe(
          format(
              " Attempting to create a MolecularDynamicsOpenMM with %d barostats: this presently only allows 0-1!",
              barostats.size()));
    } else {
      barostat = barostats.get(0);
      barostat.setActive(false);
    }

    // Update set active and inactive atoms.
    forceFieldEnergyOpenMM.setActiveAtoms();

    thermostatType = thermostat;
    integratorType = integrator;
    integratorToString(integratorType);

    // Pseudo-random number generator used to seed the OpenMM velocity generator method.
    Random random = new Random();
    if (properties.containsKey("velRandomSeed")) {
      random.setSeed(properties.getInt("velRandomSeed", 0));
    } else {
      random.setSeed(0);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling extended system
   * variables. Will throw an UnsupportedOperationException.
   */
  @Override
  public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
    throw new UnsupportedOperationException(
        " MolecularDynamicsOpenMM does not support extended system variables!");
  }

  /**
   * {@inheritDoc}
   *
   * <p>UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling extended system
   * variables. Will throw an UnsupportedOperationException.
   */
  @Override
  public void detachExtendedSystem() {
    throw new UnsupportedOperationException(
        " MolecularDynamicsOpenMM does not support extended system variables!");
  }

  /**
   * {@inheritDoc}
   *
   * <p>Start sets up context, write out file name, restart file name, sets the integrator and
   * determines whether the simulation is starting out from a previous molecular dynamics run (.dyn)
   * or if the initial velocities are determined by a Maxwell Boltzmann distribution. This method
   * then calls methods openMMUpdate and takeOpenMMSteps to run the molecular dynamics simulation.
   */
  @Override
  public void dynamic(
      long numSteps,
      double timeStep,
      double printInterval,
      double saveInterval,
      double temperature,
      boolean initVelocities,
      File dyn) {
    // Return if already running and a second thread calls the dynamic method.
    if (!done) {
      logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
      return;
    }

    init(numSteps, timeStep, printInterval, saveInterval, fileType, restartInterval,
        temperature, initVelocities, dyn);

    if (intervalSteps == 0 || intervalSteps > numSteps) {
      // Safe cast: if intervalSteps > numSteps, then numSteps must be less than Integer.MAX_VALUE.
      intervalSteps = (int) numSteps;
    }

    // Initialization, including reading a dyn file or initialization of velocity.
    preRunOps();

    // Send coordinates, velocities, etc. to OpenMM.
    setOpenMMState();

    // Retrieve starting energy values.
    getOpenMMEnergies();
    initialKinetic = currentKineticEnergy;
    initialPotential = currentPotentialEnergy;
    initialTotal = currentTotalEnergy;
    initialTemp = currentTemperature;

    // Check that our context is using correct Integrator, time step, and target temperature.
    forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature, false);

    // Pre-run operations (mostly logging) that require knowledge of system energy.
    postInitEnergies();

    // Run the MD steps.
    mainLoop(numSteps);

    // Post-run cleanup operations.
    postRun();
  }

  /** {@inheritDoc} */
  @Override
  public int getIntervalSteps() {
    return intervalSteps;
  }

  /**
   * Setter for the field <code>intervalSteps</code>.
   *
   * @param intervalSteps a int.
   */
  @Override
  public void setIntervalSteps(int intervalSteps) {
    this.intervalSteps = intervalSteps;
  }

  /** {@inheritDoc} */
  @Override
  public double getTimeStep() {
    return dt;
  }

  /** {@inheritDoc} */
  @Override
  public void init(
      long numSteps,
      double timeStep,
      double loggingInterval,
      double trajectoryInterval,
      String fileType,
      double restartInterval,
      double temperature,
      boolean initVelocities,
      File dyn) {

    super.init(
        numSteps,
        timeStep,
        loggingInterval,
        trajectoryInterval,
        fileType,
        restartInterval,
        temperature,
        initVelocities,
        dyn);

    boolean isLangevin = integratorType.equals(IntegratorEnum.STOCHASTIC);

    ForceFieldEnergyOpenMM.System system = forceFieldEnergyOpenMM.getSystem();
    if (!isLangevin && !thermostatType.equals(ThermostatEnum.ADIABATIC)) {
      // Add Andersen thermostat, or if already present update its target temperature.
      system.addAndersenThermostatForce(targetTemperature);
    }

    if (constantPressure) {
      // Add an isotropic Monte Carlo barostat, or if already present update its target temperature,
      // pressure
      // and frequency.
      double pressure = barostat.getPressure();
      int frequency = barostat.getMeanBarostatInterval();
      system.addMonteCarloBarostatForce(pressure, targetTemperature, frequency);
    }

    // For Langevin/Stochastic dynamics, center of mass motion will not be removed.
    if (!isLangevin) {
      // No action is taken if a COMMRemover is already present.
      system.addCOMMRemoverForce();
    }

    forceFieldEnergyOpenMM.setLambda(forceFieldEnergyOpenMM.getLambda());
  }

  @Override
  public void revertState() throws Exception {
    super.revertState();
    setOpenMMState();
  }

  /** {@inheritDoc} */
  @Override
  public void setFileType(String fileType) {
    this.fileType = fileType;
  }

  /**
   * Sets whether or not to obtain all variables (velocities, gradients) from OpenMM, or just
   * positions and energies.
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
  public void writeRestart() {
    if (!getAllVars) {
      // If !getAllVars, need to ensure all variables are synced before writing the restart.
      getAllOpenMMVariables();
    }
    super.writeRestart();
  }

  @Override
  protected void appendSnapshot(String[] extraLines) {
    if (!getAllVars) {
      // If !getAllVars, need to ensure coordinates are synced before writing a snapshot.
      getOpenMMEnergiesAndPositions();
    }
    super.appendSnapshot(extraLines);
  }

  /**
   * Integrate the simulation using the defined Context and Integrator.
   *
   * @param intervalSteps Number of MD steps to take.
   */
  private void takeOpenMMSteps(int intervalSteps) {
    Context context = forceFieldEnergyOpenMM.getContext();
    context.integrate(intervalSteps);
  }

  /** Load coordinates, box vectors and velocities. */
  private void setOpenMMState() {
    Context context = forceFieldEnergyOpenMM.getContext();
    context.setOpenMMPositions(x);
    context.setPeriodicBoxVectors();
    context.setOpenMMVelocities(v);
  }

  /** Get OpenMM Energies. */
  private void getOpenMMEnergies() {
    State state = forceFieldEnergyOpenMM.createState(false, true, false, false);
    currentPotentialEnergy = state.potentialEnergy;
    currentKineticEnergy = state.kineticEnergy;
    currentTotalEnergy = state.totalEnergy;
    currentTemperature = state.temperature;
    state.free();
  }

  /** Do some logging of the beginning energy values. */
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

  /** Get OpenMM Energies and Positions. */
  private void getOpenMMEnergiesAndPositions() {
    State state = forceFieldEnergyOpenMM.createState(true, true, false, false);
    currentPotentialEnergy = state.potentialEnergy;
    currentKineticEnergy = state.kineticEnergy;
    currentTotalEnergy = state.totalEnergy;
    currentTemperature = state.temperature;
    state.getPositions(x);
    state.getPeriodicBoxVectors();
    state.free();
  }

  /** Get OpenMM energies, positions, velocities, and accelerations. */
  private void getAllOpenMMVariables() {
    State state = forceFieldEnergyOpenMM.createState(true, true, true, true);
    currentPotentialEnergy = state.potentialEnergy;
    currentKineticEnergy = state.kineticEnergy;
    currentTotalEnergy = state.totalEnergy;
    currentTemperature = state.temperature;
    x = state.getPositions(x);
    state.getPeriodicBoxVectors();
    v = state.getVelocities(v);
    a = state.getAccelerations(a);
    state.free();
  }

  /**
   * updateFromOpenMM obtains the state of the simulation from OpenMM, completes some logging, and
   * saves restart files.
   *
   * @param i Number of OpenMM MD rounds.
   * @param running True if OpenMM MD rounds have begun running.
   */
  private void updateFromOpenMM(long i, boolean running) {

    double priorPE = currentPotentialEnergy;

    obtainVariables.run();

    double defaultDeltaPEThresh = 1.0E6;
    detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

    if (running) {
      if (i == 0) {
        logger.log(
            basicLogging,
            format(
                "\n  %8s %12s %12s %12s %8s %8s",
                "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
        logger.log(
            basicLogging,
            format(
                "  %8s %12s %12s %12s %8s %8s",
                "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
        logger.log(
            basicLogging,
            format(
                "  %8s %12.4f %12.4f %12.4f %8.2f",
                "",
                currentKineticEnergy,
                currentPotentialEnergy,
                currentTotalEnergy,
                currentTemperature));
      }
      time = logThermoForTime(i, time);

      if (automaticWriteouts) {
        writeFilesForStep(i, true, true);
      }
    }
  }

  /**
   * integratorToString.
   *
   * @param integrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum} object.
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
          logger.warning(
              String.format(
                  " Integrator %s incompatible with " + "OpenMM MD integration; defaulting to %s",
                  integratorType, integratorString));
          break;
      }
    }
  }
}
