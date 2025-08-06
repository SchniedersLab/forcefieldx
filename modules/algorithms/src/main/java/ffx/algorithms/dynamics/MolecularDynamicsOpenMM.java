// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.MolecularAssembly;
import ffx.potential.openmm.OpenMMContext;
import ffx.potential.openmm.OpenMMPotential;
import ffx.potential.openmm.OpenMMState;
import ffx.potential.openmm.OpenMMSystem;
import ffx.potential.UnmodifiableState;

import java.io.File;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Michael J. Schnieders
 */
public class MolecularDynamicsOpenMM extends MolecularDynamics {

  private static final Logger logger = Logger.getLogger(MolecularDynamicsOpenMM.class.getName());
  /**
   * Integrator Type.
   */
  private final IntegratorEnum integratorType;
  /**
   * Thermostat Type.
   */
  private final ThermostatEnum thermostatType;
  /**
   * The potential energy function.
   */
  private final CrystalPotential crystalPotential;
  /**
   * Features specific to OpenMM.
   */
  private final OpenMMPotential openMMPotential;
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
   * Obtain all variables with each update (i.e., include velocities, gradients).
   */
  private boolean getAllVars = true;
  /**
   * Method to run on update for obtaining variables. Will either grab everything (default) or
   * energies + positions (MC-OST).
   */
  private Runnable obtainVariables = this::getAllOpenMMVariables;

  /**
   * Constructs an MolecularDynamicsOpenMM object, to perform molecular dynamics using native OpenMM
   * routines, avoiding the cost of communicating coordinates, gradients, and energies back and forth
   * across the PCI bus.
   *
   * @param assembly   MolecularAssembly to operate on
   * @param potential  Either a ForceFieldEnergyOpenMM, or a Barostat.
   * @param listener   a {@link ffx.algorithms.AlgorithmListener} object.
   * @param thermostat May have to be slightly modified for native OpenMM routines
   * @param integrator May have to be slightly modified for native OpenMM routines
   */
  public MolecularDynamicsOpenMM(MolecularAssembly assembly, Potential potential,
                                 AlgorithmListener listener, ThermostatEnum thermostat, IntegratorEnum integrator) {
    super(assembly, potential, listener, thermostat, integrator);

    logger.info("\n Initializing OpenMM molecular dynamics.");

    // Initialization specific to MolecularDynamicsOpenMM
    running = false;
    crystalPotential = (CrystalPotential) potential;
    openMMPotential = (OpenMMPotential) potential;

    // Update the set of active and inactive atoms.
    openMMPotential.setActiveAtoms();

    thermostatType = thermostat;
    integratorType = integrator;
    integratorToString(integratorType);
  }

  /**
   * Set the barostat for this MolecularDynamicsOpenMM instance.
   * @param barostat The Barostat to set, or null to disable constant pressure.
   */
  public void setBarostat(Barostat barostat) {
    if (barostat != null) {
      this.barostat = barostat;
      this.constantPressure = true;
      barostat.setActive(false);
    } else {
      this.barostat = null;
      this.constantPressure = false;
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Execute <code>numSteps</code> of dynamics using the provided <code>timeStep</code> and
   * <code>temperature</code>. The <code>printInterval</code> and <code>saveInterval</code>
   * control logging the state of the system to the console and writing a restart file, respectively.
   * If the <code>dyn</code> File is not null, the simulation will be initialized from the contents.
   * If the <code>iniVelocities</code> is true, the velocities will be initialized from a Maxwell
   * Boltzmann distribution.
   */
  @Override
  public void dynamic(long numSteps, double timeStep, double printInterval, double saveInterval,
                      double temperature, boolean initVelocities, File dyn) {
    // Return if already running and a second thread calls the dynamic method.
    if (!done) {
      logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
      return;
    }

    // Call the init method.
    init(numSteps, timeStep, printInterval, saveInterval, fileType, restartInterval, temperature,
        initVelocities, dyn);

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

    // Store the initial state.
    initialState = new UnmodifiableState(state);

    // Check that our context is using correct Integrator, time step, and target temperature.
    openMMPotential.updateContext(integratorString, dt, targetTemperature, false);

    // Pre-run operations (mostly logging) that require knowledge of system energy.
    postInitEnergies();

    // Run the MD steps.
    mainLoop(numSteps);

    // Post-run cleanup operations.
    postRun();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getIntervalSteps() {
    return intervalSteps;
  }

  /**
   * Setter for the field <code>intervalSteps</code>.
   *
   * @param intervalSteps The number of interval steps.
   */
  @Override
  public void setIntervalSteps(int intervalSteps) {
    this.intervalSteps = intervalSteps;
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
  public void init(long numSteps, double timeStep, double loggingInterval, double trajectoryInterval,
                   String fileType, double restartInterval, double temperature, boolean initVelocities,
                   File dyn) {

    super.init(numSteps, timeStep, loggingInterval, trajectoryInterval, fileType, restartInterval,
        temperature, initVelocities, dyn);

    boolean isLangevin = IntegratorEnum.isStochastic(integratorType);

    OpenMMSystem openMMSystem = openMMPotential.getSystem();
    if (!isLangevin && !thermostatType.equals(ThermostatEnum.ADIABATIC)) {
      // Add Andersen thermostat, or if already present update its target temperature.
      openMMSystem.addAndersenThermostatForce(targetTemperature);
    }

    if (constantPressure) {
      // Add an isotropic Monte Carlo barostat.
      // If it is already present, update its target temperature, pressure and frequency.
      double pressure = barostat.getPressure();
      int frequency = barostat.getMeanBarostatInterval();
      openMMSystem.addMonteCarloBarostatForce(pressure, targetTemperature, frequency);
    }

    // For Langevin/Stochastic dynamics, center of mass motion will not be removed.
    if (!isLangevin) {
      // No action is taken if a COMMRemover is already present.
      openMMSystem.addCOMMRemoverForce();
    }

    // Set the current value of lambda.
    if (crystalPotential instanceof LambdaInterface lambdaInferface) {
      lambdaInferface.setLambda(lambdaInferface.getLambda());
    }

  }

  @Override
  public void revertState() throws Exception {
    super.revertState();
    setOpenMMState();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setFileType(String fileType) {
    this.fileType = fileType;
  }

  /**
   * Sets whether to obtain all variables (velocities, gradients) from OpenMM, or just positions and
   * energies.
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
    OpenMMContext openMMContext = openMMPotential.getContext();
    openMMContext.integrate(intervalSteps);
  }

  /**
   * Load coordinates, box vectors and velocities.
   */
  private void setOpenMMState() {
    OpenMMContext openMMContext = openMMPotential.getContext();
    openMMContext.setPositions(state.x());
    openMMContext.setPeriodicBoxVectors(crystalPotential.getCrystal());
    openMMContext.setVelocities(state.v());
  }

  /**
   * Get OpenMM Energies.
   */
  private void getOpenMMEnergies() {
    OpenMMState openMMState = openMMPotential.getOpenMMState(OpenMM_State_Energy);
    state.setKineticEnergy(openMMState.kineticEnergy);
    state.setPotentialEnergy(openMMState.potentialEnergy);
    state.setTemperature(openMMPotential.getSystem().getTemperature(openMMState.kineticEnergy));
    openMMState.destroy();
  }

  /**
   * Do some logging of the beginning energy values.
   */
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
   * Get OpenMM Energies and Positions.
   */
  private void getOpenMMEnergiesAndPositions() {
    int mask = OpenMM_State_Energy | OpenMM_State_Positions;
    OpenMMState openMMState = openMMPotential.getOpenMMState(mask);
    state.setPotentialEnergy(openMMState.potentialEnergy);
    state.setKineticEnergy(openMMState.kineticEnergy);
    state.setTemperature(openMMPotential.getSystem().getTemperature(openMMState.kineticEnergy));
    openMMState.getPositions(state.x());


    Crystal crystal = crystalPotential.getCrystal();
    if (!crystal.aperiodic()) {
      double[][] cellVectors = openMMState.getPeriodicBoxVectors();
      crystal.setCellVectors(cellVectors);
      crystalPotential.setCrystal(crystal);
    }
    openMMState.destroy();
  }

  /**
   * Get OpenMM energies, positions, velocities, and accelerations.
   */
  private void getAllOpenMMVariables() {
    int mask = OpenMM_State_Energy | OpenMM_State_Positions | OpenMM_State_Velocities | OpenMM_State_Forces;
    OpenMMState openMMState = openMMPotential.getOpenMMState(mask);
    state.setPotentialEnergy(openMMState.potentialEnergy);
    state.setKineticEnergy(openMMState.kineticEnergy);
    state.setTemperature(openMMPotential.getSystem().getTemperature(openMMState.kineticEnergy));
    openMMState.getPositions(state.x());
    Crystal crystal = crystalPotential.getCrystal();
    if (!crystal.aperiodic()) {
      double[][] cellVectors = openMMState.getPeriodicBoxVectors();
      crystal.setCellVectors(cellVectors);
      crystalPotential.setCrystal(crystal);
    }
    openMMState.getVelocities(state.v());
    openMMState.getAccelerations(state.a());
    openMMState.destroy();
  }

  /**
   * updateFromOpenMM obtains the state of the simulation from OpenMM, completes some logging, and
   * saves restart files.
   *
   * @param i       Number of OpenMM MD rounds.
   * @param running True if OpenMM MD rounds have begun running.
   */
  private void updateFromOpenMM(long i, boolean running) {

    double priorPE = state.getPotentialEnergy();

    obtainVariables.run();

    if (running) {
      if (i == 0) {
        logger.log(basicLogging,
            format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp",
                "CPU"));
        logger.log(basicLogging,
            format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K",
                "sec"));
        logger.log(basicLogging,
            format("  %8s %12.4f %12.4f %12.4f %8.2f", "", state.getKineticEnergy(),
                state.getPotentialEnergy(), state.getTotalEnergy(), state.getTemperature()));
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
      logger.info(" An integrator was not specified. Verlet will be used.");
    } else {
      switch (integratorType) {
        default -> integratorString = "VERLET";
        case STOCHASTIC, LANGEVIN -> integratorString = "LANGEVIN";
        case RESPA, MTS -> integratorString = "MTS";
        case STOCHASTIC_MTS, LANGEVIN_MTS -> integratorString = "LANGEVIN-MTS";
      }
    }
  }
}
