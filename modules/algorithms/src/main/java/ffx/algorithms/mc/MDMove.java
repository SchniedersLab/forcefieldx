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
package ffx.algorithms.mc;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MDVerbosity;
import ffx.algorithms.dynamics.MDWriteAction;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.MolecularDynamicsOpenMM;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.SystemState;
import ffx.potential.parsers.DYNFilter;
import ffx.utilities.Constants;

import javax.annotation.Nullable;
import java.io.File;
import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

/**
 * Use MD as a coordinate based MC move.
 *
 * @author Mallory R. Tollefson
 */
public class MDMove implements MCMove {

  private static final Logger logger = Logger.getLogger(MDMove.class.getName());

  /**
   * THe MD instance that executes the move.
   */
  private final MolecularDynamics molecularDynamics;
  /**
   * Number of MD steps per move.
   */
  private final long mdSteps;
  /**
   * Time step in femtoseconds.
   */
  private final double timeStep;
  /**
   * Print interval in picoseconds.
   */
  private final double printInterval;
  /**
   * Temperature in Kelvin.
   */
  private final double temperature;
  /**
   * Snapshot appending interval in psec.
   */
  private final double saveInterval;
  /**
   * Potential to operate on.
   */
  private final Potential potential;
  /**
   * Number of MD moves.
   */
  private int mdMoveCounter = 0;
  /**
   * The total energy change for the current move.
   */
  private double energyChange;
  /**
   * Absolute energy drift over all moves.
   */
  private double energyDriftAbs;
  /**
   * Net energy drift over all moves.
   */
  private double energyDriftNet;
  /**
   * Kinetic energy at the start of the move.
   */
  private double initialKinetic;
  /**
   * Potential energy at the start of the move.
   */
  private double initialPotential;
  /**
   * Total energy at the start of the move.
   */
  private double initialTotal;

  private File dynRestartFile;

  /**
   * Constructor for MDMove.
   *
   * @param assembly        a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param listener        a {@link ffx.algorithms.AlgorithmListener} object.
   * @param dynamics        CLI object containing key MD information.
   * @param stepsPerCycle   Number of MD steps per MC cycle.
   * @param dynRestartFile  File to read restart from.
   */
  public MDMove(MolecularAssembly assembly, Potential potentialEnergy, AlgorithmListener listener,
                DynamicsOptions dynamics, long stepsPerCycle, @Nullable File dynRestartFile) {
    this.potential = potentialEnergy;
    molecularDynamics = MolecularDynamics.dynamicsFactory(assembly, potentialEnergy, listener,
        dynamics.thermostat, dynamics.integrator);
    molecularDynamics.setAutomaticWriteouts(false);

    timeStep = dynamics.getDt();
    double dtPs = timeStep * Constants.FSEC_TO_PSEC;
    mdSteps = stepsPerCycle;

    molecularDynamics.setVerbosityLevel(MDVerbosity.QUIET);
    molecularDynamics.setObtainVelAcc(false);
    molecularDynamics.setRestartFrequency(dynamics.getCheckpoint());
    this.saveInterval = dynamics.getSnapshotInterval();

    double requestedPrint = dynamics.getReport();
    double maxPrintInterval = dtPs * mdSteps;
    // Log at least once/cycle.
    printInterval = min(requestedPrint, maxPrintInterval);
    this.temperature = dynamics.getTemperature();

    this.dynRestartFile = dynRestartFile;
    // Load the restart file if it exists.
    if (this.dynRestartFile != null && this.dynRestartFile.exists()) {
      molecularDynamics.setVerbosityLevel(MDVerbosity.SILENT);
      molecularDynamics.dynamic(1, timeStep, printInterval, saveInterval, temperature, false, dynRestartFile);
      collectEnergies();
      revertMove();
      molecularDynamics.setVerbosityLevel(MDVerbosity.QUIET);
    }
  }

  /**
   * Get the total energy change for the current move.
   *
   * @return Total energy change.
   */
  public double getEnergyChange() {
    return energyChange;
  }

  public double getInitialKinetic() {
    return initialKinetic;
  }

  public double getInitialPotential() {
    return initialPotential;
  }

  public double getInitialTotal() {
    return initialTotal;
  }

  /**
   * getKineticEnergy.
   *
   * @return a double.
   */
  public double getKineticEnergy() {
    return molecularDynamics.getKineticEnergy();
  }

  public MolecularDynamics getMD() {
    return molecularDynamics;
  }

  /**
   * getPotentialEnergy.
   *
   * @return a double.
   */
  public double getPotentialEnergy() {
    return molecularDynamics.getPotentialEnergy();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void move() {
    move(MDVerbosity.QUIET);
  }

  /**
   * Performs an MDMove.
   *
   * @param verbosityLevel How verbose to be.
   */
  public void move(MDVerbosity verbosityLevel) {
    MDVerbosity origLevel = molecularDynamics.getVerbosityLevel();
    molecularDynamics.setVerbosityLevel(verbosityLevel);
    if (mdMoveCounter == 0 && dynRestartFile != null && dynRestartFile.exists()) {
      molecularDynamics.dynamic(mdSteps, timeStep, printInterval, saveInterval, temperature, false, dynRestartFile);
    } else {
      molecularDynamics.dynamic(mdSteps, timeStep, printInterval, saveInterval, temperature, true, null);
    }
    mdMoveCounter++;
    collectEnergies();
    energyChange = molecularDynamics.getTotalEnergy() - initialTotal;

    if (molecularDynamics instanceof MolecularDynamicsOpenMM && logger.isLoggable(Level.FINE)) {
      energyDriftNet += energyChange;
      energyDriftAbs += abs(energyChange);
      double energyDriftAverageNet = energyDriftNet / mdMoveCounter;
      double energyDriftAverageAbs = energyDriftAbs / mdMoveCounter;
      logger.fine(format(" Mean signed/unsigned energy drift:                   %8.4f/%8.4f",
          energyDriftAverageNet, energyDriftAverageAbs));

      double dt = molecularDynamics.getTimeStep();
      int intervalSteps = molecularDynamics.getIntervalSteps();
      int nAtoms = potential.getNumberOfVariables() / 3;
      // TODO: Determine if the *1000 factor is an old artifact of MolecularDynamicsOpenMM being the
      // one thing which (used to) store dt in fsec.
      double normalizedEnergyDriftNet = (energyDriftAverageNet / (dt * intervalSteps * nAtoms)) * 1000;
      double normalizedEnergyDriftAbs = (energyDriftAverageAbs / (dt * intervalSteps * nAtoms)) * 1000;
      logger.fine(format(" Mean singed/unsigned energy drift per psec per atom: %8.4f/%8.4f\n",
          normalizedEnergyDriftNet, normalizedEnergyDriftAbs));
    }
    molecularDynamics.setVerbosityLevel(origLevel);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void revertMove() {
    try {
      molecularDynamics.revertState();
    } catch (Exception ex) {
      logger.severe(" The MD state could not be reverted.");
    }
  }

  public void setMDIntervalSteps(int intervalSteps) {
    molecularDynamics.setIntervalSteps(intervalSteps);
  }

  /**
   * Write restart and trajectory files if the provided step matches the frequency.
   *
   * @param mdStep      MD step (not MC cycle number) to write files (if any) for.
   * @param trySnapshot If false, do not write snapshot even if the time step is correct.
   * @param tryRestart  If false, do not write a restart file even if the time step is correct.
   * @return Returns the write actions.
   */
  public EnumSet<MDWriteAction> writeFilesForStep(long mdStep, boolean trySnapshot, boolean tryRestart) {
    return molecularDynamics.writeFilesForStep(mdStep, trySnapshot, tryRestart);
  }

  private void collectEnergies() {
    initialTotal = molecularDynamics.getInitialTotalEnergy();
    initialKinetic = molecularDynamics.getInitialKineticEnergy();
    initialPotential = molecularDynamics.getInitialPotentialEnergy();
    // If total energy is non-tiny, assert that kinetic + potential = total to within tolerance.
    assert abs(initialTotal) < 1.0E-3
        || abs((initialKinetic + initialPotential - initialTotal) / initialTotal) < 1.0E-7;
  }
}
