// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.cli;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.Integrator;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.WriteoutOptions;
import java.util.logging.Logger;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that calculate thermodynamics.
 *
 * @author Michael J. Schnieders
 * @author Hernan V. Bernabe
 * @since 1.0
 */
public class DynamicsOptions {

  private static final Logger logger = Logger.getLogger(DynamicsOptions.class.getName());
  /** Thermostat. */
  public ThermostatEnum thermostat;
  /** Integrator. */
  public IntegratorEnum integrator;

  /**
   * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of 2.0 is possible for
   * the RESPA integrator.
   */
  @Option(
      names = {"-d", "--dt"},
      paramLabel = "1.0",
      defaultValue = "1.0",
      description = "Time discretization step in femtoseconds.")
  private double dt;

  /**
   * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or
   * Bussi.
   */
  @Option(
      names = {"-b", "--thermostat"},
      paramLabel = "Bussi",
      defaultValue = "Bussi",
      description = "Thermostat: [Adiabatic / Berendsen / Bussi].")
  private String thermostatString;

  /**
   * -i or --integrator sets the desired integrator: current choices are Beeman, RESPA, Stochastic
   * (i.e. Langevin dynamics) or Verlet.
   */
  @Option(
      names = {"-i", "--integrator"},
      paramLabel = "Verlet",
      defaultValue = "Verlet",
      description = "Integrator: [Beeman / Respa / Stochastic / Verlet].")
  private String integratorString;

  /**
   * -r or --report sets the thermodynamics reporting frequency in picoseconds (0.1 psec default).
   */
  @Option(
      names = {"-r", "--report"},
      paramLabel = "0.25",
      defaultValue = "0.25",
      description = "Interval in psec to report thermodynamics (psec).")
  private double report;

  /** -w or --write sets snapshot save frequency in picoseconds (1.0 psec default). */
  @Option(
      names = {"-w", "--write"},
      paramLabel = "10.0",
      defaultValue = "10.0",
      description = "Interval in psec to write out coordinates (psec).")
  private double write;

  /** -t or --temperature sets the simulation temperature (Kelvin). */
  @Option(
      names = {"-t", "--temperature"},
      paramLabel = "298.15",
      defaultValue = "298.15",
      description = "Temperature (Kelvin).")
  private double temperature;

  /** -n or --steps sets the number of molecular dynamics steps (default is 1 nsec). */
  @Option(
      names = {"-n", "--numberOfSteps"},
      paramLabel = "1000000",
      defaultValue = "1000000",
      description = "Number of molecular dynamics steps.")
  private long steps;

  /** -z or --trajSteps Number of steps for each OpenMM MD cycle. */
  @Option(
      names = {"-z", "--trajSteps"},
      paramLabel = "100",
      defaultValue = "100",
      description = "Number of steps per MD cycle (--mdE = OpenMM only).")
  private int trajSteps;

  /**
   * -o or --optimize saves low-energy snapshots discovered (only for single topology simulations).
   */
  @Option(
      names = {"-o", "--optimize"},
      defaultValue = "false",
      description = "Optimize and save low-energy snapshots.")
  private boolean optimize;

  /** -k or --checkpoint sets the restart save frequency in picoseconds (1.0 psec default). */
  @Option(
      names = {"-k", "--checkpoint"},
      paramLabel = "1.0",
      defaultValue = "1.0",
      description = "Interval in psec to write out restart files (.dyn, .his, etc).")
  private double checkpoint = 1.0;

  /**
   * --mdE or --molecularDynamicsEngine over-rides the default engine choice for integrating the
   * equations of motion
   */
  @Option(
      names = {"--mdE", "--molecularDynamicsEngine"},
      paramLabel = "FFX",
      description = "Use FFX or OpenMM to integrate dynamics.")
  private String engineString = null;

  private MolecularDynamics.DynamicsEngine engine = null;

  /**
   * The restart save frequency in picoseconds (1.0 psec default).
   *
   * @return Restart file interval in picoseconds.
   */
  public double getCheckpoint() {
    return checkpoint;
  }

  public void setCheckpoint(double checkpoint) {
    this.checkpoint = checkpoint;
  }

  /**
   * The timestep in femtoseconds (default of 1.0). A value of 2.0 is possible for the RESPA
   * integrator.
   *
   * @return Timestep in femtoseconds.
   */
  public double getDt() {
    return dt;
  }

  public void setDt(double dt) {
    this.dt = dt;
  }

  /**
   * Initialize a MolecularDynamics from the parsed options.
   *
   * @param potential a {@link ffx.numerics.Potential} object.
   * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param algorithmListener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param writeoutOptions a {@link WriteoutOptions} object.
   * @return a {@link MolecularDynamics} object.
   */
  public MolecularDynamics getDynamics(
      WriteoutOptions writeoutOptions,
      Potential potential,
      MolecularAssembly activeAssembly,
      AlgorithmListener algorithmListener) {
    MolecularDynamics molDyn;
    if (engine == null) {
      molDyn =
          MolecularDynamics.dynamicsFactory(
              activeAssembly,
              potential,
              activeAssembly.getProperties(),
              algorithmListener,
              thermostat,
              integrator);
    } else {
      molDyn =
          MolecularDynamics.dynamicsFactory(
              activeAssembly,
              potential,
              activeAssembly.getProperties(),
              algorithmListener,
              thermostat,
              integrator,
              engine);
    }
    molDyn.setFileType(writeoutOptions.getFileType());
    molDyn.setRestartFrequency(checkpoint);
    molDyn.setIntervalSteps(trajSteps);

    return molDyn;
  }

  public long getNumSteps() {
    return steps;
  }

  /**
   * Getter for the field <code>optimize</code>.
   *
   * @return Whether to optimize structures.
   */
  public boolean getOptimize() {
    return optimize;
  }

  /**
   * The thermodynamics reporting frequency in picoseconds (0.1 psec default).
   *
   * @return Thermodynamics logging interval in picoseconds.
   */
  public double getReport() {
    return report;
  }

  public void setReport(double report) {
    this.report = report;
  }

  /**
   * Write/snapshot appending interval.
   *
   * @return Interval between appending snapshots in psec.
   */
  public double getSnapshotInterval() {
    return write;
  }

  /**
   * The simulation temperature (Kelvin).
   *
   * @return Temperature in Kelvins.
   */
  public double getTemperature() {
    return temperature;
  }

  public void setTemperature(double temperature) {
    this.temperature = temperature;
  }

  /** Parse the thermostat and integrator. */
  public void init() {
    thermostat = Thermostat.parseThermostat(thermostatString);
    integrator = Integrator.parseIntegrator(integratorString);
    if (engineString != null) {
      try {
        engine = MolecularDynamics.DynamicsEngine.valueOf(engineString.toUpperCase());
      } catch (Exception ex) {
        logger.warning(
            String.format(
                " Could not parse %s as a valid dynamics engine! Defaulting to the Platform-recommended engine.",
                engineString));
        engine = null;
      }
    }
  }

  /**
   * The thermostat; intended for MC-OST going to ADIABATIC.
   *
   * @param thermo Thermostat to replace the requested one with.
   */
  public void setThermostat(ThermostatEnum thermo) {
    thermostat = thermo;
  }

  /**
   * The desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
   *
   * @return Returns a String for the requested thermostat.
   */
  public String getThermostatString() {
    return thermostatString;
  }

  public void setThermostatString(String thermostatString) {
    this.thermostatString = thermostatString;
  }

  /**
   * The integrator: current choices are Beeman, RESPA, Stochastic (Langevin) or Verlet.
   *
   * @return Returns a String for the requested integrator.
   */
  public String getIntegratorString() {
    return integratorString;
  }

  public void setIntegratorString(String integratorString) {
    this.integratorString = integratorString;
  }

  /**
   * Snapshot save frequency in picoseconds (1.0 psec default).
   *
   * @return Returns the frequency to save snapshots.
   */
  public double getWrite() {
    return write;
  }

  public void setWrite(double write) {
    this.write = write;
  }

  /**
   * The number of molecular dynamics steps (default is 1,000,000).
   *
   * @return Returns the number of MD time steps.
   */
  public long getSteps() {
    return steps;
  }

  public void setSteps(long steps) {
    this.steps = steps;
  }

  /**
   * Number of steps for each OpenMM MD cycle.
   *
   * @return Returns the number of steps for OpenMM MD cycles.
   */
  public int getTrajSteps() {
    return trajSteps;
  }

  public void setTrajSteps(int trajSteps) {
    this.trajSteps = trajSteps;
  }

  /**
   * Saves low-energy snapshots discovered (only for single topology simulations).
   *
   * @return Returns true if low-energy snapshots should be saved.
   */
  public boolean isOptimize() {
    return optimize;
  }

  public void setOptimize(boolean optimize) {
    this.optimize = optimize;
  }

  /**
   * The default engine choice for integrating the equations of motion
   *
   * @return Returns a String for the requested engine.
   */
  public String getEngineString() {
    return engineString;
  }

  public void setEngineString(String engineString) {
    this.engineString = engineString;
  }
}
