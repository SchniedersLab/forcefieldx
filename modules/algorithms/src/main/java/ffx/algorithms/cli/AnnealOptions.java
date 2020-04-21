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

import static java.lang.String.format;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.anneal.AnnealingSchedule;
import ffx.algorithms.optimize.anneal.FlatEndAnnealSchedule;
import ffx.algorithms.optimize.anneal.SimulatedAnnealing;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import java.io.File;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize simulated annealing.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class AnnealOptions {

  private static final Logger logger = Logger.getLogger(AnnealOptions.class.getName());

  /** -w or --windows Number of annealing windows (10). */
  @Option(
      names = {"-W", "--windows"},
      paramLabel = "10",
      defaultValue = "10",
      description = "Number of annealing windows.")
  private int windows;

  /** -l or --low Low temperature limit in degrees Kelvin (10.0). */
  @Option(
      names = {"--tl", "--temperatureLow"},
      paramLabel = "10.0",
      defaultValue = "10.0",
      description = "Low temperature limit (Kelvin).")
  private double low;

  /** -u or --upper Upper temperature limit in degrees Kelvin (1000.0). */
  @Option(
      names = {"--tu", "--temperatureUpper"},
      paramLabel = "1000.0",
      defaultValue = "1000.0",
      description = "High temperature limit (Kelvin).")
  private double upper;

  /**
   * --rv or --reinitVelocities forces simulated annealing to re-initialize velocities to the new
   * temperature at each annealing step, rather than letting the thermostat shift temperature
   * downwards.
   */
  @Option(
      names = {"--rv", "--reinitVelocities"},
      paramLabel = "false",
      defaultValue = "false",
      description = "Re-initialize velocities before each round of annealing.")
  private boolean reinitVelocities;

  /** --tmS or --temperingSchedule sets the schedule to be used. */
  @Option(
      names = {"--tmS", "--temperingSchedule"},
      paramLabel = "EXP",
      defaultValue = "EXP",
      description = "Tempering schedule: choose between EXP (exponential) or LINEAR")
  private String temperString;

  /**
   * --tmB or --temperingBefore sets the number of annealing windows to hold flat at the high
   * temperature (in addition to normal windows).
   */
  @Option(
      names = {"--tmB", "--temperingBefore"},
      paramLabel = "0",
      defaultValue = "0",
      description = "Number of (annealing, not MD/MC) steps to remain at the high temperature")
  private int temperBefore;

  /**
   * --tmA or --temperingAfter sets the number of annealing windows to hold flat at the low
   * temperature (in addition to normal windows).
   */
  @Option(
      names = {"--tmA", "--temperingAfter"},
      paramLabel = "0",
      defaultValue = "0",
      description = "Number of (annealing, not MD/MC) steps to remain at the low temperature")
  private int temperAfter;

  /**
   * Creates a SimulatedAnnealing object.
   *
   * @param dynamicsOptions Dynamics options to use.
   * @param molecularAssembly MolecularAssembly
   * @param potential Potential
   * @param compositeConfiguration Properties
   * @param algorithmListener AlgorithmListener
   * @return SimulatedAnnealing
   */
  public SimulatedAnnealing createAnnealer(
      DynamicsOptions dynamicsOptions,
      MolecularAssembly molecularAssembly,
      Potential potential,
      CompositeConfiguration compositeConfiguration,
      AlgorithmListener algorithmListener) {
    return createAnnealer(
        dynamicsOptions,
        molecularAssembly,
        potential,
        compositeConfiguration,
        algorithmListener,
        null);
  }

  /**
   * Creates a SimulatedAnnealing object.
   *
   * @param dynamicsOptions Dynamics options to use.
   * @param molecularAssembly MolecularAssembly
   * @param potential Potential
   * @param compositeConfiguration Properties
   * @param algorithmListener AlgorithmListener
   * @param dynFile Dynamics restart file.
   * @return SimulatedAnnealing
   */
  public SimulatedAnnealing createAnnealer(
      DynamicsOptions dynamicsOptions,
      MolecularAssembly molecularAssembly,
      Potential potential,
      CompositeConfiguration compositeConfiguration,
      AlgorithmListener algorithmListener,
      File dynFile) {
    AnnealingSchedule schedule = getSchedule();
    double totNormLen = schedule.totalWindowLength();
    long totSteps = dynamicsOptions.getSteps();
    long perWindowSteps = (long) (totSteps / totNormLen);

    if (totSteps > perWindowSteps * totNormLen) {
      ++perWindowSteps;
    }

    long minWindowSteps = (long) (perWindowSteps * schedule.minWindowLength());
    long maxWindowSteps = (long) (perWindowSteps * schedule.maxWindowLength());
    int nWindows = schedule.getNumWindows();

    if (minWindowSteps == maxWindowSteps && minWindowSteps == perWindowSteps) {
      logger.info(
          format(
              " Each of %d simulated annealing windows will have %d steps each, for a "
                  + "total duration of %d timesteps",
              nWindows, perWindowSteps, perWindowSteps * nWindows));
    } else {
      logger.info(
          format(
              " Each of %d simulated annealing windows will have %d-%d steps each, "
                  + "with a \"normal\" length of %d steps, for a total duration of %d timesteps",
              nWindows,
              minWindowSteps,
              maxWindowSteps,
              perWindowSteps,
              (int) (perWindowSteps * schedule.totalWindowLength())));
    }

    if (nWindows < 201) {
      StringBuilder sb =
          new StringBuilder(
              "\n Simulated annealing windows [index,MD steps, temperature (K)]:\n [");
      for (int i = 0; i < nWindows; i++) {
        double len = schedule.windowLength(i);
        double temp = schedule.getTemperature(i);
        sb.append(format("[%d,%d,%10.4g]", (i + 1), (int) (len * perWindowSteps), temp));
        if (i == nWindows - 1) {
          sb.append("]\n");
        } else if (i % 10 == 9) {
          sb.append("\n");
        } else {
          sb.append(",");
        }
      }
      logger.info(sb.toString());
    } else {
      logger.info(
          " Skipping printout of window lengths/temperatures (max printout at 200 windows)");
    }

    return new SimulatedAnnealing(
        molecularAssembly,
        potential,
        compositeConfiguration,
        algorithmListener,
        dynamicsOptions.thermostat,
        dynamicsOptions.integrator,
        schedule,
        perWindowSteps,
        dynamicsOptions.getDt(),
        isReinitVelocities(),
        dynFile);
  }

  /**
   * Constructs an AnnealingSchedule.
   *
   * @return An AnnealingSchedule.
   */
  public AnnealingSchedule getSchedule() {
    SimulatedAnnealing.Schedules schedules = SimulatedAnnealing.Schedules.parse(getTemperString());
    AnnealingSchedule annealingSchedule = schedules.generate(getWindows(), getLow(), getUpper());
    if (getTemperBefore() > 0 || getTemperAfter() > 0) {
      annealingSchedule =
          new FlatEndAnnealSchedule(
              annealingSchedule, getLow(), getUpper(), getTemperBefore(), getTemperAfter());
    }
    return annealingSchedule;
  }

  /**
   * Number of annealing windows.
   *
   * @return Returns the number of windows.
   */
  public int getWindows() {
    return windows;
  }

  public void setWindows(int windows) {
    this.windows = windows;
  }

  /**
   * Low temperature limit in degrees Kelvin.
   *
   * @return Returns the low temperature limit.
   */
  public double getLow() {
    return low;
  }

  public void setLow(double low) {
    this.low = low;
  }

  /**
   * Upper temperature limit in degrees Kelvin.
   *
   * @return Returns the upper temperature limit.
   */
  public double getUpper() {
    return upper;
  }

  public void setUpper(double upper) {
    this.upper = upper;
  }

  /**
   * Forces simulated annealing to re-initialize velocities to the new temperature at each annealing
   * step, rather than letting the thermostat shift temperature downwards.
   *
   * @return Returns true for re-initialization of velocities.
   */
  public boolean isReinitVelocities() {
    return reinitVelocities;
  }

  public void setReinitVelocities(boolean reinitVelocities) {
    this.reinitVelocities = reinitVelocities;
  }

  /**
   * Sets the schedule to be used.
   *
   * @return Returns a String representation of the tempering schedule.
   */
  public String getTemperString() {
    return temperString;
  }

  public void setTemperString(String temperString) {
    this.temperString = temperString;
  }

  /**
   * Sets the number of annealing windows to hold flat at the high temperature (in addition to
   * normal windows).
   *
   * @return Returns the number of annealing windows to hold flat at the high temperature.
   */
  public int getTemperBefore() {
    return temperBefore;
  }

  public void setTemperBefore(int temperBefore) {
    this.temperBefore = temperBefore;
  }

  /**
   * Sets the number of annealing windows to hold flat at the low temperature (in addition to normal
   * windows).
   *
   * @return Returns the number of annealing windows to hold flat at the low temperature.
   */
  public int getTemperAfter() {
    return temperAfter;
  }

  public void setTemperAfter(int temperAfter) {
    this.temperAfter = temperAfter;
  }
}
