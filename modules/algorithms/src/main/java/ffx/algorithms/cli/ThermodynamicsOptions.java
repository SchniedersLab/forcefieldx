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
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.WriteoutOptions;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;
import picocli.CommandLine;

/**
 * Represents command line options for scripts that calculate thermodynamics.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class ThermodynamicsOptions {
  private static final Logger logger = Logger.getLogger(ThermodynamicsOptions.class.getName());

  /**
   * -Q or --equilibrate sets the number of equilibration steps prior to production OST counts
   * begin.
   */
  @CommandLine.Option(
      names = {"-Q", "--equilibrate"},
      paramLabel = "1000",
      description = "Number of equilibration steps before evaluation of thermodynamics.")
  private long nEquil = 1000;

  /**
   * -rn or --resetNumSteps, ignores steps detected in .lam lambda-restart files and thus resets the
   * histogram; use -rn false to continue from the end of any prior simulation.
   */
  @CommandLine.Option(
      names = {"--rn", "--resetNumSteps"},
      description = "Ignore prior steps logged in .lam or similar files")
  private boolean resetNumSteps = false;

  /**
   * --tA or --thermodynamicsAlgorithm specifies the algorithm to be used; currently serves as a
   * switch between OST and window-based methods.
   */
  @CommandLine.Option(
      names = {"--tA", "--thermodynamicsAlgorithm"},
      paramLabel = "OST",
      description =
          "Choice of thermodynamics algorithm. THe default is OST, while FIXED runs MD at a fixed lambda value (e.g. BAR)")
  private String thermoAlgoString = "OST";

  /**
   * Return the selected Thermodynamics algorithm as an enumerated type.
   *
   * @return Corresponding thermodynamics algorithm
   */
  public ThermodynamicsAlgorithm getAlgorithm() {
    return ThermodynamicsAlgorithm.parse(thermoAlgoString);
  }

  /**
   * getEquilSteps.
   *
   * @return a int.
   */
  public long getEquilSteps() {
    return nEquil;
  }

  /**
   * Getter for the field <code>resetNumSteps</code>.
   *
   * @return a boolean.
   */
  public boolean getResetNumSteps() {
    return resetNumSteps;
  }

  /**
   * Run an alchemical free energy window.
   *
   * @param topologies All involved MolecularAssemblies.
   * @param potential The Potential to be sampled.
   * @param dynamics DynamicsOptions.
   * @param writeOut WriteoutOptions
   * @param dyn MD restart file
   * @param aListener AlgorithmListener
   * @return The MolecularDynamics object constructed.
   */
  public MolecularDynamics runFixedAlchemy(
      MolecularAssembly[] topologies,
      CrystalPotential potential,
      DynamicsOptions dynamics,
      WriteoutOptions writeOut,
      File dyn,
      AlgorithmListener aListener) {
    dynamics.init();

    MolecularDynamics molDyn = dynamics.getDynamics(writeOut, potential, topologies[0], aListener);
    for (int i = 1; i < topologies.length; i++) {
      molDyn.addAssembly(topologies[i], topologies[i].getProperties());
    }

    boolean initVelocities = true;
    long nSteps = dynamics.steps;
    molDyn.setRestartFrequency(dynamics.getCheckpoint());
    // Start sampling.
    if (nEquil > 0) {
      logger.info("\n Beginning equilibration");
      runDynamics(molDyn, nEquil, dynamics, writeOut, true, dyn);
      logger.info(" Beginning fixed-lambda alchemical sampling");
      initVelocities = false;
    } else {
      logger.info(" Beginning fixed-lambda alchemical sampling without equilibration");
      if (!resetNumSteps) {
        /*int nEnergyCount = ttOSRW.getEnergyCount();
        if (nEnergyCount > 0) {
            nSteps -= nEnergyCount;
            logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
            initVelocities = false;
        }*/
        // Temporary workaround for being unable to pick up preexisting steps.
        initVelocities = true;
      }
    }

    if (nSteps > 0) {
      runDynamics(molDyn, nSteps, dynamics, writeOut, initVelocities, dyn);
    } else {
      logger.info(" No steps remaining for this process!");
    }
    return molDyn;
  }

  private void runDynamics(
      MolecularDynamics md,
      long nSteps,
      DynamicsOptions dynamics,
      WriteoutOptions writeOut,
      boolean initVelocities,
      File dyn) {
    md.dynamic(
        nSteps,
        dynamics.dt,
        dynamics.report,
        dynamics.write,
        dynamics.temp,
        initVelocities,
        writeOut.getFileType(),
        dynamics.getCheckpoint(),
        dyn);
  }

  /**
   * Represents categories of thermodynamics algorithms that must be handled differentially. For
   * legacy reasons, MC-OST and MD-OST are both just "OST", and the differences are handled in
   * OSTOptions and Thermodynamics.groovy. Introduced primarily to get BAR working.
   */
  public enum ThermodynamicsAlgorithm {
    // TODO: Separate MC-OST from MD-OST. Requires coupled changes elsewhere.
    // Fixed represents generation of snapshots for estimators like BAR, FEP, etc.
    OST("OST", "MC-OST", "MD-OST", "DEFAULT"),
    FIXED("BAR", "MBAR", "FEP", "WINDOWED");

    private final Set<String> aliases;

    ThermodynamicsAlgorithm(String... aliases) {
      // If, for some reason, there are 100+ aliases, might change to a HashSet.
      Set<String> names = new TreeSet<>(Arrays.asList(aliases));
      names.add(this.name());
      this.aliases = Collections.unmodifiableSet(names);
    }

    /**
     * Parse a String to a corresponding thermodynamics algorithm, recognizing aliases.
     *
     * @param name Name to parse
     * @return A ThermodynamicsAlgorithm.
     * @throws IllegalArgumentException If name did not correspond to any alias of any
     *     ThermodynamicsAlgorithm.
     */
    public static ThermodynamicsAlgorithm parse(String name) throws IllegalArgumentException {
      String ucName = name.toUpperCase();
      for (ThermodynamicsAlgorithm ta : values()) {
        if (ta.aliases.contains(ucName)) {
          return ta;
        }
      }
      throw new IllegalArgumentException(
          String.format(" Could not parse %s as a ThermodynamicsAlgorithm", name));
    }
  }
}
