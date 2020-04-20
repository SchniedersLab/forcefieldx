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
package ffx.potential.cli;

import static ffx.utilities.StringUtils.parseAtomRanges;
import static java.lang.String.format;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.PotentialsFunctions;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize alchemistry on at least one topology.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class AlchemicalOptions {

  /** The logger for this class. */
  public static final Logger logger = Logger.getLogger(AlchemicalOptions.class.getName());

  /** A regular expression used to parse ranges of atoms. */
  public static final Pattern rangeRegEx = Pattern.compile("([0-9]+)-?([0-9]+)?");

  /** -l or --lambda sets the initial lambda value. */
  @Option(
      names = {"-l", "--lambda"},
      paramLabel = "-1",
      description = "Initial lambda value.")
  double initialLambda = -1.0;

  /** --ac or --alchemicalAtoms Specify alchemical atoms [ALL, NONE, Range(s): 1-3,6-N]." */
  @Option(
      names = {"--ac", "--alchemicalAtoms"},
      paramLabel = "<selection>",
      defaultValue = "",
      description = "Specify alchemical atoms [ALL, NONE, Range(s): 1-3,6-N].")
  String alchemicalAtoms;

  /**
   * --uc or --unchargedAtoms Specify atoms without electrostatics [ALL, NONE, Range(s): 1-3,6-N]."
   */
  @Option(
      names = {"--uc", "--unchargedAtoms"},
      paramLabel = "<selection>",
      defaultValue = "",
      description = "Specify atoms without electrostatics [ALL, NONE, Range(s): 1-3,6-N].")
  String unchargedAtoms;

  /**
   * Sets the alchemical atoms for a MolecularAssembly.
   *
   * @param assembly Assembly to which the atoms belong.
   * @param alchemicalAtoms Alchemical atoms selection string.
   */
  public static void setAlchemicalAtoms(MolecularAssembly assembly, String alchemicalAtoms) {
    if (alchemicalAtoms == null || alchemicalAtoms.equalsIgnoreCase("")) {
      // Empty or null string -- no changes to alchemical atoms.
      logger.info(" Empty alchemical selection.");
      return;
    }

    Atom[] atoms = assembly.getAtomArray();

    // No alchemical are atoms.
    if (alchemicalAtoms.equalsIgnoreCase("NONE")) {
      for (Atom atom : atoms) {
        atom.setApplyLambda(false);
      }
      logger.info(" No atoms are alchemical \n");
      return;
    }

    // All atoms are alchemical.
    if (alchemicalAtoms.equalsIgnoreCase("ALL")) {
      for (Atom atom : atoms) {
        atom.setApplyLambda(true);
      }
      logger.info(" All atoms are alchemical \n");
      return;
    }

    // A range(s) of atoms are alchemical.
    int nAtoms = atoms.length;
    for (Atom atom : atoms) {
      atom.setApplyLambda(false);
    }
    List<Integer> alchemicalAtomRanges =
        parseAtomRanges(" Alchemical atoms", alchemicalAtoms, nAtoms);
    for (int i : alchemicalAtomRanges) {
      atoms[i].setApplyLambda(true);
    }
    logger.info(" Alchemical atoms set to: " + alchemicalAtoms + "\n");
  }

  /**
   * Sets the alchemical atoms for a MolecularAssembly.
   *
   * @param assembly Assembly to which the atoms belong.
   * @param unchargedAtoms Uncharged atoms selection string.
   */
  public static void setUnchargedAtoms(MolecularAssembly assembly, String unchargedAtoms) {
    if (unchargedAtoms == null || unchargedAtoms.equalsIgnoreCase("")) {
      // Empty or null string -- no changes to uncharged atoms.
      return;
    }

    Atom[] atoms = assembly.getAtomArray();

    // No alchemical are atoms.
    if (unchargedAtoms.equalsIgnoreCase("NONE")) {
      for (Atom atom : atoms) {
        atom.setElectrostatics(false);
      }
      logger.info(" No atoms are uncharged.\n");
      return;
    }

    // All atoms are alchemical.
    if (unchargedAtoms.equalsIgnoreCase("ALL")) {
      for (Atom atom : atoms) {
        atom.setElectrostatics(true);
      }
      logger.info(" All atoms are uncharged.\n");
      return;
    }

    // A range(s) of atoms are uncharged.
    int nAtoms = atoms.length;
    for (Atom atom : atoms) {
      atom.setElectrostatics(true);
    }
    List<Integer> unchargedAtomRanges = parseAtomRanges(" Uncharged atoms", unchargedAtoms, nAtoms);
    for (int i : unchargedAtomRanges) {
      atoms[i].setElectrostatics(false);
    }
    logger.info(" Uncharged atoms set to: " + unchargedAtoms + "\n");
  }

  /**
   * Gets the initial value of lambda.
   *
   * @return Initial lambda.
   */
  public double getInitialLambda() {
    return getInitialLambda(false);
  }

  /**
   * Gets the initial value of lambda.
   *
   * @param quiet No logging if quiet.
   * @return Initial lambda.
   */
  public double getInitialLambda(boolean quiet) {
    return getInitialLambda(1, 0, quiet);
  }

  /**
   * Gets the initial value of lambda.
   *
   * @param size The number of processes.
   * @param rank THe rank of this process.
   * @return Initial lambda.
   */
  public double getInitialLambda(int size, int rank) {
    return getInitialLambda(size, rank, false);
  }

  /**
   * Gets the initial value of lambda.
   *
   * @param size The number of processes.
   * @param rank THe rank of this process.
   * @param quiet No logging if quiet.
   * @return Initial lambda.
   */
  public double getInitialLambda(int size, int rank, boolean quiet) {
    if (initialLambda < 0.0 || initialLambda > 1.0) {
      if (rank == 0 || size < 2) {
        initialLambda = 0.0;
        if (!quiet) logger.info(format(" Setting lambda to %5.3f", initialLambda));
      } else if (rank == size - 1) {
        initialLambda = 1.0;
        if (!quiet) logger.info(format(" Setting lambda to %5.3f", initialLambda));
      } else {
        double dL = 1.0 / (size - 1);
        initialLambda = dL * rank;
        if (!quiet) logger.info(format(" Setting lambda to %5.3f.", initialLambda));
      }
    }
    return initialLambda;
  }

  /**
   * If any softcore Atoms have been detected.
   *
   * @return Presence of softcore Atoms.
   */
  public boolean hasSoftcore() {
    return (alchemicalAtoms != null
        && !alchemicalAtoms.equalsIgnoreCase("NONE")
        && !alchemicalAtoms.equalsIgnoreCase(""));
  }

  /**
   * Opens a File to a MolecularAssembly for alchemistry.
   *
   * @param potentialFunctions A utility object for opening Files into MolecularAssemblies.
   * @param topologyOptions TopologyOptions in case a dual-topology or greater is to be used.
   * @param threadsPer Number of threads to be used for this MolecularAssembly.
   * @param toOpen The name of the File to be opened.
   * @param topNum The index of this topology.
   * @return The processed MolecularAssembly.
   */
  public MolecularAssembly openFile(
      PotentialsFunctions potentialFunctions,
      TopologyOptions topologyOptions,
      int threadsPer,
      String toOpen,
      int topNum) {
    potentialFunctions.openAll(toOpen, threadsPer);
    MolecularAssembly molecularAssembly = potentialFunctions.getActiveAssembly();
    return processFile(topologyOptions, molecularAssembly, topNum);
  }

  /**
   * Performs processing on a MolecularAssembly for alchemistry.
   *
   * @param topologyOptions TopologyOptions in case a dual-topology or greater is to be used.
   * @param molecularAssembly The MolecularAssembly to be processed.
   * @param topNum The index of this topology, 0-indexed.
   * @return The processed MolecularAssembly.
   */
  public MolecularAssembly processFile(
      TopologyOptions topologyOptions, MolecularAssembly molecularAssembly, int topNum) {

    int remainder = (topNum % 2) + 1;
    switch (remainder) {
      case 1:
        setFirstSystemAlchemistry(molecularAssembly);
        setFirstSystemUnchargedAtoms(molecularAssembly);
        break;
      case 2:
        if (topologyOptions == null) {
          throw new IllegalArgumentException(
              " For >= 2 systems, topologyOptions must not be empty!");
        }
        topologyOptions.setSecondSystemAlchemistry(molecularAssembly);
        topologyOptions.setSecondSystemUnchargedAtoms(molecularAssembly);
        break;
    }

    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    ForceFieldEnergy energy = molecularAssembly.getPotentialEnergy();
    energy.getCrystal().setSpecialPositionCutoff(0.0);

    return molecularAssembly;
  }

  /**
   * Set the alchemical atoms for this molecularAssembly.
   *
   * @param topology a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setFirstSystemAlchemistry(MolecularAssembly topology) {
    setAlchemicalAtoms(topology, alchemicalAtoms);
  }

  /**
   * Set uncharged atoms for this molecularAssembly.
   *
   * @param topology a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setFirstSystemUnchargedAtoms(MolecularAssembly topology) {
    setUnchargedAtoms(topology, unchargedAtoms);
  }
}
