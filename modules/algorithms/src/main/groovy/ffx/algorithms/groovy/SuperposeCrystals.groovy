//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.math.RunningStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.ProgressiveAlignmentOfCrystals
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FilenameUtils.concat
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.io.FilenameUtils.getFullPath

/**
 * Quantifies the similarity of input crystals based on progressive alignments.
 * This script is based off of the PACCOM code created by Okimasa Okada.
 *
 * @author Okimasa OKADA
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * @author Aaron J. Nessler and Michael J. Schnieders
 * ported to FFX by Aaron Nessler and Micheal Schnieders 2020
 * revised by Aaron Nessler and Michael Schnieders 2021
 * <br>
 * Usage:
 * <br>
 * ffxc test.SuperposeCrystals &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Determine the RMSD for crystal polymorphs using the Progressive Alignment of Crystals (PAC) algorithm.",
    name = "ffxc SuperposeCrystals")
class SuperposeCrystals extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * --na or --numAU AUs in the RMSD.
   */
  @Option(names = ['--na', '--numAU'], paramLabel = '20', defaultValue = '20',
      description = 'AUs in the RMSD.')
  private int numAU

  /**
   * --ni or --numInflatedAU AUs in the expanded crystal.
   */
  @Option(names = ['--ni', '--numInflatedAU'], paramLabel = '500', defaultValue = '500',
      description = 'AUs in the expanded crystal.')
  private int numInflatedAU

  /**
   * --ns or --numSearch Crystal 1 AUs to compare for unique conformations.
   */
  @Option(names = ['--ns', '--numSearch'], paramLabel = '1', defaultValue = '1',
      description = 'Crystal 1 AUs to compare for unique conformations.')
  private int numSearch

  /**
   * --ns2 or --numSearch2 Crystal 2 AUs to compare for unique conformations.
   */
  @Option(names = ['--ns2', '--numSearch2'], paramLabel = '1', defaultValue = '1',
      description = 'Crystal 2 AUs to compare for unique conformations.')
  private int numSearch2

  /**
   * --zp or --zPrime Z' for crystal 1 (-1 to autodetect).
   */
  @Option(names = ['--zp', '--zPrime'], paramLabel = '-1', defaultValue = '-1',
      description = "Z'' for crystal 1 (-1 to autodetect).")
  private int zPrime

  /**
   * --zp2 or --zPrime2 Z' for crystal 2 (-1 to autodetect).
   */
  @Option(names = ['--zp2', '--zPrime2'], paramLabel = '-1', defaultValue = '-1',
      description = "Z'' for crystal 2 (-1 to autodetect).")
  private int zPrime2

  /**
   * -w or --write Write out the PAC RMSD matrix.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write out the PAC RMSD matrix.')
  private static boolean write

  /**
   * -r or --restart Restart from a previously written RMSD matrix (if one exists).
   */
  @Option(names = ['-r', '--restart'], paramLabel = "false", defaultValue = "false",
      description = 'Restart from a previously written RMSD matrix (if one exists).')
  private static boolean restart

  /**
   * --sp or --savePDB Save PDB files for the superposed crystals.
   */
  @Option(names = ['--sp', '--savePDB'], paramLabel = "false", defaultValue = "false",
      description = 'Save PDB files for the superposed crystals.')
  private static boolean savePDB

  /**
   * --ex or --exhaustive Perform an exhaustive comparison to handle multiple conformations (more expensive, but may find lower RMSD).
   */
  @Option(names = ['--ex', '--exhaustive'], paramLabel = "false", defaultValue = "false",
      description = 'Perform an exhaustive comparison to handle multiple conformations (more expensive, but may find lower RMSD).')
  private static boolean exhaustive

  /**
   * --ac or --alphaCarbons Consider only alpha carbons for proteins.
   */
  @Option(names = ['--ac', '--alphaCarbons'], paramLabel = "false", defaultValue = "false",
      description = 'Consider only alpha carbons for proteins.')
  private static boolean alphaCarbons

  /**
   * --nh or --noHydrogen Ignore hydrogen atoms.
   */
  @Option(names = ['--nh', '--noHydrogen'], paramLabel = "false", defaultValue = "false",
      description = 'Ignore hydrogen atoms.')
  private static boolean noHydrogen

  /**
   * --sm or --saveMachineLearning Save out PDB and CSV for machine learning.
   */
  @Option(names = ['--sm', '--saveMachineLearning'], paramLabel = "false", defaultValue = "false",
      description = 'Final structures for each comparison will be written out with RMSD in a CSV.')
  private static boolean machineLearning

  /**
   * --mw or --massWeighted Use mass-weighted atomic coordinates for alignment.
   */
  @Option(names = ['--mw', '--massWeighted'], paramLabel = "false", defaultValue = "false",
      description = 'Use mass-weighted atomic coordinates for alignment.')
  private static boolean massWeighted

  /**
   * -l or --linkage Single (0), Average (1), or Complete (2) coordinate linkage for molecule prioritization.
   */
  @Option(names = ['-l', '--linkage'], paramLabel = '0', defaultValue = '0',
      description = 'Single (0), Average (1), or Complete (2) coordinate linkage for molecule prioritization.')
  private int linkage

  /**
   * --fo or --fileOrder Prioritize matching molecules of the first command line crystal (rather than using a density criteria).
   */
  @Option(names = ['--fo', '--fileOrder'], paramLabel = "false", defaultValue = "false",
      description = 'Prioritize matching molecules of the first command line crystal (rather than using a density criteria).')
  private static boolean fileOrder

  /**
   * --ld or --lowDensity Prioritize matching molecules of the lower density crystal (default uses higher density).
   */
  @Option(names = ['--ld', '--lowDensity'], paramLabel = "false", defaultValue = "false",
      description = 'Prioritize matching molecules of the lower density crystal (default uses higher density).')
  private static boolean lowDensity

  /**
   * --pc or --prioritizeCrystals Prioritize the crystals being compared based on high density (0), low density (1), or file order (2).
   */
  @Option(names = ['--pc', '--prioritizeCrystals'], paramLabel = '0', defaultValue = '0',
          description = 'High density (0), low density (1), or file order (2) prioritization of submitted crystals.')
  private int crystalPriority

  /**
   * The final argument(s) should be two or more filenames (same file twice if comparing same structures).
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = 'Atomic coordinate file(s) to compare in XYZ format.')
  List<String> filenames = null

  /**
   * CrystalSuperpose Test requires a public variable containing observables to test.
   */
  public RunningStatistics runningStatistics

  /**
   * CrystalSuperpose Constructor.
   */
  SuperposeCrystals() {
    this(new Binding())
  }

  /**
   * CrystalSuperpose Constructor.
   * @param binding Groovy Binding to use.
   */
  SuperposeCrystals(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SuperposeCrystals run() {

    // Turn off non-bonded interactions for speed.
    System.setProperty("vdwterm", "false")

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Ensure file exists.
    if (filenames == null) {
      logger.info(helpString())
      return this
    }

    // SystemFilter containing structures stored in file 0.
    SystemFilter baseFilter
    // SystemFilter containing structures stored in file 1 (or file 0 if file 1 does not exist).
    SystemFilter targetFilter

    algorithmFunctions.openAll(filenames.get(0))
    baseFilter = algorithmFunctions.getFilter()
    // Example atoms to determine single molecule characteristics (e.g. number of atoms, hydrogen, etc.)
    MolecularAssembly activeAssembly = baseFilter.getActiveMolecularSystem()

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    // Number of files to read in.
    boolean isSymmetric = false
    int numFiles = filenames.size()
    if (numFiles == 1) {
      logger.info(
          "\n PAC will be applied between all pairs of conformations within the supplied file.\n")
      isSymmetric = true
      // If only one file is supplied, compare all structures in that file to each other.
      algorithmFunctions.openAll(filenames.get(0))
      targetFilter = algorithmFunctions.getFilter()
    } else {
      // Otherwise, compare structures from first file those in the second.
      logger.info(
          "\n PAC will compare all conformations in the first file to all those in the second file.\n")
      algorithmFunctions.openAll(filenames.get(1))
      targetFilter = algorithmFunctions.getFilter()
    }
    atomSelectionOptions.setActiveAtoms(targetFilter.getActiveMolecularSystem())

    // Compare structures in baseFilter and targetFilter.
    ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
        isSymmetric)

    // Define the filename to use for the RMSD values.
    String filename = filenames.get(0)
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".txt")

    // To save in ARES format a PDB must be written out.
    if (machineLearning) {
      savePDB = true
    }

    if (linkage == 0) {
      logger.finer(" Single linkage will be used.")
    } else if (linkage == 2) {
      logger.finer(" Complete linkage will be used.")
    } else if (linkage == 1) {
      logger.finer(" Average linkage will be used.")
    } else {
      logger.warning(
          "Prioritization method specified incorrectly (--pm {0, 1, 2}). Using default of average linkage.")
      linkage = 1
    }

    runningStatistics =
        pac.comparisons(numAU, numInflatedAU, numSearch, numSearch2, zPrime, zPrime2, alphaCarbons,
            noHydrogen, massWeighted, crystalPriority, exhaustive, savePDB,
            restart, write,
            machineLearning, linkage, pacFilename)

    return this
  }
}