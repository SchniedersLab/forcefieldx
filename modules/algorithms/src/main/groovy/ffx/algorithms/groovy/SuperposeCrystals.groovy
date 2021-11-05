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
   * --na or --numAU Number of asymmetric units to include from each crystal in RMSD comparison.
   */
  @Option(names = ['--na', '--numAU'], paramLabel = '20', defaultValue = '20',
      description = 'Set the number of asymmetric units to include in final RMSD.')
  int numAU

  /**
   * --ni or --numInflatedAU Number of asymmetric units in the inflated sphere.
   */
  @Option(names = ['--ni', '--numInflatedAU'], paramLabel = '500', defaultValue = '500',
      description = 'Specifies the number asymmetric units in the expanded crystal.')
  int numInflatedAU

  /**
   * --ns or --numSearch Number of asymmetric units for each handedness in the first crystal.
   */
  @Option(names = ['--ns', '--numSearch'], paramLabel = '1', defaultValue = '1',
      description = 'Set the number of asymmetric units to search in the 1st crystal to check for additional conformations.')
  int numSearch

  /**
   * --ns2 or --numSearch2 Number asymmetric units for each handedness in the second crystal.
   */
  @Option(names = ['--ns2', '--numSearch2'], paramLabel = '1', defaultValue = '1',
      description = 'Set the number of asymmetric units to search in the 2nd crystal to check for additional conformations.')
  int numSearch2

  /**
   * --zp or --zPrime Overrides number of species in the asymmetric unit.
   */
  @Option(names = ['--zp', '--zPrime'], paramLabel = '-1', defaultValue = '-1',
          description = 'Number of species in asymmetric unit of first crystal.')
  int zPrime

  /**
   * --zp2 or --zPrime2 Overrides number of species in the asymmetric unit.
   */
  @Option(names = ['--zp2', '--zPrime2'], paramLabel = '-1', defaultValue = '-1',
          description = 'Number of species in asymmetric unit of second crystal.')
  int zPrime2

  /**
   * -w or --write Write out the RMSD matrix.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write out the RMSD matrix.')
  private static boolean write

  /**
   * -r or --restart Attempt to restart from a previously written RMSD matrix.
   */
  @Option(names = ['-r', '--restart'], paramLabel = "false", defaultValue = "false",
      description = 'Attempt to restart from a previously written RMSD matrix.')
  private static boolean restart

  /**
   * --sp or --savePDB Save out a PDB.
   */
  @Option(names = ['--sp', '--savePDB'], paramLabel = "false", defaultValue = "false",
      description = 'Save a PDB file for the superposed crystal.')
  private static boolean savePDB

  /**
   * --ex or --exhaustive Perform an exhaustive comparison to handle multiple conformations (more expensive, but may find lower RMSD).
   */
  @Option(names = ['--ex', '--exhaustive'], paramLabel = "false", defaultValue = "false",
      description = 'Perform an exhaustive comparison to handle multiple conformations (more expensive, but may find lower RMSD).')
  private static boolean exhaustive

  /**
   * --ac or --alphaCarbons Protein RMSD will only include alpha carbons.
   */
  @Option(names = ['--ac', '--alphaCarbons'], paramLabel = "false", defaultValue = "false",
          description = 'Protein RMSD will only include alpha carbons for comparison.')
  private static boolean alphaCarbons

  /**
   * --nh or --noHydrogen RMSD will not include hydrogen atoms.
   */
  @Option(names = ['--nh', '--noHydrogen'], paramLabel = "false", defaultValue = "false",
      description = 'RMSD will not include hydrogen atoms.')
  private static boolean noHydrogen

  /**
   * --sa or --saveAres Save out PDB in ARES input format.
   */
  @Option(names = ['--sa', '--saveAres'], paramLabel = "false", defaultValue = "false",
          description = 'Final structures for each comparison will be written out in ARES input format.')
  private static boolean ares

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

    // Compare structures in baseFilter and targetFilter.
    ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter, isSymmetric)

    // Define the filename to use for the RMSD values.
    String filename = filenames.get(0)
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".txt")

    // To save in ARES format a PDB must be written out.
    if(ares){
      savePDB = true
    }

    runningStatistics = pac.comparisons(numAU, numInflatedAU, numSearch, numSearch2, zPrime, zPrime2, alphaCarbons,
        noHydrogen, exhaustive, savePDB, restart, write, ares, pacFilename)

    return this
  }
}