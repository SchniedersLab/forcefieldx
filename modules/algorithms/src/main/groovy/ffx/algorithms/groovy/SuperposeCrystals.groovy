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
import ffx.potential.bonded.Atom
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.ProgressiveAlignmentOfCrystals
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.concat
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.io.FilenameUtils.getFullPath

/**
 * Quantifies the similarity of input crystals based on progressive alignments.
 * This script is based off of the PACCOM code created by Okimasa Okada.
 *
 * @author Okimasa OKADA
 * @author Aaron J. Nessler and Michael J. Schnieders
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler and Micheal Schnieders 2020
 * revised by Aaron Nessler and Michael Schnieders 2021
 * <br>
 * Usage:
 * <br>
 * ffxc test.SuperposeCrystals &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Compare crystal polymorphs based on a RMSD of aligned crystal coordinates.", name = "ffxc SuperposeCrystals")
class SuperposeCrystals extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * --na or --numAU Number of asymmetric units to include from each crystal in RMSD comparison.
   */
  @Option(names = ['--na', '--numAU'], paramLabel = '20', defaultValue = '20',
      description = 'Set the number of asymmetric units to include in final PAC RMSD.')
  int numAU

  /**
   * --ni or --numInflatedAU Number of asymmetric units in the inflated sphere.
   */
  @Option(names = ['--ni', '--numInflatedAU'], paramLabel = '500', defaultValue = '500',
      description = 'Specifies the number asymmetric units in the expanded crystal.')
  int numInflatedAU

  /**
   * --ns or --numSearch Number of molecules for each handedness in the first crystal.
   */
  @Option(names = ['--ns', '--numSearch'], paramLabel = '3', defaultValue = '3',
      description = 'Set the number of asymmetric units to search in the 1st crystal to check for mirrored conformations.')
  int numSearch

  /**
   * --ns2 or --numSearch2 Number of molecules for each handedness in the second crystal.
   */
  @Option(names = ['--ns2', '--numSearch2'], paramLabel = '3', defaultValue = '3',
      description = 'Set the number of asymmetric units to search in the 2nd crystal to check for mirrored conformations.')
  int numSearch2

  /**
   * --zp or --zPrime Overrides number of molecules in the asymmetric unit.
   */
  @Option(names = ['--zp', '--zPrime'], paramLabel = '-1', defaultValue = '-1',
          description = 'Number of species in asymmetric unit should be detected by default (Z\'). This flag should only be used if the default detection fails.')
  int zPrime

  /**
   * --nms or --noMirrorSearch Do not loop over asymmetric units to check for mirrored conformations.
   */
  @Option(names = ['--nms', '--noMirrorSearch'], paramLabel = "false", defaultValue = "false",
      description = 'Do not loop over asymmetric units to check for mirrored conformations.')
  private static boolean noMirrorSearch

  /**
   * -w or --write Write out the PAC RMSD matrix.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write out the PAC RMSD matrix.')
  private static boolean write

  /**
   * -r or --restart Attempt to restart from a previously written PAC RMSD matrix.
   */
  @Option(names = ['-r', '--restart'], paramLabel = "false", defaultValue = "false",
      description = 'Attempt to restart from a previously written PAC RMSD matrix.')
  private static boolean restart

  /**
   * --sp or --savePDB Save out a PDB.
   */
  @Option(names = ['--sp', '--savePDB'], paramLabel = "false", defaultValue = "false",
      description = 'Save a PDB file for the superposed crystal.')
  private static boolean savePDB

  /**
   * -f or --full Perform the full comparison (takes considerably longer, but more information returned).
   */
  @Option(names = ['-f', '--full'], paramLabel = "false", defaultValue = "false",
      description = 'Use the full number of comparisons (tradeoff between speed and information returned).')
  private static boolean full

  /**
   * --ac or --alphaCarbons PAC RMSD will only include alpha carbons.
   */
  @Option(names = ['--ac', '--alphaCarbons'], paramLabel = "false", defaultValue = "false",
          description = 'PAC RMSD will only include alpha carbons.')
  private static boolean alphaCarbons

  /**
   * --nh or --noHydrogen PAC RMSD will not include hydrogen atoms.
   */
  @Option(names = ['--nh', '--noHydrogen'], paramLabel = "false", defaultValue = "false",
      description = 'PAC RMSD will not include hydrogen atoms.')
  private static boolean noHydrogen

  // TODO finish implementing symmetric flag.
  /**
   * --sym or --symmetric Enforce generation of a symmetric PAC RMSD matrix.
   */
  @Option(names = ['--sym', '--symmetric'], paramLabel = "false", defaultValue = "false",
      description = 'Enforce generation of a symmetric PAC RMSD matrix.')
  private static boolean symmetric

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

    // Atom array from the 1st assembly.
    Atom[] baseAtoms = activeAssembly.getAtomArray()
    int nAtoms = baseAtoms.size()

    // Collect selected atoms.
    ArrayList<Integer> atomList = new ArrayList<>()
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = baseAtoms[i]
      if (atom.isActive()) {
        String atomName = atom.getName()
        int atomAtNum = atom.getAtomicNumber()
        boolean proteinCheck = atomName == "CA" && atomAtNum == 6
        boolean aminoCheck = (atomName == "N1" || atomName == "N9") && atomAtNum == 7
        if(alphaCarbons){
          if(proteinCheck||aminoCheck){
            atomList.add(i)
          }
        }else if (!noHydrogen || !atom.isHydrogen()) {
          atomList.add(i)
        }
      }
      // Reset all atoms to active once the selection is recorded.
      atom.setActive(true)
    }

    if (atomList.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD.")
      return this
    }

    // Number of atoms included in the PAC RMSD.
    int nPACAtoms = atomList.size()
    logger.info(format("\n %d atoms will be used for the PAC RMSD out of %d.\n", nPACAtoms, nAtoms))

    // If search is desired ensure inner loop will be used. Else use single comparison.
    if (noMirrorSearch) {
      numSearch = 1
      numSearch2 = 1
    }

    // Compare structures in baseFilter and targetFilter.
    ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter, isSymmetric)

    // Define the filename to use for the PAC RMSD values.
    String filename = filenames.get(0)
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".txt")

    // To save in ARES format a PDB must be written out.
    if(ares){
      savePDB = true
    }

    runningStatistics = pac.comparisons(atomList, numAU, numInflatedAU,
        numSearch, numSearch2, zPrime, full, savePDB, restart, write, ares, pacFilename)

    return this
  }
}