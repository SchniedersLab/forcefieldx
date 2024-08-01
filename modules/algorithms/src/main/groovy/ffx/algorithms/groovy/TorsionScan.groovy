//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.TorsionSearch
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FileUtils
import org.apache.commons.io.FilenameUtils
import org.apache.commons.math3.util.FastMath
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The TorsionSearch command enumerates conformations of a molecule using torsional scans around rotatable bonds.
 *
 * @author Aaron J. Nessler
 * @author Matthew J. Speranza
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc TorsionSearch &lt;filename&gt;
 */
@Command(description = " The TorsionSearch command enumerates conformations of a molecule using torsional scans around rotatable bonds.",
    name = "TorsionSearch")
class TorsionScan extends AlgorithmsScript {

  /**
   * --th or --theta Step size for bond rotations.
   */
  @Option(names = ['--th', '--theta'], paramLabel = '60.0', defaultValue = '60.0',
      description = "Step size for bond rotations (in Degrees).")
  private double increment = 60.0

  /**
   * --saveNumStates
   */
  @Option(names = ['--saveNumStates', "--sns"], paramLabel = '10', defaultValue = '10',
      description = 'Save this many of the lowest energy states per worker. This is the default.')
  private int saveNumStates = 10

  /**
   * --elimMax
   */
  @Option(names = ['--elimMax', "--em"], paramLabel = '0', defaultValue = '0',
      description = 'Eliminate bonds where one torsion causes high energies. Reduces the complexity of the search.')
  private int elimMax = 0

  /**
   * --saveAll
   */
  @Option(names = ['--saveAll'], paramLabel = 'false', defaultValue = 'false',
      description = 'Save out all states. Not recommended for large systems.')
  private boolean saveAll = false

  /**
   * --sc or --staticComparison Hold angles fixed.
   */
  @Option(names = ['--sc', '--staticComparison'], paramLabel = "false", defaultValue = "false",
      description = 'If set, each bond is rotated independently (faster, but fewer permutations).')
  private static boolean staticCompare

  /**
   * --eliminationThreshold
   */
  @Option(names = ['--eliminationThreshold', '--et'], paramLabel = '50.0', defaultValue = '50.0',
      description = "Remove bonds that cause > this energy change during static analysis (kcal/mol).")
  private double eliminationThreshold = 50.0

  /**
   * --startIndex
   */
  @Option(names = ['--startIndex'], paramLabel = '0', defaultValue = '0',
      description = "Start at this hilbert index.")
  private long startIndex = 0

  /**
   * --endIndex
   */
  @Option(names = ['--endIndex'], paramLabel = '0', defaultValue = '0',
      description = "End at this hilbert index.")
  private long endIndex = 0

  /**
   * The final argument should be the structure filename.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = 'Atomic coordinate file(s) to permute in XYZ format.')
  String filename = null

  /**
   * CrystalSuperpose Constructor.
   */
  TorsionScan() {
    this(new Binding())
  }

  /**
   * CrystalSuperpose Constructor.
   * @param binding Groovy Binding to use.
   */
  TorsionScan(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  TorsionScan run() {
    // Direct is sufficient to scan for atom clashes.
    System.setProperty("polarization", "direct")

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Ensure file exists.
    if (filename == null) {
      logger.info(helpString())
      return this
    }

    // Setup
    activeAssembly = getActiveAssembly(filename)
    SystemFilter systemFilter = algorithmFunctions.getFilter()
    MolecularAssembly ma = systemFilter.getActiveMolecularSystem()
    logger.info("")
    ForceFieldEnergy forceFieldEnergy = activeAssembly.potentialEnergy
    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    Comm world = Comm.world()
    int rank = world.rank()
    int size = world.size()
    if (saveAll) {
      saveNumStates = -1
    }

    // Run the torsion search.
    int numTorsions = FastMath.ceil(360 / increment) as int
    TorsionSearch torsionSearch = new TorsionSearch(ma, ma.getMoleculeArray()[0], numTorsions, saveNumStates)
    // Don't repeat static analysis if we're running multiple workers.
    if (size > 1 && rank == 0) {
      torsionSearch.staticAnalysis(elimMax, eliminationThreshold)
    } else if (size == 1) {
      torsionSearch.staticAnalysis(elimMax, eliminationThreshold)
    }

    if (world.size() > 1 && !staticCompare) {
      torsionSearch.buildWorker(rank, size)
      torsionSearch.runWorker()
    } else if (!staticCompare) {
      endIndex = endIndex == 0 ? torsionSearch.getEnd() : endIndex
      torsionSearch.spinTorsions(startIndex, endIndex)
    }

    // Get states, energies, and hilbert indices.
    List<AssemblyState> states = torsionSearch.getStates()
    List<Double> energies = torsionSearch.getEnergies()
    List<Long> hilbertIndices = torsionSearch.getHilbertIndices()

    // Save out states.
    int count = 0
    logger.info("\n Saving " + states.size() + " states.")
    String extension = "_rot.arc"
    if (world.size() > 1) {
      extension = "_rank" + world.rank() + extension
    }
    File saveLocation = new File(FilenameUtils.removeExtension(filename) + extension)
    logger.info(" Logging structures into: " + saveLocation)
    XYZFilter xyzFilter = new XYZFilter(saveLocation,
        activeAssembly,
        activeAssembly.getForceField(),
        activeAssembly.properties)
    while (!states.empty) {
      AssemblyState assembly = states.get(0)
      states.remove(0)
      assembly.revertState()
      forceFieldEnergy.getCoordinates(x)
      double e = energies.get(0)
      energies.remove(0)
      long hilbertIndex = hilbertIndices.get(0)
      hilbertIndices.remove(0)
      logger.info(format(" Writing to file. Configuration #%-6d energy: %-12.5f Hilbert index: %-15d",
          (count + 1),
          e,
          hilbertIndex))
      xyzFilter.writeFile(saveLocation, true)
      count++
    }
    logger.info("\n -1 indices come from static torsion scan.")

    // Create properties/key file
    File key = new File(FilenameUtils.removeExtension(filename) + ".key")
    File properties = new File(FilenameUtils.removeExtension(filename) + ".properties")
    try {
      if (key.exists()) {
        File keyComparison = new File(FilenameUtils.removeExtension(filename) + "_rot.key")
        if (keyComparison.createNewFile()) {
          FileUtils.copyFile(key, keyComparison)
        }
      } else if (properties.exists()) {
        File propertiesComparison = new File(FilenameUtils.removeExtension(filename) + "_rot.properties")
        if (propertiesComparison.createNewFile()) {
          FileUtils.copyFile(properties, propertiesComparison)
        }
      } else {
        logger.info(" No key or properties file found.")
      }
    } catch (Exception e) {
      e.printStackTrace()
    }

    if (world.rank() == 0 && world.size() > 1) {
      // Combine all the rank's files into filename_rot.arc
      logger.info("\n Combining all rank files into one file.")
      File combined = new File(FilenameUtils.removeExtension(filename) + "_rot.arc")
      try {
        if (combined.createNewFile()) {
          FileOutputStream fos = new FileOutputStream(combined)
          for (int i = 0; i < world.size(); i++) {
            File rankFile = new File(FilenameUtils.removeExtension(filename) + "_rank" + i + "_rot.arc")
            if (rankFile.exists()) {
              FileInputStream fis = new FileInputStream(rankFile)
              byte[] buffer = new byte[1024]
              int length
              while ((length = fis.read(buffer)) > 0) {
                fos.write(buffer, 0, length)
              }
              fis.close()
            }
          }
          fos.close()
        }
      } catch (Exception e) {
        e.printStackTrace()
      }
    }
    return this
  }
}