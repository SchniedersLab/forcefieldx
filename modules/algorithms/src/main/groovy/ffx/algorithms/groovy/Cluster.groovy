//******************************************************************************
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
//******************************************************************************

package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.Crystal
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.utils.Clustering.Conformation
import org.apache.commons.math3.ml.clustering.CentroidCluster
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.potential.parsers.DistanceMatrixFilter.readDistanceMatrix
import static ffx.potential.utils.Clustering.*
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.floorDiv

/**
 * The Cluster script clusters structures utilizing RMSD.
 *
 * @author Aaron J. Nessler
 * @author Mallory R. Tollefson
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc Cluster [options] &lt;filename&gt;
 */
@Command(description = " Cluster structures using an RMSD matrix.", name = "ffxc Cluster")
class Cluster extends AlgorithmsScript {

  /**
   * -a or --algorithm Clustering algorithm to use.
   * Algorithm: Multi-K-Means++ (0) or Hierarchical (1)
   */
  @Option(names = ['-a', '--algorithm'], paramLabel = "0", defaultValue = "0",
      description = "Algorithm: Multi-K-Means++ (0) or Hierarchical (1)")
  private int algorithm

  /**
   * -t or --trials Number of trials for Multi-K-Means.
   */
  @Option(names = ['-t', '--trials'], paramLabel = "1", defaultValue = "1",
      description = "Number of trials for Multi-K-Means.")
  private int trials

  /**
   * -k or --clusters Maximum number of kmeans clusters.
   */
  @Option(names = ['-k', '--clusters'], paramLabel = "0", defaultValue = "0",
      description = "Maximum number of kmeans clusters.")
  private int maxClusters

  /**
   * --rs or --randomSeed Set the random seed for deterministic clustering.
   */
  @Option(names = ['--rs', '--randomSeed'], paramLabel = "-1", defaultValue = "-1",
      description = 'Set the random seed for deterministic clustering (-1 uses the current time).')
  private long randomSeed

  /**
   * --td or --threshold
   * RMSD value for dividing clusters from a hierarchical tree.
   */
  @Option(names = ['--td', '--threshold'], paramLabel = "2.0", defaultValue = "2.0",
      description = "The distance value where a hierarchical tree flattened into clusters.")
  private double threshold

  /**
   * -w or --write Write out an archive of a representative structure from each cluster.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write an archive containing a representative from each cluster.')
  private boolean write

  /**
   * The final argument(s) should be one or two filenames.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = 'The RMSD distance matrix and optionally the corresponding structure archive (required for --write).')
  private List<String> filenames = null

  /**
   * Generated list of clusters.
   */
  private List<CentroidCluster<Conformation>> clusterList

  /**
   * Cluster constructor.
   */
  Cluster() {
    this(new Binding())
  }

  /**
   * Cluster constructor.
   * @param binding The Groovy Binding to use.
   */
  Cluster(Binding binding) {
    super(binding)
  }

  /**
   * Return the Clusters.
   * @return Returns the generated clusters.
   */
  List<CentroidCluster<Conformation>> getClusterList() {
    return clusterList
  }

  /**
   * Execute the script.
   */
  @Override
  Cluster run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    if (filenames == null || filenames.isEmpty()) {
      logger.info(helpString())
      return this
    }

    List<double[]> distMatrix = new ArrayList<double[]>()

    String filename = filenames.get(0)
    if (!readDistanceMatrix(filename, distMatrix)) {
      logger.info(format(" The distance matrix %s could not be read in.", filename))
      return this
    }

    // Either use Multi-K-Means++ or Hierarchical agglomerative clustering.
    switch (algorithm) {
      case 1:
        // Hierarchical clustering.
        logger.info(" Performing Hierarchical Clustering")
        logger.info(format("  Cluster separation threshold: %6.4f A.\n", threshold))
        clusterList = hierarchicalClustering(distMatrix, threshold)
        break
      default:
        // K-Means++ and Multi K-Means++.
        logger.info(" Performing K-Means++ Clustering")

        // Ensure cluster are within bounds
        int dim = distMatrix.size()
        if (maxClusters <= 0 || maxClusters >= dim) {
          int newSize = floorDiv(dim, 2)
          maxClusters = newSize
        }
        logger.info(format("  Number of clusters: %d", maxClusters))

        // Array representing indices of structures that represent their cluster.
        if (randomSeed == -1) {
          randomSeed = System.nanoTime()
        } else {
          logger.info(format("  Random seed:        %d", randomSeed))
        }
        logger.info(format("  Number of trials:   %d\n", trials))
        clusterList = kMeansClustering(distMatrix, maxClusters, trials, randomSeed)
    }

    boolean verbose = true
    List<Integer> representativeConformer = new ArrayList<>()
    analyzeClusters(clusterList, representativeConformer, verbose)

    if (write) {
      writeStructures(representativeConformer)
    }

    return this
  }

  /**
   * Write out structures corresponding to the representative of each cluster.
   * @param distMatrix Distances between structures (metric to determine clusters).
   * @param repStructs Array list for index of representative structures.
   */
  private void writeStructures(List<Integer> repStructs) {
    if (filenames.size() < 2) {
      logger.info(" Please supply an ARC or PDB file corresponding to the distance matrix.")
      return
    }

    // Turn off nonbonded terms.
    System.setProperty("vdwterm", "false")

    String saveName = filenames.get(1)
    File saveFile = new File(algorithmFunctions.versionFile(saveName))
    logger.info(format("\n Saving representative cluster members to %s.", saveFile.toString()))

    MolecularAssembly[] molecularAssemblies = algorithmFunctions.openAll(saveName)
    activeAssembly = molecularAssemblies[0]

    SystemFilter systemFilter = algorithmFunctions.getFilter()
    PDBFilter pdbFilter = null
    XYZFilter xyzFilter = null
    if (systemFilter instanceof PDBFilter) {
      pdbFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
    } else if (systemFilter instanceof XYZFilter) {
      xyzFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
    } else {
      logger.info(" Unknown file type.")
      return
    }

    logger.info("\n Conformations")
    int structNum = 0
    do {
      if (repStructs.contains(structNum)) {
        int clusterNum = repStructs.indexOf(structNum)
        // Make sure the crystal is is updated.
        Crystal crystal = activeAssembly.getCrystal()
        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
        forceFieldEnergy.setCrystal(crystal)
        if (crystal.aperiodic()) {
          logger.info(format(" Conformation %d, Cluster %d)", structNum + 1, clusterNum + 1))
        } else {
          logger.info(format(" Conformation %d, Cluster %d (%s)", structNum + 1, clusterNum + 1,
              crystal.getUnitCell().toShortString()))
        }
        if (systemFilter instanceof PDBFilter) {
          pdbFilter.writeFile(saveFile, true, false, false)
        } else if (systemFilter instanceof XYZFilter) {
          String[] message = new String[1]
          message[0] = format(" Conformation %d, Cluster %d ", structNum + 1, clusterNum + 1)
          xyzFilter.writeFile(saveFile, true, message)
        }
      }
      structNum++
    } while (systemFilter.readNext(false, false))

    if (systemFilter instanceof PDBFilter) {
      saveFile.append("END\n")
    }

    systemFilter.closeReader()
  }
}


