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

package ffx.potential.groovy

import com.apporiented.algorithm.clustering.ClusteringAlgorithm
import com.apporiented.algorithm.clustering.CompleteLinkageStrategy
import com.apporiented.algorithm.clustering.DefaultClusteringAlgorithm
import com.apporiented.algorithm.clustering.visualization.DendrogramPanel
import ffx.crystal.Crystal
import ffx.potential.ForceFieldEnergy
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.swing.*
import java.awt.*
import java.util.List
import java.util.logging.Level

import static ffx.potential.utils.Cluster.kMeansCluster
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.floorDiv

/**
 * The Cluster script clusters structures utilizing RMSD.
 * TODO: Create a unit test for the Cluster script.
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
class Cluster extends PotentialScript {

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
   * -r or --readInDistMat Read in an NxN distance matrix.
   */
  @Option(names = ['-r', '--readInDistMat'], paramLabel = "false", defaultValue = "false",
      description = "Read in an NxN distance matrix.")
  private boolean readIn

  /**
   * -s or --start Atom number where RMSD calculation of structure will begin.
   */
  @Option(names = ['-s', '--start'], paramLabel = "1", defaultValue = "1",
      description = 'Starting atom to include in the RMSD calculation.')
  private String start

  /**
   * -f or --final Atom number where RMSD calculation of structure will end.
   */
  @Option(names = ['-f', '--final'], paramLabel = "nAtoms",
      description = 'Final atom to include in the RMSD calculation.')
  private String finish = Integer.MAX_VALUE.toString()

  /**
   * --rs or --randomSeed Set the random seed for deterministic clustering.
   */
  @Option(names = ['--rs', '--randomSeed'], paramLabel = "-1", defaultValue = "-1",
      description = 'Set the random seed for deterministic clustering (-1 uses the current time).')
  private long randomSeed

  /**
   * --td or --treeDistance
   * Distance value for dividing clusters from hierarchical tree.
   * The Dill Group at Stony Brook University uses a value of 2.0.
   */
  @Option(names = ['--td', '--treeDistance'], paramLabel = "2.0", defaultValue = "2.0",
      description = "The distance value where a hierarchical tree should be divided into clusters.")
  private double treeDistance

  /**
   * -w or --write Write out an archive of a representative structure from each cluster.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write an archive containing a representative from each cluster (algorithm=0).')
  private boolean write

  /**
   * The final argument(s) should be one or two filenames.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = 'The RMSD distance matrix (and an arc of the structures if using write (-w) flag.')
  private List<String> filenames = null

  /**
   * Generated list of clusters.
   */
  private final List<List<String>> clusterList = new ArrayList<>()

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
  List<List<String>> getClusterList() {
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

    ArrayList<double[]> distMatrix = new ArrayList<double[]>()

    // Either read in the distance matrix or calculate the distance matrix on the fly.
    if (readIn) {
      distMatrix = readInDistanceMatrix(distMatrix)
    } else {
      File file = null
      if (filenames != null && filenames.size() > 0) {
        activeAssembly = potentialFunctions.open(filenames.get(0))
        file = new File(filenames.get(0))
      } else if (activeAssembly == null) {
        logger.info(helpString())
        return null
      }

      distMatrix = calcDistanceMatrix(distMatrix, file)
    }

    // Either use Multi-K-Means++ or Hierarchical agglomerative clustering.
    if (algorithm == 0) {
      int dim = distMatrix.size()
      // Ensure cluster are within bounds
      if (maxClusters <= 0 || maxClusters >= dim) {
        int newSize = floorDiv(dim, 2)
        logger.warning(format(
            " Cluster of size %3d is out of bounds, using max cluster size of %3d.", maxClusters,
            newSize))
        maxClusters = newSize
      }
      // Array representing indices of structures that represent their cluster.
      int[] repStructs = new int[maxClusters]
      if (randomSeed == -1) {
        randomSeed = System.nanoTime()
      }
      kMeansCluster(distMatrix, maxClusters, trials, repStructs, randomSeed, true)
      if (write) {
        writeStructures(repStructs)
      }
    } else if (algorithm == 1) {
      hierarchicalAgglomerativeCluster(distMatrix)
    } else {
      logger.severe("Clustering algorithm has not been set.")
    }

    return this
  }

  /**
   * Write out structures corresponding to the representative of each cluster.
   * @param distMatrix Distances between structures (metric to determine clusters).
   * @param repStructs Array list for index of representative structures.
   */
  private void writeStructures(int[] repStructs) {
    String coordFileName = filenames.get(1)
    potentialFunctions.openAll(coordFileName)
    SystemFilter systemFilter = potentialFunctions.getFilter()
    activeAssembly = systemFilter.getActiveMolecularSystem()
    String fileName = FilenameUtils.getName(coordFileName)
    String ext = FilenameUtils.getExtension(fileName)
    fileName = FilenameUtils.removeExtension(fileName)
    File saveFile
    SystemFilter writeFilter
    if (ext.toUpperCase().contains("XYZ")) {
      saveFile = new File(fileName + ".xyz")
      writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
    } else if (ext.toUpperCase().contains("ARC")) {
      saveFile = new File(fileName + ".arc")
      saveFile = potentialFunctions.versionFile(saveFile)
      writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      logger.info("SaveFile: " + saveFile.getAbsolutePath())
      saveFile.createNewFile()
    } else {
      saveFile = new File(fileName + ".pdb")
      saveFile = potentialFunctions.versionFile(saveFile)
      writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
      int numModels = systemFilter.countNumModels()
      if (numModels > 1) {
        writeFilter.setModelNumbering(0)
      }
      writeFilter.writeFile(saveFile, true, false, false)
    }

    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
      int structNum = 0
      do {
        if (repStructs.contains(structNum++)) {
          Crystal crystal = activeAssembly.getCrystal()
          ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
          forceFieldEnergy.setCrystal(crystal)
          if (systemFilter instanceof PDBFilter) {
            saveFile.append("ENDMDL\n")
            PDBFilter pdbFilter = (PDBFilter) systemFilter
            pdbFilter.writeFile(saveFile, true, false, false)
          } else if (systemFilter instanceof XYZFilter) {
            writeFilter.writeFile(saveFile, true)
          }
        }
      } while (systemFilter.readNext())

      if (systemFilter instanceof PDBFilter) {
        saveFile.append("END\n")
      }
    }
  }

  /**
   * This method performs hierarchical clustering on a distance matrix. If the system isn't headless, a dendrogram
   * is printed of the clustered results. A PDB file for the centroid of each cluster is saved.
   *
   * @param distMatrix An ArrayList<double[]> that holds the distance matrix.
   */
  private void hierarchicalAgglomerativeCluster(ArrayList<double[]> distMatrix) {
    //Convert the distance matrix to a double[][] for the clustering algorithm.
    int distMatrixLength = distMatrix.size()
    double[][] distMatrixArray = new double[distMatrixLength][distMatrixLength]
    String[] names = new String[distMatrixLength]
    for (int i = 0; i < distMatrixLength; i++) {
      distMatrixArray[i] = distMatrix.get(i)
      //Set names of the clustered elements equal to the model number in the arc/pdb
      // by creating string of sequential numbers.
      names[i] = i.toString()
    }

    //Cluster the data. Note that the "cluster" object is actually the root node for the tree.
    ClusteringAlgorithm clusteringAlgorithm = new DefaultClusteringAlgorithm()
    com.apporiented.algorithm.clustering.Cluster rootNode =
        clusteringAlgorithm.performClustering(distMatrixArray,
            names, new CompleteLinkageStrategy())

    //Separate clusters based on the user-supplied treeDistance and fill the clusterList with lists holding
    //the model numbers that belong to a particular cluster.
    parseClusters(rootNode)

    //If the system is headless, skip all graphical components. Otherwise print the dendrogram.
    printDendrogram(rootNode)

    //Find the index for the centroid of each cluster in the clusterList.
    ArrayList<Integer> indicesOfCentroids = findCentroids(distMatrixArray)
    HashMap<Integer, Integer> pdbsToWrite = new HashMap<Integer, Integer>(indicesOfCentroids.size())

    //Get and store the index of each centroid in context of ALL models, not just the index relative to the cluster
    //the centroid belongs to.
    int counter = 0
    for (Integer centroidIndex : indicesOfCentroids) {
      ArrayList<String> cluster = clusterList.get(counter)
      pdbsToWrite.put(Integer.valueOf(cluster.get(centroidIndex)), counter)
      counter++
    }

    // Print out size of each cluster.
    logger.info(" ========== Cluster Sizes ========== ")
    counter = 0
    for (List cluster : clusterList) {
      logger.info(" Cluster " + counter + " Size: " + cluster.size())
      counter++
    }

    //Write out PDB files for each centroid structure of each cluster.
    final Map<Integer, Integer> sortedIds = pdbsToWrite.toSorted()
    SystemFilter systemFilter = potentialFunctions.getFilter()
    if (!sortedIds.isEmpty() && sortedIds.containsKey(0)) {
      String fileName = "centroid" + sortedIds.get(0).toString()
      potentialFunctions.saveAsPDB(activeAssembly, new File(fileName + ".pdb"))
      sortedIds.remove(0)
    }
    while ((!sortedIds.isEmpty()) && systemFilter.readNext()) {
      int modelNumber = systemFilter.getSnapshot() - 1
      if (sortedIds.containsKey(modelNumber)) {
        String fileName = "centroid" + sortedIds.get(modelNumber).toString()
        potentialFunctions.saveAsPDB(activeAssembly, new File(fileName + ".pdb"))
        sortedIds.remove(modelNumber)
      }
    }

    if (!sortedIds.isEmpty()) {
      logger.info(" Some models from clustering not found while parsing: "
          + Arrays.asList(sortedIds))
    }
  }

  /**
   * This method finds the centroid for each cluster in the clusterList. The index for the location of
   * the centroid of each cluster is returned.
   *
   * @param distMatrixArray The all vs. all distance matrix of RMSD values for each model/node in the hierarchical tree.
   * @return An ArrayList<Integer> containing the index for the centroid of each cluster in the clusterList.
   */
  private ArrayList<Integer> findCentroids(double[][] distMatrixArray) {
    // Find the centroid of each cluster.
    ArrayList<Integer> indicesOfCentroids = new ArrayList<Integer>()
    // Loop through every cluster.
    for (ArrayList<String> clusterNodes : clusterList) {
      ArrayList<Double> rmsds = new ArrayList<>()
      // Loop through every node in a cluster.
      for (String node1 : clusterNodes) {
        double rmsd = 0
        int counter = 0
        // Loop through every node in a cluster again for comparison.
        for (String node2 : clusterNodes) {
          if (node1 != node2) {
            // Find the rmsd of the two nodes from the all vs. all distance matrix of rmsds
            // and add it to the total.
            rmsd += distMatrixArray[node1.toInteger()][node2.toInteger()]
            counter++
          }
        }
        //Calculate the average rmsd for a node to all other nodes in a cluster.
        rmsd = rmsd / counter
        rmsds.add(rmsd.toDouble())
      }
      //Find minimum average rmsd a node has to all other nodes in a cluster. This node is the centroid of the cluster.
      Double minimum = Collections.min(rmsds)
      indicesOfCentroids.add(rmsds.indexOf(minimum).toInteger())
    }
    return indicesOfCentroids
  }

  /**
   * This method prints a dendrogram to the screen if the system is not headless.
   *
   * @param cluster The root cluster of the hierarchical tree.
   */
  private static void printDendrogram(com.apporiented.algorithm.clustering.Cluster cluster) {
    String headless = System.getProperty("java.awt.headless")
    if (!headless) {
      JFrame frame = new JFrame()
      frame.setSize(400, 300)
      frame.setLocation(400, 300)
      frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
      JPanel content = new JPanel()
      DendrogramPanel dp = new DendrogramPanel()
      frame.setContentPane(content)
      content.setBackground(Color.red)
      content.setLayout(new BorderLayout())
      content.add(dp, BorderLayout.CENTER)
      dp.setBackground(Color.WHITE)
      dp.setLineColor(Color.BLACK)
      dp.setScaleValueDecimals(0)
      dp.setScaleValueInterval(1)
      dp.setShowDistances(false)
      dp.setModel(cluster)
      frame.setVisible(true)
    }
  }

  /**
   * This method parses clusters in the hierarchical tree and can print the model that belongs to each cluster.
   *
   * @param root The root node of the hierarchical tree.
   */
  private void parseClusters(final com.apporiented.algorithm.clustering.Cluster root) {
    populateChildren(root, false)
    // Print out 0 indexed model numbers that belong to each cluster.
    /*
        for (final List<String> curCluster : clusterList) {
        System.out.println("==========Iterating over next cluster==========")
        for (final String element : curCluster) {
            System.out.println(element)
        }
    }
    */
  }

  /**
   * This method adds node names to a cluster list to indicate which nodes belong to a particular cluster. Clusters
   * from the hierarchical tree are determined based on the treeDistance cutoff.
   *
   * @param rootNode The cluster object being iterated over.
   * @param curCluster A boolean indicating that a cluster has been identified based on the treeDistance cutoff.
   * The curCluster is set to true once for each cluster.
   */
  private void populateChildren(final com.apporiented.algorithm.clustering.Cluster rootNode,
      boolean curCluster) {
    final double distance = rootNode.getDistanceValue()
    final List<com.apporiented.algorithm.clustering.Cluster> children = rootNode.getChildren()
    if (!curCluster && (distance <= treeDistance)) {
      curCluster = true
      List<String> clusterSubList = new ArrayList<>()
      clusterList.add(clusterSubList)
      populateChildren(rootNode, curCluster)
    } else if (!children.empty) {
      for (final com.apporiented.algorithm.clustering.Cluster child : children) {
        populateChildren(child, curCluster)
      }
    } else {
      final int clusterListSize = clusterList.size()
      if (clusterListSize != 0) {
        clusterList.get(clusterListSize - 1).add(rootNode.getName())
      } else {
        logger.severe(" SEVERE: A node cannot be added to the tree.")
      }
    }
  }

  /**
   * This method reads in the distance matrix from an input file.
   *
   * @param distMatrix An empty ArrayList<double[]> to hold the distance matrix values.
   * @return ArrayList < double [ ] >   that holds all values for the read in distance matrix.
   */
  private ArrayList<double[]> readInDistanceMatrix(ArrayList<double[]> distMatrix) {
    File file = new File(filenames.get(0))
    int nDim = 0

    // Read in the RMSD matrix.
    try (FileReader fr = new FileReader(file)
        BufferedReader br = new BufferedReader(fr)) {
      String data = br.readLine()
      // Check for blank lines at the top of the file
      while (data != null && data.trim() == "") {
        data = br.readLine()
      }
      if (data == null) {
        logger.severe("No data in RMSD file.")
      }
      String[] tokens = data.trim().split(" +")
      // Expect a n x n matrix of distance values.
      nDim = tokens.size()
      for (int i = 0; i < nDim; i++) {
        double[] tokens2 = new double[nDim]
        for (int j = 0; j < nDim; j++) {
          tokens2[j] = tokens[j].toDouble()
        }
        distMatrix.add(tokens2)
        data = br.readLine()
        if (data != null) {
          tokens = data.trim().split(" +")
        }
      }
    } catch (IOException e) {
      logger.severe(e.toString())
    }

    if (distMatrix == null) {
      logger.severe("Input read attempt failed.")
    }

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(format(" Original Distance Matrix:\n"))
      String tempString = ""
      for (double[] i : distMatrix) {
        for (int j = 0; j < nDim; j++) {
          tempString += format("%f\t", i[j])
        }
        tempString += "\n"
      }
      logger.finest(tempString)
    }

    return distMatrix
  }

  /**
   * This method calculates the distance matrix of all molecular assemblies in an arc/multiple model file.
   *
   * @param distMatrix An empty ArrayList<double[]> to hold the distance matrix values.
   * @return ArrayList < double [ ] >     that holds all values for the read in distance matrix.
   */
  private ArrayList<double[]> calcDistanceMatrix(ArrayList<double[]> distMatrix, File file) {
    //Prepare the superpose object and binding.
    Binding binding = new Binding()
    Superpose superpose = new Superpose()
    superpose.setBinding(binding)

    // Set-up the input arguments for the Superpose script.
    String[] args = ["--aS", "2", "-A", "-s", start, "-f", finish, "--store", "-v", file]
    binding.setVariable("args", args)

    // Evaluate the superpose script to get the distance matrix of RMSD values.
    superpose.run()
    double[][] tempDistMatrix = superpose.getDistanceMatrix()
    int matrixLength = tempDistMatrix.length

    for (int i = 0; i < matrixLength; i++) {
      distMatrix.add(tempDistMatrix[i])
    }

    return distMatrix
  }
}


