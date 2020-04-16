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
import ffx.potential.cli.PotentialScript
import ffx.potential.groovy.Superpose
import ffx.potential.parsers.SystemFilter
import org.apache.commons.math3.ml.clustering.CentroidCluster
import org.apache.commons.math3.ml.clustering.Clusterable
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.swing.*
import java.awt.*
import java.util.List
import java.util.logging.Level

import static org.apache.commons.math3.util.FastMath.pow
import static org.apache.commons.math3.util.FastMath.sqrt

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
class Cluster extends PotentialScript {

  /**
   * -a or --algorithm Clustering algorithm to use.
   * Choices are kmeans (0), multikmeans (1), and hierarchical (2).
   */
  @Option(names = ['-a', '--algorithm'], paramLabel = "0", defaultValue = "0",
      description = "Algorithm to be used during clustering: kmeans (0), multikmeans (1), hierarchical (2)")
  int algorithm = 0

  /**
   * -k or --clusters Clustering algorithm to use.
   */
  @Option(names = ['-k', '--clusters'], paramLabel = "0", defaultValue = "0",
      description = "Maximum number of kmeans clusters for the input data.")
  private int maxClusters = 0

  /**
   * -r or --readInDistMat The algorithm should read in a provided distance matrix rather than the matrix being generated on the fly.
   */
  @Option(names = ['-r', '--readInDistMat'], paramLabel = "false", defaultValue = "false",
      description = "Tells algorithm to read in the distance matrix from an input file.")
  Boolean readIn = false

  /**
   * -numI or --numIterations The number of times each kmeans cluster should be determined.
   */
  @Option(names = ['--numI', '--numIterations'], paramLabel = "1", defaultValue = "1",
      description = "Number of repetitions for each kmeans cluster.")
  private int numIterations = 1

  /**
   * -s or --start Atom number where RMSD calculation of structure will begin.
   */
  @Option(names = ['-s', '--start'], paramLabel = "1", defaultValue = "1",
      description = 'Starting atom to include in the RMSD calculation.')
  private String start = "1"

  /**
   * -f or --final Atom number where RMSD calculation of structure will end.
   */
  @Option(names = ['-f', '--final'], paramLabel = "nAtoms",
      description = 'Final atom to include in the RMSD calculation.')
  private String finish = Integer.MAX_VALUE.toString()

  /**
   * --clusDis or --clusterDistance Distance value for dividing clusters from hierarchical tree.
   * The Dill Group at Stony Brook University uses a value of 2.0.
   */
  @Option(names = ['--treeDis', '--treeDistance'], paramLabel = "2.0", defaultValue = "2.0",
      description = "The distance value where a hierarchical tree should be divided into clusters.")
  private double treeDistance = 2.0

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = 'The RMSD matrix.')
  List<String> filenames = null

  private File baseDir = null

  private final List<List<String>> clusterList = new ArrayList<>()

  void setBaseDir(File baseDir) {
    this.baseDir = baseDir
  }

  /**
   * Execute the script.
   */
  @Override
  Cluster run() {
    if (!init()) {
      return null
    }

    if (filenames == null || filenames.isEmpty()) {
      logger.info(helpString())
      return null
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

    // Either use kmeans clustering or hierarchical agglomerative clustering.
    if (algorithm == 0 || algorithm == 1) {
      kmeansCluster(distMatrix)
    } else if (algorithm == 2) {
      hierarchicalAgglomerativeCluster(distMatrix)
    } else {
      logger.severe("Clustering algorithm has not been set.")
    }

    return this
  }

  void kmeansCluster(ArrayList<double[]> distMatrix) {
    // Input the RMSD matrix to the clustering algorithm
    // Use the org.apache.commons.math3.ml.clustering package.
    int minClusters = 1

    if (maxClusters <= 0 || maxClusters > distMatrix.size() - 1) {
      maxClusters = distMatrix.size() - 1
    }

    double[] twss = new double[distMatrix.size() - 1]
    if (algorithm == 0) {
      minClusters = maxClusters
      twss = new double[1]
    }

    for (int i = 0; i < twss.size(); i++) {
      twss[i] = Double.MAX_VALUE
    }
    for (int clusters = minClusters; clusters <= maxClusters; clusters++) {
      logger.info(String.format("%d Clusters:", clusters))
      for (int k = 0; k < numIterations; k++) {
        double currentTWSS = 0
        KMeansPlusPlusClusterer<ClusterWrapper> kClust1 = new KMeansPlusPlusClusterer<ClusterWrapper>(
            clusters, 10000)
        List<ClusterWrapper> myClusterables = new ArrayList<ClusterWrapper>()
        int id = 0
        for (double[] i : distMatrix) {
          myClusterables.add(new ClusterWrapper(i, id))
          id++
        }
        List<CentroidCluster<ClusterWrapper>> kClusters = kClust1.cluster(myClusterables)

        if (algorithm == 1) {
          MultiKMeansPlusPlusClusterer<ClusterWrapper> kClust2 = new MultiKMeansPlusPlusClusterer<>(
              kClust1, 10000)
          kClusters = kClust2.cluster(myClusterables)
        }

        // TODO: Output the clusters in a useful way.
        //Temp output method prints to screen
        for (int i = 0; i < kClusters.size(); i++) {
          double wss = 0 // Reset cluster within distance
          double[] sum = new double[kClusters.get(0).getPoints()[0].getPoint().size()]
          for (ClusterWrapper clusterWrapper : kClusters.get(i).getPoints()) {
            double[] distArray = clusterWrapper.getPoint()
            // Implement WSS
            for (int j = 0; j < sum.size(); j++) {
              wss += pow(distArray[j] - kClusters.get(i).getCenter().getPoint()[j], 2)
            }
          }
          wss = sqrt(wss)
          currentTWSS += wss
        }
        if (algorithm == 1 && currentTWSS < twss[clusters - 1]) {
          twss[clusters - 1] = currentTWSS
        }
      }
      logger.info(String.format("Total WSS: %f\n", twss[clusters - 1]))
    }
    if (algorithm == 0) {
      double[] d2TWSS = new double[twss.size() - 2]
      for (int i = 1; i < twss.size() - 1; i++) {
        d2TWSS[i - 1] = twss[i + 1] - 2 * twss[i] + twss[i - 1]
        logger.info(String.format("Second Derivative: %f", d2TWSS[i - 1]))
      }
    }
  }

  /**
   * This method performs hierarchical clustering on a distance matrix. If the system isn't headless, a dendrogram
   * is printed of the clustered results. A PDB file for the centroid of each cluster is saved.
   *
   * @param distMatrix An ArrayList<double[]> that holds the distance matrix.
   */
  void hierarchicalAgglomerativeCluster(ArrayList<double[]> distMatrix) {
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

    //Print out size of each cluster.
    System.out.println("==========Cluster Sizes==========")
    counter = 0
    for (List cluster : clusterList) {
      System.out.println(" Cluster " + counter + " Size: " + cluster.size())
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
      System.out.println(" Some models from clustering not found while parsing: " + Arrays.asList(
          sortedIds))
    }
  }

  /**
   * This method finds the centroid for each cluster in the clusterList. The index for the location of
   * the centroid of each cluster is returned.
   *
   * @param distMatrixArray The all vs. all distance matrix of RMSD values for each model/node in the hierarchical tree.
   * @return An ArrayList<Integer> containing the index for the centroid of each cluster in the clusterList.
   */
  ArrayList<Integer> findCentroids(double[][] distMatrixArray) {
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
  static void printDendrogram(com.apporiented.algorithm.clustering.Cluster cluster) {
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
  void parseClusters(final com.apporiented.algorithm.clustering.Cluster root) {
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
  void populateChildren(final com.apporiented.algorithm.clustering.Cluster rootNode,
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
  ArrayList<double[]> readInDistanceMatrix(ArrayList<double[]> distMatrix) {
    File file = new File(filenames.get(0))
    int nDim = 0

    // Read in the RMSD matrix.
    try {
      FileReader fr = new FileReader(file)
      BufferedReader br = new BufferedReader(fr)
      String data = br.readLine()
      // Check for blank lines at the top of the file
      while (data != null && data.trim() == "") {
        data = br.readLine()
      }
      if (data == null) {
        logger.severe("No data in RMSD file.")
      }
      String[] tokens = data.trim().split("\t")
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
          tokens = data.trim().split("\t")
        }
      }
      br.close()
      fr.close()
    } catch (IOException e) {
      logger.severe(e.toString())
    }
    if (distMatrix == null) {
      logger.severe("Input read attempt failed.")
    }
    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(String.format("Original Distance Matrix:\n"))
      String tempString = ""
      for (double[] i : distMatrix) {
        for (int j = 0; j < nDim; j++) {
          tempString += String.format("%f\t", i[j])
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
  ArrayList<double[]> calcDistanceMatrix(ArrayList<double[]> distMatrix, File file) {
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

  class ClusterWrapper implements Clusterable {

    private double[] point
    private final int UUID

    ClusterWrapper(double[] distances, int ID) {
      this.point = distances
      UUID = ID
    }

    double[] getPoint() {
      return point
    }

    int getUUID() {
      return UUID
    }

    /* Min-Max normalization of distances (not important if all inputs are on the same scale)
          double minimumDist=0;
          double maximumDist=0;
          for (double[] distArray in distMatrix){
              for (double dist in distArray){
                  if( minimumDist>dist){
                      minimumDist = dist;
                  }
                  if(maximumDist<dist){
                      maximumDist = dist;
                  }
              }
          }
          for(int i = 0; i<distMatrix.size(); i++) {
              for (int j = 0; j < distMatrix.get(i).size(); j++) {
                  distMatrix.get(i)[j] = (distMatrix.get(i)[j] - minimumDist) / (maximumDist - minimumDist);
              }
          }

          if (logger.isLoggable(Level.FINEST)) {
              logger.finest(String.format("\nNormalized Matrix:\n"));
              String tempString2 = "";
              for (double[] i : distMatrix) {
                  for (int j = 0; j < nDim; j++) {
                      tempString2 += String.format("%f\t", i[j]);
                  }
                  tempString2 += "\n";
              }
              logger.finest(tempString2);
          } */
  }
}


