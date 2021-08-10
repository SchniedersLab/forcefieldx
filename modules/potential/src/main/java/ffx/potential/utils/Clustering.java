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
package ffx.potential.utils;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

import com.apporiented.algorithm.clustering.Cluster;
import com.apporiented.algorithm.clustering.ClusteringAlgorithm;
import com.apporiented.algorithm.clustering.DefaultClusteringAlgorithm;
import com.apporiented.algorithm.clustering.SingleLinkageStrategy;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.evaluation.SumOfClusterVariances;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Cluster contains methods utilized in the Cluster.groovy file.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 */
public class Clustering {

  private static final Logger logger = Logger.getLogger(PotentialsUtils.class.getName());

  /**
   * Number of iterations for k-means clustering
   */
  private static final int NUM_ITERATIONS = 10000;

  /**
   * Perform a k-means clustering for a specified number of clusters.
   *
   * @param distMatrix Coordinate input serves as the data points.
   * @param maxClusters Number of clusters to use (k).
   * @param numTrials Number of trials for the Multi K-Means++ algorithm.
   * @param seed The seed to use for clustering (-1 uses the current system time).
   * @return The clusters.
   */
  public static List<CentroidCluster<Conformation>> kMeansClustering(List<double[]> distMatrix,
      int maxClusters, int numTrials, long seed) {
    // Square distance matrix size (dim x dim).
    int dim = distMatrix.size();
    List<Conformation> conformationList = new ArrayList<>();
    for (int i = 0; i < dim; i++) {
      double[] row = distMatrix.get(i);
      // Check that the input data is appropriate.
      if (row.length != dim) {
        logger.severe(format(" Row %d of the distance matrix (%d x %d) has %d columns.",
            i, dim, dim, row.length));
      }
      conformationList.add(new Conformation(row, i));
    }

    if (maxClusters < 1 || maxClusters >= dim) {
      maxClusters = dim / 2;
    }

    // Input the RMSD matrix to the clustering algorithm
    // Use the org.apache.commons.math3.ml.clustering package.
    KMeansPlusPlusClusterer<Conformation> kMeansPlusPlusClusterer = new KMeansPlusPlusClusterer<>(
        maxClusters, NUM_ITERATIONS);
    // Set the random seed for deterministic clustering.
    RandomGenerator randomGenerator = kMeansPlusPlusClusterer.getRandomGenerator();
    randomGenerator.setSeed(seed);
    // Create a MultiKMeansPlusPlusClusterer
    MultiKMeansPlusPlusClusterer<Conformation> multiKMeansPlusPlusClusterer =
        new MultiKMeansPlusPlusClusterer<>(kMeansPlusPlusClusterer, numTrials);

    // Perform the clustering.
    return multiKMeansPlusPlusClusterer.cluster(conformationList);
  }

  /**
   * This method performs hierarchical clustering on a distance matrix. If the system isn't headless,
   * a dendrogram is printed of the clustered results. A PDB file for the centroid of each cluster is
   * saved.
   *
   * @param distanceMatrix An ArrayList<double[]> that holds the distance matrix.
   * @param threshold the distance used to separate clusters.
   * @return Return a list of CentroidClusters.
   */
  public static List<CentroidCluster<Conformation>> hierarchicalClustering(
      List<double[]> distanceMatrix, double threshold) {

    // Convert the distanceMatrix to a double[][] for the clustering algorithm.
    int dim = distanceMatrix.size();

    double[][] distMatrixArray = new double[dim][];
    String[] names = new String[dim];
    for (int i = 0; i < dim; i++) {
      distMatrixArray[i] = distanceMatrix.get(i);
      names[i] = Integer.toString(i);
    }

    // Cluster the data.
    ClusteringAlgorithm clusteringAlgorithm = new DefaultClusteringAlgorithm();
    List<Cluster> clusterList = clusteringAlgorithm.performFlatClustering(distMatrixArray, names,
        new SingleLinkageStrategy(), threshold);

    // Convert the cluster format to CentroidClusters
    return clustersToCentroidClusters(distMatrixArray, clusterList);
  }

  /**
   * Analyze a list of CentroidClusters.
   *
   * @param clusters The List of CentroidClusters to analyze.
   * @param repStructs Store a representative conformation for each cluster.
   * @param verbose If true, use verbose printing.
   */
  public static void analyzeClusters(List<CentroidCluster<Conformation>> clusters, List<Integer> repStructs,
      boolean verbose) {
    // Number of clusters.
    int nClusters = clusters.size();
    double meanClusterRMSD = 0.0;

    // Loop over clusters
    for (int i = 0; i < nClusters; i++) {
      CentroidCluster<Conformation> clusterI = clusters.get(i);

      List<Conformation> conformations = clusterI.getPoints();
      int nConformers = conformations.size();
      StringBuilder sb = new StringBuilder(
          format(" Cluster %d with %d conformations\n  Conformations:", i + 1, nConformers));

      double minRMS = Double.MAX_VALUE;
      int minID = -1;
      for (Conformation conformation : conformations) {
        int row = conformation.index;
        sb.append(format(" %d", row + 1));
        if (nConformers > 1) {
          // Get a row of the RMSD matrix.
          double[] rmsd = conformation.getPoint();
          // Calculate the sum of squares.
          double wrms = 0;
          for (Conformation conformation2 : conformations) {
            int col = conformation2.index;
            if (col == row) {
              continue;
            }
            wrms += rmsd[col] * rmsd[col];
          }
          // Calculate the root mean sum of squares distance within the cluster.
          wrms = sqrt(wrms / (nConformers - 1));
          if (wrms < minRMS) {
            minRMS = wrms;
            minID = conformation.getIndex();
          }
        } else {
          // Only 1 conformer in this cluster.
          minID = row;
          minRMS = 0.0;
        }
      }

      // Calculate the RMSD within the cluster.
      if (nConformers > 1) {
        double clusterRMSD = clusterRMSD(conformations);
        meanClusterRMSD += clusterRMSD;
        sb.append(format("\n  RMSD within the cluster:\t %6.4f A.\n", clusterRMSD));

        // minID contains the index for the representative conformer.
        sb.append(format("  Minimum RMSD conformer %d:\t %6.4f A.\n", minID + 1, minRMS));
      } else {
        sb.append("\n");
      }

      repStructs.add(minID);

      if (verbose) {
        logger.info(sb.toString());
      }
    }

    if (verbose) {
      logger.info(
          format(" Mean RMSD within clusters: \t %6.4f A.", meanClusterRMSD / nClusters));
      double sumOfClusterVariances = sumOfClusterVariances(clusters);
      logger.info(
          format(" Sum of cluster variances:  \t %6.4f A.\n", sumOfClusterVariances));
    }

  }

  /**
   * Compute the RMSD of one cluster.
   *
   * @param conformations Conformers for this cluster.
   * @return The RMSD for the cluster.
   */
  private static double clusterRMSD(List<Conformation> conformations) {
    int nConformers = conformations.size();

    // If there is only 1 conformer, the RMSD is 0.
    if (nConformers == 1) {
      return 0.0;
    }

    // Calculate the RMSD within the cluster.
    double sum = 0.0;
    int count = 0;
    for (int j = 0; j < nConformers; j++) {
      Conformation conformation = conformations.get(j);
      double[] rmsd = conformation.rmsd;
      for (int k = j + 1; k < nConformers; k++) {
        Conformation conformation2 = conformations.get(k);
        int col = conformation2.index;
        sum += rmsd[col] * rmsd[col];
        count++;
      }
    }

    return sqrt(sum / count);
  }

  /**
   * Compute the Sum of Cluster Variances.
   *
   * @param clusters The cluster to operate on.
   * @return The sum of cluster variances.
   */
  private static double sumOfClusterVariances(List<CentroidCluster<Conformation>> clusters) {
    SumOfClusterVariances<Conformation> sumOfClusterVariances = new SumOfClusterVariances<>(
        new EuclideanDistance());
    return sumOfClusterVariances.score(clusters);
  }

  /**
   * Convert clusters defined by a List of Strings to Apache Math style CentroidClusters.
   *
   * @param distMatrixArray Distance matrix.
   * @param clusterList Input List of Clusters.
   * @return Return a List of CentroidClusters.
   */
  private static List<CentroidCluster<Conformation>> clustersToCentroidClusters(
      double[][] distMatrixArray, List<Cluster> clusterList) {

    List<CentroidCluster<Conformation>> centroidClusters = new ArrayList<>();
    int dim = distMatrixArray.length;

    // Loop over clusters defined by lists of Strings.
    for (Cluster cluster : clusterList) {

      List<String> names = new ArrayList<>();
      collectNames(cluster, names);

      // Collect conformations for this cluster.
      List<Conformation> conformations = new ArrayList<>();
      for (String name : names) {
        int index = Integer.parseInt(name);
        conformations.add(new Conformation(distMatrixArray[index], index));
      }

      // Compute its centroid.
      Conformation centroid = centroidOf(conformations, dim);

      // Create a new CentroidCluster and add the conformations.
      CentroidCluster<Conformation> centroidCluster = new CentroidCluster<>(centroid);
      centroidClusters.add(centroidCluster);
      for (Conformation conformation : conformations) {
        centroidCluster.addPoint(conformation);
      }
    }
    return centroidClusters;
  }

  /**
   * Collect the names of each leaf in Cluster.
   *
   * @param cluster The cluster to operate on.
   * @param names A List of leaf names.
   */
  private static void collectNames(Cluster cluster, List<String> names) {
    if (cluster.isLeaf()) {
      names.add(cluster.getName());
    } else {
      for (Cluster c : cluster.getChildren()) {
        collectNames(c, names);
      }
    }
  }


  /**
   * Computes the centroid for a set of points.
   *
   * @param points the set of points
   * @param dimension the point dimension
   * @return the computed centroid for the set of points
   */
  private static Conformation centroidOf(final List<Conformation> points, final int dimension) {
    final double[] centroid = new double[dimension];
    for (final Conformation p : points) {
      final double[] point = p.getPoint();
      for (int i = 0; i < centroid.length; i++) {
        centroid[i] += point[i];
      }
    }
    for (int i = 0; i < centroid.length; i++) {
      centroid[i] /= points.size();
    }
    return new Conformation(centroid, 0);
  }

  /**
   * Class for cluster objects.
   */
  public static class Conformation implements Clusterable {

    private final double[] rmsd;
    private final int index;

    Conformation(double[] rmsd, int index) {
      this.rmsd = rmsd;
      this.index = index;
    }

    public double[] getPoint() {
      return rmsd;
    }

    int getIndex() {
      return index;
    }
  }
}
