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
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
public class Cluster {

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
   * @param repStructs Identity of representative structures to write to archive.
   * @param seed The seed to use for clustering (-1 uses the current system time).
   * @param verbose If true, log cluster results.
   * @return The clusters.
   */
  public static List<CentroidCluster<Conformation>> kMeansCluster(ArrayList<double[]> distMatrix,
      int maxClusters, int numTrials, int[] repStructs, long seed, boolean verbose) {
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
    List<CentroidCluster<Conformation>> clusters = multiKMeansPlusPlusClusterer.cluster(
        conformationList);
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
      double clusterRMSD = clusterRMSD(conformations);
      meanClusterRMSD += clusterRMSD;
      sb.append(format("\n  RMSD within the cluster:\t %6.4f A.\n", clusterRMSD));

      // minID contains the index for the representative conformer.
      // TODO: Use energy/density as an alternative than centroid distance.
      repStructs[i] = minID;
      sb.append(format("  Minimum RMSD conformer %d:\t %6.4f A.\n", minID + 1, minRMS));

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

    return clusters;
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
    SumOfClusterVariances sumOfClusterVariances = new SumOfClusterVariances(new EuclideanDistance());
    return sumOfClusterVariances.score(clusters);
  }

  /**
   * Perform k-means clustering with multiple cluster sizes (k).
   *
   * @param distMatrix Coordinate input serves as data.
   * @param maxClusters Most clusters to include.
   */
  public static void kMeansCluster(ArrayList<double[]> distMatrix, int maxClusters) {
    // Input the RMSD matrix to the clustering algorithm
    // Use the org.apache.commons.math3.ml.clustering package.
    int minClusters = 1;
    int distMatrixLength = distMatrix.size();

    int dim = distMatrix.size();
    if (maxClusters < 1 || maxClusters >= dim) {
      maxClusters = dim / 2;
    }

    double[] twss = new double[distMatrixLength - 1];
    int twssSize = twss.length;
    for (int i = 0; i < twssSize; i++) {
      twss[i] = Double.MAX_VALUE;
    }

    for (int clusters = minClusters; clusters <= maxClusters; clusters++) {
      logger.info(format(" %d Cluster(s):", clusters));
      // algorithm 1 is multi k-means clustering
      // Displays Total Within Sum of Squares to give quantification of cluster quality
      // Cutoff should be made where d2TWSS is small.
      // Total Within Sum of Squares
      double currentTWSS = 0;
      KMeansPlusPlusClusterer<Conformation> kClust1 = new KMeansPlusPlusClusterer<>(
          clusters, NUM_ITERATIONS);
      List<Conformation> myClusterables = new ArrayList<>();
      int id = 0;
      for (double[] i : distMatrix) {
        myClusterables.add(new Conformation(i, id));
        id++;
      }

      MultiKMeansPlusPlusClusterer<Conformation> kClust2 = new MultiKMeansPlusPlusClusterer<>(
          kClust1, NUM_ITERATIONS);
      List<CentroidCluster<Conformation>> kClusters = kClust2.cluster(myClusterables);

      int size = kClusters.get(0).getPoints().get(0).getPoint().length;
      for (CentroidCluster<Conformation> centroid : kClusters) {
        // Reset cluster within distance
        double wss = 0;
        for (Conformation conformation : centroid.getPoints()) {
          double[] distArray = conformation.getPoint();
          // calculate WSS
          for (int j = 0; j < size; j++) {
            wss += pow(distArray[j] - centroid.getCenter().getPoint()[j], 2);
          }
        }
        wss = sqrt(wss);
        currentTWSS += wss;
      }
      // Identify cluster with the minimal WSS.
      if (currentTWSS < twss[clusters - 1]) {
        twss[clusters - 1] = currentTWSS;
      }
      // TODO: Make cluster output more useful.
      if (clusters == 1) {
        logger.info(format(" TWSS: %8.4f\n", twss[clusters - 1]));
      } else if (clusters == 2) {
        if (twssSize > 1) {
          logger.info(format(" TWSS: %8.4f\tdTWSS: %8.4f\n", twss[clusters - 1],
              twss[clusters - 2] - twss[clusters - 1]));
        } else {
          logger.warning(" Distance matrix may be too small for number of clusters.");
        }
      } else {
        // Finite difference approximation of d2TWSS
        if (twssSize > 2) {
          double d2TWSS = twss[clusters - 1] - 2 * twss[clusters - 2] + twss[clusters - 3];
          logger.info(format(" TWSS: %8.4f\tdTWSS: %8.4f\td2TWSS: %8.4f\n",
              twss[clusters - 1], twss[clusters - 2] - twss[clusters - 1], d2TWSS));
        } else {
          logger.warning(" Distance Matrix size is less than three.");
        }
      }
    }
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
