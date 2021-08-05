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

import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
    private static final int numIterations = 10000;

    /**
     * Perform a k-means clustering for a specified number of clusters.
     * @param distMatrix Coordinate input serves as the data points.
     * @param maxClusters Number of clusters to use (k).
     * @param repStructs Identity of representative structures to write to archive.
     */
    public static void kMeansCluster(ArrayList<double[]> distMatrix, int maxClusters, int[] repStructs) {
        // Input the RMSD matrix to the clustering algorithm
        // Use the org.apache.commons.math3.ml.clustering package.
        KMeansPlusPlusClusterer<ClusterWrapper> kClust1 = new KMeansPlusPlusClusterer<ClusterWrapper>(
                maxClusters, numIterations);
        List<ClusterWrapper> myClusterables = new ArrayList<ClusterWrapper>();
        int id = 0;
        for (double[] i : distMatrix) {
            myClusterables.add(new ClusterWrapper(i, id));
            id++;
        }
        List<CentroidCluster<ClusterWrapper>> kClusters = kClust1.cluster(myClusterables);
        // Number of observations:
        int size = kClusters.get(0).getPoints().get(0).getPoint().length;
        // Loop over clusters:
        for (int i = 0 ; i<kClusters.size(); i++) {
            logger.info(format(" Cluster: %2d", i));
            // Reset cluster within distance
            double minWSS = Double.MAX_VALUE;
            int minID = -1;
            double wss = 0;
            for (ClusterWrapper datum : kClusters.get(i).getPoints()) {
                logger.info(format("   structure: %2d", datum.getUUID()));
                double[] distArray = datum.getPoint();
                // calculate WSS
                for (int j = 0; j < size; j++) {
                    wss += pow(distArray[j] - kClusters.get(i).getCenter().getPoint()[j], 2);
                }
                wss = sqrt(wss);
                if (wss < minWSS) {
                    minWSS = wss;
                    minID = datum.getUUID();
                }
            }
            // minID should contain the integer for the representative structure.
            // TODO use energy/density rather than centroid distance
            repStructs[i] = minID;
            logger.info(format("Structure %3d represents cluster %2d", minID, i + 1));
        }
    }

    /**
     * Perform k-means clustering with multiple cluster sizes (k).
     * @param distMatrix Coordinate input serves as data.
     * @param maxClusters Most clusters to include.
     */
    public static void kMeansCluster(ArrayList<double[]> distMatrix, int maxClusters) {
        // Input the RMSD matrix to the clustering algorithm
        // Use the org.apache.commons.math3.ml.clustering package.
        int minClusters = 1;
        int distMatrixLength = distMatrix.size();

        if (maxClusters <= 0 || maxClusters > distMatrixLength - 1) {
            maxClusters = distMatrixLength - 1;
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
            KMeansPlusPlusClusterer<ClusterWrapper> kClust1 = new KMeansPlusPlusClusterer<ClusterWrapper>(
                    clusters, numIterations);
            List<ClusterWrapper> myClusterables = new ArrayList<ClusterWrapper>();
            int id = 0;
            for (double[] i : distMatrix) {
                myClusterables.add(new ClusterWrapper(i, id));
                id++;
            }
            List<CentroidCluster<ClusterWrapper>> kClusters = kClust1.cluster(myClusterables);
            MultiKMeansPlusPlusClusterer<ClusterWrapper> kClust2 = new MultiKMeansPlusPlusClusterer<>(
                    kClust1, numIterations);
            kClusters = kClust2.cluster(myClusterables);

            int size = kClusters.get(0).getPoints().get(0).getPoint().length;
            for (CentroidCluster<ClusterWrapper> centroid : kClusters) {
                // Reset cluster within distance
                double wss = 0;
                for (ClusterWrapper clusterWrapper : centroid.getPoints()) {
                    double[] distArray = clusterWrapper.getPoint();
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
                    logger.info(format(" TWSS: %8.4f\tdTWSS: %8.4f\n", twss[clusters - 1], twss[clusters - 2] - twss[clusters - 1]));
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
    public static class ClusterWrapper implements Clusterable {

        private final double[] point;
        private final int UUID;

        ClusterWrapper(double[] distances, int ID) {
            this.point = distances;
            UUID = ID;
        }

        public double[] getPoint() {
            return point;
        }

        int getUUID() {
            return UUID;
        }

        // Code to be used if coordinates differ in unit.
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
