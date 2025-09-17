// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
// ******************************************************************************
package ffx.numerics.clustering;

import java.util.ArrayList;
import java.util.List;

/**
 * Clustering algorithm that consumes a condensed (pdist-style) upper-triangular
 * distance array to produce hierarchical agglomerative clusters. Supports flat
 * clustering by threshold; weighted inputs delegate to unweighted behavior.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PDistClusteringAlgorithm implements ClusteringAlgorithm {

  /**
   * Performs hierarchical agglomerative clustering using a condensed pdist-like matrix.
   *
   * @param distances       a 1 x M array holding the upper-triangular distances in row-major pdist order
   * @param clusterNames    names corresponding to the N items (M = N*(N-1)/2)
   * @param linkageStrategy linkage criterion used during agglomeration
   * @return root Cluster of the resulting hierarchy
   */
  @Override
  public Cluster performClustering(double[][] distances,
                                   String[] clusterNames, LinkageStrategy linkageStrategy) {

    /* Argument checks */
    if (distances == null || distances.length == 0) {
      throw new IllegalArgumentException("Invalid distance matrix");
    }
    if (distances[0].length != clusterNames.length
        * (clusterNames.length - 1) / 2) {
      throw new IllegalArgumentException("Invalid cluster name array");
    }
    if (linkageStrategy == null) {
      throw new IllegalArgumentException("Undefined linkage strategy");
    }

    /* Setup model */
    List<Cluster> clusters = createClusters(clusterNames);
    DistanceMap linkages = createLinkages(distances, clusters);

    /* Process */
    HierarchyBuilder builder = new HierarchyBuilder(clusters, linkages);
    while (!builder.isTreeComplete()) {
      builder.agglomerate(linkageStrategy);
    }

    return builder.getRootCluster();
  }

  /**
   * Produces a flat clustering from a condensed distance matrix by agglomerating until the
   * threshold is exceeded.
   *
   * @param distances       a 1 x M condensed distance array (pdist order)
   * @param clusterNames    names of the N items
   * @param linkageStrategy linkage criterion used during agglomeration
   * @param threshold       maximum allowed linkage distance for merging
   * @return list of clusters at the chosen cut
   */
  @Override
  public List<Cluster> performFlatClustering(double[][] distances,
                                             String[] clusterNames, LinkageStrategy linkageStrategy, Double threshold) {

    /* Argument checks */
    if (distances == null || distances.length == 0) {
      throw new IllegalArgumentException("Invalid distance matrix");
    }
    if (distances[0].length != clusterNames.length
        * (clusterNames.length - 1) / 2) {
      throw new IllegalArgumentException("Invalid cluster name array");
    }
    if (linkageStrategy == null) {
      throw new IllegalArgumentException("Undefined linkage strategy");
    }

    /* Setup model */
    List<Cluster> clusters = createClusters(clusterNames);
    DistanceMap linkages = createLinkages(distances, clusters);

    /* Process */
    HierarchyBuilder builder = new HierarchyBuilder(clusters, linkages);
    return builder.flatAgg(linkageStrategy, threshold);
  }

  /**
   * Weighted variant for condensed inputs; currently delegates to unweighted clustering as
   * weights are not applied with condensed input in this implementation.
   *
   * @param distances       a 1 x M condensed distance array (pdist order)
   * @param clusterNames    names of the N items
   * @param weights         weights for the N items (unused)
   * @param linkageStrategy linkage criterion used during agglomeration
   * @return root Cluster of the resulting hierarchy
   */
  @Override
  public Cluster performWeightedClustering(double[][] distances, String[] clusterNames,
                                           double[] weights, LinkageStrategy linkageStrategy) {
    return performClustering(distances, clusterNames, linkageStrategy);
  }

  /**
   * Builds initial linkages from a condensed upper-triangular distance array.
   *
   * @param distances a 1 x M condensed distance array (pdist order)
   * @param clusters  list of initial singleton clusters
   * @return a DistanceMap containing inter-cluster linkages
   */
  private DistanceMap createLinkages(double[][] distances,
                                     List<Cluster> clusters) {
    DistanceMap linkages = new DistanceMap();
    for (int col = 0; col < clusters.size(); col++) {
      Cluster cluster_col = clusters.get(col);
      for (int row = col + 1; row < clusters.size(); row++) {
        ClusterPair link = new ClusterPair();
        Double d = distances[0][accessFunction(row, col,
            clusters.size())];
        link.setLinkageDistance(d);
        link.setlCluster(cluster_col);
        link.setrCluster(clusters.get(row));
        linkages.add(link);
      }
    }
    return linkages;
  }

  /**
   * Creates initial singleton clusters, one per input name.
   *
   * @param clusterNames names for singleton clusters
   * @return list of newly created singleton clusters
   */
  private List<Cluster> createClusters(String[] clusterNames) {
    List<Cluster> clusters = new ArrayList<>();
    for (String clusterName : clusterNames) {
      Cluster cluster = new Cluster(clusterName);
      cluster.addLeafName(clusterName);
      clusters.add(cluster);
    }
    return clusters;
  }

  // Credit to this function goes to
  // http://stackoverflow.com/questions/13079563/how-does-condensed-distance-matrix-work-pdist

  /**
   * Indexing function mapping (i,j) with i>j to the position in the condensed pdist array.
   *
   * @param i row index (0-based)
   * @param j column index (0-based), with i > j
   * @param n the total number of items (matrix dimension)
   * @return index into the condensed upper-triangular storage
   */
  private static int accessFunction(int i, int j, int n) {
    return n * j - j * (j + 1) / 2 + i - 1 - j;
  }

}
