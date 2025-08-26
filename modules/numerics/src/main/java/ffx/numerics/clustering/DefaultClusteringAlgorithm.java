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
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class DefaultClusteringAlgorithm implements ClusteringAlgorithm {

  @Override
  public Cluster performClustering(double[][] distances,
                                   String[] clusterNames, LinkageStrategy linkageStrategy) {

    checkArguments(distances, clusterNames, linkageStrategy);
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

  @Override
  public List<Cluster> performFlatClustering(double[][] distances,
                                             String[] clusterNames, LinkageStrategy linkageStrategy, Double threshold) {

    checkArguments(distances, clusterNames, linkageStrategy);
    /* Setup model */
    List<Cluster> clusters = createClusters(clusterNames);
    DistanceMap linkages = createLinkages(distances, clusters);

    /* Process */
    HierarchyBuilder builder = new HierarchyBuilder(clusters, linkages);
    return builder.flatAgg(linkageStrategy, threshold);
  }

  private void checkArguments(double[][] distances, String[] clusterNames,
                              LinkageStrategy linkageStrategy) {
    if (distances == null || distances.length == 0
        || distances[0].length != distances.length) {
      throw new IllegalArgumentException("Invalid distance matrix");
    }
    if (distances.length != clusterNames.length) {
      throw new IllegalArgumentException("Invalid cluster name array");
    }
    if (linkageStrategy == null) {
      throw new IllegalArgumentException("Undefined linkage strategy");
    }
    int uniqueCount = new HashSet<String>(Arrays.asList(clusterNames)).size();
    if (uniqueCount != clusterNames.length) {
      throw new IllegalArgumentException("Duplicate names");
    }
  }

  @Override
  public Cluster performWeightedClustering(double[][] distances, String[] clusterNames,
                                           double[] weights, LinkageStrategy linkageStrategy) {

    checkArguments(distances, clusterNames, linkageStrategy);

    if (weights.length != clusterNames.length) {
      throw new IllegalArgumentException("Invalid weights array");
    }

    /* Setup model */
    List<Cluster> clusters = createClusters(clusterNames, weights);
    DistanceMap linkages = createLinkages(distances, clusters);

    /* Process */
    HierarchyBuilder builder = new HierarchyBuilder(clusters, linkages);
    while (!builder.isTreeComplete()) {
      builder.agglomerate(linkageStrategy);
    }

    return builder.getRootCluster();
  }

  private DistanceMap createLinkages(double[][] distances,
                                     List<Cluster> clusters) {
    DistanceMap linkages = new DistanceMap();
    for (int col = 0; col < clusters.size(); col++) {
      for (int row = col + 1; row < clusters.size(); row++) {
        ClusterPair link = new ClusterPair();
        Cluster lCluster = clusters.get(col);
        Cluster rCluster = clusters.get(row);
        link.setLinkageDistance(distances[col][row]);
        link.setlCluster(lCluster);
        link.setrCluster(rCluster);
        linkages.add(link);
      }
    }
    return linkages;
  }

  private List<Cluster> createClusters(String[] clusterNames) {
    List<Cluster> clusters = new ArrayList<Cluster>();
    for (String clusterName : clusterNames) {
      Cluster cluster = new Cluster(clusterName);
      cluster.addLeafName(clusterName);
      clusters.add(cluster);
    }
    return clusters;
  }

  private List<Cluster> createClusters(String[] clusterNames, double[] weights) {
    List<Cluster> clusters = new ArrayList<Cluster>();
    for (int i = 0; i < weights.length; i++) {
      Cluster cluster = new Cluster(clusterNames[i]);
      cluster.setDistance(new Distance(0.0, weights[i]));
      clusters.add(cluster);
    }
    return clusters;
  }

}
