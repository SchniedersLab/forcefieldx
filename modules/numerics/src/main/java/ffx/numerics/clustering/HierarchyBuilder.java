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
import java.util.Collection;
import java.util.List;

/**
 * Performs agglomerative steps to build a clustering hierarchy from an initial set
 * of singleton clusters and a map of pairwise distances.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class HierarchyBuilder {

  private final DistanceMap distances;
  private final List<Cluster> clusters;
  private int globalClusterIndex = 0;

  /**
   * Gets the DistanceMap used to track inter-cluster distances during agglomeration.
   *
   * @return the DistanceMap backing this builder
   */
  public DistanceMap getDistances() {
    return distances;
  }

  /**
   * Returns the current working list of clusters (not necessarily a single root).
   *
   * @return the list of current clusters
   */
  public List<Cluster> getClusters() {
    return clusters;
  }

  /**
   * Constructs a HierarchyBuilder with an initial set of clusters and inter-cluster distances.
   *
   * @param clusters  initial clusters (typically singletons)
   * @param distances map of inter-cluster distances
   */
  public HierarchyBuilder(List<Cluster> clusters, DistanceMap distances) {
    this.clusters = clusters;
    this.distances = distances;
  }

  /**
   * Performs agglomeration until the minimal inter-cluster distance exceeds the threshold,
   * and returns the remaining clusters (flat clustering at that cut).
   *
   * @param linkageStrategy linkage strategy to compute inter-cluster distances
   * @param threshold       maximum allowed linkage distance for merging
   * @return flat list of clusters remaining at the specified threshold
   */
  public List<Cluster> flatAgg(LinkageStrategy linkageStrategy, Double threshold) {
    while ((!isTreeComplete()) && (distances.minDist() != null) && (distances.minDist() <= threshold)) {
      agglomerate(linkageStrategy);
    }
    return clusters;
  }

  /**
   * Performs one agglomerative step by merging the two closest clusters and updating linkages.
   *
   * @param linkageStrategy strategy to compute new distances to the merged cluster
   */
  public void agglomerate(LinkageStrategy linkageStrategy) {
    ClusterPair minDistLink = distances.removeFirst();
    if (minDistLink != null) {
      clusters.remove(minDistLink.getrCluster());
      clusters.remove(minDistLink.getlCluster());

      Cluster oldClusterL = minDistLink.getlCluster();
      Cluster oldClusterR = minDistLink.getrCluster();
      Cluster newCluster = minDistLink.agglomerate(++globalClusterIndex);

      for (Cluster iClust : clusters) {
        ClusterPair link1 = findByClusters(iClust, oldClusterL);
        ClusterPair link2 = findByClusters(iClust, oldClusterR);
        ClusterPair newLinkage = new ClusterPair();
        newLinkage.setlCluster(iClust);
        newLinkage.setrCluster(newCluster);
        Collection<Distance> distanceValues = new ArrayList<>();

        if (link1 != null) {
          Double distVal = link1.getLinkageDistance();
          Double weightVal = link1.getOtherCluster(iClust).getWeightValue();
          distanceValues.add(new Distance(distVal, weightVal));
          distances.remove(link1);
        }
        if (link2 != null) {
          Double distVal = link2.getLinkageDistance();
          Double weightVal = link2.getOtherCluster(iClust).getWeightValue();
          distanceValues.add(new Distance(distVal, weightVal));
          distances.remove(link2);
        }

        Distance newDistance = linkageStrategy.calculateDistance(distanceValues);

        newLinkage.setLinkageDistance(newDistance.getDistance());
        distances.add(newLinkage);
      }
      clusters.add(newCluster);
    }
  }

  private ClusterPair findByClusters(Cluster c1, Cluster c2) {
    return distances.findByCodePair(c1, c2);
  }

  /**
   * Returns true if only a single cluster remains (i.e., the hierarchy has a root).
   *
   * @return true if a single root cluster remains; false otherwise
   */
  public boolean isTreeComplete() {
    return clusters.size() == 1;
  }

  /**
   * Returns the root cluster if the hierarchy is complete.
   *
   * @return the single remaining root cluster
   * @throws RuntimeException if the tree is not complete
   */
  public Cluster getRootCluster() {
    if (!isTreeComplete()) {
      throw new RuntimeException("No root available");
    }
    return clusters.getFirst();
  }

}
