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

public class HierarchyBuilder {

  private DistanceMap distances;
  private List<Cluster> clusters;
  private int globalClusterIndex = 0;

  public DistanceMap getDistances() {
    return distances;
  }

  public List<Cluster> getClusters() {
    return clusters;
  }

  public HierarchyBuilder(List<Cluster> clusters, DistanceMap distances) {
    this.clusters = clusters;
    this.distances = distances;
  }

  /**
   * Returns Flattened clusters, i.e. clusters that are at least apart by a given threshold
   *
   * @param linkageStrategy
   * @param threshold
   * @return flat list of clusters
   */
  public List<Cluster> flatAgg(LinkageStrategy linkageStrategy, Double threshold) {
    while ((!isTreeComplete()) && (distances.minDist() != null) && (distances.minDist() <= threshold)) {
      //System.out.println("Cluster Distances: " + distances.toString());
      //System.out.println("Cluster Size: " + clusters.size());
      agglomerate(linkageStrategy);
    }

    //System.out.println("Final MinDistance: " + distances.minDist());
    //System.out.println("Tree complete: " + isTreeComplete());
    return clusters;
  }

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
        Collection<Distance> distanceValues = new ArrayList<Distance>();

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

  public boolean isTreeComplete() {
    return clusters.size() == 1;
  }

  public Cluster getRootCluster() {
    if (!isTreeComplete()) {
      throw new RuntimeException("No root available");
    }
    return clusters.get(0);
  }

}
