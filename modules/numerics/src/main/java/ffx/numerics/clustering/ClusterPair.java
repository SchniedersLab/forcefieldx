// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

/**
 * Immutable-like holder describing a pair of clusters and the linkage distance
 * between them; used as entries within the DistanceMap during agglomeration.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ClusterPair implements Comparable<ClusterPair> {

  private Cluster lCluster;
  private Cluster rCluster;
  private Double linkageDistance;

  /**
   * Creates an empty ClusterPair.
   */
  public ClusterPair() {
  }

  /**
   * Creates a ClusterPair linking two clusters at a given distance.
   *
   * @param left     the left cluster
   * @param right    the right cluster
   * @param distance the linkage distance between clusters
   */
  public ClusterPair(Cluster left, Cluster right, Double distance) {
    lCluster = left;
    rCluster = right;
    linkageDistance = distance;
  }

  /**
   * Returns the opposite cluster of the provided one in this pair.
   *
   * @param c one cluster in this pair
   * @return the other cluster in the pair
   */
  public Cluster getOtherCluster(Cluster c) {
    return lCluster == c ? rCluster : lCluster;
  }

  /**
   * Gets the left cluster.
   *
   * @return the left cluster
   */
  public Cluster getlCluster() {
    return lCluster;
  }

  /**
   * Sets the left cluster.
   *
   * @param lCluster the left cluster
   */
  public void setlCluster(Cluster lCluster) {
    this.lCluster = lCluster;
  }

  /**
   * Gets the right cluster.
   *
   * @return the right cluster
   */
  public Cluster getrCluster() {
    return rCluster;
  }

  /**
   * Sets the right cluster.
   *
   * @param rCluster the right cluster
   */
  public void setrCluster(Cluster rCluster) {
    this.rCluster = rCluster;
  }

  /**
   * Gets the linkage distance.
   *
   * @return the linkage distance
   */
  public Double getLinkageDistance() {
    return linkageDistance;
  }

  /**
   * Sets the linkage distance.
   *
   * @param distance the linkage distance to set
   */
  public void setLinkageDistance(Double distance) {
    this.linkageDistance = distance;
  }

  /**
   * Creates a new ClusterPair with left and right clusters swapped.
   *
   * @return a new ClusterPair with the two left/right inverted
   */
  public ClusterPair reverse() {
    return new ClusterPair(getrCluster(), getlCluster(), getLinkageDistance());
  }


  /**
   * {@inheritDoc}
   */
  @Override
  public int compareTo(ClusterPair o) {
    int result;
    if (o == null || o.getLinkageDistance() == null) {
      result = -1;
    } else if (getLinkageDistance() == null) {
      result = 1;
    } else {
      result = getLinkageDistance().compareTo(o.getLinkageDistance());
    }

    return result;
  }

  /**
   * Agglomerates left and right clusters under an auto-generated name.
   *
   * @param clusterIdx index appended to the generated cluster name
   * @return the new parent Cluster
   */
  public Cluster agglomerate(int clusterIdx) {
    return agglomerate("clstr#" + clusterIdx);
  }

  /**
   * Agglomerates left and right clusters into a new parent with the given name.
   *
   * @param name name of the new parent cluster
   * @return the new parent Cluster
   */
  public Cluster agglomerate(String name) {
    Cluster cluster = new Cluster(name);
    cluster.setDistance(new Distance(getLinkageDistance()));
    //New clusters will track their children's leaf names; i.e. each cluster knows what part of the original data it contains
    cluster.appendLeafNames(lCluster.getLeafNames());
    cluster.appendLeafNames(rCluster.getLeafNames());
    cluster.addChild(lCluster);
    cluster.addChild(rCluster);
    lCluster.setParent(cluster);
    rCluster.setParent(cluster);

    Double lWeight = lCluster.getWeightValue();
    Double rWeight = rCluster.getWeightValue();
    double weight = lWeight + rWeight;
    cluster.getDistance().setWeight(weight);

    return cluster;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    if (lCluster != null) {
      sb.append(lCluster.getName());
    }
    if (rCluster != null) {
      if (!sb.isEmpty()) {
        sb.append(" + ");
      }
      sb.append(rCluster.getName());
    }
    sb.append(" : ").append(linkageDistance);
    return sb.toString();
  }

}
