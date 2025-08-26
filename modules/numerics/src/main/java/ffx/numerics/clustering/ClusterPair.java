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

public class ClusterPair implements Comparable<ClusterPair> {

  private Cluster lCluster;
  private Cluster rCluster;
  private Double linkageDistance;

  public ClusterPair() {
  }

  public ClusterPair(Cluster left, Cluster right, Double distance) {
    lCluster = left;
    rCluster = right;
    linkageDistance = distance;
  }

  public Cluster getOtherCluster(Cluster c) {
    return lCluster == c ? rCluster : lCluster;
  }

  public Cluster getlCluster() {
    return lCluster;
  }

  public void setlCluster(Cluster lCluster) {
    this.lCluster = lCluster;
  }

  public Cluster getrCluster() {
    return rCluster;
  }

  public void setrCluster(Cluster rCluster) {
    this.rCluster = rCluster;
  }

  public Double getLinkageDistance() {
    return linkageDistance;
  }

  public void setLinkageDistance(Double distance) {
    this.linkageDistance = distance;
  }

  /**
   * @return a new ClusterPair with the two left/right inverted
   */
  public ClusterPair reverse() {
    return new ClusterPair(getrCluster(), getlCluster(), getLinkageDistance());
  }


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

  public Cluster agglomerate(int clusterIdx) {
    return agglomerate("clstr#" + clusterIdx);
  }

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

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    if (lCluster != null) {
      sb.append(lCluster.getName());
    }
    if (rCluster != null) {
      if (sb.length() > 0) {
        sb.append(" + ");
      }
      sb.append(rCluster.getName());
    }
    sb.append(" : ").append(linkageDistance);
    return sb.toString();
  }

}
