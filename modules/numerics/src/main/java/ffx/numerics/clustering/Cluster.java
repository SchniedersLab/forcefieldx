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

public class Cluster {

  private String name;

  private Cluster parent;

  private List<Cluster> children;

  private List<String> leafNames;

  private Distance distance = new Distance();


  public Cluster(String name) {
    this.name = name;
    leafNames = new ArrayList<String>();
  }

  public Distance getDistance() {
    return distance;
  }

  public Double getWeightValue() {
    return distance.getWeight();
  }

  public Double getDistanceValue() {
    return distance.getDistance();
  }

  public void setDistance(Distance distance) {
    this.distance = distance;
  }

  public List<Cluster> getChildren() {
    if (children == null) {
      children = new ArrayList<Cluster>();
    }

    return children;
  }

  public void addLeafName(String lname) {
    leafNames.add(lname);
  }

  public void appendLeafNames(List<String> lnames) {
    leafNames.addAll(lnames);
  }

  public List<String> getLeafNames() {
    return leafNames;
  }

  public void setChildren(List<Cluster> children) {
    this.children = children;
  }

  public Cluster getParent() {
    return parent;
  }

  public void setParent(Cluster parent) {
    this.parent = parent;
  }


  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public void addChild(Cluster cluster) {
    getChildren().add(cluster);

  }

  public boolean contains(Cluster cluster) {
    return getChildren().contains(cluster);
  }

  @Override
  public String toString() {
    return "Cluster " + name;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    Cluster other = (Cluster) obj;
    if (name == null) {
      if (other.name != null) {
        return false;
      }
    } else if (!name.equals(other.name)) {
      return false;
    }
    return true;
  }

  @Override
  public int hashCode() {
    return (name == null) ? 0 : name.hashCode();
  }

  public boolean isLeaf() {
    return getChildren().size() == 0;
  }

  public int countLeafs() {
    return countLeafs(this, 0);
  }

  public int countLeafs(Cluster node, int count) {
    if (node.isLeaf()) count++;
    for (Cluster child : node.getChildren()) {
      count += child.countLeafs();
    }
    return count;
  }

  public void toConsole(int indent) {
    for (int i = 0; i < indent; i++) {
      System.out.print("  ");

    }
    String name = getName() + (isLeaf() ? " (leaf)" : "") + (distance != null ? "  distance: " + distance : "");
    System.out.println(name);
    for (Cluster child : getChildren()) {
      child.toConsole(indent + 1);
    }
  }

  public String toNewickString(int indent) {
    String cdtString = "";
    if (!isLeaf()) cdtString += "(";

    for (int i = 0; i < indent; i++) cdtString += " ";


    if (isLeaf()) {
      cdtString += getName();
    }

    List<Cluster> children = getChildren();

    boolean firstChild = true;
    for (Cluster child : children) {
      cdtString += child.toNewickString(indent);
      String distanceString = distance.getDistance().toString().replace(",", ".");
      String weightString = distance.getWeight().toString().replace(",", ".");
      if (firstChild) cdtString += ":" + distanceString + ",";
      else cdtString += ":" + weightString;

      firstChild = false;
    }

    for (int i = 0; i < indent; i++)
      cdtString += " ";

    if (!isLeaf()) cdtString += ")";

    return cdtString;
  }

  public double getTotalDistance() {
    Double dist = getDistance() == null ? 0 : getDistance().getDistance();
    if (getChildren().size() > 0) {
      dist += children.get(0).getTotalDistance();
    }
    return dist;

  }

}
