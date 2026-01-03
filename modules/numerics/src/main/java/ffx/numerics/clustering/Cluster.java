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

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a node in a hierarchical clustering tree (dendrogram).
 * A Cluster may be a leaf (no children) or an internal node with two children.
 * It tracks its name, parent, children, the list of leaf names contained beneath it,
 * and an associated Distance used to store linkage distance and aggregate weight.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Cluster {

  private String name;

  private Cluster parent;

  private List<Cluster> children;

  private final List<String> leafNames;

  private Distance distance = new Distance();


  /**
   * Creates a new Cluster with the provided name.
   *
   * @param name the cluster name
   */
  public Cluster(String name) {
    this.name = name;
    leafNames = new ArrayList<>();
  }

  /**
   * Gets the Distance metadata for this cluster (linkage distance and weight).
   *
   * @return the Distance object
   */
  public Distance getDistance() {
    return distance;
  }

  /**
   * Convenience accessor for this cluster's aggregate weight.
   *
   * @return the weight value stored in the Distance
   */
  public Double getWeightValue() {
    return distance.getWeight();
  }

  /**
   * Convenience accessor for this cluster's linkage distance value.
   *
   * @return the linkage distance, or null if unset
   */
  public Double getDistanceValue() {
    return distance.getDistance();
  }

  /**
   * Sets the Distance metadata for this cluster.
   *
   * @param distance the Distance to set
   */
  public void setDistance(Distance distance) {
    this.distance = distance;
  }

  /**
   * Returns the list of child clusters, creating it lazily if needed.
   *
   * @return mutable list of child clusters (possibly empty)
   */
  public List<Cluster> getChildren() {
    if (children == null) {
      children = new ArrayList<>();
    }

    return children;
  }

  /**
   * Adds a single leaf name contained within this cluster's subtree.
   *
   * @param lname the leaf name to add
   */
  public void addLeafName(String lname) {
    leafNames.add(lname);
  }

  /**
   * Appends a list of leaf names to this cluster's leaf name list.
   *
   * @param lnames list of leaf names to append
   */
  public void appendLeafNames(List<String> lnames) {
    leafNames.addAll(lnames);
  }

  /**
   * Returns the list of leaf names beneath this cluster.
   *
   * @return list of leaf names
   */
  public List<String> getLeafNames() {
    return leafNames;
  }

  /**
   * Sets the list of child clusters for this node.
   *
   * @param children the children to set
   */
  public void setChildren(List<Cluster> children) {
    this.children = children;
  }

  /**
   * Gets the parent cluster of this node, or null if this is the root.
   *
   * @return the parent cluster, or null
   */
  public Cluster getParent() {
    return parent;
  }

  /**
   * Sets the parent cluster of this node.
   *
   * @param parent the parent cluster to set
   */
  public void setParent(Cluster parent) {
    this.parent = parent;
  }


  /**
   * Gets the name of this cluster.
   *
   * @return the cluster name
   */
  public String getName() {
    return name;
  }

  /**
   * Sets the name of this cluster.
   *
   * @param name the name to set
   */
  public void setName(String name) {
    this.name = name;
  }

  /**
   * Adds a child cluster to this node.
   *
   * @param cluster the child to add
   */
  public void addChild(Cluster cluster) {
    getChildren().add(cluster);

  }

  /**
   * Tests whether this cluster has the specified child.
   *
   * @param cluster the child to look for
   * @return true if present; false otherwise
   */
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
      return other.name == null;
    } else return name.equals(other.name);
  }

  @Override
  public int hashCode() {
    return (name == null) ? 0 : name.hashCode();
  }

  /**
   * Returns true if this cluster is a leaf (has no children).
   *
   * @return true if leaf; false otherwise
   */
  public boolean isLeaf() {
    return getChildren().isEmpty();
  }

  /**
   * Counts the number of leaf descendants beneath this cluster.
   *
   * @return number of leaf nodes
   */
  public int countLeafs() {
    return countLeafs(this, 0);
  }

  /**
   * Recursive helper to count leaves under the specified node.
   *
   * @param node  the node to inspect
   * @param count running count of leaves
   * @return updated count including leaves beneath node
   */
  public int countLeafs(Cluster node, int count) {
    if (node.isLeaf()) count++;
    for (Cluster child : node.getChildren()) {
      count += child.countLeafs();
    }
    return count;
  }

  /**
   * Prints this cluster and its subtree to the console with indentation.
   *
   * @param indent number of indentation levels for this node
   */
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

  /**
   * Serializes this cluster subtree into a simple Newick-like string.
   * The first child is annotated with distance, the second with weight.
   *
   * @param indent indentation spaces to include (for readability)
   * @return a Newick-like representation of this subtree
   */
  public String toNewickString(int indent) {
    StringBuilder cdtString = new StringBuilder();
    if (!isLeaf()) cdtString.append("(");

    cdtString.append(" ".repeat(Math.max(0, indent)));


    if (isLeaf()) {
      cdtString.append(getName());
    }

    List<Cluster> children = getChildren();

    boolean firstChild = true;
    for (Cluster child : children) {
      cdtString.append(child.toNewickString(indent));
      String distanceString = distance.getDistance().toString().replace(",", ".");
      String weightString = distance.getWeight().toString().replace(",", ".");
      if (firstChild) cdtString.append(":").append(distanceString).append(",");
      else cdtString.append(":").append(weightString);

      firstChild = false;
    }

    cdtString.append(" ".repeat(Math.max(0, indent)));

    if (!isLeaf()) cdtString.append(")");

    return cdtString.toString();
  }

  /**
   * Computes the cumulative distance down the leftmost branch of this cluster.
   *
   * @return total distance along the first-child path
   */
  public double getTotalDistance() {
    double dist = getDistance() == null ? 0 : getDistance().getDistance();
    if (!getChildren().isEmpty()) {
      dist += children.getFirst().getTotalDistance();
    }
    return dist;

  }

}
