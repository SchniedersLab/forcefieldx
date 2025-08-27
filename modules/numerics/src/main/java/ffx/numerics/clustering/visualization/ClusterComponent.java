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
package ffx.numerics.clustering.visualization;

import ffx.numerics.clustering.Cluster;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

/**
 * Visual component representing a single cluster node within a dendrogram.
 * Responsible for rendering the node, its label, and the connector lines to
 * its parent/children based on virtual coordinates.
 *
 * <p>Used by DendrogramPanel to draw hierarchical clustering results.</p>
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ClusterComponent implements Paintable {

  private Cluster cluster;
  private VCoord linkPoint;
  private VCoord initPoint;
  private boolean printName;
  private int dotRadius = 2;
  private int namePadding = 6;

  private List<ClusterComponent> children;

  /**
   * Returns the child visual components corresponding to child clusters.
   *
   * @return list of child ClusterComponents (lazy-initialized)
   */
  public List<ClusterComponent> getChildren() {
    if (children == null) {
      children = new ArrayList<>();
    }
    return children;
  }

  /**
   * Gets the pixel padding between a leaf node and its name text.
   *
   * @return name padding in pixels
   */
  public int getNamePadding() {
    return namePadding;
  }

  /**
   * Sets the pixel padding between a leaf node and its name text.
   *
   * @param namePadding name padding in pixels
   */
  public void setNamePadding(int namePadding) {
    this.namePadding = namePadding;
  }

  /**
   * Gets the radius of node dots in pixels.
   *
   * @return dot radius in pixels
   */
  public int getDotRadius() {
    return dotRadius;
  }

  /**
   * Sets the radius of node dots in pixels.
   *
   * @param dotRadius dot radius in pixels
   */
  public void setDotRadius(int dotRadius) {
    this.dotRadius = dotRadius;
  }

  /**
   * Sets the list of child visual components.
   *
   * @param children list of child components
   */
  public void setChildren(List<ClusterComponent> children) {
    this.children = children;
  }

  /**
   * Gets the virtual coordinate where this node connects to its parent.
   *
   * @return link point coordinate
   */
  public VCoord getLinkPoint() {
    return linkPoint;
  }

  /**
   * Sets the virtual coordinate where this node connects to its parent.
   *
   * @param linkPoint link point coordinate
   */
  public void setLinkPoint(VCoord linkPoint) {
    this.linkPoint = linkPoint;
  }

  /**
   * Gets the virtual coordinate at which this node is drawn.
   *
   * @return initial coordinate for this node
   */
  public VCoord getInitPoint() {
    return initPoint;
  }

  /**
   * Sets the virtual coordinate at which this node is drawn.
   *
   * @param initPoint initial coordinate for this node
   */
  public void setInitPoint(VCoord initPoint) {
    this.initPoint = initPoint;
  }

  /**
   * Gets the Cluster model represented by this component.
   *
   * @return the associated Cluster
   */
  public Cluster getCluster() {
    return cluster;
  }

  /**
   * Sets the Cluster model represented by this component.
   *
   * @param cluster the Cluster to associate with this component
   */
  public void setCluster(Cluster cluster) {
    this.cluster = cluster;
  }

  /**
   * Returns whether the node name should be drawn.
   *
   * @return true if the name should be drawn
   */
  public boolean isPrintName() {
    return printName;
  }

  /**
   * Sets whether the node name should be drawn.
   *
   * @param printName true to draw the node name
   */
  public void setPrintName(boolean printName) {
    this.printName = printName;
  }

  /**
   * Constructs a visual node component for a Cluster.
   *
   * @param cluster   the cluster represented by this component
   * @param printName whether to render the cluster name
   * @param initPoint the initial virtual coordinate of this node
   */
  public ClusterComponent(Cluster cluster, boolean printName, VCoord initPoint) {
    this.printName = printName;
    this.cluster = cluster;
    this.initPoint = initPoint;
    this.linkPoint = initPoint;
  }

  /**
   * {@inheritDoc}
   * Draws the node dot, horizontal and vertical connectors, and optionally labels.
   */
  @Override
  public void paint(Graphics2D g, int xDisplayOffset, int yDisplayOffset, double xDisplayFactor, double yDisplayFactor, boolean decorated) {
    int x1, y1, x2, y2;
    FontMetrics fontMetrics = g.getFontMetrics();
    x1 = (int) (initPoint.x() * xDisplayFactor + xDisplayOffset);
    y1 = (int) (initPoint.y() * yDisplayFactor + yDisplayOffset);
    x2 = (int) (linkPoint.x() * xDisplayFactor + xDisplayOffset);
    y2 = y1;
    g.fillOval(x1 - dotRadius, y1 - dotRadius, dotRadius * 2, dotRadius * 2);
    g.drawLine(x1, y1, x2, y2);

    if (cluster.isLeaf()) {
      g.drawString(cluster.getName(), x1 + namePadding, y1 + (fontMetrics.getHeight() / 2) - 2);
    }
    if (decorated && cluster.getDistance() != null && !cluster.getDistance().isNaN() && cluster.getDistance().getDistance() > 0) {
      String s = String.format("%.2f", cluster.getDistance().getDistance());
      Rectangle2D rect = fontMetrics.getStringBounds(s, g);
      g.drawString(s, x1 - (int) rect.getWidth(), y1 - 2);
    }

    x1 = x2;
    y1 = y2;
    y2 = (int) (linkPoint.y() * yDisplayFactor + yDisplayOffset);
    g.drawLine(x1, y1, x2, y2);


    for (ClusterComponent child : children) {
      child.paint(g, xDisplayOffset, yDisplayOffset, xDisplayFactor, yDisplayFactor, decorated);
    }
  }

  /**
   * Computes the minimal X value of this component and its children in model space.
   *
   * @return minimal X coordinate across subtree
   */
  public double getRectMinX() {

    // TODO Better use closure / callback here
    assert initPoint != null && linkPoint != null;
    double val = Math.min(initPoint.x(), linkPoint.x());
    for (ClusterComponent child : getChildren()) {
      val = Math.min(val, child.getRectMinX());
    }
    return val;
  }

  /**
   * Computes the minimal Y value of this component and its children in model space.
   *
   * @return minimal Y coordinate across subtree
   */
  public double getRectMinY() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.min(initPoint.y(), linkPoint.y());
    for (ClusterComponent child : getChildren()) {
      val = Math.min(val, child.getRectMinY());
    }
    return val;
  }

  /**
   * Computes the maximal X value of this component and its children in model space.
   *
   * @return maximal X coordinate across subtree
   */
  public double getRectMaxX() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.max(initPoint.x(), linkPoint.x());
    for (ClusterComponent child : getChildren()) {
      val = Math.max(val, child.getRectMaxX());
    }
    return val;
  }

  /**
   * Computes the maximal Y value of this component and its children in model space.
   *
   * @return maximal Y coordinate across subtree
   */
  public double getRectMaxY() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.max(initPoint.y(), linkPoint.y());
    for (ClusterComponent child : getChildren()) {
      val = Math.max(val, child.getRectMaxY());
    }
    return val;
  }

  /**
   * Computes the width in pixels of this node's name label.
   *
   * @param g               graphics context used for font metrics
   * @param includeNonLeafs if true, include internal nodes; otherwise only leaf names
   * @return width in pixels of the label (0 if not drawn)
   */
  public int getNameWidth(Graphics2D g, boolean includeNonLeafs) {
    int width = 0;
    if (includeNonLeafs || cluster.isLeaf()) {
      Rectangle2D rect = g.getFontMetrics().getStringBounds(cluster.getName(), g);
      width = (int) rect.getWidth();
    }
    return width;
  }

  /**
   * Recursively computes the maximal name width across this node and its children.
   *
   * @param g               graphics context used for font metrics
   * @param includeNonLeafs if true, include internal nodes; otherwise only leaf names
   * @return maximum label width in pixels across the subtree
   */
  public int getMaxNameWidth(Graphics2D g, boolean includeNonLeafs) {
    int width = getNameWidth(g, includeNonLeafs);
    for (ClusterComponent comp : getChildren()) {
      int childWidth = comp.getMaxNameWidth(g, includeNonLeafs);
      if (childWidth > width) {
        width = childWidth;
      }
    }
    return width;
  }
}
