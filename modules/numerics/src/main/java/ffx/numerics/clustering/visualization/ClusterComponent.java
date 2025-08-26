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

public class ClusterComponent implements Paintable {

  private Cluster cluster;
  private VCoord linkPoint;
  private VCoord initPoint;
  private boolean printName;
  private int dotRadius = 2;
  private int namePadding = 6;

  private List<ClusterComponent> children;

  public List<ClusterComponent> getChildren() {
    if (children == null) {
      children = new ArrayList<ClusterComponent>();
    }
    return children;
  }

  public int getNamePadding() {
    return namePadding;
  }

  public void setNamePadding(int namePadding) {
    this.namePadding = namePadding;
  }

  public int getDotRadius() {
    return dotRadius;
  }

  public void setDotRadius(int dotRadius) {
    this.dotRadius = dotRadius;
  }

  public void setChildren(List<ClusterComponent> children) {
    this.children = children;
  }

  public VCoord getLinkPoint() {
    return linkPoint;
  }

  public void setLinkPoint(VCoord linkPoint) {
    this.linkPoint = linkPoint;
  }

  public VCoord getInitPoint() {
    return initPoint;
  }

  public void setInitPoint(VCoord initPoint) {
    this.initPoint = initPoint;
  }

  public Cluster getCluster() {
    return cluster;
  }

  public void setCluster(Cluster cluster) {
    this.cluster = cluster;
  }

  public boolean isPrintName() {
    return printName;
  }

  public void setPrintName(boolean printName) {
    this.printName = printName;
  }

  public ClusterComponent(Cluster cluster, boolean printName, VCoord initPoint) {
    this.printName = printName;
    this.cluster = cluster;
    this.initPoint = initPoint;
    this.linkPoint = initPoint;
  }

  @Override
  public void paint(Graphics2D g, int xDisplayOffset, int yDisplayOffset, double xDisplayFactor, double yDisplayFactor, boolean decorated) {
    int x1, y1, x2, y2;
    FontMetrics fontMetrics = g.getFontMetrics();
    x1 = (int) (initPoint.getX() * xDisplayFactor + xDisplayOffset);
    y1 = (int) (initPoint.getY() * yDisplayFactor + yDisplayOffset);
    x2 = (int) (linkPoint.getX() * xDisplayFactor + xDisplayOffset);
    y2 = y1;
    g.fillOval(x1 - dotRadius, y1 - dotRadius, dotRadius * 2, dotRadius * 2);
    g.drawLine(x1, y1, x2, y2);

    if (cluster.isLeaf()) {
      g.drawString(cluster.getName(), x1 + namePadding, y1 + (fontMetrics.getHeight() / 2) - 2);
    }
    if (decorated && cluster.getDistance() != null && !cluster.getDistance().isNaN() && cluster.getDistance().getDistance() > 0) {
      String s = String.format("%.2f", cluster.getDistance());
      Rectangle2D rect = fontMetrics.getStringBounds(s, g);
      g.drawString(s, x1 - (int) rect.getWidth(), y1 - 2);
    }

    x1 = x2;
    y1 = y2;
    y2 = (int) (linkPoint.getY() * yDisplayFactor + yDisplayOffset);
    g.drawLine(x1, y1, x2, y2);


    for (ClusterComponent child : children) {
      child.paint(g, xDisplayOffset, yDisplayOffset, xDisplayFactor, yDisplayFactor, decorated);
    }
  }

  public double getRectMinX() {

    // TODO Better use closure / callback here
    assert initPoint != null && linkPoint != null;
    double val = Math.min(initPoint.getX(), linkPoint.getX());
    for (ClusterComponent child : getChildren()) {
      val = Math.min(val, child.getRectMinX());
    }
    return val;
  }

  public double getRectMinY() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.min(initPoint.getY(), linkPoint.getY());
    for (ClusterComponent child : getChildren()) {
      val = Math.min(val, child.getRectMinY());
    }
    return val;
  }

  public double getRectMaxX() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.max(initPoint.getX(), linkPoint.getX());
    for (ClusterComponent child : getChildren()) {
      val = Math.max(val, child.getRectMaxX());
    }
    return val;
  }

  public double getRectMaxY() {

    // TODO Better use closure here
    assert initPoint != null && linkPoint != null;
    double val = Math.max(initPoint.getY(), linkPoint.getY());
    for (ClusterComponent child : getChildren()) {
      val = Math.max(val, child.getRectMaxY());
    }
    return val;
  }

  public int getNameWidth(Graphics2D g, boolean includeNonLeafs) {
    int width = 0;
    if (includeNonLeafs || cluster.isLeaf()) {
      Rectangle2D rect = g.getFontMetrics().getStringBounds(cluster.getName(), g);
      width = (int) rect.getWidth();
    }
    return width;
  }

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
