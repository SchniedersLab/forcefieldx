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

import javax.swing.JPanel;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;

public class DendrogramPanel extends JPanel {

  private static final long serialVersionUID = 1L;

  private final static BasicStroke SOLID_STROKE =
      new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);

  private Cluster model;
  private ClusterComponent component;
  private Color lineColor = Color.BLACK;
  private boolean showDistanceValues = false;
  private boolean showScale = true;
  private int borderTop = 20;
  private int borderLeft = 20;
  private int borderRight = 20;
  private int borderBottom = 20;
  private int scalePadding = 10;
  private int scaleTickLength = 4;
  private int scaleTickLabelPadding = 4;
  private double scaleValueInterval = 0;
  private int scaleValueDecimals = 0;

  private double xModelOrigin = 0.0;
  private double yModelOrigin = 0.0;
  private double wModel = 0.0;
  private double hModel = 0.0;

  public boolean isShowDistanceValues() {
    return showDistanceValues;
  }

  public void setShowDistances(boolean showDistanceValues) {
    this.showDistanceValues = showDistanceValues;
  }

  public boolean isShowScale() {
    return showScale;
  }

  public void setShowScale(boolean showScale) {
    this.showScale = showScale;
  }

  public int getScalePadding() {
    return scalePadding;
  }

  public void setScalePadding(int scalePadding) {
    this.scalePadding = scalePadding;
  }

  public int getScaleTickLength() {
    return scaleTickLength;
  }

  public void setScaleTickLength(int scaleTickLength) {
    this.scaleTickLength = scaleTickLength;
  }

  public double getScaleValueInterval() {
    return scaleValueInterval;
  }

  public void setScaleValueInterval(double scaleTickInterval) {
    this.scaleValueInterval = scaleTickInterval;
  }

  public int getScaleValueDecimals() {
    return scaleValueDecimals;
  }

  public void setScaleValueDecimals(int scaleValueDecimals) {
    this.scaleValueDecimals = scaleValueDecimals;
  }

  public int getBorderTop() {
    return borderTop;
  }

  public void setBorderTop(int borderTop) {
    this.borderTop = borderTop;
  }

  public int getBorderLeft() {
    return borderLeft;
  }

  public void setBorderLeft(int borderLeft) {
    this.borderLeft = borderLeft;
  }

  public int getBorderRight() {
    return borderRight;
  }

  public void setBorderRight(int borderRight) {
    this.borderRight = borderRight;
  }

  public int getBorderBottom() {
    return borderBottom;
  }

  public void setBorderBottom(int borderBottom) {
    this.borderBottom = borderBottom;
  }

  public Color getLineColor() {
    return lineColor;
  }

  public void setLineColor(Color lineColor) {
    this.lineColor = lineColor;
  }

  public Cluster getModel() {
    return model;
  }

  public void setModel(Cluster model) {
    this.model = model;
    component = createComponent(model);
    updateModelMetrics();
  }

  private void updateModelMetrics() {
    double minX = component.getRectMinX();
    double maxX = component.getRectMaxX();
    double minY = component.getRectMinY();
    double maxY = component.getRectMaxY();

    xModelOrigin = minX;
    yModelOrigin = minY;
    wModel = maxX - minX;
    hModel = maxY - minY;
  }

  private ClusterComponent createComponent(Cluster cluster, VCoord initCoord, double clusterHeight) {

    ClusterComponent comp = null;
    if (cluster != null) {
      comp = new ClusterComponent(cluster, cluster.isLeaf(), initCoord);
      double leafHeight = clusterHeight / cluster.countLeafs();
      double yChild = initCoord.getY() - (clusterHeight / 2);
      double distance = cluster.getDistanceValue() == null ? 0 : cluster.getDistanceValue();
      for (Cluster child : cluster.getChildren()) {
        int childLeafCount = child.countLeafs();
        double childHeight = childLeafCount * leafHeight;
        double childDistance = child.getDistanceValue() == null ? 0 : child.getDistanceValue();
        VCoord childInitCoord = new VCoord(
            initCoord.getX() + (distance - childDistance),
            yChild + childHeight / 2.0);
        yChild += childHeight;

        /* Traverse cluster node tree */
        ClusterComponent childComp = createComponent(child, childInitCoord, childHeight);

        childComp.setLinkPoint(initCoord);
        comp.getChildren().add(childComp);
      }
    }
    return comp;

  }

  private ClusterComponent createComponent(Cluster model) {

    double virtualModelHeight = 1;
    VCoord initCoord = new VCoord(0, virtualModelHeight / 2);

    ClusterComponent comp = createComponent(model, initCoord, virtualModelHeight);
    comp.setLinkPoint(initCoord);
    return comp;
  }

  @Override
  public void paint(Graphics g) {
    super.paint(g);
    Graphics2D g2 = (Graphics2D) g;
    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    g2.setColor(lineColor);
    g2.setStroke(SOLID_STROKE);

    int wDisplay = getWidth() - borderLeft - borderRight;
    int hDisplay = getHeight() - borderTop - borderBottom;
    int xDisplayOrigin = borderLeft;
    int yDisplayOrigin = borderBottom;

    if (component != null) {

      int nameGutterWidth = component.getMaxNameWidth(g2, false) + component.getNamePadding();
      wDisplay -= nameGutterWidth;

      if (showScale) {
        Rectangle2D rect = g2.getFontMetrics().getStringBounds("0", g2);
        int scaleHeight = (int) rect.getHeight() + scalePadding + scaleTickLength + scaleTickLabelPadding;
        hDisplay -= scaleHeight;
        yDisplayOrigin += scaleHeight;
      }

      /* Calculate conversion factor and offset for display */
      double xFactor = wDisplay / wModel;
      double yFactor = hDisplay / hModel;
      int xOffset = (int) (xDisplayOrigin - xModelOrigin * xFactor);
      int yOffset = (int) (yDisplayOrigin - yModelOrigin * yFactor);
      component.paint(g2, xOffset, yOffset, xFactor, yFactor, showDistanceValues);

      if (showScale) {
        int x1 = xDisplayOrigin;
        int y1 = yDisplayOrigin - scalePadding;
        int x2 = x1 + wDisplay;
        int y2 = y1;
        g2.drawLine(x1, y1, x2, y2);

        double totalDistance = component.getCluster().getTotalDistance();
        double xModelInterval;
        if (scaleValueInterval <= 0) {
          xModelInterval = totalDistance / 10.0;
        } else {
          xModelInterval = scaleValueInterval;
        }

        int xTick = xDisplayOrigin + wDisplay;
        y1 = yDisplayOrigin - scalePadding;
        y2 = yDisplayOrigin - scalePadding - scaleTickLength;
        double distanceValue = 0;
        double xDisplayInterval = xModelInterval * xFactor;
        while (xTick >= xDisplayOrigin) {
          g2.drawLine(xTick, y1, xTick, y2);

          String distanceValueStr = String.format("%." + scaleValueDecimals + "f", distanceValue);
          Rectangle2D rect = g2.getFontMetrics().getStringBounds(distanceValueStr, g2);
          g2.drawString(distanceValueStr, (int) (xTick - (rect.getWidth() / 2)), y2 - scaleTickLabelPadding);
          xTick -= xDisplayInterval;
          distanceValue += xModelInterval;
        }

      }
    } else {

      /* No data available */
      String str = "No data";
      Rectangle2D rect = g2.getFontMetrics().getStringBounds(str, g2);
      int xt = (int) (wDisplay / 2.0 - rect.getWidth() / 2.0);
      int yt = (int) (hDisplay / 2.0 - rect.getHeight() / 2.0);
      g2.drawString(str, xt, yt);
    }
  }
}
