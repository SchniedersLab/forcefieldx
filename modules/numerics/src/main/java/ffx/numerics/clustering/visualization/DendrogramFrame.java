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
import ffx.numerics.clustering.ClusteringAlgorithm;
import ffx.numerics.clustering.CompleteLinkageStrategy;
import ffx.numerics.clustering.LinkageStrategy;
import ffx.numerics.clustering.PDistClusteringAlgorithm;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Frame;

/**
 * Simple Swing JFrame that hosts a DendrogramPanel to visualize a clustering
 * result. Provides a demo main method to render example dendrograms.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DendrogramFrame extends JFrame {

  /**
   * Creates a frame displaying a dendrogram for the provided clustering result.
   *
   * @param cluster the root Cluster to visualize
   */
  public DendrogramFrame(Cluster cluster) {
    setSize(500, 400);
    setLocation(100, 200);
    setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

    JPanel content = new JPanel();
    DendrogramPanel dp = new DendrogramPanel();

    setContentPane(content);
    content.setBackground(Color.red);
    content.setLayout(new BorderLayout());
    content.add(dp, BorderLayout.CENTER);
    dp.setBackground(Color.WHITE);
    dp.setLineColor(Color.BLACK);
    dp.setScaleValueDecimals(0);
    dp.setScaleValueInterval(1);
    dp.setShowDistances(false);

    dp.setModel(cluster);
    setVisible(true);
  }

  /**
   * Demo entry point that creates two frames with example dendrograms.
   *
   * @param args CLI arguments (unused)
   */
  public static void main(String[] args) {
    LinkageStrategy strategy = new CompleteLinkageStrategy();
    Frame f1 = new DendrogramFrame(createSampleCluster(strategy));
    f1.setSize(500, 400);
    f1.setLocation(100, 200);
    Frame f2 = new DendrogramFrame(createSampleCluster2(strategy));
    f2.setSize(500, 400);
    f2.setLocation(600, 200);
  }

  /**
   * Creates a small sample Cluster for demonstration purposes.
   *
   * @param strategy the LinkageStrategy used by the clustering algorithm
   * @return a sample root Cluster
   */
  private static Cluster createSampleCluster(LinkageStrategy strategy) {
    double[][] distances = new double[][]{
        {1, 9, 7, 11, 14, 4, 3, 8, 10, 9, 2, 8, 6, 13, 10}
    };
    String[] names = new String[]{"O1", "O2", "O3", "O4", "O5", "O6"};
    ClusteringAlgorithm alg = new PDistClusteringAlgorithm();
    Cluster cluster = alg.performClustering(distances, names, strategy);
    cluster.toConsole(0);
    return cluster;
  }

  /**
   * Creates a second sample Cluster for demonstration purposes.
   *
   * @param strategy the LinkageStrategy used by the clustering algorithm
   * @return a sample root Cluster
   */
  private static Cluster createSampleCluster2(LinkageStrategy strategy) {
    double[][] distances = new double[][]{
        {1, 9, 7, 11, 14, 12, 4, 3, 8, 10, 12, 9, 2, 8, 9, 6, 13, 11, 10, 7, 2}
    };
    String[] names = new String[]{"O1", "O2", "O3", "O4", "O5", "O6", "07"};
    ClusteringAlgorithm alg = new PDistClusteringAlgorithm();
    Cluster cluster = alg.performClustering(distances, names, strategy);
    cluster.toConsole(0);
    return cluster;
  }
}
