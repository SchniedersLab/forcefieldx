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

import ffx.numerics.clustering.visualization.DendrogramPanel;
import ffx.utilities.FFXTest;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Lars Behnke (lars.behnke@bruker.com)
 */
public class CoordTest extends FFXTest {

    private static Cluster importCluster() throws IOException {
        List<Coord> coords = readCoordinates();

        double[][] distances = new double[coords.size()][coords.size()];
        String[] names = new String[coords.size()];
        for (int row = 0; row < coords.size(); row++) {
            Coord coord1 = coords.get(row);
            for (int col = row+1; col < coords.size(); col++) {
                Coord coord2 = coords.get(col);
                double d = Math.sqrt(Math.pow(coord2.getX()-coord1.getX(), 2)+ Math.pow(coord2.getY()-coord1.getY(), 2));
                distances[row][col] = d;
                distances[col][row] = d;
            }
            names[row] = ""+row;
        }
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new AverageLinkageStrategy());
        return cluster;
    }


    public static void main(String[] args) throws Exception {
        JFrame frame = new JFrame();
        frame.setSize(1024, 768);
        frame.setLocation(400, 300);
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JPanel content = new JPanel();
        DendrogramPanel dp = new DendrogramPanel();

        frame.setContentPane(content);
        content.setBackground(Color.red);
        content.setLayout(new BorderLayout());
        content.add(dp, BorderLayout.CENTER);
        dp.setBackground(Color.WHITE);
        dp.setLineColor(Color.BLACK);
        dp.setScaleValueDecimals(0);
        dp.setScaleValueInterval(1);
        dp.setShowDistances(false);

        Cluster cluster = importCluster();
        dp.setModel(cluster);
        frame.setVisible(true);
    }

    private static List<Coord> readCoordinates() throws IOException {
        List<Coord> coordList = new ArrayList<Coord>();
        BufferedReader br = new BufferedReader(new InputStreamReader(CoordTest.class.getResourceAsStream("/testData1.txt")));
        String line;

        while ((line = br.readLine()) != null) {
            String[] elems = line.split(" ");
            if (elems.length != 2) {
                continue;
            }
            int x;
            int y;

            try {
                x = Integer.parseInt(elems[0]);
                y = Integer.parseInt(elems[1]);
            } catch (Exception e) {
                continue;
            }
            coordList.add(new Coord(x, y));
        }
        return coordList;
    }

    public static class Coord {
        private double x;
        private double y;
        public Coord(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double getX() {
            return x;
        }

        public double getY() {
            return y;
        }
    }
}
