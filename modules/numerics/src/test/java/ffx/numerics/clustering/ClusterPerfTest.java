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

import ffx.utilities.FFXTest;
import org.junit.Test;

import java.util.Random;

import static java.lang.String.format;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class ClusterPerfTest extends FFXTest {

    Cluster randomCluster(int n) {
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();

        return alg.performClustering(randomDataDist(n), randomDataNames(n),
                new AverageLinkageStrategy());
    }

    private double[][] randomDataDist(int n) {
        Random rnd = new Random();
        double[][] mat = new double[n][n];
        for (int i = 0; i < n; i++) {
            mat[i][i] = 0;
            for (int j = i + 1; j < n; j++) {
                double r = Math.floor(rnd.nextDouble() * 100) * 0.1;
                mat[i][j] = r;
                mat[j][i] = r;
            }
        }

        return mat;
    }

    private String[] randomDataNames(int n) {
        String[] ret = new String[n];
        for (int i = 0; i < n; i++) {
            ret[i] = "" + i;
        }
        return ret;
    }

    @Test
    public void testRandomDataDist() throws Exception {
        double[][] dist = randomDataDist(4);
        assertEquals(dist.length, 4);
    }

    @Test
    public void testRandomDataNames() {
        String[] names = randomDataNames(4);
        assertEquals(names.length, 4);
        String[] exp = {"0", "1", "2", "3"};
        assertArrayEquals(names, exp);
    }

    private Long timeN(int n) {
        Long t0 = System.currentTimeMillis();
        Cluster cluster = randomCluster(n);
        return System.currentTimeMillis() - t0;
    }

    @Test
    public void testn() {
        for (int n = 1; n <= 1024; n = n * 2) {
            Long t = timeN(n);
            logger.info(format("%3d nodes -> %5d ms", n, t));
        }
    }
}
