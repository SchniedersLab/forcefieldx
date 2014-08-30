/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.autoparm;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

import ffx.numerics.LBFGS;

/**
 * Test the Potential2 class.
 *
 * Things I want to test: Number of points in the grid Multipole Vals Avg.
 * Potential Before (and other stuff before like avg. unsigned diff) Avg.
 * Potential After (and other stuff after like avg. unsigned diff) Test
 * Phenobarbital and 12-ethanediol for now
 *
 * @author Gaurav Chattree
 *
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class Potential2Test {

    private static final Logger logger = Logger.getLogger(Potential2Test.class.getName());
    private final double tolerance = 1.0e-5;

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                "12-Ethanediol Benchmark",
                "ffx/autoparm/structures/12-ethanediol/12-ethanediol.xyz",
                .1,
                7471,
                new double[]{0.15787932025820384, 0.0, 0.04763753450168774, 0.08730356716095355, -0.08438672145153087,
                    -0.0029168457094226885, 0.0, -0.14288948649965377, 0.0, 0.0070311494580525165, 0.0, -0.010208713719432605,
                    -0.01893339298739448, -0.029403617130171186, 0.048337010117565665, 0.0, 0.015200446852008324, 0.0,
                    0.08836068958950538, 0.0, 0.19309961014679017, 0.039029095472027965, -0.2620911893802569, 0.22306209390822893,
                    0.0, -0.13352971381310716, 0.0, 0.0, 0.0, -0.00769802568195431, -0.01582751526374996, -0.007987655242155419,
                    0.02381367050590538, 0.0, 0.0, 0.0},
                new double[]{5.08619, 5.31849, -0.0774428, 0.881366, 1.09829},
                new double[]{5.30394, 5.31849, -0.0252355, 0.0467198, 0.0714136}
            }
        });
    }
    private int npoints;
    private double mpoles[];
    private double potstats_before[];
    private double potstats_after[];
    private double stats_before[];
    private double stats_after[];
    private Potential2 p;
    private String info;

    @SuppressWarnings("LeakingThisInConstructor")
    public Potential2Test(String info, String xyz_filename, double eps, int npoints, double[] mpoles, double[] potstats_before, double[] potstats_after) {
        //gen_pot_grid - gives npoints
        //output_stats
        //Run the optimizer
        //Get mpole data
        //get output_stats data

        this.npoints = npoints;
        this.mpoles = mpoles;
        this.potstats_before = potstats_before;
        this.potstats_after = potstats_after;
        this.info = info;

        ClassLoader cl = this.getClass().getClassLoader();
        String xyzfname = cl.getResource(xyz_filename).getPath();
        System.setProperty("polarization", "mutual");
        try {

            p = new Potential2(0, xyzfname, null, null);
            p.target_grid = p.gen_pot_grid(p.structure_cube, p.atoms, 0);
            p.pme.set_target_grid(p.target_grid);
            p.output_stats();
            stats_before = p.stats.clone();
            int m = 7;
            p.nvars = p.pme.getNumberOfVariables();
            p.x = new double[p.nvars];
            p.grad = new double[p.nvars];
            p.scaling = new double[p.nvars];
            p.pme.getCoordinates(p.x);
            double e = p.pme.energyAndGradient(p.x, p.grad);
            long at = System.nanoTime();
            int status = LBFGS.minimize(p.nvars, m, p.x, e, p.grad, eps, p.pme, p);
            long bt = System.nanoTime();
            System.out.println("TIME FOR LBFGS: " + (bt - at) * 1e-9 + " seconds\n");
            switch (status) {
                case 0:
                    logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", p.grms));
                    break;
                case 1:
                    logger.info(String.format("\n Optimization terminated at step %d.\n", p.nSteps));
                    break;
                default:
                    logger.warning("\n Optimization failed.\n");
            }
            p.pme.init_prms();
            p.output_stats();
            stats_after = p.stats.clone();
        } catch (IOException ex) {
            Logger.getLogger(Potential2Test.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    @Test
    public void testPotential() {

        assertEquals(info, npoints, p.target_grid[0].length);
        for (int i = 0; i < stats_before.length; i++) {
            assertEquals(info, stats_before[i], potstats_before[i], tolerance);
        }
        for (int i = 0; i < stats_after.length; i++) {
            assertEquals(info, stats_after[i], potstats_after[i], tolerance);
        }
        double[] m = p.pme.getallmpoles(p.x);
        for (int i = 0; i < mpoles.length; i++) {
            assertEquals(info, mpoles[i], m[i], tolerance);
        }
    }
}
