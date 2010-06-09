/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.xray;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;
import java.util.Arrays;

/**
 * This class implements a spatial decomposition based on partitioning a
 * grid into octants. The over-ridden "selectAtoms" method selects atoms
 * that are not in the asymmetric unit, but are within the supplied cutoff
 * radius.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class BulkSolventDensityRegion extends SpatialDensityRegion {

    /**
     * Atoms within cutoff angstroms of an asymmetric unit atom will
     * be selected.
     */
    private double cutoff2;

    public BulkSolventDensityRegion(int gX, int gY, int gZ, double grid[],
                                    int basisSize, int nSymm, int minWork,
                                    int threadCount, Crystal crystal,
                                    Atom atoms[], double coordinates[][][],
                                    double cutoff) {
        super(gX, gY, gZ, grid, basisSize, nSymm, minWork,
              threadCount, crystal, atoms, coordinates);

        this.cutoff2 = cutoff * cutoff;
        // Asymmetric unit atoms are always not selected by this class.
        Arrays.fill(select[0], false);
    }

    public BulkSolventDensityRegion(int gX, int gY, int gZ, float grid[],
                                    int basisSize, int nSymm, int minWork,
                                    int threadCount, Crystal crystal,
                                    Atom atoms[], double coordinates[][][],
                                    double cutoff) {
        super(gX, gY, gZ, grid, basisSize, nSymm, minWork,
              threadCount, crystal, atoms, coordinates);

        this.cutoff2 = cutoff * cutoff;
        // Asymmetric unit atoms are always not selected by this class.
        Arrays.fill(select[0], false);
    }

       @Override
    public void run() {
        int ti = getThreadIndex();
        int actualWork1 = actualWork - 1;
        SpatialDensityLoop loop = spatialDensityLoop[ti];
        try {
            execute(0, actualWork1, loop.setOctant(0));
            // Fractional chunks along the C-axis.
            if (nC > 1) {
                execute(0, actualWork1, loop.setOctant(1));
                // Fractional chunks along the B-axis.
                if (nB > 1) {
                    execute(0, actualWork1, loop.setOctant(2));
                    execute(0, actualWork1, loop.setOctant(3));
                    // Fractional chunks along the A-axis.
                    if (nA > 1) {
                        execute(0, actualWork1, loop.setOctant(4));
                        execute(0, actualWork1, loop.setOctant(5));
                        execute(0, actualWork1, loop.setOctant(6));
                        execute(0, actualWork1, loop.setOctant(7));
                    }
                }
            }
        } catch (Exception e) {
            logger.severe(e.toString());
        }
    }

    @Override
    public void selectAtoms() {
        final double xyz[][] = coordinates[0];
        final double x[] = xyz[0];
        final double y[] = xyz[1];
        final double z[] = xyz[2];
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
            boolean selected[] = select[iSymm];
            Arrays.fill(selected, false);
            final double xyzs[][] = coordinates[iSymm];
            final double xs[] = xyzs[0];
            final double ys[] = xyzs[1];
            final double zs[] = xyzs[2];
            // loop over symmetry mate atoms
            for (int i = 0; i < nAtoms; i++) {
                double xsi = xs[i];
                double ysi = ys[i];
                double zsi = zs[i];
                /**
                 * Loop over asymmetric unit atoms looking for an interaction
                 * distance within the cutoff.
                 */
                for (int j=0; j < nAtoms; j++) {
                    double xr = xsi - x[j];
                    double yr = ysi - y[j];
                    double zr = zsi - z[j];
                    double r2 = crystal.image(xr,yr,zr);
                    if (r2 < cutoff2) {
                        selected[i] = true;
                        break;
                    }
                }
            }
        }
    }
}
