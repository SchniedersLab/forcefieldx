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
package ffx.potential.nonbonded;

import edu.rit.pj.IntegerForLoop;
import java.util.logging.Logger;

/**
 * Loop over a list of atoms and assign their density to a grid.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public abstract class SpatialDensityLoop extends IntegerForLoop {

    private static final Logger logger = Logger.getLogger(SpatialDensityLoop.class.getName());

    private int nSymm;

    private final SpatialDensityRegion region;
    private int octant = 0;

    public SpatialDensityLoop(SpatialDensityRegion region, int Nsymm) {
        this.region = region;
        this.nSymm = nSymm;
        assert (nSymm <= region.nSymm);
    }

    public void setNsymm(int nSymm) {
        this.nSymm = nSymm;
        assert(nSymm <= region.nSymm);
    }

    public SpatialDensityLoop setOctant(int octant) {
        this.octant = octant;
        return this;
    }

    @Override
    public void run(int lb, int ub) {
        // Loop over work cells
        for (int icell = lb; icell <= ub; icell++) {
            int ia = region.workA[icell];
            int ib = region.workB[icell];
            int ic = region.workC[icell];
            switch (octant) {
                // Case 0 -> In place.
                case 0:
                    gridCell(ia, ib, ic);
                    break;
                // Case 1: Step along the C-axis.
                case 1:
                    gridCell(ia, ib, ic + 1);
                    break;
                // Case 2 & 3: Step along the B-axis.
                case 2:
                    gridCell(ia, ib + 1, ic);
                    break;
                case 3:
                    gridCell(ia, ib + 1, ic + 1);
                    break;
                // Case 4-7: Step along the A-axis.
                case 4:
                    gridCell(ia + 1, ib, ic);
                    break;
                case 5:
                    gridCell(ia + 1, ib, ic + 1);
                    break;
                case 6:
                    gridCell(ia + 1, ib + 1, ic);
                    break;
                case 7:
                    gridCell(ia + 1, ib + 1, ic + 1);
                    break;
                default:
                    String message = "Programming error in PermanentDensityLoop.\n";
                    logger.severe(message);
            }
        }
    }

    private void gridCell(int ia, int ib, int ic) {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int pairList[] = region.cellList[iSymm];
            final int index = region.index(ia, ib, ic);
            final int start = region.cellStart[iSymm][index];
            final int stop = start + region.cellCount[iSymm][index];
            for (int i = start; i < stop; i++) {
                int n = pairList[i];
                gridDensity(iSymm, n);
            }
        }
    }

    public abstract void gridDensity(int iSymm, int iAtom);

}
