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
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;

/**
 * @author Armin Avdic
 */
public abstract class RowLoop extends IntegerForLoop {

    /**
     * Constant <code>logger</code>
     */
    private static final Logger logger = Logger.getLogger(RowLoop.class.getName());
    private static final double toSeconds = 1.0e-9;
    int nAtoms;
    int nSymm;
    RowRegion rowRegion;

    public RowLoop(int nAtoms, int nSymm, RowRegion rowRegion) {
        this.nAtoms = nAtoms;
        this.nSymm = nSymm;
        this.rowRegion = rowRegion;
    }

    /**
     * <p>
     * setNsymm</p>
     *
     * @param nSymm a int.
     */
    public void setNsymm(int nSymm) {
        this.nSymm = nSymm;
        assert (nSymm <= rowRegion.nSymm);
    }

    @Override
    public void run(int lb, int ub) throws Exception {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                if (rowRegion.select[iSymm][iAtom]) {
                    //logger.info(String.format(" SymOp %d Atom %d", n, i));
                    gridDensity(iSymm, iAtom, lb, ub);
                }
            }
        }
    }

    /**
     * Apply electron density "as normal", but check that the z index is greater
     * than or equal to lb and less than or equal to ub.
     *
     * @param iSymm the SymOp to apply.
     * @param iAtom the index of the Atom to put onto the grid.
     * @param lb the lower bound along the z-axis.
     * @param ub the upper bound along the z-axis.
     */
    public abstract void gridDensity(int iSymm, int iAtom, int lb, int ub);
}
