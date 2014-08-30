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
public abstract class SliceLoop extends IntegerForLoop {

    /**
     * Constant <code>logger</code>
     */
    private static final Logger logger = Logger.getLogger(SliceLoop.class.getName());
    private static final double toSeconds = 1.0e-9;
    int nAtoms;
    int nSymm;
    SliceRegion sliceRegion;
    double sliceLoopTime;

    public SliceLoop(int nAtoms, int nSymm, SliceRegion sliceRegion) {
        this.nAtoms = nAtoms;
        this.nSymm = nSymm;
        this.sliceRegion = sliceRegion;
    }

    @Override
    public void start() {
        sliceLoopTime -= System.nanoTime();
    }

    /**
     * <p>
     * setNsymm</p>
     *
     * @param nSymm a int.
     */
    public void setNsymm(int nSymm) {
        this.nSymm = nSymm;
        assert (nSymm <= sliceRegion.nSymm);
    }

    @Override
    public void run(int lb, int ub) throws Exception {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                if (sliceRegion.select[iSymm][iAtom]) {
                    //logger.info(String.format(" SymOp %d Atom %d", n, i));
                    gridDensity(iSymm, iAtom, lb, ub);
                }
            }
        }
    }

    @Override
    public void finish() {
        sliceLoopTime += System.nanoTime();

//        logger.info(String.format("Slice Loop Time: %7.4f (sec)",sliceLoopTime * toSeconds));
    }

    /**
     * Apply electron density "as normal", but check that the z index is greater
     * than or equal to lb and less than or equal to ub.
     *
     * @param iSymm
     * @param iAtom
     * @param lb
     * @param ub
     */
    public abstract void gridDensity(int iSymm, int iAtom, int lb, int ub);
}
