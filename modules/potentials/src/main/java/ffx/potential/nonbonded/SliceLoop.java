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

import edu.rit.pj.IntegerForLoop;

/**
 * @author Armin Avdic
 */
public abstract class SliceLoop extends IntegerForLoop {

    int nAtoms;
    int nSymm;

    public SliceLoop(int nAtoms, int nSymm) {
        this.nAtoms = nAtoms;
        this.nSymm = nSymm;
    }

    @Override
    public void run(int lb, int ub) throws Exception {
        for (int n = 0; n < nSymm; n++) {
            for (int i = 0; i < nAtoms; i++) {
                gridDensity(n, i, lb, ub);
            }
        }
    }

    /**
     * Apply electron density "as normal", but check that the z index is lb <= z
     * <= ub.
     *
     * @param atom
     * @param iSymm
     * @param lb
     * @param ub
     */
    abstract void gridDensity(int atom, int iSymm, int lb, int ub);
}
