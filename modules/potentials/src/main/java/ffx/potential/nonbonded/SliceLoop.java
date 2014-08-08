/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.nonbonded;

import edu.rit.pj.IntegerForLoop;

/**
 *
 * @author avdic
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

        /**
         * Todo: pre-compute lists of atoms that can contribute to this thread's
         * slices.
         */
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
