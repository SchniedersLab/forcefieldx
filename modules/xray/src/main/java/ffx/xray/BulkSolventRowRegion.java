/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.util.Arrays.fill;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.RowRegion;

public class BulkSolventRowRegion extends RowRegion {

    /**
     * Constant <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(BulkSolventRowRegion.class.getName());

    private final BulkSolventList bulkSolventList;
    private final int gZ;
    private final int gY;

    /**
     * <p>
     * Constructor for BulkSolventDensityRegion.</p>
     *
     * @param gX a int.
     * @param gY a int.
     * @param gZ a int.
     * @param grid an array of double.
     * @param basisSize a int.
     * @param nSymm a int.
     * @param threadCount a int.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param coordinates an array of double.
     * @param cutoff a double.
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
     */
    public BulkSolventRowRegion(int gX, int gY, int gZ, double grid[],
            int basisSize, int nSymm, int threadCount, Crystal crystal,
            Atom atoms[], double coordinates[][][],
            double cutoff, ParallelTeam parallelTeam) {
        super(gX, gY, gZ, grid, basisSize, nSymm,
                threadCount, crystal, atoms, coordinates);

        this.gZ = gZ;
        this.gY = gY;
        // Asymmetric unit atoms never selected by this class.
        fill(select[0], false);
        bulkSolventList = new BulkSolventList(crystal, atoms, cutoff, parallelTeam);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        try {
            execute(0, (gZ * gY) - 1, rowLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = " Exception in BulkSolventRowRegion.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void selectAtoms() {
        bulkSolventList.buildList(coordinates, select, false);
    }
}
