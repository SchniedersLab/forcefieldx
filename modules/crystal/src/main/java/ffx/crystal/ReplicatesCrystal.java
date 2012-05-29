/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
package ffx.crystal;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The ReplicatesCrystal class extends Crystal to generate additional symmetry
 * operators needed to describe a "replicated" crystal. The replicated crystal
 * cell edges are of length (l*a, m*b, n*c) where l, m and n are integers
 * and a, b and c are the original unit cell edge lengths. The replicates
 * integers l, m and n are chosen large enough so that the ReplicatesCrystal
 * is consistent with application of the minimum image convention. In other
 * words, each replicated cell edge is more than twice the cutoff.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class ReplicatesCrystal extends Crystal {

    private static final Logger logger = Logger.getLogger(ReplicatesCrystal.class.getName());
    private final Crystal unitCell;
    private final int l;
    private final int m;
    private final int n;

    /**
     * Constructor for a ReplicatesCrystal.
     *
     * @param unitCell The base unit cell.
     * @param l Number of replicates along the a-axis.
     * @param m Number of replicates along the b-axis.
     * @param n Number of replicates along the c-axis.
     * @since 1.0
     */
    public ReplicatesCrystal(Crystal unitCell, int l, int m, int n) {
        super(unitCell.a * l, unitCell.b * m, unitCell.c * n, unitCell.alpha,
                unitCell.beta, unitCell.gamma, unitCell.spaceGroup.shortName);
        this.unitCell = unitCell;

        assert (l >= 1);
        assert (m >= 1);
        assert (n >= 1);

        this.l = l;
        this.m = m;
        this.n = n;

        /**
         * At this point, the ReplicatesCrystal is ready to go, except
         * that it references a SpaceGroup instance whose symmetry operators
         * are inconsistent. We need to generate symmetry operators that fill
         * up the ReplicatesCrystal based on the asymmetric unit.
         */
        List<SymOp> symOps = spaceGroup.symOps;
        /**
         * First, we remove the existing symmetry operators.
         */
        symOps.clear();
        
        /**
         * Now create symmetry operators for each replicate. Note that
         * the first symOp is still equivalent to the asymmetric unit and the
         * first set of symOps are still equivalent to the unit cell.
         */
        double dX = 1.0 / (double) l;
        double dY = 1.0 / (double) m;
        double dZ = 1.0 / (double) n;
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < n; k++) {
                    int ii = 0;
                    for (SymOp symOp : unitCell.spaceGroup.symOps) {
                        double repTrans[] = new double[3];
                        repTrans[0] = (symOp.tr[0] + i) * dX;
                        repTrans[1] = (symOp.tr[1] + j) * dY;
                        repTrans[2] = (symOp.tr[2] + k) * dZ;
                        SymOp repSymOp = new SymOp(symOp.rot, repTrans);
                        symOps.add(repSymOp);
                        if (logger.isLoggable(Level.FINEST)) {
                            logger.finest(String.format("\n SymOp (%2d,%2d,%2d): %d", i, j, k, ii));
                            logger.finest(repSymOp.toString());
                        }
                        ii++;
                    }
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     *
     * Returns the unit cell for this ReplicatesCrystal. This is useful for
     * the reciprocal space portion of PME that operates on the unit cell
     * even though the real space cutoff requires a ReplicatesCrystal.
     */
    @Override
    public Crystal getUnitCell() {
        return unitCell;
    }

    /**
     * {@inheritDoc}
     *
     * Include information about the base unit cell and replicates cell.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(unitCell.toString());

        sb.append("\n\n Replicates Cell\n");
        sb.append(String.format(" Dimension:                     (%3d x%3d x%3d)\n", l,m,n));
        sb.append(String.format(" A-axis:                               %8.3f\n", a));
        sb.append(String.format(" B-axis:                               %8.3f\n", b));
        sb.append(String.format(" C-axis:                               %8.3f\n", c));
        sb.append(String.format(" Alpha:                                %8.3f\n", alpha));
        sb.append(String.format(" Beta:                                 %8.3f\n", beta));
        sb.append(String.format(" Gamma:                                %8.3f\n", gamma));
        sb.append(String.format(" Total Symmetry Operators:             %8d\n", spaceGroup.getNumberOfSymOps()));
        return sb.toString();
    }
}
