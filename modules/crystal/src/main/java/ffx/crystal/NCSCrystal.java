/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.crystal;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * The NCSCrystal class extends Crystal to support non-crystallographic symmetry (NCS).
 * The NCS operators can be obtained from the BIOMT records of a PDB file, and are
 * permuted with the space group symmetry operators.
 *
 * @author Aaron J. Nessler & Michael J. Schnieders
 *
 * @since 1.0
 *
 * @see Crystal
 */
public class NCSCrystal extends Crystal {

    /**
     * The logger.
     */
    private static final Logger logger = Logger.getLogger(NCSCrystal.class.getName());
    /**
     * The base unit cell for the system being simulated.
     */
    private final Crystal unitCell;
    /**
     * The non-crystallographic symmetry operators needed to be applied to the unit cell.
     */
    private final List<SymOp> NCSsymOps;

    /**
     * Constructor for a NCSCrystal.
     *
     * @param unitCell The base unit cell.
     * @param NCSsymOps Non-crystallographic symmetry operators applied to the unit cell.
     * @since 1.0
     */
    private NCSCrystal(Crystal unitCell, List<SymOp> NCSsymOps) {
        super(unitCell.a, unitCell.b, unitCell.c, unitCell.alpha,
                unitCell.beta, unitCell.gamma, unitCell.spaceGroup.shortName);
        this.unitCell = unitCell;
        this.NCSsymOps = NCSsymOps;

        /**
         * At this point, the NCSCrystal references a SpaceGroup instance
         * that is lacking symmetry operators. This is corrected by generating symmetry
         * operators to fill up the non-crystallographic, which is added to the unit 
         * cell based on the asymmetric unit.
         */
        updateNCSOperators();
    }

    /**
     * Update the list of symmetry operators to contain the non-crystallographic operations as well.
     */
    private void updateNCSOperators() {
        List<SymOp> symOps = spaceGroup.symOps;
        /**
         * First, we remove the existing symmetry operators.
         */
        symOps.clear();

        /**
         * Symmetry operators are produced that are the combination of the non-crystallographic and the unit cell symmetry operations.
         * Note that the first symOp is still equivalent to the asymmetric unit and the first set of
         * symOps are still equivalent to the unit cell.
         */
        int ii = 0;
        int soRemoved = 0;
        for (SymOp NCSsymOp : NCSsymOps) {
            // All UC symOps seem to have collision issues... Use 1st (identity) for now...
            for (SymOp symOp : unitCell.spaceGroup.symOps) {
//            SymOp symOp=unitCell.spaceGroup.symOps.get(0); //Identity matrix is standard first matrix. (P1)
                double NCSTrans[] = new double[3];
                RealMatrix NCS = MatrixUtils.createRealMatrix(NCSsymOp.rot);
                RealMatrix UC = MatrixUtils.createRealMatrix(symOp.rot);
                RealMatrix result = NCS.multiply(UC); //Abelian groups order doesnt matter...
                double NCSRot[][] = result.getData();
                NCSTrans[0] = symOp.tr[0] + NCSsymOp.tr[0];
                NCSTrans[1] = symOp.tr[1] + NCSsymOp.tr[1];
                NCSTrans[2] = symOp.tr[2] + NCSsymOp.tr[2];
                SymOp NCSSymOp = new SymOp(NCSRot, NCSTrans);
                if (!symOps.contains(NCSSymOp)) {
                    symOps.add(NCSSymOp);
                } else {
                    soRemoved++;
                }
                if (logger.isLoggable(Level.FINEST)) {
                    logger.finest(format("\n SymOp: %d", ii));
                    logger.finest(NCSSymOp.toString());
                }
                ii++;
            }
        }
        if (soRemoved != 0) {
            logger.warning(format("\n NCS Replicated SymOps Removed: %d", soRemoved));
        }
    }

    /**
     * Change the cell parameters for the base unit cell, which is followed by
     * an update of the ReplicateCrystal parameters and possibly the number of
     * replicated cells.
     *
     * @param a The length of the a-axis for the base unit cell (in Angstroms).
     * @param b The length of the b-axis for the base unit cell (in Angstroms).
     * @param c The length of the c-axis for the base unit cell (in Angstroms).
     * @param alpha The angle between the b-axis and c-axis (in Degrees).
     * @param beta The angle between the a-axis and c-axis (in Degrees).
     * @param gamma The angle between the a-axis and b-axis (in Degrees).
     * @return True is returned if the unit cell and replicates cell are updated
     * successfully.
     */

    @Override
    public boolean changeUnitCellParameters(double a, double b, double c,
                                            double alpha, double beta, double gamma) {
        /**
         * First, update the parameters of the unit cell.
         */
        if (unitCell.changeUnitCellParameters(a, b, c, alpha, beta, gamma)) {
            /**
             * Then, update the parameters of the NCSCrystal.
             */
            if (super.changeUnitCellParameters(a, b, c, alpha, beta, gamma)) {
                // Finally, update the NCS operators.
                updateNCSOperators();
                return true;
            }
        }
        return false;
    }

    /**
     * Two crystals are equal only if all unit cell parameters are exactly the
     * same.
     */
    @Override
    public boolean strictEquals(Object obj) {
        if (obj == null) {
            return false;
        }

        if (!(obj instanceof NCSCrystal)) {
            return false;
        }

        if (this == obj) {
            return true;
        }

        NCSCrystal other = (NCSCrystal) obj;
        return this.unitCell.strictEquals(other.unitCell) && (a == other.a && b == other.b && c == other.c
                && alpha == other.alpha && beta == other.beta && gamma == other.gamma
                && spaceGroup.number == other.spaceGroup.number
                && spaceGroup.symOps.size() == other.spaceGroup.symOps.size());
    }

    /**
     * {@inheritDoc}
     *
     * Returns the unit cell for this NCSCrystal. This is useful for the
     * reciprocal space portion of PME that operates on the unit cell.
     */
    @Override
    public Crystal getUnitCell() {
        return unitCell;
    }

    /**
     * {@inheritDoc}
     *
     * Include information about the base unit cell and NCS cell.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(unitCell.toString());

        sb.append("\n\n Non-Crystallographic Cell\n");
        sb.append(String.format("  A-axis:                              %8.3f\n", a));
        sb.append(String.format("  B-axis:                              %8.3f\n", b));
        sb.append(String.format("  C-axis:                              %8.3f\n", c));
        sb.append(String.format("  Alpha:                               %8.3f\n", alpha));
        sb.append(String.format("  Beta:                                %8.3f\n", beta));
        sb.append(String.format("  Gamma:                               %8.3f\n", gamma));
        sb.append(String.format("  UnitCell Symmetry Operators:         %8d\n", unitCell.spaceGroup.symOps.size()));
        sb.append(String.format("  NCS Symmetry Operators:              %8d\n", NCSsymOps.size()));
        sb.append(String.format("  Total Symmetry Operators:            %8d", spaceGroup.getNumberOfSymOps()));

        return sb.toString();
    }

    /**
     * Returns an NCSCrystal by expanding the orignal unit cell with the symmetry operators provided by the BIOMT
     * records in the PDB files. See REMARK 350.
     *
     * @param unitCell The unit cell of the crystal.
     * @param symOps Symmetry operators for non-crystallographic symmetry
     * @return A Crystal or NCSCrystal
     */

    public static Crystal NCSCrystalFactory(Crystal unitCell, List<SymOp> symOps) {

        if (unitCell == null || unitCell.aperiodic()) {
            return unitCell;
        } else {
            return new NCSCrystal(unitCell, symOps);
        }
    }
}
