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

import static java.lang.Math.floor;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedBooleanArray;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

/**
 * The BulkSolventList class builds a list of atoms in symmetry mates that are
 * within a cutoff distance of an atom in the asymmetric unit. This is done
 * in parallel via a spatial decomposition.
 * <ol>
 * <li>
 * The unit cell is partitioned into <code>nA * nB * nC</code> smaller
 * axis-aligned cells, where nA, nB and nC are chosen as large as possible
 * subject to the criteria that the length of each side of a sub-volume
 * (rCellA, rCellB, rCellC) multiplied by (nEdgeA, nEdgeB, nEdgeC),
 * respectively, must be greater than the cutoff distance <code>Rcut</code>
 * plus a buffer distance <code>delta</code>:
 * <center><code>rCellA * nEdgeA >= (Rcut + delta)</code></center>
 * <center><code>rCellB * nEdgeB >= (Rcut + delta)</code></center>
 * <center><code>rCellC * nEdgeC >= (Rcut + delta)</code></center>
 * All neighbors of an atom are in a block of
 * (2*nEdgeA+1)(2*nEdgeB+1)(2*nEdgeC+1)
 * neighborCells.
 * </li>
 * </ol>
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class BulkSolventList extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(BulkSolventList.class.getName());
    /**
     * The crystal object defines the unit cell dimensions and spacegroup.
     */
    private final Crystal crystal;
    /**
     * The number of asymmetric units in the unit cell.
     */
    private final int nSymm;
    /**
     * The array of atoms in the asymmetric unit.
     */
    private final Atom atoms[];
    /**
     * The number of atoms in the asymmetric unit.
     */
    private final int nAtoms;
    /**
     * Reduced coordinates for each symmetry copy. [nsymm][3][natom]
     */
    private double coordinates[][][];
    /**
     * The selected atoms. [nsymm][natom]
     */
    private final SharedBooleanArray[] sharedSelect;
    private boolean selected[][];
    /**
     * Number of selected atoms.
     */
    private int nSelected;
    private final int nEdgeA, nEdgeB, nEdgeC;
    /**
     * The number of divisions along the A-axis.
     */
    private int nA;
    /**
     * The number of divisions along the B-axis.
     */
    private int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    private int nC;
    /**
     * The number of cells in one plane (nDivisions^2).
     */
    private int nAB;
    /**
     * The number of cells (nDivisions^3).
     */
    private final int nCells;
    /**
     * A temporary array that holds the index of the cell each atom is assigned
     * to.
     */
    private final int cellIndex[][];
    /**
     * The cell indices of each atom along a A-axis.
     */
    private final int cellA[];
    /**
     * The cell indices of each atom along a B-axis.
     */
    private final int cellB[];
    /**
     * The cell indices of each atom along a C-axis.
     */
    private final int cellC[];
    /**
     * The list of atoms in each cell. [nsymm][natom] = atom index
     */
    private final int cellList[][];
    /**
     * The offset of each atom from the start of the cell. The first atom atom
     * in the cell has 0 offset. [nsymm][natom] = offset of the atom
     */
    private final int cellOffset[][];
    /**
     * The number of atoms in each cell. [nsymm][ncell]
     */
    private final int cellCount[][];
    /**
     * The index of the first atom in each cell. [nsymm][ncell]
     */
    private final int cellStart[][];
    /**
     * The cutoff beyound which the pairwise energy is zero.
     */
    private final double cutoff;
    private final double minLengthA, minLengthB, minLengthC;
    /**
     * For distance comparisons without taking a sqrt.
     */
    private final double cutoff2;
    /**
     * The array of fractional "a", "b", and "c" coordinates.
     */
    private final double frac[][];
    /***************************************************************************
     * Parallel variables.
     */
    /**
     * The ParallelTeam coordinates use of threads and their schedules.
     */
    private final ParallelTeam parallelTeam;
    /**
     * Number of threads used by the parallelTeam.
     */
    private final int threadCount;
    /**
     * A Verlet list loop for each thread.
     */
    private final SelectionListLoop selectionListLoop[];
    private long time;
    private long cellTime, selectTime, totalTime;
    private double toSeconds = 1.0e-9;

    /**
     * Constructor for the NeighborList class.
     *
     * @param crystal Definition of the unit cell and space group.
     * @param atoms The atoms to generate Verlet lists for.
     * @param cutoff The cutoff distance.
     * @param parallelTeam Specifies the parallel environment.
     *
     * @since 1.0
     */
    public BulkSolventList(Crystal crystal, Atom atoms[], double cutoff,
            ParallelTeam parallelTeam) {
        this.crystal = crystal;
        this.atoms = atoms;
        this.cutoff = cutoff;
        this.parallelTeam = parallelTeam;
        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.symOps.size();

        cutoff2 = cutoff * cutoff;
        final double side = min(min(crystal.a, crystal.b), crystal.c);

        assert (side > 2.0 * cutoff);

        /**
         * nEdgeA, nEdgeB and nEdgeC must be >= 1.
         */
        nEdgeA = 2;
        nEdgeB = nEdgeA;
        nEdgeC = nEdgeA;
        minLengthA = cutoff / (double) nEdgeA;
        minLengthB = cutoff / (double) nEdgeB;
        minLengthC = cutoff / (double) nEdgeC;

        nA = (int) floor(crystal.a / minLengthA);
        nB = (int) floor(crystal.b / minLengthB);
        nC = (int) floor(crystal.c / minLengthC);

        if (nA < nEdgeA * 2 + 1) {
            nA = 1;
        }
        if (nB < nEdgeB * 2 + 1) {
            nB = 1;
        }
        if (nC < nEdgeC * 2 + 1) {
            nC = 1;
        }

        nAB = nA * nB;
        nCells = nAB * nC;

        cellList = new int[nSymm][nAtoms];
        cellIndex = new int[nSymm][nAtoms];
        cellOffset = new int[nSymm][nAtoms];
        cellStart = new int[nSymm][nCells];
        cellCount = new int[nSymm][nCells];
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
        frac = new double[3][nAtoms];
        // Parallel constructs.
        threadCount = parallelTeam.getThreadCount();
        selectionListLoop = new SelectionListLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            selectionListLoop[i] = new SelectionListLoop();
        }
        sharedSelect = new SharedBooleanArray[nSymm];
        for (int i = 0; i < nSymm; i++) {
            sharedSelect[i] = new SharedBooleanArray(nAtoms);
        }
    }

    /**
     * This method can be called as necessary to build/rebuild the neighbor
     * lists.
     *
     * @param coordinates The coordinates of each atom [nSymm][nAtoms*3].
     * @param selected The list of selected atoms [nSymm][nAtoms].
     *
     * @since 1.0
     */
    public void buildList(final double coordinates[][][], final boolean selected[][],
            boolean log) {
        this.coordinates = coordinates;
        this.selected = selected;

        cellTime = -System.nanoTime();
        assignAtomsToCells();
        cellTime += System.nanoTime();

        selectTime = -System.nanoTime();
        selectAtoms();
        selectTime += System.nanoTime();
        totalTime = cellTime + selectTime;

        if (log) {
            log();
        }
    }

    private void log() {
        StringBuilder sb = new StringBuilder(format(" The cutoff is %5.2f angstroms.\n", cutoff));
        sb.append(format(" Assignment to cells: %8.3f\n", cellTime * toSeconds));
        sb.append(format(" Selection:           %8.3f\n", selectTime * toSeconds));
        sb.append(format(" Total:               %8.3f (sec)\n", totalTime * toSeconds));
        if (nSymm > 1) {
            double speedup = (double) (nAtoms * nSymm) / (double) (nAtoms + nSelected);
            sb.append(format(" Atoms in the asymmetric unit: %12d\n", nAtoms));
            sb.append(format(" Selected symmetry mate atoms: %12d\n", nSelected));
            sb.append(format(" Babinet speed up factor:      %12.3f\n", speedup));
        }
        sb.append("\n");
        logger.info(sb.toString());
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells. This is very fast;
     * there is little to be gained from parallelizing it at this point.
     *
     * @since 1.0
     */
    private void assignAtomsToCells() {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int cellIndexs[] = cellIndex[iSymm];
            final int cellCounts[] = cellCount[iSymm];
            final int cellStarts[] = cellStart[iSymm];
            final int cellLists[] = cellList[iSymm];
            final int cellOffsets[] = cellOffset[iSymm];
            // Zero out the cell counts.
            for (int i = 0; i < nCells; i++) {
                cellCounts[i] = 0;
            }
            // Convert to fractional coordinates.
            final double xyz[][] = coordinates[iSymm];


            crystal.toFractionalCoordinates(nAtoms, xyz[0], xyz[1], xyz[2], frac[0], frac[1], frac[2]);

            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i < nAtoms; i++) {
                double xu = frac[0][i];
                double yu = frac[1][i];
                double zu = frac[2][i];
                // Move the atom into the range 0.0 <= x < 1.0
                while (xu >= 1.0) {
                    xu -= 1.0;
                }
                while (xu < 0.0) {
                    xu += 1.0;
                }
                while (yu >= 1.0) {
                    yu -= 1.0;
                }
                while (yu < 0.0) {
                    yu += 1.0;
                }
                while (zu >= 1.0) {
                    zu -= 1.0;
                }
                while (zu < 0.0) {
                    zu += 1.0;
                }
                // The cell indices of this atom.
                final int a = (int) floor(xu * nA);
                final int b = (int) floor(yu * nB);
                final int c = (int) floor(zu * nC);
                if (iSymm == 0) {
                    cellA[i] = a;
                    cellB[i] = b;
                    cellC[i] = c;
                }
                // The cell index of this atom.
                final int index = a + b * nA + c * nAB;
                cellIndexs[i] = index;
                // The offset of this atom from the beginning of the cell.
                cellOffsets[i] = cellCounts[index]++;
            }
            // Define the starting indices.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
            // Move atom locations into a list ordered by cell.
            for (int i = 0; i < nAtoms; i++) {
                final int index = cellIndexs[i];
                cellLists[cellStarts[index]++] = i;
            }
            // Define the starting indices again.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
        }
    }

    /**
     * Execute the parallel atom selection builder.
     *
     * @since 1.0
     */
    private void selectAtoms() {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            if (selected[iSymm] == null) {
                selected[iSymm] = new boolean[nAtoms];
            }
            for (int i = 0; i < nAtoms; i++) {
                sharedSelect[iSymm].set(i, false);
            }
        }
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list.\n";
            logger.log(Level.SEVERE, message, e);
        }
        nSelected = 0;
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
            boolean select[] = selected[iSymm];
            SharedBooleanArray shared = sharedSelect[iSymm];
            for (int i = 0; i < nAtoms; i++) {
                select[i] = shared.get(i);
                if (select[i]) {
                    nSelected++;
                }
            }
        }
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void start() {
        time = System.nanoTime();
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, selectionListLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list in thread: " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * since 0.1
     */
    @Override
    public void finish() {
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format("Parallel Neighbor List: %10.3f seconds",
                    (System.nanoTime() - time) * 1e-9));
        }
    }

    /**
     * The SelectionListLoop class encapsulates thread local variables and methods
     * for selecting atoms based on a spatial decomposition of the unit cell.
     *
     * @author Michael J. Schnieders
     *
     * @since 1.0
     */
    private class SelectionListLoop extends IntegerForLoop {

        private int iSymm;
        private int atomIndex;
        private double xyz[][];
        private SharedBooleanArray select;
        private final IntegerSchedule schedule;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        public SelectionListLoop() {
            super();
            schedule = IntegerSchedule.dynamic(10);
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            xyz = coordinates[0];
        }

        @Override
        public void run(final int lb, final int ub) {
            for (iSymm = 1; iSymm < nSymm; iSymm++) {
                select = sharedSelect[iSymm];
                // Loop over all atoms.
                for (atomIndex = lb; atomIndex <= ub; atomIndex++) {
                    final int a = cellA[atomIndex];
                    final int b = cellB[atomIndex];
                    final int c = cellC[atomIndex];

                    int aStart = a - nEdgeA;
                    int aStop = a + nEdgeA;
                    int bStart = b - nEdgeB;
                    int bStop = b + nEdgeB;
                    int cStart = c - nEdgeC;
                    int cStop = c + nEdgeC;

                    /**
                     * If the number of divisions is 1 in any direction
                     * then set the loop limits to the current cell
                     * value.
                     */
                    if (nA == 1) {
                        aStart = a;
                        aStop = a;
                    }
                    if (nB == 1) {
                        bStart = b;
                        bStop = b;
                    }
                    if (nC == 1) {
                        cStart = c;
                        cStop = c;
                    }

                    /**
                     * Check atoms in symmetry mate cells.
                     */
                    for (int ai = aStart; ai <= aStop; ai++) {
                        for (int bi = bStart; bi <= bStop; bi++) {
                            for (int ci = cStart; ci <= cStop; ci++) {
                                selectAsymmetricAtoms(image(ai, bi, ci));
                            }
                        }
                    }

                }
            }
        }

        /**
         * If the index is >= to nX, it is mapped back into the periodic unit
         * cell by subtracting nX. If the index is < 0, it is mapped into the
         * periodic unit cell by adding nX. The Neighbor list algorithm never
         * requires multiple additions or subtractions of nX.
         *
         * @param i The index along the a-axis.
         * @param j The index along the b-axis.
         * @param k The index along the c-axis.
         * @return The pointer into the 1D cell array.
         */
        private int image(int i, int j, int k) {
            if (i >= nA) {
                i -= nA;
            } else if (i < 0) {
                i += nA;
            }
            if (j >= nB) {
                j -= nB;
            } else if (j < 0) {
                j += nB;
            }
            if (k >= nC) {
                k -= nC;
            } else if (k < 0) {
                k += nC;
            }
            return i + j * nA + k * nAB;
        }

        private void selectAsymmetricAtoms(final int pairCellIndex) {
            final double xi = xyz[0][atomIndex];
            final double yi = xyz[1][atomIndex];
            final double zi = xyz[2][atomIndex];
            final int pairList[] = cellList[iSymm];
            int start = cellStart[iSymm][pairCellIndex];
            final int pairStop = start + cellCount[iSymm][pairCellIndex];
            final double pair[][] = coordinates[iSymm];
            // Loop over atoms in the "pair" cell.
            for (int j = start; j < pairStop; j++) {
                final int aj = pairList[j];
                // If this atom is already selected, continue.
                if (select.get(aj)) {
                    continue;
                }
                final double xj = pair[0][aj];
                final double yj = pair[1][aj];
                final double zj = pair[2][aj];
                final double xr = xi - xj;
                final double yr = yi - yj;
                final double zr = zi - zj;
                final double d2 = crystal.image(xr, yr, zr);
                if (d2 <= cutoff2) {
                    select.set(aj, true);
                }
            }
        }
    }
    private final static int XX = 0;
    private final static int YY = 1;
    private final static int ZZ = 2;
}
