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

import static java.lang.Math.floor;
import static java.lang.Math.min;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

/**
 * The NeighborList class builds Verlet lists in parallel via a
 * spatial decomposition.
 * <ol>
 * <li>
 * The unit cell is partitioned into <code>nA * nB * nC</code> smaller
 * isomorphous cells, where nA, nB and nC are chosen as large as possible
 * subject to the criteria that the length of the shortest side of a subvolume
 * <code>Rcell</code> must be at least half the cutoff distance
 * <code>Rcut</code> plus a buffer distance <code>delta</code>:
 * <center><code>Rcell >= Rcut + delta</code></center>
 * All atoms are then assigned to a cell.
 * </li>
 * <p>
 * <li>
 * Based on this spatial decomposition, neighbors of an atom must be in its cell
 * or the 124 surrounding cells. Interactions between an atom and neighbors in
 * the asymmetric unit require only half the neighboring cells (62) to be
 * searched, rather than 124. This avoids including the same interaction twice.
 * However, enumeration of interactions between an atom in the asymmetric unit
 * and its neighbors in a symmetry mate require all 125 cells to be searched.
 * </li>
 * <p>
 * <li>
 * Verlet lists from the search are stored, which reduces the number of
 * neigbors whose distances must be calculated by a factor of approximately:
 * <center><code>(4/3*Pi*Rcut^3)/(125*Vcell)</code></center>
 * About 1/3 as many interactions are contained in the Verlet lists as in the
 * 125 neighboring cells.
 * </li>
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class NeighborList extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(NeighborList.class.getName());

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
     * The maksing rules to apply when building the neighbor list.
     */
    private final MaskingInterface maskingRules;
    /**
     * Reduced coordinates for each symmetry copy. [nsymm][3][natom]
     */
    private double coordinates[][];
    /**
     * The reduced coordinates of the asymmetric unit when the list was last
     * rebuilt.
     */
    private final double previous[];
    /**
     * The Verlet lists. [nsymm][natom][ neighbors... ]
     */
    private int lists[][][];
    /**
     * Number of interactions between atoms in the asymmetric unit.
     */
    private int asymmetricUnitCount;
    /**
     * Number of interactions between atoms in the asymmetric unit and atoms in
     * symmetry mates.
     */
    private int symmetryMateCount;
    /**
     * The minimum number of cells along each axis.
     */
    private int minDivisions;
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
    /**
     * A buffer, which is added to the cutoff distance, such that the Verlet
     * lists do not need to be calculated for all coordinate changes.
     */
    private final double buffer;
    /**
     * The maximum squared displacement allowed before list rebuild.
     */
    private final double motion2;
    /**
     * The sum of the cutoff + buffer.
     */
    private final double total;
    private final double halftotal;
    /**
     * Total^2 for distance comparisons without taking a sqrt.
     */
    private final double total2;
    /**
     * The array of fractional "a", "b", and "c" coordinates.
     */
    private final double frac[];
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
    private final VerletListLoop verletListLoop[];
    private long time;
    private int len = 1000;

    /**
     * Constructor for the NeighborList class.
     *
     * @param maskingRules This parameter may be null.
     * @param crystal Definition of the unit cell and space group.
     * @param atoms The atoms to generate Verlet lists for.
     * @param cutoff The cutoff distance.
     * @param buffer The buffer distance.
     * @param parallelTeam Specifies the parallel environment.
     *
     * @since 1.0
     */
    public NeighborList(MaskingInterface maskingRules, Crystal crystal,
            Atom atoms[], double cutoff, double buffer,
            ParallelTeam parallelTeam) {
        this.maskingRules = maskingRules;
        this.crystal = crystal;
        this.atoms = atoms;
        this.cutoff = cutoff;
        this.buffer = buffer;
        this.parallelTeam = parallelTeam;
        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.symOps.size();
        total = cutoff + buffer;
        halftotal = total / 2.0;
        total2 = total * total;
        motion2 = (buffer / 2.0) * (buffer / 2.0);
        previous = new double[nAtoms * 3];
        final double side = min(min(crystal.a, crystal.b), crystal.c);
        assert (side > 2.0 * total);
        nA = (int) floor(crystal.a / halftotal);
        nB = (int) floor(crystal.b / halftotal);
        nC = (int) floor(crystal.c / halftotal);
        minDivisions = min(min(nA, nB), nC);
        /**
         * If the number of divisions is greater than >=3, we'll use a spatial
         * decomposition. Otherwise, the unit cell is too small relative to the
         * cutoff plus buffer distance use the logic implemented here.
         */
        StringBuffer sb = new StringBuffer("\n NEIGHBOR LIST BUILDER\n");
        if (minDivisions >= 5) {
            nAB = nA * nB;
            nCells = nAB * nC;
            sb.append(" The unit cell is partitioned into " + nCells + " volumes.");
        } else {
            minDivisions = 1;
            nA = 1;
            nB = 1;
            nC = 1;
            nAB = 1;
            nCells = 1;
            sb.append(" The unit cell will not be partitioned into cells.");
        }
        cellList = new int[nSymm][nAtoms];
        cellIndex = new int[nSymm][nAtoms];
        cellOffset = new int[nSymm][nAtoms];
        cellStart = new int[nSymm][nCells];
        cellCount = new int[nSymm][nCells];
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
        frac = new double[3 * nAtoms];
        // Parallel constructs.
        threadCount = parallelTeam.getThreadCount();
        verletListLoop = new VerletListLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            verletListLoop[i] = new VerletListLoop();
        }
        logger.info(sb.toString());
    }

    /**
     * This method can be called as necessary to build/rebuild the neighbor
     * lists.
     *
     * @param coordinates The coordinates of each atom [nSymm][x/y/z][nAtoms].
     * @param lists The neighbor lists [nSymm][nAtoms][nPairs].
     * @param rebuild If true, the list is rebuilt even if no atom has moved
     *      half the buffer size.
     *
     * @since 1.0
     */
    public void buildList(final double coordinates[][], final int lists[][][],
            boolean rebuild, boolean print) {
        this.coordinates = coordinates;
        this.lists = lists;
        if (rebuild || motion()) {
            /**
             * Save the current coordinates.
             */
            double current[] = coordinates[0];
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                int iX = i3 + XX;
                int iY = i3 + YY;
                int iZ = i3 + ZZ;
                previous[iX] = current[iX];
                previous[iY] = current[iY];
                previous[iZ] = current[iZ];
            }

            final long startTime = System.nanoTime();
            assignAtomsToCells();
            long currentTime = System.nanoTime();
            final long cellTime = currentTime - startTime;
            createVerletLists();
            currentTime = System.nanoTime();
            final long totalTime = currentTime - startTime;
            final long verletTime = totalTime - cellTime;
            if (print) {
                StringBuffer sb = new StringBuffer(" The cutoff is " + cutoff + " angstroms.\n");
                final double toSeconds = 0.000000001;
                sb.append(String.format(" Assignment to cells: %8.3f\n", cellTime * toSeconds)
                        + String.format(" Verlet lists:        %8.3f\n", verletTime * toSeconds)
                        + String.format(" Total:               %8.3f (sec)\n", totalTime * toSeconds));
                sb.append(String.format(" Neighbors in the asymmetric unit: %12d\n", asymmetricUnitCount));
                if (nSymm > 1) {
                    int num = (int) (asymmetricUnitCount * nSymm + symmetryMateCount * (nSymm * 0.5));
                    double speedup = ((double) num) / (asymmetricUnitCount + symmetryMateCount);
                    sb.append(String.format(" Neighbors in symmetry mates:      %12d\n", symmetryMateCount));
                    sb.append(String.format(" Neighbors in the unit cell:       %12d\n", num));
                    sb.append(String.format(" Space group speed up factor:      %12.3f\n", speedup));
                }
                logger.info(sb.toString() + "\n");
            }
            int atomCount[] = new int[nAtoms];
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int i = 0; i < nAtoms; i++) {
                    atomCount[i] += lists[iSymm][i].length;
                }
            }
            /*
            sb = new StringBuffer("\n");
            for (int i = 0; i < nAtoms; i++) {
            if (i < 200) {
            sb.append(String.format("%d: %d\n", i, atomCount[i]));
            }
            }
            logger.info(sb.toString());
            System.exit(1); */
        }
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells. This is very fast;
     * there is little to be gained from parallelizing it at this point.
     *
     * @since 1.0
     */
    private void assignAtomsToCells() {
        // Zero out the cell counts.
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int cellIndexs[] = cellIndex[iSymm];
            final int cellCounts[] = cellCount[iSymm];
            final int cellStarts[] = cellStart[iSymm];
            final int cellLists[] = cellList[iSymm];
            final int cellOffsets[] = cellOffset[iSymm];
            for (int i = 0; i < nCells; i++) {
                cellCounts[i] = 0;
            }
            // Convert to fractional coordinates.
            final double xyz[] = coordinates[iSymm];
            crystal.toFractionalCoordinates(nAtoms, xyz, frac);
            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                double xu = frac[i3 + XX];
                double yu = frac[i3 + YY];
                double zu = frac[i3 + ZZ];
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
                final int a = (int) Math.floor(xu * nA);
                final int b = (int) Math.floor(yu * nB);
                final int c = (int) Math.floor(zu * nC);
                if (iSymm == 0) {
                    cellA[i] = a;
                    cellB[i] = b;
                    cellC[i] = c;
                }
                /*
                if (iSymm == 0 && i < 100) {
                System.out.println(String.format("%5d (%10.5f,%10.5f,%10.5f)/(%5d,%5d,%5d)",i+1,xu,yu,zu,a,b,c));
                }
                 */
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
     * Execute the parallel Verlet list builder.
     *
     * @since 1.0
     */
    private void createVerletLists() {
        asymmetricUnitCount = 0;
        symmetryMateCount = 0;
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            if (lists[iSymm] == null) {
                lists[iSymm] = new int[nAtoms][];
            }
        }
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
        int list[][] = lists[0];
        for (int i = 0; i < nAtoms; i++) {
            asymmetricUnitCount += list[i].length;
        }
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
            list = lists[iSymm];
            for (int i = 0; i < nAtoms; i++) {
                symmetryMateCount += list[i].length;
            }
        }
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 0.l
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
            execute(0, nAtoms - 1, verletListLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list in thread: " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
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
                    (System.nanoTime() - time) * 0.000000001));
        }
    }

    /**
     * Detect if any atom has moved 1/2 the buffer size.
     *
     * @return True if an atom has moved 1/2 the buffer size.
     *
     * @since 1.0
     */
    private boolean motion() {
        boolean motion = false;
        double current[] = coordinates[0];
        for (int i = 0; i < nAtoms; i++) {
            int i3 = i * 3;
            int iX = i3 + XX;
            int iY = i3 + YY;
            int iZ = i3 + ZZ;
            double dx = previous[iX] - current[iX];
            double dy = previous[iY] - current[iY];
            double dz = previous[iZ] - current[iZ];
            double dr2 = crystal.image(dx, dy, dz);
            if (dr2 > motion2) {
                motion = true;
                break;
            }
        }
        return motion;
    }

    /**
     * The VerletListLoop class encapsulates thread local variables and methods
     * for building Verlet lists based on a spatial decomposition of the unit
     * cell.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class VerletListLoop extends IntegerForLoop {

        private int n;
        private int iSymm;
        private int atomIndex;
        private double xyz[];
        private final int pairs[];
        private final double mask[];
        private final int asymmetricIndex[];
        private final IntegerSchedule schedule;

        public VerletListLoop() {
            super();
            pairs = new int[len];
            asymmetricIndex = cellIndex[0];
            mask = new double[nAtoms];
            for (int i = 0; i < nAtoms; i++) {
                mask[i] = 1.0;
            }
            schedule = IntegerSchedule.guided(10);
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
            for (iSymm = 0; iSymm < nSymm; iSymm++) {
                int list[][] = lists[iSymm];
                // Loop over all atoms.
                for (atomIndex = lb; atomIndex <= ub; atomIndex++) {
                    n = 0;
                    if (minDivisions == 1) {
                        atomCellPairs(0);
                    } else {
                        final int a = cellA[atomIndex];
                        final int b = cellB[atomIndex];
                        final int c = cellC[atomIndex];
                        final int index = a + b * nA + c * nAB;
                        final int a1 = a + 1;
                        final int a2 = a + 2;
                        final int b1 = b + 1;
                        final int b2 = b + 2;
                        final int c1 = c + 1;
                        final int c2 = c + 2;
                        final int bm2 = b - 2;
                        final int cm2 = c - 2;
                        if (iSymm == 0) {
                            // Interactions within the "self-volume".
                            atomCellPairs(index);
                            // Half (62) of the neighboring volumes (124) are
                            // searched to avoid double counting.
                            // 2 cells (a, b+1..2, c)
                            atomCellPairs(image(a, b1, c));
                            atomCellPairs(image(a, b2, c));
                            // 10 cells (a, b-2..2, c+1..2)
                            for (int bi = bm2; bi <= b2; bi++) {
                                atomCellPairs(image(a, bi, c1));
                                atomCellPairs(image(a, bi, c2));
                            }
                            // 50 cells (a+1..2, b-2..2, c-2..2)
                            for (int bi = bm2; bi <= b2; bi++) {
                                for (int ci = cm2; ci <= c2; ci++) {
                                    atomCellPairs(image(a1, bi, ci));
                                    atomCellPairs(image(a2, bi, ci));
                                }
                            }
                        } else {
                            // Interactions with the 125 adjacent symmetry mate
                            // cells.
                            for (int ai = a - 2; ai <= a2; ai++) {
                                for (int bi = bm2; bi <= b2; bi++) {
                                    for (int ci = cm2; ci <= c2; ci++) {
                                        atomCellPairs(image(ai, bi, ci));
                                    }
                                }
                            }
                        }
                    }
                    list[atomIndex] = new int[n];
                    System.arraycopy(pairs, 0, list[atomIndex], 0, n);
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

        private void atomCellPairs(final int pairCellIndex) {
            final int atomCellIndex = asymmetricIndex[atomIndex];
            final int i3 = atomIndex * 3;
            final double xi = xyz[i3 + XX];
            final double yi = xyz[i3 + YY];
            final double zi = xyz[i3 + ZZ];
            final int pairList[] = cellList[iSymm];
            int start = cellStart[iSymm][pairCellIndex];
            final int pairStop = start + cellCount[iSymm][pairCellIndex];
            final double pair[] = coordinates[iSymm];
            // Check if this pair search is over atoms in the asymmetric unit.
            if (iSymm == 0) {
                // Interactions between atoms in the asymmetrc unit may be
                // masked.
                if (maskingRules != null) {
                    maskingRules.applyMask(mask, atomIndex);
                }
                // If the self-volume is being searched for pairs, we must avoid
                // double counting.
                if (atomCellIndex == pairCellIndex) {
                    start += cellOffset[iSymm][atomIndex] + 1;
                }
            }
            // Loop over atoms in the "pair" cell.
            for (int j = start; j < pairStop; j++) {
                final int aj = pairList[j];
                if (mask[aj] > 0.0) {
                    int aj3 = aj * 3;
                    final double xj = pair[aj3 + XX];
                    final double yj = pair[aj3 + YY];
                    final double zj = pair[aj3 + ZZ];
                    final double xr = xi - xj;
                    final double yr = yi - yj;
                    final double zr = zi - zj;
                    final double d2 = crystal.image(xr, yr, zr);
                    /*
                    if (atomIndex == 0 && aj < 100) {
                    System.out.print(String.format("%5d %5d %10.5f\n",atomIndex+1,aj+1,Math.sqrt(d2)));
                    }
                     */
                    if (d2 <= total2) {
                        pairs[n++] = aj;
                    }
                }
            }
            if (iSymm == 0 && maskingRules != null) {
                maskingRules.removeMask(mask, atomIndex);
            }
        }
    }
    private final static int XX = 0;
    private final static int YY = 1;
    private final static int ZZ = 2;
}
