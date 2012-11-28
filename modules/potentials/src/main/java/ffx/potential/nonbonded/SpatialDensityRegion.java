/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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

import java.nio.DoubleBuffer;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.floor;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

/**
 * This class implements a spatial decomposition based on partitioning a grid
 * into octants.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class SpatialDensityRegion extends ParallelRegion {

    /**
     * Constant
     * <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(SpatialDensityRegion.class.getName());
    /**
     * The number of divisions along the A-axis.
     */
    protected int nA;
    /**
     * The number of divisions along the B-axis.
     */
    protected int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    protected int nC;
    /**
     * The number of cells in one plane (nDivisions^2).
     */
    private int nAB;
    /**
     * The number of cells (nDivisions^3).
     */
    private final int nCells;
    /**
     * Number of octant work cells.
     */
    private final int nWork;
    /**
     * Number of octant work cells with at least one atom (actualWork is less
     * than or equal to nWork).
     */
    protected int actualWork;
    public int actualCount[];
    /**
     * A temporary array that holds the index of the cell each atom is assigned
     * to.
     */
    private final int cellIndex[][];
    /**
     * The A index of each octant (0..nA - 1) that may not have any atoms.
     */
    protected final int workA[];
    /**
     * The B index of each octant (0..nB - 1) that may not have any atoms.
     */
    protected final int workB[];
    /**
     * The C index of each octant (0..nC - 1) that may not have any atoms.
     */
    protected final int workC[];
    /**
     * The A index of each octant (0..nA - 1) that may not have any atoms.
     */
    protected final int actualA[];
    /**
     * The B index of each octant (0..nB - 1) that may not have any atoms.
     */
    protected final int actualB[];
    /**
     * The C index of each octant (0..nC - 1) that may not have any atoms.
     */
    protected final int actualC[];
    /**
     * The list of atoms in each cell. [nsymm][natom] = atom index
     */
    protected final int cellList[][];
    /**
     * The offset of each atom from the start of the cell. The first atom atom
     * in the cell has 0 offset. [nsymm][natom] = offset of the atom
     */
    protected final int cellOffset[][];
    /**
     * The number of atoms in each cell. [nsymm][ncell]
     */
    protected final int cellCount[][];
    /**
     * The index of the first atom in each cell. [nsymm][ncell]
     */
    protected final int cellStart[][];
    protected final int nSymm;
    protected final double coordinates[][][];
    protected final boolean select[][];
    private final double xf[];
    private final double yf[];
    private final double zf[];
    protected final Crystal crystal;
    public final int nAtoms;
    public final int nThreads;
    private final int gridSize;
    private double grid[] = null;
    private DoubleBuffer gridBuffer;
    private double initValue = 0.0;
    protected SpatialDensityLoop spatialDensityLoop[];
    private GridInitLoop gridInitLoop;

    /**
     * <p>Constructor for SpatialDensityRegion.</p>
     *
     * @param gX a int.
     * @param gY a int.
     * @param gZ a int.
     * @param grid an array of double.
     * @param basisSize a int.
     * @param nSymm a int.
     * @param minWork a int.
     * @param threadCount a int.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param coordinates an array of double.
     */
    public SpatialDensityRegion(int gX, int gY, int gZ, double grid[],
            int basisSize, int nSymm, int minWork,
            int threadCount, Crystal crystal,
            Atom atoms[], double coordinates[][][]) {
        this(gX, gY, gZ, basisSize, nSymm, minWork, threadCount, crystal, atoms, coordinates);
        this.grid = grid;
        if (grid != null) {
            gridBuffer = DoubleBuffer.wrap(grid);
        }
    }

    private SpatialDensityRegion(int gX, int gY, int gZ,
            int basisSize, int nSymm, int minWork,
            int threadCount, Crystal crystal,
            Atom atoms[], double coordinates[][][]) {
        /**
         * Chop up the 3D unit cell domain into fractional coordinate chunks to
         * allow multiple threads to put charge density onto the grid without
         * needing the same grid point. First, we partition the X-axis, then the
         * Y-axis, and finally the Z-axis if necessary.
         */
        this.crystal = crystal;
        this.coordinates = coordinates;
        this.nSymm = nSymm;
        this.nAtoms = atoms.length;
        this.nThreads = threadCount;

        gridInitLoop = new GridInitLoop();
        gridSize = gX * gY * gZ * 2;

        xf = new double[nAtoms];
        yf = new double[nAtoms];
        zf = new double[nAtoms];

        int nX = gX / basisSize;
        int nY = gY / basisSize;
        int nZ = gZ / basisSize;
        int div = 1;
        int currentWork = 0;
        if (threadCount > 1 && nZ > 1) {
            if (nZ % 2 != 0) {
                nZ--;
            }
            nC = nZ;
            div = 2;
            currentWork = nC / div / threadCount;
            // If we have enough work per thread, stop dividing the domain.
            if (currentWork >= minWork || nY < 2) {
                nA = 1;
                nB = 1;
                // Reduce the number of divisions along the Z-axis if possible
                while (currentWork >= 2 * minWork) {
                    nC -= 2;
                    currentWork = nC / div / threadCount;
                }

            } else {
                if (nY % 2 != 0) {
                    nY--;
                }
                nB = nY;
                div = 4;
                currentWork = nB * nC / div / threadCount;
                // If we have 4 * threadCount * minWork chunks, stop dividing the domain.
                if (currentWork >= minWork || nX < 2) {
                    nA = 1;
                    while (currentWork >= 2 * minWork) {
                        nB -= 2;
                        currentWork = nB * nC / div / threadCount;
                    }
                } else {
                    if (nX % 2 != 0) {
                        nX--;
                    }
                    nA = nX;
                    div = 8;
                    currentWork = nA * nB * nC / div / threadCount;
                    while (currentWork >= 2 * minWork) {
                        nA -= 2;
                        currentWork = nA * nB * nC / div / threadCount;
                    }
                }
            }
            nAB = nA * nB;
            nCells = nAB * nC;
            nWork = nCells / div;
        } else {
            nA = 1;
            nB = 1;
            nC = 1;
            nAB = 1;
            nCells = 1;
            nWork = 1;
        }

        logger.fine(String.format(" nA %d nB %d nC %d nWork %d", nA, nB, nC, nWork));

        workA = new int[nWork];
        workB = new int[nWork];
        workC = new int[nWork];
        actualA = new int[nWork];
        actualB = new int[nWork];
        actualC = new int[nWork];
        actualCount = new int[nWork];
        int index = 0;
        for (int h = 0; h < nA; h += 2) {
            for (int k = 0; k < nB; k += 2) {
                for (int l = 0; l < nC; l += 2) {
                    workA[index] = h;
                    workB[index] = k;
                    workC[index++] = l;
                }
            }
        }
        select = new boolean[nSymm][nAtoms];
        for (int i = 0; i < nSymm; i++) {
            Arrays.fill(select[i], true);
        }
        cellList = new int[nSymm][nAtoms];
        cellIndex = new int[nSymm][nAtoms];
        cellOffset = new int[nSymm][nAtoms];
        cellStart = new int[nSymm][nCells];
        cellCount = new int[nSymm][nCells];
    }

    /**
     * <p>Getter for the field
     * <code>grid</code>.</p>
     *
     * @return an array of double.
     */
    public double[] getGrid() {
        return grid;
    }

    public void setGridBuffer(DoubleBuffer grid) {
        gridBuffer = grid;
    }

    /**
     * <p>getNsymm</p>
     *
     * @return a int.
     */
    public int getNsymm() {
        return nSymm;
    }

    /**
     * <p>setDensityLoop</p>
     *
     * @param loops an array of
     * {@link ffx.potential.nonbonded.SpatialDensityLoop} objects.
     */
    public void setDensityLoop(SpatialDensityLoop loops[]) {
        spatialDensityLoop = loops;
    }

    private class GridInitLoop extends IntegerForLoop {

        private final IntegerSchedule schedule = IntegerSchedule.fixed();
        // Extra padding to avert cache interference.
        long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        long pad8, pad9, pada, padb, padc, padd, pade, padf;

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void run(int lb, int ub) {
            if (gridBuffer != null) {
                //if (grid != null) {
                for (int i = lb; i <= ub; i++) {
                    //grid[i] = initValue;
                    gridBuffer.put(i, initValue);
                }
            }
        }
    }

    /**
     * <p>Setter for the field
     * <code>initValue</code>.</p>
     *
     * @param initValue a double.
     */
    public void setInitValue(double initValue) {
        this.initValue = initValue;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        int ti = getThreadIndex();
        int actualWork1 = actualWork - 1;
        SpatialDensityLoop loop = spatialDensityLoop[ti];

        /**
         * This lets the same SpatialDensityLoops be used with different
         * SpatialDensityRegions.
         */
        loop.setNsymm(nSymm);

        try {
            execute(0, gridSize - 1, gridInitLoop);
            execute(0, actualWork1, loop.setOctant(0));
            // Fractional chunks along the C-axis.
            if (nC > 1) {
                execute(0, actualWork1, loop.setOctant(1));
                // Fractional chunks along the B-axis.
                if (nB > 1) {
                    execute(0, actualWork1, loop.setOctant(2));
                    execute(0, actualWork1, loop.setOctant(3));
                    // Fractional chunks along the A-axis.
                    if (nA > 1) {
                        execute(0, actualWork1, loop.setOctant(4));
                        execute(0, actualWork1, loop.setOctant(5));
                        execute(0, actualWork1, loop.setOctant(6));
                        execute(0, actualWork1, loop.setOctant(7));
                    }
                }
            }
        } catch (Exception e) {
            String message = " Exception in SpatialDensityRegion.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Select atoms that should be assigned to cells. The default is to include
     * all atoms, which is set up in the constructor. This function should be
     * over-ridden by subclasses that want finer control.
     */
    public void selectAtoms() {
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells. This is very fast;
     * there is little to be gained from parallelizing it at this point.
     */
    public void assignAtomsToCells() {
        // Call the selectAtoms method of subclasses.
        selectAtoms();

        // Zero out the cell counts.
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int cellIndexs[] = cellIndex[iSymm];
            final int cellCounts[] = cellCount[iSymm];
            final int cellStarts[] = cellStart[iSymm];
            final int cellLists[] = cellList[iSymm];
            final int cellOffsets[] = cellOffset[iSymm];
            final boolean selected[] = select[iSymm];
            for (int i = 0; i < nCells; i++) {
                cellCounts[i] = 0;
            }

            // Convert to fractional coordinates.
            final double xyz[][] = coordinates[iSymm];
            final double x[] = xyz[0];
            final double y[] = xyz[1];
            final double z[] = xyz[2];
            crystal.toFractionalCoordinates(nAtoms, x, y, z, xf, yf, zf);

            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i < nAtoms; i++) {
                if (!selected[i]) {
                    continue;
                }

                double xu = xf[i];
                double yu = yf[i];
                double zu = zf[i];

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
                int a = (int) floor(xu * nA);
                int b = (int) floor(yu * nB);
                int c = (int) floor(zu * nC);

                // Check to make sure a, b and c are less than nA, nB and nC, respectively.
                if (a >= nA) {
                    a = nA - 1;
                }
                if (b >= nB) {
                    b = nB - 1;
                }
                if (c >= nC) {
                    c = nC - 1;
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
                if (!selected[i]) {
                    continue;
                }
                final int index = cellIndexs[i];
                cellLists[cellStarts[index]++] = i;
            }

            // Redefine the starting indices again.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }

            // Loop over work chunks and get rid of empty chunks.
            actualWork = 0;
            for (int icell = 0; icell < nWork; icell++) {
                int ia = workA[icell];
                int ib = workB[icell];
                int ic = workC[icell];
                int ii = count(ia, ib, ic);
                // Fractional chunks along the C-axis.
                if (nC > 1) {
                    ii += count(ia, ib, ic + 1);
                    // Fractional chunks along the B-axis.
                    if (nB > 1) {
                        ii += count(ia, ib + 1, ic);
                        ii += count(ia, ib + 1, ic + 1);
                        // Fractional chunks along the A-axis.
                        if (nA > 1) {
                            ii += count(ia + 1, ib, ic);
                            ii += count(ia + 1, ib, ic + 1);
                            ii += count(ia + 1, ib + 1, ic);
                            ii += count(ia + 1, ib + 1, ic + 1);
                        }
                    }
                }

                // If there is work in this chunk, include it.
                if (ii > 0) {
                    actualA[actualWork] = ia;
                    actualB[actualWork] = ib;
                    actualC[actualWork] = ic;
                    actualCount[actualWork++] = ii;
                }
            }

            if (logger.isLoggable(Level.FINEST)) {
                logger.fine(String.format(" Empty chunks: %d out of %d.",
                        nWork - actualWork, nWork));
            }
        }
    }

    private int count(int ia, int ib, int ic) {
        int count = 0;
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int index = index(ia, ib, ic);
            final int start = cellStart[iSymm][index];
            final int stop = start + cellCount[iSymm][index];
            count += (stop - start);
        }
        return count;
    }

    /**
     * <p>index</p>
     *
     * @param ia a int.
     * @param ib a int.
     * @param ic a int.
     * @return a int.
     */
    public int index(int ia, int ib, int ic) {
        return ia + ib * nA + ic * nAB;
    }
}
