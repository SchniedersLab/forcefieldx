/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

import static java.lang.Math.round;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.Resolution;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Real3DParallel;
import ffx.potential.bonded.Atom;

/**
 */
public class CrystalReciprocalSpace {

    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private final Crystal crystal;
    private final Resolution resolution;
    private final int nSymm;
    private final Atom atoms[];
    private final double xf[];
    private final double yf[];
    private final double zf[];
    private final int nAtoms;
    private final int fftX, fftY, fftZ;
    private final int ngridx, ngridy, ngridz;
    private static final double arad = 2.0;
    private final int nfftTotal;
    private final double densityGrid[];
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
    private final int nWork;
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
     * The cell indices of each atom along a A-axis.
     */
    private final int workA[];
    /**
     * The cell indices of each atom along a B-axis.
     */
    private final int workB[];
    /**
     * The cell indices of each atom along a C-axis.
     */
    private final int workC[];
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
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final PermanentDensityRegion permanentDensity;
    private final Real3DParallel realFFT3D;

    /**
     * Crystal Reciprocal Space.
     */
    public CrystalReciprocalSpace(Crystal crystal, Resolution resolution, Atom atoms[],
            int nAtoms, ParallelTeam parallelTeam) {
        this.crystal = crystal;
        this.resolution = resolution;
        this.atoms = atoms;
        this.nAtoms = nAtoms;
        this.parallelTeam = parallelTeam;
        threadCount = parallelTeam.getThreadCount();

        double density = 2.0 * resolution.sampling_limit();
        double res = resolution.res_limit();

        // Set default FFT grid size from unit cell dimensions.
        int nX = (int) Math.floor(crystal.a * density / res) + 1;
        int nY = (int) Math.floor(crystal.b * density / res) + 1;
        int nZ = (int) Math.floor(crystal.c * density / res) + 1;
        if (nX % 2 != 0) {
            nX += 1;
        }
        // Use preferred dimensions.
        while (!Complex.preferredDimension(nX)) {
            nX += 2;
        }
        while (!Complex.preferredDimension(nY)) {
            nY += 1;
        }
        while (!Complex.preferredDimension(nZ)) {
            nZ += 1;
        }

        fftX = nX;
        fftY = nY;
        fftZ = nZ;
        // number of grid points to sample density on
        ngridx = (int) Math.floor(arad * nX / crystal.a) + 1;
        ngridy = (int) Math.floor(arad * nY / crystal.b) + 1;
        ngridz = (int) Math.floor(arad * nZ / crystal.c) + 1;
        nfftTotal = (fftX + 2) * fftY * fftZ;
        // this.nSymm = crystal.spaceGroup.symOps.size();
        this.nSymm = 1;
        densityGrid = new double[nfftTotal];
        /**
         * Chop up the 3D unit cell domain into fractional coordinate chunks to
         * allow multiple threads to put charge density onto the grid without
         * needing the same grid point. First, we partition the X-axis, then
         * the Y-axis, and finally the Z-axis if necesary.
         */
        nX = 1;
        nY = 1;
        nZ = 1;
        int div = 1;
        if (threadCount > 1) {
            nZ = fftZ / (ngridz * 2);
            if (nZ % 2 != 0) {
                nZ--;
            }
            nC = nZ;
            div *= 2;
            // If we have 2 * threadCount chunks, stop dividing the domain.
            if (nC / threadCount >= div) {
                nA = 1;
                nB = 1;
            } else {
                nY = fftY / (ngridy * 2);
                if (nY % 2 != 0) {
                    nY--;
                }
                nB = nY;
                div *= 2;
                // If we have 4 * threadCount chunks, stop dividing the domain.
                if (nB * nC / threadCount >= div) {
                    nA = 1;
                } else {
                    nX = fftX / (ngridx * 2);
                    if (nX % 2 != 0) {
                        nX--;
                    }
                    nA = nX;
                    div *= 2;
                }
            }
            nAB = nA * nB;
            nCells = nAB * nC;
            nWork = nA * nB * nC / div;
        } else {
            nA = 1;
            nB = 1;
            nC = 1;
            nAB = 1;
            nCells = 1;
            nWork = 1;
        }
        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append(String.format(" Form Factor Width:         (%d,%d,%d)\n", ngridx, ngridy, ngridz));
            sb.append(String.format(" Grid density:               %8.3f\n", density));
            sb.append(String.format(" Grid dimensions:           (%d,%d,%d)\n", fftX, fftY, fftZ));
            sb.append(String.format(" Grid chunks per thread:    (%dx%dx%d)/%d = %d\n",
                    nA, nB, nC, div, nWork));
            logger.info(sb.toString());
        }

        workA = new int[nWork];
        workB = new int[nWork];
        workC = new int[nWork];
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
        cellList = new int[nSymm][nAtoms];
        cellIndex = new int[nSymm][nAtoms];
        cellOffset = new int[nSymm][nAtoms];
        cellStart = new int[nSymm][nCells];
        cellCount = new int[nSymm][nCells];
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
        xf = new double[nAtoms];
        yf = new double[nAtoms];
        zf = new double[nAtoms];

        realFFT3D = new Real3DParallel(fftX, fftY, fftZ, parallelTeam);
        permanentDensity = new PermanentDensityRegion();
    }

    /**
     * Note that the Java function "signum" and the FORTRAN version have
     * different definitions for an argument of zero.
     * <p>
     * JAVA: Math.signum(0.0) == 0.0
     * <p>
     * FORTRAN: signum(0.0) .eq. 1.0
     */
    public void permanent() {
        assignAtomsToCells();
        try {
            parallelTeam.execute(permanentDensity);
            realFFT3D.fft(densityGrid);
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }
    private static double toSeconds = 0.000000001;

    public double getNfftX() {
        return fftX;
    }

    public double getNfftY() {
        return fftY;
    }

    public double getNfftZ() {
        return fftZ;
    }

    private class PermanentDensityRegion extends ParallelRegion {

        private final GridInitLoop gridInitLoop;
        private final PermanentDensityLoop permanentDensityLoop[];

        public PermanentDensityRegion() {
            gridInitLoop = new GridInitLoop();
            permanentDensityLoop = new PermanentDensityLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                permanentDensityLoop[i] = new PermanentDensityLoop();
            }
        }

        @Override
        public void run() {
            int ti = getThreadIndex();
            int work1 = nWork - 1;
            PermanentDensityLoop thisLoop = permanentDensityLoop[ti];
            try {
                execute(0, nfftTotal - 1, gridInitLoop);
                execute(0, work1, thisLoop.setOctant(0));
                // Fractional chunks along the C-axis.
                if (nC > 1) {
                    execute(0, work1, thisLoop.setOctant(1));
                    // Fractional chunks along the B-axis.
                    if (nB > 1) {
                        execute(0, work1, thisLoop.setOctant(2));
                        execute(0, work1, thisLoop.setOctant(3));
                        // Fractional chunks along the A-axis.
                        if (nA > 1) {
                            execute(0, work1, thisLoop.setOctant(4));
                            execute(0, work1, thisLoop.setOctant(5));
                            execute(0, work1, thisLoop.setOctant(6));
                            execute(0, work1, thisLoop.setOctant(7));
                        }
                    }
                }
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class GridInitLoop extends IntegerForLoop {

            private final IntegerSchedule schedule = IntegerSchedule.fixed();

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    densityGrid[i] = 0.0;
                }
            }
        }

        private class PermanentDensityLoop extends IntegerForLoop {

            private final double r00;
            private final double r01;
            private final double r02;
            private final double r10;
            private final double r11;
            private final double r12;
            private final double r20;
            private final double r21;
            private final double r22;
            private int octant = 0;
            private final IntegerSchedule schedule = IntegerSchedule.fixed();
            // 128 bytes of extra padding to avert cache interference.
            private long p0, p1, p2, p3, p4, p5, p6, p7;
            private long p8, p9, pa, pb, pc, pd, pe, pf;

            public PermanentDensityLoop() {
                r00 = crystal.recip[0][0];
                r01 = crystal.recip[0][1];
                r02 = crystal.recip[0][2];
                r10 = crystal.recip[1][0];
                r11 = crystal.recip[1][1];
                r12 = crystal.recip[1][2];
                r20 = crystal.recip[2][0];
                r21 = crystal.recip[2][1];
                r22 = crystal.recip[2][2];
            }

            public PermanentDensityLoop setOctant(int octant) {
                this.octant = octant;
                return this;
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(int lb, int ub) {
                // Loop over work cells
                for (int icell = lb; icell <= ub; icell++) {
                    int ia = workA[icell];
                    int ib = workB[icell];
                    int ic = workC[icell];
                    switch (octant) {
                        // Case 0 -> In place.
                        case 0:
                            gridCell(ia, ib, ic);
                            break;
                        // Case 1: Step along the C-axis.
                        case 1:
                            gridCell(ia, ib, ic + 1);
                            break;
                        // Case 2 & 3: Step along the B-axis.
                        case 2:
                            gridCell(ia, ib + 1, ic);
                            break;
                        case 3:
                            gridCell(ia, ib + 1, ic + 1);
                            break;
                        // Case 4-7: Step along the A-axis.
                        case 4:
                            gridCell(ia + 1, ib, ic);
                            break;
                        case 5:
                            gridCell(ia + 1, ib, ic + 1);
                            break;
                        case 6:
                            gridCell(ia + 1, ib + 1, ic);
                            break;
                        case 7:
                            gridCell(ia + 1, ib + 1, ic + 1);
                            break;
                        default:
                            String message = "Programming error in PermanentDensityLoop.\n";
                            logger.severe(message);
                            System.exit(-1);
                    }
                }
            }

            private void gridCell(int ia, int ib, int ic) {
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    final int pairList[] = cellList[iSymm];
                    final int index = ia + ib * nA + ic * nAB;
                    final int start = cellStart[iSymm][index];
                    final int stop = start + cellCount[iSymm][index];
                    for (int i = start; i < stop; i++) {
                        int n = pairList[i];
                        gridPermanent(iSymm, n);
                    }
                }
            }

            private void gridPermanent(int iSymm, int n) {

                FormFactor atomff = new FormFactor(atoms[n]);
                final double xyz[] = atoms[n].getXYZ();
                final double xi = xyz[0];
                final double yi = xyz[1];
                final double zi = xyz[2];

                // Logic to loop within the cutoff box.

                final double wx = xi * r00 + yi * r10 + zi * r20;
                final double ux = wx - round(wx) + 0.5;
                final double frx = fftX * ux;
                final int ifrx = (int) frx;

                final double wy = xi * r01 + yi * r11 + zi * r21;
                final double uy = wy - round(wy) + 0.5;
                final double fry = fftY * uy;
                final int ifry = (int) fry;

                final double wz = xi * r02 + yi * r12 + zi * r22;
                final double uz = wz - round(wz) + 0.5;
                final double frz = fftZ * uz;
                final int ifrz = (int) frz;

                for (int ix = ifrx - ngridx; ix <= ifrx + ngridx; ix++) {
                    for (int iy = ifry - ngridy; iy <= ifry + ngridy; iy++) {
                        for (int iz = ifrz - ngridz; iz <= ifrz + ngridz; iz++) {
                            // apply minimum image
                            int gix = ix;
                            if (gix < 0) {
                                gix = fftX + ix;
                            }
                            if (gix >= fftX) {
                                gix = ix - fftX;
                            }
                            int giy = iy;
                            if (giy < 0) {
                                giy = fftY + iy;
                            }
                            if (giy >= fftY) {
                                giy = iy - fftY;
                            }
                            int giz = iz;
                            if (giz < 0) {
                                giz = fftZ + iz;
                            }
                            if (giz >= fftZ) {
                                giz = iz - fftZ;
                            }

                            double xc[] = new double[3];
                            double xf[] = new double[3];
                            xf[0] = gix / fftX;
                            xf[1] = giy / fftY;
                            xf[2] = giz / fftZ;
                            crystal.toCartesianCoordinates(xf, xc);

                            int index = giz * (fftY * fftX * 2) + giy * (fftX * 2) + gix * 2;
                            densityGrid[index] += atomff.rho(xc);
                        }
                    }
                }
            }
        }
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells. This is very fast;
     * there is little to be gained from parallelizing it at this point.
     */
    private void assignAtomsToCells() {
        // Zero out the cell counts.
        for (int iSymm = 0; iSymm
                < nSymm; iSymm++) {
            final int cellIndexs[] = cellIndex[iSymm];
            final int cellCounts[] = cellCount[iSymm];
            final int cellStarts[] = cellStart[iSymm];
            final int cellLists[] = cellList[iSymm];
            final int cellOffsets[] = cellOffset[iSymm];

            for (int i = 0; i
                    < nCells; i++) {
                cellCounts[i] = 0;
            }

            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i
                    < nAtoms; i++) {
                final double xc[] = atoms[i].getXYZ();
                crystal.toFractionalCoordinates(xc[0], xc[1], xc[2],
                        xf[i], yf[i], zf[i]);

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
                final int a = (int) Math.floor(xu * nA);
                final int b = (int) Math.floor(yu * nB);
                final int c = (int) Math.floor(zu * nC);

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
            } // Define the starting indices.
            cellStarts[0] = 0;

            for (int i = 1; i
                    < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
            // Move atom locations into a list ordered by cell.
            for (int i = 0; i
                    < nAtoms; i++) {
                final int index = cellIndexs[i];
                cellLists[cellStarts[index]++] = i;
            } // Redefine the starting indices again.
            cellStarts[0] = 0;

            for (int i = 1; i
                    < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
        }
    }
}
