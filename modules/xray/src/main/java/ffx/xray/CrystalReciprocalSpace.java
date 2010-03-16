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

import static java.lang.Math.exp;
import static ffx.numerics.fft.Complex3D.iComplex3D;

import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Vector;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SymOp;
import ffx.numerics.ComplexNumber;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DParallel;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;

/**
 * Structure factor calculation (including bulk solvent structure factors)
 * 
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.
 * 
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank">
 * J. A. Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank">
 * J. S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.
 */
public class CrystalReciprocalSpace {

    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private static double toSeconds = 0.000000001;
    private static double badd = 2.0;
    private static final double arad = 2.0;
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;
    private boolean solvent = false;
    private double solvent_a = 11.5;
    private double solvent_sd = 0.75;
    private final int nSymm;
    private final Atom atoms[];
    private final int nAtoms;
    private final double coordinates[][][];
    private final int fftX, fftY, fftZ;
    private final int complexFFT3DSpace;
    private final int halfFFTX, halfFFTY, halfFFTZ;
    private final int aradgrid;
    private final double densityGrid[];
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final SpatialDensityRegion spatialDensityRegion;
    private final AtomicDensityLoop atomicDensityLoops[];
    private final SolventDensityLoop solventDensityLoops[];
    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;
    private final IntegerSchedule recipSchedule;

    /**
     * Crystal Reciprocal Space.
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, false);
    }

    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.resolution = reflectionlist.resolution;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;
        this.solvent = solventmask;
        threadCount = parallelTeam.getThreadCount();

        double density = 2.0 / resolution.sampling_limit();
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
        halfFFTX = nX / 2;
        halfFFTY = nY / 2;
        halfFFTZ = nZ / 2;
        // number of grid points to sample density on
        aradgrid = (int) Math.floor(arad * nX / crystal.a) + 1;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;

        /*
         * need to rework this - currently, since all the non-masked points
         * are set to 1.0, then non-ASU points are also set to 1.0 when they
         * should be 0.0. So, we have to go over all the symops to fill in all
         * the solvent mask.
         * so we need to either work with ASUs or come up with a better idea,
         * since this will slow things down.
         */
        if (solvent) {
            this.nSymm = crystal.spaceGroup.getNumberOfSymOps();
        } else {
            this.nSymm = 1;
        }

        coordinates = new double[nSymm][3][nAtoms];
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        for (int i = 0; i < nSymm; i++) {
            double x[] = coordinates[i][0];
            double y[] = coordinates[i][1];
            double z[] = coordinates[i][2];
            for (int j = 0; j < nAtoms; j++) {
                Atom aj = atoms[j];
                double xyz[] = new double[3];
                crystal.applySymOp(atoms[j].getXYZ(), xyz, symops.get(i));
                x[j] = xyz[0];
                y[j] = xyz[1];
                z[j] = xyz[2];
            }
        }

        densityGrid = new double[complexFFT3DSpace];

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append(String.format(" Form Factor Grid Radius        %d\n",
                    aradgrid));
            sb.append(String.format(" Grid density:               %8.3f\n",
                    density));
            sb.append(String.format(" Grid dimensions:           (%d,%d,%d)\n",
                    fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        spatialDensityRegion =
                new SpatialDensityRegion(fftX, fftY, fftZ,
                densityGrid, aradgrid * 2, nSymm,
                threadCount, crystal, atoms, coordinates);

        if (solvent) {
            recipSchedule = IntegerSchedule.fixed();
            solventDensityLoops = new SolventDensityLoop[threadCount];
            atomicDensityLoops = null;
            for (int i = 0; i < threadCount; i++) {
                solventDensityLoops[i] =
                        new SolventDensityLoop(spatialDensityRegion);
            }
        } else {
            // for Fc, only fill the ASU - better to use dynamic scheduler
            recipSchedule = IntegerSchedule.dynamic();
            atomicDensityLoops = new AtomicDensityLoop[threadCount];
            solventDensityLoops = null;
            for (int i = 0; i < threadCount; i++) {
                atomicDensityLoops[i] =
                        new AtomicDensityLoop(spatialDensityRegion);
            }
        }
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
        // complexFFT3D.setRecip(createReciprocalLattice());
    }

    public void setSolventA(double a) {
        this.solvent_a = a;
    }

    public void setSolventsd(double s) {
        this.solvent_sd = s;
    }

    public void computeDensity(double hkldata[][]) {
        if (solvent) {
            computeSolventDensity(hkldata);
        } else {
            computeAtomicDensity(hkldata);
        }
    }

    public void computeAtomicDensity(double hkldata[][]) {
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        spatialDensityRegion.setDensityLoop(atomicDensityLoops);
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            long permanentDensityTime = System.nanoTime() - startTime;
            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" Grid Atomic Density: %8.3f", permanentDensityTime * toSeconds));
            }

            startTime = System.nanoTime();
            complexFFT3D.fft(densityGrid);
            long fftTime = System.nanoTime() - startTime;
            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" FFT: %8.3f", fftTime * toSeconds));
            }

            // extract structure factors
            long symtime = -System.nanoTime();
            int nsym = crystal.spaceGroup.symOps.size();
            Vector<SymOp> symops = crystal.spaceGroup.symOps;
            for (HKL ih : reflectionlist.hkllist) {
                double fc[] = hkldata[ih.index()];

                // apply symmetry
                for (int j = 0; j < nsym; j++) {
                    HKL ij = new HKL();
                    crystal.applyTransSymRot(ih, ij, symops.get(j));
                    double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                    int h = Crystal.mod(ij.h(), fftX);
                    int k = Crystal.mod(ij.k(), fftY);
                    int l = Crystal.mod(ij.l(), fftZ);

                    if (h < halfFFTX + 1) {
                        final int ii = iComplex3D(h, k, l, fftX, fftY, fftZ);
                        ComplexNumber c = new ComplexNumber(densityGrid[ii],
                                densityGrid[ii + 1]);
                        fc[0] += c.phase_shift(shift).re();
                        fc[1] += c.phase_shift(shift).im();
                    } else {
                        h = (fftX - h) % fftX;
                        k = (fftY - k) % fftY;
                        l = (fftZ - l) % fftZ;
                        final int ii = iComplex3D(h, k, l, fftX, fftY, fftZ);
                        ComplexNumber c = new ComplexNumber(densityGrid[ii],
                                -densityGrid[ii + 1]);
                        fc[0] += c.phase_shift(shift).re();
                        fc[1] += c.phase_shift(shift).im();
                    }
                }
            }

            // scale
            double scale = crystal.volume / (fftX * fftY * fftZ);
            for (HKL ih : reflectionlist.hkllist) {
                double fc[] = hkldata[ih.index()];
                ComplexNumber c = new ComplexNumber(fc[0], fc[1]);
                // remove Badd
                double s = Crystal.invressq(crystal, ih);
                c = c.times(scale * exp(0.25 * badd * s));

                fc[0] = c.conjugate().re();
                fc[1] = c.conjugate().im();
            }
            symtime += System.nanoTime();

            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" Symmetry/scaling: %8.3f", symtime * toSeconds));
            }

        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void computeSolventDensity(double hkldata[][]) {
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        spatialDensityRegion.setDensityLoop(solventDensityLoops);
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            long permanentDensityTime = System.nanoTime() - startTime;
            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" Grid Solvent Density: %8.3f", permanentDensityTime * toSeconds));
            }

            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY, fftZ);
                        densityGrid[ii] = exp(-solvent_a * densityGrid[ii]);
                    }
                }
            }

            startTime = System.nanoTime();
            complexFFT3D.fft(densityGrid);
            long fftTime = System.nanoTime() - startTime;
            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" FFT: %8.3f", fftTime * toSeconds));
            }

            // extract structure factors
            long symtime = -System.nanoTime();
            int nsym = 1;
            Vector<SymOp> symops = crystal.spaceGroup.symOps;
            for (HKL ih : reflectionlist.hkllist) {
                double fc[] = hkldata[ih.index()];

                // apply symmetry
                for (int j = 0; j < nsym; j++) {
                    HKL ij = new HKL();
                    crystal.applyTransSymRot(ih, ij, symops.get(j));
                    double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                    int h = Crystal.mod(ij.h(), fftX);
                    int k = Crystal.mod(ij.k(), fftY);
                    int l = Crystal.mod(ij.l(), fftZ);

                    if (h < halfFFTX + 1) {
                        final int ii = iComplex3D(h, k, l, fftX, fftY, fftZ);
                        ComplexNumber c = new ComplexNumber(densityGrid[ii],
                                densityGrid[ii + 1]);
                        fc[0] += c.phase_shift(shift).re();
                        fc[1] += c.phase_shift(shift).im();
                    } else {
                        h = (fftX - h) % fftX;
                        k = (fftY - k) % fftY;
                        l = (fftZ - l) % fftZ;
                        final int ii = iComplex3D(h, k, l, fftX, fftY, fftZ);
                        ComplexNumber c = new ComplexNumber(densityGrid[ii],
                                -densityGrid[ii + 1]);
                        fc[0] += c.phase_shift(shift).re();
                        fc[1] += c.phase_shift(shift).im();
                    }
                }
            }

            // scale
            double scale = crystal.volume / (fftX * fftY * fftZ);
            for (HKL ih : reflectionlist.hkllist) {
                double fc[] = hkldata[ih.index()];
                ComplexNumber c = new ComplexNumber(fc[0], fc[1]);
                c = c.times(scale);

                fc[0] = c.conjugate().re();
                fc[1] = c.conjugate().im();
            }
            symtime += System.nanoTime();

            if (logger.isLoggable(Level.INFO)) {
                logger.info(String.format(" Symmetry/scaling: %8.3f", symtime * toSeconds));
            }

        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class AtomicDensityLoop extends SpatialDensityLoop {

        @Override
        public IntegerSchedule schedule() {
            return recipSchedule;
        }

        public AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.nSymm);
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            Vector<SymOp> symops = crystal.spaceGroup.symOps;
            double xyz[] = {coordinates[iSymm][0][n],
                coordinates[iSymm][1][n],
                coordinates[iSymm][2][n]};
            FormFactor atomff = new FormFactor(atoms[n], badd, xyz);
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;

            for (int ix = ifrx - aradgrid; ix <= ifrx + aradgrid; ix++) {
                int gix = Crystal.mod(ix, fftX);
                for (int iy = ifry - aradgrid; iy <= ifry + aradgrid; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    for (int iz = ifrz - aradgrid; iz <= ifrz + aradgrid; iz++) {
                        int giz = Crystal.mod(iz, fftZ);

                        double xc[] = new double[3];
                        double xf[] = new double[3];
                        xf[0] = ix / (double) fftX;
                        xf[1] = iy / (double) fftY;
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY, fftZ);
                        densityGrid[ii] += atomff.rho(xc);
                    }
                }
            }
        }
    }

    private class SolventDensityLoop extends SpatialDensityLoop {

        @Override
        public IntegerSchedule schedule() {
            return recipSchedule;
        }

        public SolventDensityLoop(SpatialDensityRegion region) {
            super(region, region.nSymm);
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            Vector<SymOp> symops = crystal.spaceGroup.symOps;
            double xyz[] = {coordinates[iSymm][0][n],
                coordinates[iSymm][1][n],
                coordinates[iSymm][2][n]};
            FormFactor atomff = new FormFactor(atoms[n], 0.0, xyz);
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;

            for (int ix = ifrx - aradgrid; ix <= ifrx + aradgrid; ix++) {
                int gix = Crystal.mod(ix, fftX);
                for (int iy = ifry - aradgrid; iy <= ifry + aradgrid; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    for (int iz = ifrz - aradgrid; iz <= ifrz + aradgrid; iz++) {
                        int giz = Crystal.mod(iz, fftZ);

                        double xc[] = new double[3];
                        double xf[] = new double[3];
                        xf[0] = ix / (double) fftX;
                        xf[1] = iy / (double) fftY;
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY, fftZ);
                        densityGrid[ii] += atomff.rho_gauss(xc, solvent_sd);
                    }
                }
            }
        }
    }

    public double getXDim() {
        return fftX;
    }

    public double getYDim() {
        return fftY;
    }

    public double getZDim() {
        return fftZ;
    }
}
