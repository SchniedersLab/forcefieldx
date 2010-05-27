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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Vector;

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
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank">
 * J. A. Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank">
 * J. S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.
 */
public class CrystalReciprocalSpace {

    public static interface SolventModel {

        public static final int NONE = 1;
        public static final int BINARY = 2;
        public static final int POLYNOMIAL = 3;
        public static final int GAUSSIAN = 4;
    }
    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private static double toSeconds = 0.000000001;
    private static double badd = 2.0;
    private static double htwopi = Math.sqrt(2.0 * Math.PI) / 2.0;
    private double arad = 4.5;
    private int aradgrid;
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;
    protected int solventmodel;
    private boolean use_3g = true;
    protected double solvent_a;
    protected double solvent_b;
    private final int nSymm;
    private final Atom atoms[];
    private final int nAtoms;
    // not final for purposes of finite differences
    private double coordinates[][][];
    private final int fftX, fftY, fftZ;
    private final double fftscale;
    private final int complexFFT3DSpace;
    private final int halfFFTX, halfFFTY, halfFFTZ;
    private final double densityGrid[];
    protected double solventGrid[];
    protected double dfGrid[];
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final SpatialDensityRegion spatialDensityRegion;
    private final AtomicDensityLoop atomicDensityLoops[];
    private final AtomicGradientRegion atomicGradientRegion;
    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;

    /**
     * Crystal Reciprocal Space.
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, SolventModel.POLYNOMIAL);
    }

    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam, int solventmodel) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.resolution = reflectionlist.resolution;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;
        this.solventmodel = solventmodel;
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
        fftscale = crystal.volume;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;

        // number of grid points to sample density on
        aradgrid = (int) Math.floor(arad * nX / crystal.a) + 1;
        densityGrid = new double[complexFFT3DSpace];
        String solventmodeltype = "none";
        switch (solventmodel) {
            case SolventModel.NONE:
                solventGrid = null;
                dfGrid = null;
                break;
            case SolventModel.BINARY:
                solventmodeltype = "binary";
                solvent_a = 1.0;
                solvent_b = 1.0;
                solventGrid = new double[complexFFT3DSpace];
                dfGrid = null;
                break;
            case SolventModel.POLYNOMIAL:
                solventmodeltype = "polynomial";
                solvent_b = 0.8;
                solventGrid = new double[complexFFT3DSpace];
                dfGrid = new double[complexFFT3DSpace];
                break;
            case SolventModel.GAUSSIAN:
                solventmodeltype = "Gaussian";
                solvent_a = 11.5;
                solventGrid = new double[complexFFT3DSpace];
                dfGrid = new double[complexFFT3DSpace];
                break;
            default:
                break;
        }

        // symmetry is applied in reciprocal space, so only go over identity symop
        this.nSymm = 1;
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

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format(" bulk solvent type:      %s\n",
                    solventmodeltype));
            sb.append(String.format(" Form Factor Grid Radius        %d\n",
                    aradgrid));
            sb.append(String.format(" Grid density:               %8.3f\n",
                    density));
            sb.append(String.format(" Grid dimensions:           (%d,%d,%d)\n",
                    fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        /**
         * Create nSymm pieces of work per thread; the empty pieces will
         * be removed leaving 1 piece of work per thread.
         */
        int minWork = crystal.spaceGroup.getNumberOfSymOps();
        spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ,
                densityGrid, (aradgrid + 2) * 2, nSymm, minWork,
                threadCount, crystal, atoms, coordinates);
        spatialDensityRegion.setInitValue(0.0);
        atomicDensityLoops = new AtomicDensityLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            atomicDensityLoops[i] =
                    new AtomicDensityLoop(spatialDensityRegion);
        }
        atomicGradientRegion = new AtomicGradientRegion();
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
        // complexFFT3D.setRecip(createReciprocalLattice());
    }

    public void setSolventA(double a) {
        this.solvent_a = a;
    }

    public void setSolventB(double b) {
        this.solvent_b = b;
    }

    public void setUse3G(boolean use_3g) {
        this.use_3g = use_3g;
    }

    public void deltaX(int n, double delta) {
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[0] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][0][n] = symxyz[0];
        }
    }

    public void deltaY(int n, double delta) {
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[1] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][1][n] = symxyz[1];
        }
    }

    public void deltaZ(int n, double delta) {
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[2] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][2][n] = symxyz[2];
        }
    }

    public void computeDensity(double fcdata[][], double fsdata[][]) {
        computeAtomicDensity(fcdata, fsdata);
    }

    public void computeSolventMask(double fsdata[][]) {
        if (solventmodel != SolventModel.GAUSSIAN
                && solventmodel != SolventModel.POLYNOMIAL) {
            return;
        }
        // zero out the density
        Arrays.fill(dfGrid, 0.0);

        StringBuilder sb = new StringBuilder();
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        // int nsym = 1;
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        for (HKL ih : reflectionlist.hkllist) {
            double fs[] = fsdata[ih.index()];
            if (Double.isNaN(fs[0])) {
                continue;
            }

            c.re(fs[0]);
            c.im(fs[1]);
            // scale
            c = c.times(2.0 / fftscale);

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                HKL ij = new HKL();
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                // why = and not += for dfgrid?
                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    dfGrid[ii] = c.phase_shift(shift).re();
                    dfGrid[ii + 1] = -c.phase_shift(shift).im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    dfGrid[ii] = c.phase_shift(shift).re();
                    dfGrid[ii + 1] = c.phase_shift(shift).im();
                }
            }
        }
        symtime += System.nanoTime();
        sb.append(String.format(" solvent mask symmetry: %8.3f\n", symtime * toSeconds));

        long startTime = System.nanoTime();
        complexFFT3D.ifft(dfGrid);
        long fftTime = System.nanoTime() - startTime;
        sb.append(String.format(" inverse solvent mask FFT: %8.3f\n", fftTime * toSeconds));

        for (int i = 0; i < dfGrid.length; i++) {
            dfGrid[i] *= 1.001;
            dfGrid[i] = Math.max(Math.min(htwopi + dfGrid[i], 1.0), 0.0);
        }

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(dfGrid);

        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }
    }

    public void computeAtomicGradients(double fcdata[][], double fsdata[][],
            double fsmaskdata[][], int freer[], int flag) {

        // compute bulk solvent mask, needed for derivatives
        if (solventmodel == SolventModel.GAUSSIAN
                || solventmodel == SolventModel.POLYNOMIAL) {
            if (fsmaskdata == null) {
                logger.severe("bulk solvent mask not provided for gradients");
            }
            computeSolventMask(fsmaskdata);
        }

        // zero out the density
        Arrays.fill(densityGrid, 0.0);
        if (solventmodel != SolventModel.NONE) {
            Arrays.fill(solventGrid, 0.0);
        }

        int nfree = 0;
        StringBuilder sb = new StringBuilder();
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        // int nsym = 1;
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        ComplexNumber s = new ComplexNumber();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = fcdata[ih.index()];
            if (Double.isNaN(fc[0])) {
                continue;
            }
            double fs[] = null;
            if (solventmodel != SolventModel.NONE) {
                fs = fsdata[ih.index()];
            }

            // cross validation check!!!
            if (freer != null) {
                if (freer[ih.index()] == flag) {
                    nfree++;
                    continue;
                }
            }

            c.re(fc[0]);
            c.im(fc[1]);
            // scale
            c = c.times(2.0 / fftscale);
            if (solventmodel != SolventModel.NONE) {
                s.re(fs[0]);
                s.im(fs[1]);
                // scale
                s = s.times(2.0 / fftscale);
            }

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                HKL ij = new HKL();
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    densityGrid[ii] += c.phase_shift(shift).re();
                    densityGrid[ii + 1] += -c.phase_shift(shift).im();
                    if (solventmodel != SolventModel.NONE) {
                        solventGrid[ii] += s.phase_shift(shift).re();
                        solventGrid[ii + 1] += -s.phase_shift(shift).im();
                    }
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    densityGrid[ii] += c.phase_shift(shift).re();
                    densityGrid[ii + 1] += c.phase_shift(shift).im();
                    if (solventmodel != SolventModel.NONE) {
                        solventGrid[ii] += s.phase_shift(shift).re();
                        solventGrid[ii + 1] += -s.phase_shift(shift).im();
                    }
                }
            }
        }
        symtime += System.nanoTime();
        sb.append(String.format(" symmetry: %8.3f\n", symtime * toSeconds));

        long startTime = System.nanoTime();
        complexFFT3D.ifft(densityGrid);
        long fcfftTime = System.nanoTime() - startTime;
        sb.append(String.format(" inverse atomic FFT: %8.3f\n", fcfftTime * toSeconds));

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);

        if (solventmodel != SolventModel.NONE) {
            startTime = System.nanoTime();
            complexFFT3D.ifft(solventGrid);
            long fsfftTime = System.nanoTime() - startTime;
            sb.append(String.format(" inverse solvent FFT: %8.3f\n", fsfftTime * toSeconds));
        }

        startTime = System.nanoTime();
        try {
            parallelTeam.execute(atomicGradientRegion);
            long permanentDensityTime = System.nanoTime() - startTime;
            sb.append(String.format(" Grid Atomic Gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)", nfree));
        } catch (Exception e) {
            String message = "Exception computing atomic gradients.";
            logger.log(Level.SEVERE, message, e);
        }

        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }
    }

    public void computeAtomicDensity(double fcdata[][], double fsdata[][]) {
        StringBuilder sb = new StringBuilder();
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            fcdata[i][0] = fcdata[i][1] = 0.0;
        }
        if (fsdata != null) {
            for (int i = 0; i < n; i++) {
                fsdata[i][0] = fsdata[i][1] = 0.0;
            }
        }

        // set the solventdensity
        if (solventmodel == SolventModel.BINARY
                || solventmodel == SolventModel.POLYNOMIAL) {
            Arrays.fill(solventGrid, 1.0);
        } else if (solventmodel == SolventModel.GAUSSIAN) {
            Arrays.fill(solventGrid, 0.0);
        }

        spatialDensityRegion.assignAtomsToCells();
        spatialDensityRegion.setDensityLoop(atomicDensityLoops);
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            long permanentDensityTime = System.nanoTime() - startTime;
            sb.append(String.format(" Grid Atomic Density: %8.3f\n", permanentDensityTime * toSeconds));
        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        // for the binary model, need to shrink mask
        long startTime = System.nanoTime();
        if (solventmodel == SolventModel.BINARY) {
            int nmap = solventGrid.length;
            double tmpGrid[] = new double[nmap];
            System.arraycopy(solventGrid, 0, tmpGrid, 0, nmap);

            // Logic to loop within the cutoff box.
            double xyz[] = {solvent_b, solvent_b, solvent_b};
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;

            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);

                        if (solventGrid[ii] == 1.0) {
                            continue;
                        }

                        for (int iz = k - ifrz; iz <= k + ifrz; iz++) {
                            int giz = Crystal.mod(iz, fftZ);
                            for (int iy = j - ifry; iy <= j + ifry; iy++) {
                                int giy = Crystal.mod(iy, fftY);
                                for (int ix = i - ifrx; ix <= i + ifrx; ix++) {
                                    int gix = Crystal.mod(ix, fftX);

                                    final int jj = iComplex3D(gix, giy, giz, fftX, fftY);

                                    if (solventGrid[jj] == 1.0) {
                                        tmpGrid[ii] = 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // copy back
            System.arraycopy(tmpGrid, 0, solventGrid, 0, nmap);
        }

        // generate the characteristic function
        if (solventmodel != SolventModel.NONE) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        if (solventmodel == SolventModel.BINARY
                                || solventmodel == SolventModel.POLYNOMIAL) {
                            solventGrid[ii] = 1.0 - solventGrid[ii];
                        } else if (solventmodel == SolventModel.GAUSSIAN) {
                            solventGrid[ii] = 1.0 - exp(-solvent_a * densityGrid[ii]);
                        }
                    }
                }
            }
            long expTime = System.nanoTime() - startTime;
            sb.append(String.format(" bulk solvent exponentiation: %8.3f\n",
                    expTime * toSeconds));
        }

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo_unique.map");
        // mapout.write(solventGrid);

        // FFTs
        startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fcfftTime = System.nanoTime() - startTime;
        sb.append(String.format(" atomic FFT: %8.3f\n", fcfftTime * toSeconds));

        if (solventmodel != SolventModel.NONE) {
            startTime = System.nanoTime();
            complexFFT3D.fft(solventGrid);
            long fsfftTime = System.nanoTime() - startTime;
            sb.append(String.format(" solvent FFT: %8.3f\n", fsfftTime * toSeconds));
        }

        // extract structure factors
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        HKL ij = new HKL();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = fcdata[ih.index()];
            double fs[] = null;
            if (solventmodel != SolventModel.NONE
                    && fsdata != null) {
                fs = fsdata[ih.index()];
            }

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(densityGrid[ii + 1]);
                    fc[0] += c.phase_shift(shift).re();
                    fc[1] += c.phase_shift(shift).im();

                    if (solventmodel != SolventModel.NONE
                            && fsdata != null) {
                        c.re(solventGrid[ii]);
                        c.im(solventGrid[ii + 1]);
                        fs[0] += c.phase_shift(shift).re();
                        fs[1] += c.phase_shift(shift).im();
                    }
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(-densityGrid[ii + 1]);
                    fc[0] += c.phase_shift(shift).re();
                    fc[1] += c.phase_shift(shift).im();

                    if (solventmodel != SolventModel.NONE
                            && fsdata != null) {
                        c.re(solventGrid[ii]);
                        c.im(-solventGrid[ii + 1]);
                        fs[0] += c.phase_shift(shift).re();
                        fs[1] += c.phase_shift(shift).im();
                    }
                }
            }
        }

        // scale
        double scale = (fftscale) / (fftX * fftY * fftZ);
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = fcdata[ih.index()];
            c.re(fc[0]);
            c.im(fc[1]);
            // remove Badd
            double s = Crystal.invressq(crystal, ih);
            c = c.times(scale * exp(0.25 * badd * s)).conjugate();

            fc[0] = c.re();
            fc[1] = c.im();

            if (solventmodel != SolventModel.NONE
                    && fsdata != null) {
                fc = fsdata[ih.index()];
                c.re(fc[0]);
                c.im(fc[1]);
                c = c.times(scale).conjugate();
                fc[0] = -c.re();
                fc[1] = -c.im();
            }
        }
        symtime += System.nanoTime();

        if (logger.isLoggable(Level.INFO)) {
            sb.append(String.format(" Symmetry/scaling: %8.3f\n", symtime * toSeconds));
            logger.info(sb.toString());
        }
    }

    /*
    public void computeSolventDensity(double hkldata[][]) {
    StringBuffer sb = new StringBuffer();
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
    sb.append(String.format(" Grid Solvent Density: %8.3f\n", permanentDensityTime * toSeconds));
    } catch (Exception e) {
    String message = "Fatal exception evaluating solvent electron density.";
    logger.log(Level.SEVERE, message, e);
    }

    long startTime = System.nanoTime();
    if (solventmodel == SolventModel.BINARY) {
    // need to shrink mask
    int nmap = densityGrid.length;
    System.arraycopy(densityGrid, 0, df_map, 0, nmap);

    // Logic to loop within the cutoff box.
    double xyz[] = {solvent_sd, solvent_sd, solvent_sd};
    double uvw[] = new double[3];
    crystal.toFractionalCoordinates(xyz, uvw);
    final double frx = fftX * uvw[0];
    final int ifrx = (int) frx;

    final double fry = fftY * uvw[1];
    final int ifry = (int) fry;

    final double frz = fftZ * uvw[2];
    final int ifrz = (int) frz;

    for (int k = 0; k < fftZ; k++) {
    for (int j = 0; j < fftY; j++) {
    for (int i = 0; i < fftX; i++) {
    final int ii = iComplex3D(i, j, k, fftX, fftY);

    if (densityGrid[ii] == 1.0) {
    continue;
    }

    for (int iz = k - ifrz; iz <= k + ifrz; iz++) {
    int giz = Crystal.mod(iz, fftZ);
    for (int iy = j - ifry; iy <= j + ifry; iy++) {
    int giy = Crystal.mod(iy, fftY);
    for (int ix = i - ifrx; ix <= i + ifrx; ix++) {
    int gix = Crystal.mod(ix, fftX);

    final int jj = iComplex3D(gix, giy, giz, fftX, fftY);

    if (densityGrid[jj] == 1.0) {
    df_map[ii] = 1.0;
    }
    }
    }
    }
    }
    }
    }
    // copy back
    System.arraycopy(df_map, 0, densityGrid, 0, nmap);
    }

    double ngrid = 0.0;
    double bulkfrc = 0.0;
    double min, max, mean;
    for (int k = 0; k < fftZ; k++) {
    for (int j = 0; j < fftY; j++) {
    for (int i = 0; i < fftX; i++) {
    final int ii = iComplex3D(i, j, k, fftX, fftY);
    if (solventmodel == SolventModel.GAUSSIAN) {
    densityGrid[ii] = exp(-solvent_a * densityGrid[ii]);
    }
    ngrid += 1.0;
    if (densityGrid[ii] > 0.5) {
    bulkfrc += 1.0;
    }
    }
    }
    }
    bulkfrc /= ngrid;
    long expTime = System.nanoTime() - startTime;
    sb.append(String.format(" bulk solvent exponentiation: %8.3f\n",
    expTime * toSeconds));
    sb.append(String.format(" percent bulk solvent: %6.2f\n",
    bulkfrc * 100.0));

    // need to copy density for gradient calculation
    assert (df_map != null);
    int nmap = densityGrid.length;
    System.arraycopy(densityGrid, 0, df_map, 0, nmap);

    // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
    // mapout.write(densityGrid);

    startTime = System.nanoTime();
    complexFFT3D.fft(densityGrid);
    long fftTime = System.nanoTime() - startTime;
    sb.append(String.format(" FFT: %8.3f\n", fftTime * toSeconds));

    // extract structure factors
    long symtime = -System.nanoTime();
    // symmetry expansion done in density loop for the bulk solvent
    int nsym = 1;
    Vector<SymOp> symops = crystal.spaceGroup.symOps;
    ComplexNumber c = new ComplexNumber();
    HKL ij = new HKL();
    for (HKL ih : reflectionlist.hkllist) {
    double fc[] = hkldata[ih.index()];

    // apply symmetry
    for (int j = 0; j < nsym; j++) {
    crystal.applyTransSymRot(ih, ij, symops.get(j));
    double shift = Crystal.sym_phase_shift(ih, symops.get(j));

    int h = Crystal.mod(ij.h(), fftX);
    int k = Crystal.mod(ij.k(), fftY);
    int l = Crystal.mod(ij.l(), fftZ);

    if (h < halfFFTX + 1) {
    final int ii = iComplex3D(h, k, l, fftX, fftY);
    c.re(densityGrid[ii]);
    c.im(densityGrid[ii + 1]);
    fc[0] += c.phase_shift(shift).re();
    fc[1] += c.phase_shift(shift).im();
    } else {
    h = (fftX - h) % fftX;
    k = (fftY - k) % fftY;
    l = (fftZ - l) % fftZ;
    final int ii = iComplex3D(h, k, l, fftX, fftY);
    c.re(densityGrid[ii]);
    c.im(-densityGrid[ii + 1]);
    fc[0] += c.phase_shift(shift).re();
    fc[1] += c.phase_shift(shift).im();
    }
    }
    }

    // scale
    double scale = (fftscale) / (fftX * fftY * fftZ);
    for (HKL ih : reflectionlist.hkllist) {
    double fc[] = hkldata[ih.index()];
    c.re(fc[0]);
    c.im(fc[1]);
    c = c.times(scale).conjugate();

    fc[0] = c.re();
    fc[1] = c.im();
    }
    symtime += System.nanoTime();

    if (logger.isLoggable(Level.INFO)) {
    sb.append(String.format(" Symmetry/scaling: %8.3f\n", symtime * toSeconds));
    logger.info(sb.toString());
    }

    }
     */
    private class AtomicDensityLoop extends SpatialDensityLoop {

        public AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.nSymm, region.actualCount);
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            double xyz[] = {coordinates[iSymm][0][n],
                coordinates[iSymm][1][n],
                coordinates[iSymm][2][n]};
            FormFactor atomff = new FormFactor(atoms[n], use_3g, badd, xyz);
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);
            // int frad = aradgrid;
            int frad = Math.min(aradgrid,
                    (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

            double vdwr = atoms[n].getVDWType().radius * 0.5;
            boolean independentsolvent = false;
            if (solventmodel == SolventModel.BINARY) {
                independentsolvent = true;
                frad = Math.min(aradgrid,
                        (int) Math.floor(4.0 * fftX / crystal.a) + 1);
            } else if (solventmodel == SolventModel.POLYNOMIAL) {
                independentsolvent = true;
                frad = Math.min(aradgrid,
                        (int) Math.floor(3.5 * fftX / crystal.a) + 1);
            }

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;

            double xc[] = new double[3];
            double xf[] = new double[3];
            for (int ix = ifrx - frad; ix <= ifrx + frad; ix++) {
                int gix = Crystal.mod(ix, fftX);
                for (int iy = ifry - frad; iy <= ifry + frad; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    for (int iz = ifrz - frad; iz <= ifrz + frad; iz++) {
                        int giz = Crystal.mod(iz, fftZ);

                        xf[0] = ix / (double) fftX;
                        xf[1] = iy / (double) fftY;
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        densityGrid[ii] += atomff.rho(xc);

                        if (independentsolvent) {
                            if (solventmodel == SolventModel.BINARY) {
                                solventGrid[ii] *= atomff.rho_binary(xc, vdwr + solvent_a);
                            } else if (solventmodel == SolventModel.POLYNOMIAL) {
                                solventGrid[ii] *= atomff.rho_poly(xc, vdwr, solvent_b);
                            }
                        }
                    }
                }
            }
        }
    }

    /*
    private class SolventDensityLoop extends SpatialDensityLoop {

    @Override
    public IntegerSchedule schedule() {
    return recipSchedule;
    }

    public SolventDensityLoop(SpatialDensityRegion region) {
    super(region, region.nSymm, region.actualCount);
    }

    @Override
    public void gridDensity(int iSymm, int n) {
    Vector<SymOp> symops = crystal.spaceGroup.symOps;
    double xyz[] = {coordinates[iSymm][0][n],
    coordinates[iSymm][1][n],
    coordinates[iSymm][2][n]};
    FormFactor atomff = new FormFactor(atoms[n], use_3g, 0.0, xyz);
    double vdwr = atoms[n].getVDWType().radius * 0.5;
    double uvw[] = new double[3];
    crystal.toFractionalCoordinates(xyz, uvw);

    // Logic to loop within the cutoff box.
    final double frx = fftX * uvw[0];
    final int ifrx = (int) frx;

    final double fry = fftY * uvw[1];
    final int ifry = (int) fry;

    final double frz = fftZ * uvw[2];
    final int ifrz = (int) frz;

    double xc[] = new double[3];
    double xf[] = new double[3];
    for (int ix = ifrx - aradgrid; ix <= ifrx + aradgrid; ix++) {
    int gix = Crystal.mod(ix, fftX);
    for (int iy = ifry - aradgrid; iy <= ifry + aradgrid; iy++) {
    int giy = Crystal.mod(iy, fftY);
    for (int iz = ifrz - aradgrid; iz <= ifrz + aradgrid; iz++) {
    int giz = Crystal.mod(iz, fftZ);
    
    final int ii = iComplex3D(gix, giy, giz, fftX, fftY);

    xf[0] = ix / (double) fftX;
    xf[1] = iy / (double) fftY;
    xf[2] = iz / (double) fftZ;
    crystal.toCartesianCoordinates(xf, xc);

    if (solventmodel == SolventModel.BINARY) {
    densityGrid[ii] *= atomff.rho_binary(xc, vdwr + solvent_a);
    } else if (solventmodel == SolventModel.GAUSSIAN) {
    densityGrid[ii] += atomff.rho_gauss(xc, solvent_sd);
    } else if (solventmodel == SolventModel.POLYNOMIAL) {
    densityGrid[ii] *= atomff.rho_poly(xc, vdwr, solvent_sd);
    }
    }
    }
    }
    }
    }
     */
    private class AtomicGradientRegion extends ParallelRegion {

        private final AtomicGradientLoop atomicGradientLoop[];

        public AtomicGradientRegion() {
            atomicGradientLoop = new AtomicGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                atomicGradientLoop[i] = new AtomicGradientLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, atomicGradientLoop[getThreadIndex()]);
            } catch (Exception e) {
                e.printStackTrace();
                logger.severe(e.toString());
            }
        }

        private class AtomicGradientLoop extends IntegerForLoop {

            private final IntegerSchedule schedule = IntegerSchedule.dynamic(1);
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                double uvw[] = new double[3];
                double xc[] = new double[3];
                double xf[] = new double[3];
                for (int n = lb; n <= ub; n++) {
                    double xyz[] = {coordinates[0][0][n],
                        coordinates[0][1][n],
                        coordinates[0][2][n]};
                    FormFactor atomff = new FormFactor(atoms[n], use_3g, 0.0, xyz);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    // int dfrad = aradgrid;
                    int dfrad = Math.min(aradgrid,
                            (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);
                    if (solventmodel == SolventModel.GAUSSIAN
                            || solventmodel == SolventModel.POLYNOMIAL) {
                        dfrad = Math.min(aradgrid,
                                (int) Math.floor(3.5 * fftX / crystal.a) + 1);
                    }

                    double vdwr = atoms[n].getVDWType().radius * 0.5;

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;

                    for (int ix = ifrx - dfrad; ix <= ifrx + dfrad; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
                                int giz = Crystal.mod(iz, fftZ);

                                xf[0] = ix / (double) fftX;
                                xf[1] = iy / (double) fftY;
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                if (solventmodel == SolventModel.GAUSSIAN) {
                                    atomff.rho_grad(xc, densityGrid[ii], solventGrid[ii] * -solvent_a * dfGrid[ii]);
                                } else {
                                    atomff.rho_grad(xc, densityGrid[ii]);
                                    if (solventmodel == SolventModel.POLYNOMIAL) {
                                        atomff.rho_poly_grad(xc, vdwr, solvent_b, solventGrid[ii] * dfGrid[ii]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*
    private class SolventGradientRegion extends ParallelRegion {

    private final SolventGradientLoop solventGradientLoop[];

    public SolventGradientRegion() {
    solventGradientLoop = new SolventGradientLoop[threadCount];
    for (int i = 0; i < threadCount; i++) {
    solventGradientLoop[i] = new SolventGradientLoop();
    }
    }

    @Override
    public void run() {
    try {
    execute(0, nAtoms - 1, solventGradientLoop[getThreadIndex()]);
    } catch (Exception e) {
    logger.severe(e.toString());
    }
    }

    private class SolventGradientLoop extends IntegerForLoop {

    private final IntegerSchedule schedule = IntegerSchedule.dynamic(1);
    // Extra padding to avert cache interference.
    long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    long pad8, pad9, pada, padb, padc, padd, pade, padf;

    @Override
    public IntegerSchedule schedule() {
    return schedule;
    }

    @Override
    public void run(final int lb, final int ub) {
    final int dfrad = aradgrid;
    double uvw[] = new double[3];
    double xc[] = new double[3];
    double xf[] = new double[3];
    for (int n = lb; n <= ub; n++) {
    double xyz[] = {coordinates[0][0][n],
    coordinates[0][1][n],
    coordinates[0][2][n]};
    FormFactor atomff = new FormFactor(atoms[n], use_3g, 0.0, xyz);
    double vdwr = atoms[n].getVDWType().radius * 0.5;
    crystal.toFractionalCoordinates(xyz, uvw);

    // Logic to loop within the cutoff box.
    final double frx = fftX * uvw[0];
    final int ifrx = (int) frx;

    final double fry = fftY * uvw[1];
    final int ifry = (int) fry;

    final double frz = fftZ * uvw[2];
    final int ifrz = (int) frz;

    for (int ix = ifrx - dfrad; ix <= ifrx + dfrad; ix++) {
    int gix = Crystal.mod(ix, fftX);
    for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
    int giy = Crystal.mod(iy, fftY);
    for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
    int giz = Crystal.mod(iz, fftZ);

    xf[0] = ix / (double) fftX;
    xf[1] = iy / (double) fftY;
    xf[2] = iz / (double) fftZ;
    crystal.toCartesianCoordinates(xf, xc);

    final int ii = iComplex3D(gix, giy, giz, fftX, fftY);

    if (solventmodel == SolventModel.GAUSSIAN) {
    atomff.rho_gauss_grad(xc, solvent_sd,
    densityGrid[ii] * -solvent_a * df_map[ii]);
    } else if (solventmodel == SolventModel.POLYNOMIAL) {
    atomff.rho_poly_grad(xc, vdwr, solvent_sd,
    densityGrid[ii] * df_map[ii]);
    }
    }
    }
    }
    }
    }
    }
    }
     */
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
