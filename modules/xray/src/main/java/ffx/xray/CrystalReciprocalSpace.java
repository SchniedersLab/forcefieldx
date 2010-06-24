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
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

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
import ffx.xray.RefinementMinimize.RefinementMode;

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
        public static final int GAUSSIAN = 3;
        public static final int POLYNOMIAL = 4;
    }
    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private static double toSeconds = 0.000000001;
    private final double badd;
    private double arad;
    private int aradgrid;
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;
    private boolean solvent = false;
    protected int solventmodel;
    private boolean use_3g = true;
    protected double solvent_a;
    protected double solvent_b;
    private final int nSymm;
    private final Atom atoms[];
    private final int nAtoms;
    // not final for purposes of finite differences
    private double coordinates[][][];
    private final FormFactor ffactors[];
    private final int fftX, fftY, fftZ;
    private final double fftscale;
    private final int complexFFT3DSpace;
    private final int halfFFTX, halfFFTY, halfFFTZ;
    private final double densityGrid[];
    protected double solventGrid[];
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final SpatialDensityRegion spatialDensityRegion;
    private final SpatialDensityRegion bulkSpatialDensityRegion;
    private final AtomicDensityLoop atomicDensityLoops[];
    private final AtomicGradientRegion atomicGradientRegion;
    private final SolventDensityLoop solventDensityLoops[];
    private final SolventDensityLoop bulkSolventDensityLoops[];
    private final SolventGradientRegion solventGradientRegion;
    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;

    /**
     * Crystal Reciprocal Space.
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, false, SolventModel.POLYNOMIAL);
    }

    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, solventmask, SolventModel.POLYNOMIAL);
    }

    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask, int solventmodel) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.resolution = reflectionlist.resolution;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;
        this.solvent = solventmask;
        this.solventmodel = solventmodel;
        this.nSymm = 1;
        threadCount = parallelTeam.getThreadCount();
        // necssary for the bulksolvent expansion!
        int bulknsym = crystal.spaceGroup.symOps.size();

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
        densityGrid = new double[complexFFT3DSpace];

        String solventname = "none";
        if (solvent) {
            switch (solventmodel) {
                case SolventModel.BINARY:
                    solventname = "binary";
                    solvent_a = 1.0;
                    solvent_b = 1.0;
                    solventGrid = new double[complexFFT3DSpace];
                    break;
                case SolventModel.GAUSSIAN:
                    solventname = "Gaussian";
                    solvent_a = 11.5;
                    solvent_b = 0.55;
                    solventGrid = new double[complexFFT3DSpace];
                    break;
                case SolventModel.POLYNOMIAL:
                    solventname = "polynomial switch";
                    solvent_a = 0.0;
                    solvent_b = 0.8;
                    solventGrid = new double[complexFFT3DSpace];
                    break;
                default:
                    break;
            }
        }

        // set up form factors
        if (solvent) {
            badd = 0.0;
        } else {
            badd = 2.0;
        }
        ffactors = new FormFactor[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            ffactors[i] = new FormFactor(atoms[i], use_3g, badd);
        }

        // determine number of grid points to sample density on
        arad = -1.0;
        for (Atom a : atoms) {
            double vdwr = a.getVDWType().radius * 0.5;
            if (!solvent) {
                arad = Math.max(arad, a.getFormFactorWidth());
            } else {
                switch (solventmodel) {
                    case SolventModel.BINARY:
                        arad = Math.max(arad, vdwr + solvent_a + 0.2);
                        break;
                    case SolventModel.GAUSSIAN:
                        arad = Math.max(arad, vdwr * solvent_b + 2.0);
                        break;
                    case SolventModel.POLYNOMIAL:
                        arad = Math.max(arad, vdwr + solvent_b + 0.2);
                        break;
                }
            }
        }
        aradgrid = (int) Math.floor(arad * nX / crystal.a) + 1;

        // local copy of coordinates - note use of bulknsym
        coordinates = new double[bulknsym][3][nAtoms];
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        for (int i = 0; i < bulknsym; i++) {
            double x[] = coordinates[i][0];
            double y[] = coordinates[i][1];
            double z[] = coordinates[i][2];
            for (int j = 0; j < nAtoms; j++) {
                Atom aj = atoms[j];
                crystal.applySymOp(aj.getXYZ(), xyz, symops.get(i));
                x[j] = xyz[0];
                y[j] = xyz[1];
                z[j] = xyz[2];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            if (solvent) {
                sb.append(String.format("\n Bulk solvent grid settings:\n"));
                sb.append(String.format(" Bulk solvent model type:       %s\n",
                        solventname));
            } else {
                sb.append(String.format("\n Atomic grid settings:\n"));
            }
            sb.append(String.format(" Form factor grid radius (radius):  %d (%8.3f)\n",
                    aradgrid, arad));
            sb.append(String.format(" Grid density:               %8.3f\n",
                    density));
            sb.append(String.format(" Grid dimensions:           (%d,%d,%d)\n",
                    fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        if (solvent) {
            int minWork = nSymm;
            spatialDensityRegion =
                    new SpatialDensityRegion(fftX, fftY, fftZ,
                    densityGrid, (aradgrid + 2) * 2, nSymm, minWork,
                    threadCount, crystal, atoms, coordinates);
            if (solventmodel != SolventModel.NONE) {
                bulkSpatialDensityRegion =
                        new BulkSolventDensityRegion(fftX, fftY, fftZ,
                        solventGrid, (aradgrid + 2) * 2, bulknsym, minWork,
                        threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
            } else {
                bulkSpatialDensityRegion = null;
            }
            if (solventmodel == SolventModel.GAUSSIAN) {
                spatialDensityRegion.setInitValue(0.0);
            } else {
                spatialDensityRegion.setInitValue(1.0);
            }

            atomicDensityLoops = null;
            solventDensityLoops = new SolventDensityLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                solventDensityLoops[i] =
                        new SolventDensityLoop(spatialDensityRegion);
            }
            spatialDensityRegion.setDensityLoop(solventDensityLoops);

            if (solventmodel != SolventModel.NONE) {
                bulkSolventDensityLoops = new SolventDensityLoop[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    bulkSolventDensityLoops[i] =
                            new SolventDensityLoop(bulkSpatialDensityRegion);
                }
                bulkSpatialDensityRegion.setDensityLoop(bulkSolventDensityLoops);
            } else {
                bulkSolventDensityLoops = null;
            }

            atomicGradientRegion = null;
            solventGradientRegion = new SolventGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS);
        } else {
            /**
             * Create nSymm pieces of work per thread; the empty pieces will
             * be removed leaving 1 piece of work per thread.
             */
            int minWork = nSymm;
            spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ,
                    densityGrid, (aradgrid + 2) * 2, nSymm, minWork,
                    threadCount, crystal, atoms, coordinates);
            bulkSpatialDensityRegion = null;
            spatialDensityRegion.setInitValue(0.0);

            atomicDensityLoops = new AtomicDensityLoop[threadCount];
            solventDensityLoops = null;
            bulkSolventDensityLoops = null;
            for (int i = 0; i < threadCount; i++) {
                atomicDensityLoops[i] =
                        new AtomicDensityLoop(spatialDensityRegion);
            }
            spatialDensityRegion.setDensityLoop(atomicDensityLoops);

            atomicGradientRegion = new AtomicGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS);
            solventGradientRegion = null;
        }
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
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

    public void setCoordinates(double coords[]) {
        assert (coords != null);
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];

        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            xyz[0] = coords[index++];
            xyz[1] = coords[index++];
            xyz[2] = coords[index++];
            for (int j = 0; j < nSymm; j++) {
                crystal.applySymOp(xyz, symxyz, symops.get(j));
                coordinates[j][0][i] = symxyz[0];
                coordinates[j][1][i] = symxyz[1];
                coordinates[j][2][i] = symxyz[2];
            }
        }
    }

    public void computeDensity(double hkldata[][]) {
        computeDensity(hkldata, false);
    }

    public void computeDensity(double hkldata[][], boolean print) {
        if (solvent) {
            if (solventmodel != SolventModel.NONE) {
                computeSolventDensity(hkldata, print);
            }
        } else {
            computeAtomicDensity(hkldata, print);
        }
    }

    public void computeAtomicGradients(double hkldata[][],
            int freer[], int flag, RefinementMode refinementmode) {
        computeAtomicGradients(hkldata, freer, flag, refinementmode, false);
    }

    public void computeAtomicGradients(double hkldata[][],
            int freer[], int flag, RefinementMode refinementmode,
            boolean print) {

        if (solvent && solventmodel == SolventModel.NONE) {
            return;
        }

        // zero out the density
        for (int i = 0; i < complexFFT3DSpace; i++) {
            densityGrid[i] = 0.0;
        }

        int nfree = 0;
        StringBuilder sb = new StringBuilder();
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        // int nsym = 1;
        Vector<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];
            if (Double.isNaN(fc[0])) {
                continue;
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
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    densityGrid[ii] += c.phase_shift(shift).re();
                    densityGrid[ii + 1] += c.phase_shift(shift).im();
                }
            }
        }
        symtime += System.nanoTime();

        long startTime = System.nanoTime();
        complexFFT3D.ifft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        startTime = System.nanoTime();
        long permanentDensityTime = 0;
        try {
            if (solvent) {
                solventGradientRegion.setRefinementMode(refinementmode);
                parallelTeam.execute(solventGradientRegion);
            } else {
                atomicGradientRegion.setRefinementMode(refinementmode);
                parallelTeam.execute(atomicGradientRegion);
            }
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Exception computing atomic gradients.";
            logger.log(Level.SEVERE, message, e);
        }

        if (solvent) {
            sb.append(String.format(" Solvent symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Solvent inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid solvent gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        } else {
            sb.append(String.format(" Atomic symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Atomic inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid atomic gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        }

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    public void computeAtomicDensity(double hkldata[][]) {
        computeAtomicDensity(hkldata, false);
    }

    public void computeAtomicDensity(double hkldata[][], boolean print) {
        StringBuilder sb = new StringBuilder();
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        long permanentDensityTime = 0;
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        long startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        // extract structure factors
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
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
            // remove Badd
            double s = Crystal.invressq(crystal, ih);
            c = c.times(scale * exp(0.25 * badd * s)).conjugate();

            fc[0] = c.re();
            fc[1] = c.im();
        }
        symtime += System.nanoTime();

        sb.append(String.format(" Atomic grid density: %8.3f\n", permanentDensityTime * toSeconds));
        sb.append(String.format(" Atomic FFT: %8.3f\n", fftTime * toSeconds));
        sb.append(String.format(" Atomic symmetry/scaling: %8.3f\n", symtime * toSeconds));

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    public void computeSolventDensity(double hkldata[][]) {
        computeSolventDensity(hkldata, false);
    }

    public void computeSolventDensity(double hkldata[][], boolean print) {
        StringBuilder sb = new StringBuilder();
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        long permanentDensityTime = 0;
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        long startTime = System.nanoTime();
        if (solventmodel == SolventModel.BINARY) {
            /*
             * need to shrink mask
             * first, copy densitygrid to df_map
             * temporarily to store original map
             */
            int nmap = densityGrid.length;
            System.arraycopy(densityGrid, 0, solventGrid, 0, nmap);

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
                                        solventGrid[ii] = 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // copy back
            System.arraycopy(solventGrid, 0, densityGrid, 0, nmap);
        }

        // copy stuff over for derivatives
        System.arraycopy(densityGrid, 0, solventGrid, 0, densityGrid.length);

        // babinet principle
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    final int ii = iComplex3D(i, j, k, fftX, fftY);
                    if (solventmodel == SolventModel.BINARY
                            || solventmodel == SolventModel.POLYNOMIAL) {
                        densityGrid[ii] = 1.0 - densityGrid[ii];
                    } else if (solventmodel == SolventModel.GAUSSIAN) {
                        densityGrid[ii] = 1.0 - exp(-solvent_a * densityGrid[ii]);
                    }
                }
            }
        }
        long expTime = System.nanoTime() - startTime;

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);

        startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        // extract structure factors
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
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

            // negative: babinet
            fc[0] = -c.re();
            fc[1] = -c.im();
        }
        symtime += System.nanoTime();

        // now expand derivative mask so it includes nearby local density
        bulkSpatialDensityRegion.assignAtomsToCells();
        long solventDensityTime = 0;
        try {
            startTime = System.nanoTime();
            parallelTeam.execute(bulkSpatialDensityRegion);
            solventDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        long solventExpTime = 0;
        if (solventmodel == SolventModel.GAUSSIAN) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        solventGrid[ii] = exp(-solvent_a * solventGrid[ii]);
                    }
                }
            }
        }
        solventExpTime = System.nanoTime() - startTime;

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(solventGrid);

        sb.append(String.format(" Solvent grid density: %8.3f\n", permanentDensityTime * toSeconds));
        sb.append(String.format(" bulk solvent exponentiation: %8.3f\n", expTime * toSeconds));
        sb.append(String.format(" Solvent FFT: %8.3f\n", fftTime * toSeconds));
        sb.append(String.format(" Solvent symmetry/scaling: %8.3f\n", symtime * toSeconds));
        sb.append(String.format(" Solvent grid expansion: %8.3f\n", solventDensityTime * toSeconds));
        sb.append(String.format(" Solvent grid expansion exponentiation: %8.3f\n", solventExpTime * toSeconds));

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    private class AtomicDensityLoop extends SpatialDensityLoop {

        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            double xyz[] = {coordinates[iSymm][0][n],
                coordinates[iSymm][1][n],
                coordinates[iSymm][2][n]};
            FormFactor atomff = ffactors[n];
            atomff.setXYZ(xyz);
            atomff.updateB(badd);
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = Math.min(aradgrid,
                    (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                int gix = Crystal.mod(ix, fftX);
                xf[0] = ix / (double) fftX;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                        int giz = Crystal.mod(iz, fftZ);
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] += atomff.rho(xc);
                    }
                }
            }
        }
    }

    private class SolventDensityLoop extends SpatialDensityLoop {

        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public SolventDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            double xyz[] = {coordinates[iSymm][0][n],
                coordinates[iSymm][1][n],
                coordinates[iSymm][2][n]};
            FormFactor atomff = ffactors[n];
            atomff.setXYZ(xyz);
            atomff.updateB(badd);
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[n].getVDWType().radius * 0.5;
            int frad = aradgrid;
            switch (solventmodel) {
                case SolventModel.BINARY:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr + solvent_a + 0.2) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.GAUSSIAN:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr * solvent_b + 2.0) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.POLYNOMIAL:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr + solvent_b + 0.2) * fftX / crystal.a) + 1);
                    break;
            }

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                int gix = Crystal.mod(ix, fftX);
                xf[0] = ix / (double) fftX;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                        int giz = Crystal.mod(iz, fftZ);
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);

                        if (solventmodel == SolventModel.BINARY) {
                            grid[ii] *= atomff.rho_binary(xc, vdwr + solvent_a);
                        } else if (solventmodel == SolventModel.GAUSSIAN) {
                            grid[ii] += atomff.rho_gauss(xc, vdwr * solvent_b);
                        } else if (solventmodel == SolventModel.POLYNOMIAL) {
                            grid[ii] *= atomff.rho_poly(xc, vdwr + solvent_a, solvent_b);
                        }
                    }
                }
            }
        }
    }

    private class AtomicGradientRegion extends ParallelRegion {

        private final AtomicGradientLoop atomicGradientLoop[];
        private RefinementMode refinementmode;

        public AtomicGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            atomicGradientLoop = new AtomicGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                atomicGradientLoop[i] = new AtomicGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
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

            final double uvw[] = new double[3];
            final double xc[] = new double[3];
            final double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    double xyz[] = {coordinates[0][0][n],
                        coordinates[0][1][n],
                        coordinates[0][2][n]};
                    FormFactor atomff = ffactors[n];
                    atomff.setXYZ(xyz);
                    atomff.updateB(0.0);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    final int dfrad = Math.min(aradgrid,
                            (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;

                    for (int ix = ifrx - dfrad; ix <= ifrx + dfrad; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                atomff.rho_grad(xc, densityGrid[ii], refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class SolventGradientRegion extends ParallelRegion {

        private final SolventGradientLoop solventGradientLoop[];
        private RefinementMode refinementmode;

        public SolventGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            solventGradientLoop = new SolventGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                solventGradientLoop[i] = new SolventGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
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

            double uvw[] = new double[3];
            double xc[] = new double[3];
            double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    double xyz[] = {coordinates[0][0][n],
                        coordinates[0][1][n],
                        coordinates[0][2][n]};
                    FormFactor atomff = ffactors[n];
                    atomff.setXYZ(xyz);
                    atomff.updateB(0.0);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    double vdwr = atoms[n].getVDWType().radius * 0.5;
                    int dfrad = aradgrid;
                    switch (solventmodel) {
                        case SolventModel.BINARY:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr + solvent_a + 0.2) * fftX / crystal.a) + 1);
                            break;
                        case SolventModel.GAUSSIAN:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr * solvent_b + 2.0) * fftX / crystal.a) + 1);
                            break;
                        case SolventModel.POLYNOMIAL:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr + solvent_b + 0.2) * fftX / crystal.a) + 1);
                            break;
                    }

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;
                    final int ifrxu = ifrx + dfrad;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;
                    final int ifryu = ifry + dfrad;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;
                    final int ifrzu = ifrz + dfrad;

                    for (int ix = ifrx - dfrad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifryu; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrzu; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);

                                if (solventmodel == SolventModel.GAUSSIAN) {
                                    atomff.rho_gauss_grad(xc, vdwr * solvent_b,
                                            densityGrid[ii] * solvent_a * solventGrid[ii],
                                            refinementmode);
                                    /*
                                    atomff.rho_gauss_grad(xc, vdwr * solvent_b,
                                    densityGrid[ii] * -solvent_a * solventGrid[ii],
                                    refinementmode);
                                     */
                                } else if (solventmodel == SolventModel.POLYNOMIAL) {
                                    atomff.rho_poly_grad(xc, vdwr + solvent_a, solvent_b,
                                            densityGrid[ii] * solventGrid[ii],
                                            refinementmode);
                                }
                            }
                        }
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
