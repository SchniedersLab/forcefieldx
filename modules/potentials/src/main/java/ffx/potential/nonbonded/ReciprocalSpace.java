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

import static java.lang.Math.*;
import static java.lang.String.format;

import static ffx.crystal.Crystal.mod;
import static ffx.numerics.UniformBSpline.*;
import static ffx.numerics.fft.Complex3D.iComplex3D;
import static ffx.potential.parameters.MultipoleType.*;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.numerics.TensorRecursion;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DCuda;
import ffx.numerics.fft.Complex3DParallel;
import ffx.numerics.fft.Real3DParallel;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;

/**
 * The Reciprocal Space class computes the reciprocal space contribution to
 * {@link ParticleMeshEwald} for the AMOEBA force field.
 *
 * <ol>
 * <li>
 * Assignment of polarizable multipole charge density to the 3D grid,
 * via b-Splines, is parallelized using a spatial decomposition.
 * </li>
 * <p>
 * <li>
 * The convolution depends on methods of the {@link Real3DParallel} and
 * {@link Complex3DParallel} classes.
 * </li>
 * <p>
 * <li>
 * Finally, the electric potential and its gradients are collected,
 * in parallel, off the grid using b-Splines.
 * </li>
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ReciprocalSpace {

    private static final Logger logger = Logger.getLogger(ReciprocalSpace.class.getName());
    private final Atom atoms[];
    private final int nAtoms;
    private final double coordinates[][][];
    private final Crystal crystal;
    private final int nSymm;
    private final double fractionalMultipole[][][];
    private final double fractionalDipole[][][];
    private final double fractionalDipolep[][][];
    private final double fractionalMultipolePhi[][];
    private final double fractionalInducedDipolePhi[][];
    private final double fractionalInducedDipolepPhi[][];
    private final int fftX, fftY, fftZ;
    private final int complexFFT3DSpace;
    private final double aewald;
    private final int bSplineOrder;
    private final int derivOrder = 3;
    private final double densityGrid[];
    private final float floatGrid[];
    private final float floatRecip[];
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final BSplineRegion bSplineRegion;
    private final SpatialDensityRegion spatialDensityRegion;
    private final PermanentDensityLoop permanentDensityLoops[];
    private final PolarizationDensityLoop polarizationDensityLoops[];
    private final PermanentPhiRegion permanentPhi;
    private final PolarizationPhiRegion polarizationPhi;
    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;
    private final boolean cudaFFT;
    private final Thread cudaThread;
    private final Complex3DCuda cudaFFT3D;

    /**
     * Reciprocal Space PME contribution.
     */
    public ReciprocalSpace(Crystal crystal, ForceField forceField,
                           double coordinates[][][], Atom atoms[], double aewald,
                           ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this.crystal = crystal;
        this.coordinates = coordinates;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.aewald = aewald;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;

        threadCount = parallelTeam.getThreadCount();
        bSplineOrder = forceField.getInteger(ForceFieldInteger.PME_ORDER, 5);
        double density = forceField.getDouble(ForceFieldDouble.PME_SPACING, 1.2);
        cudaFFT = forceField.getBoolean(ForceField.ForceFieldBoolean.CUDAFFT, false);

        // Set default FFT grid size from unit cell dimensions.
        int nX = (int) Math.floor(crystal.a * density) + 1;
        int nY = (int) Math.floor(crystal.b * density) + 1;
        int nZ = (int) Math.floor(crystal.c * density) + 1;
        if (nX % 2 != 0) {
            nX += 1;
        }
        if (nY % 2 != 0) {
            nY += 1;
        }
        if (nZ % 2 != 0) {
            nZ += 1;
        }

        while (!Complex.preferredDimension(nX)) {
            nX += 2;
        }
        while (!Complex.preferredDimension(nY)) {
            nY += 2;
        }
        while (!Complex.preferredDimension(nZ)) {
            nZ += 2;
        }

        fftX = nX;
        fftY = nY;
        fftZ = nZ;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;
        a = new double[3][3];
        nSymm = crystal.spaceGroup.getNumberOfSymOps();

        if (cudaFFT) {
            densityGrid = null;
            floatGrid = new float[complexFFT3DSpace];
            floatRecip = new float[complexFFT3DSpace / 2];
        } else {
            densityGrid = new double[complexFFT3DSpace];
            floatGrid = null;
            floatRecip = null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format(" PME B-Spline order:                   %8d\n", bSplineOrder));
            sb.append(format(" PME Grid density:                     %8.3f\n", density));
            sb.append(format(" PME Grid dimensions:             (%3d,%3d,%3d)", fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        fractionalMultipole = new double[nSymm][nAtoms][10];
        fractionalDipole = new double[nSymm][nAtoms][3];
        fractionalDipolep = new double[nSymm][nAtoms][3];
        fractionalMultipolePhi = new double[nAtoms][tensorCount];
        fractionalInducedDipolePhi = new double[nAtoms][tensorCount];
        fractionalInducedDipolepPhi = new double[nAtoms][tensorCount];
        transformMultipoleMatrix(tmm);
        transformFieldMatrix(tfm);
        bSplineRegion = new BSplineRegion();
        if (cudaFFT) {
            spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ, floatGrid, bSplineOrder, nSymm,
                                                            1, threadCount, crystal, atoms, coordinates);
        } else {
            spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ, densityGrid, bSplineOrder, nSymm,
                                                            1, threadCount, crystal, atoms, coordinates);
        }
        permanentDensityLoops = new PermanentDensityLoop[threadCount];
        polarizationDensityLoops = new PolarizationDensityLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            permanentDensityLoops[i] = new PermanentDensityLoop(spatialDensityRegion, bSplineRegion);
            polarizationDensityLoops[i] = new PolarizationDensityLoop(spatialDensityRegion, bSplineRegion);
        }
        permanentPhi = new PermanentPhiRegion(bSplineRegion);
        polarizationPhi = new PolarizationPhiRegion(bSplineRegion);

        boolean available = false;
        String recipStrategy = null;
        try {
            recipStrategy = forceField.getString(ForceField.ForceFieldString.RECIP_SCHEDULE);
            IntegerSchedule.parse(recipStrategy);
            available = true;
        } catch (Exception e) {
            available = false;
        }
        IntegerSchedule recipSchedule;
        if (available) {
            recipSchedule = IntegerSchedule.parse(recipStrategy);
            logger.info(" Convolution schedule " + recipStrategy);
        } else {
            recipSchedule = IntegerSchedule.fixed();
        }
        if (!cudaFFT) {
            complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam, recipSchedule);
            complexFFT3D.setRecip(generalizedInfluenceFunction());
            cudaFFT3D = null;
            cudaThread = null;
        } else {
            complexFFT3D = null;
            double temp[] = generalizedInfluenceFunction();
            for (int i = 0; i < complexFFT3DSpace / 2; i++) {
                floatRecip[i] = (float) temp[i];
            }
            cudaFFT3D = new Complex3DCuda(fftX, fftY, fftZ, floatGrid, floatRecip);
            cudaThread = new Thread(cudaFFT3D);
            cudaThread.setPriority(Thread.MAX_PRIORITY);
            cudaThread.start();
        }
    }

    public void computeBSplines() {
        try {
            long time = -System.nanoTime();
            parallelTeam.execute(bSplineRegion);
            time += System.nanoTime();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Compute B-Splines:      %8.3f", time * toSeconds));
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating b-Splines.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void computePermanentDensity(double globalMultipoles[][][]) {
        spatialDensityRegion.assignAtomsToCells();
        spatialDensityRegion.setDensityLoop(permanentDensityLoops);
        for (int i = 0; i < threadCount; i++) {
            permanentDensityLoops[i].setPermanent(globalMultipoles);
        }
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            long permanentDensityTime = System.nanoTime() - startTime;
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Grid Permanent Density: %8.3f", permanentDensityTime * toSeconds));
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating permanent multipole density.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void computePermanentConvolution() {
        try {
            if (cudaFFT) {
                long fftTime = -System.nanoTime();
                cudaFFT3D.convolution(floatGrid);
                fftTime += System.nanoTime();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" CUDA Convolution:       %8.3f", fftTime * toSeconds));
                }
            } else {
                long convolutionTime = -System.nanoTime();
                complexFFT3D.convolution(densityGrid);
                convolutionTime += System.nanoTime();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Java Convolution:       %8.3f", convolutionTime * toSeconds));
                }
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating permanent convolution.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void computePermanentPhi(double cartesianMultipolePhi[][]) {
        try {
            long time = -System.nanoTime();
            parallelTeam.execute(permanentPhi);
            time += System.nanoTime();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Compute Phi:            %8.3f (sec)", time * toSeconds));
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating permanent reciprocal space potential.";
            logger.log(Level.SEVERE, message, e);
        }
        fractionalToCartesianPhi(fractionalMultipolePhi, cartesianMultipolePhi);
    }

    public void computeInducedDensity(double inducedDipole[][][],
                                      double inducedDipolep[][][]) {
        for (int i = 0; i < 3; i++) {
            a[0][i] = fftX * crystal.A[i][0];
            a[1][i] = fftY * crystal.A[i][1];
            a[2][i] = fftZ * crystal.A[i][2];
        }
        spatialDensityRegion.setDensityLoop(polarizationDensityLoops);
        for (int i = 0; i < threadCount; i++) {
            polarizationDensityLoops[i].setPolarization(inducedDipole, inducedDipolep);
        }
        try {
            long time = -System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            time += System.nanoTime();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Induced Density:        %8.3f", time * toSeconds));
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating induced density.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Note that the Java function "signum" and the FORTRAN version have
     * different definitions for an argument of zero.
     * <p>
     * JAVA: Math.signum(0.0) == 0.0
     * <p>
     * FORTRAN: signum(0.0) .eq. 1.0
     */
    public void computeInducedConvolution() {
        try {
            if (cudaFFT) {
                long time = -System.nanoTime();
                cudaFFT3D.convolution(floatGrid);
                time += System.nanoTime();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" CUDA Convolution:       %8.3f", time * toSeconds));
                }
            } else {
                long time = -System.nanoTime();
                complexFFT3D.convolution(densityGrid);
                time += System.nanoTime();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Java Convolution:       %8.3f", time * toSeconds));
                }
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating induced convolution.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void computeInducedPhi(double cartesianPolarizationPhi[][],
                                  double cartesianChainRulePhi[][]) {
        try {
            long time = -System.nanoTime();
            parallelTeam.execute(polarizationPhi);
            time += System.nanoTime();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Compute Induced Phi:    %8.3f (sec)", time * toSeconds));
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating induced reciprocal space potential.";
            logger.log(Level.SEVERE, message, e);
        }
        fractionalToCartesianPhi(fractionalInducedDipolePhi, cartesianPolarizationPhi);
        fractionalToCartesianPhi(fractionalInducedDipolepPhi, cartesianChainRulePhi);
    }

    public void cartesianToFractionalDipoles(double inducedDipole[][][],
                                             double inducedDipolep[][][]) {
        for (int i = 0; i < 3; i++) {
            a[0][i] = fftX * crystal.A[i][0];
            a[1][i] = fftY * crystal.A[i][1];
            a[2][i] = fftZ * crystal.A[i][2];
        }
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            for (int i = 0; i < nAtoms; i++) {
                double in[] = inducedDipole[iSymm][i];
                double out[] = fractionalDipole[iSymm][i];
                out[0] = a[0][0] * in[0] + a[0][1] * in[1] + a[0][2] * in[2];
                out[1] = a[1][0] * in[0] + a[1][1] * in[1] + a[1][2] * in[2];
                out[2] = a[2][0] * in[0] + a[2][1] * in[1] + a[2][2] * in[2];
                in = inducedDipolep[iSymm][i];
                out = fractionalDipolep[iSymm][i];
                out[0] = a[0][0] * in[0] + a[0][1] * in[1] + a[0][2] * in[2];
                out[1] = a[1][0] * in[0] + a[1][1] * in[1] + a[1][2] * in[2];
                out[2] = a[2][0] * in[0] + a[2][1] * in[1] + a[2][2] * in[2];
            }
        }
    }

    public double[][] getFractionalMultipolePhi() {
        return fractionalMultipolePhi;
    }

    public double[][] getFractionalMultipoles() {
        return fractionalMultipole[0];
    }

    public double[][] getFractionalInducedDipolePhi() {
        return fractionalInducedDipolePhi;
    }

    public double[][] getFractionalInducedDipoles() {
        return this.fractionalDipole[0];
    }

    public double[][] getFractionalInducedDipolepPhi() {
        return fractionalInducedDipolepPhi;
    }

    public double[][] getFractionalInducedDipolesp() {
        return this.fractionalDipolep[0];
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

    private class BSplineRegion extends ParallelRegion {

        private final BSplineFillLoop bSplineFillLoop[];
        public final double splineX[][][][];
        public final double splineY[][][][];
        public final double splineZ[][][][];
        public final int initGrid[][][];

        public BSplineRegion() {
            bSplineFillLoop = new BSplineFillLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                bSplineFillLoop[i] = new BSplineFillLoop();
            }
            initGrid = new int[nSymm][nAtoms][3];
            splineX = new double[nSymm][nAtoms][bSplineOrder][derivOrder + 1];
            splineY = new double[nSymm][nAtoms][bSplineOrder][derivOrder + 1];
            splineZ = new double[nSymm][nAtoms][bSplineOrder][derivOrder + 1];
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, bSplineFillLoop[getThreadIndex()]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class BSplineFillLoop extends IntegerForLoop {

            private final double r00;
            private final double r01;
            private final double r02;
            private final double r10;
            private final double r11;
            private final double r12;
            private final double r20;
            private final double r21;
            private final double r22;
            private final double bSplineWork[][];
            private final IntegerSchedule schedule = IntegerSchedule.dynamic(10);
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            public BSplineFillLoop() {
                super();
                r00 = crystal.A[0][0];
                r01 = crystal.A[0][1];
                r02 = crystal.A[0][2];
                r10 = crystal.A[1][0];
                r11 = crystal.A[1][1];
                r12 = crystal.A[1][2];
                r20 = crystal.A[2][0];
                r21 = crystal.A[2][1];
                r22 = crystal.A[2][2];
                bSplineWork = new double[bSplineOrder][bSplineOrder];
            }

            @Override
            public void run(int lb, int ub) {
                for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                    final double x[] = coordinates[iSymOp][0];
                    final double y[] = coordinates[iSymOp][1];
                    final double z[] = coordinates[iSymOp][2];
                    final int initgridi[][] = initGrid[iSymOp];
                    final double splineXi[][][] = splineX[iSymOp];
                    final double splineYi[][][] = splineY[iSymOp];
                    final double splineZi[][][] = splineZ[iSymOp];
                    for (int i = lb; i <= ub; i++) {
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final int[] grd = initgridi[i];
                        final double wx = xi * r00 + yi * r10 + zi * r20;
                        final double ux = wx - round(wx) + 0.5;
                        final double frx = fftX * ux;
                        final int ifrx = (int) frx;
                        final double bx = frx - ifrx;
                        grd[0] = ifrx - bSplineOrder;
                        bSplineDerivatives(bx, bSplineOrder, derivOrder, splineXi[i],
                                           bSplineWork);
                        final double wy = xi * r01 + yi * r11 + zi * r21;
                        final double uy = wy - round(wy) + 0.5;
                        final double fry = fftY * uy;
                        final int ifry = (int) fry;
                        final double by = fry - ifry;
                        grd[1] = ifry - bSplineOrder;
                        bSplineDerivatives(by, bSplineOrder, derivOrder, splineYi[i],
                                           bSplineWork);
                        final double wz = xi * r02 + yi * r12 + zi * r22;
                        final double uz = wz - round(wz) + 0.5;
                        final double frz = fftZ * uz;
                        final int ifrz = (int) frz;
                        final double bz = frz - ifrz;
                        grd[2] = ifrz - bSplineOrder;
                        bSplineDerivatives(bz, bSplineOrder, derivOrder, splineZi[i],
                                           bSplineWork);
                    }
                }
            }
        }
    }

    private class PermanentDensityLoop extends SpatialDensityLoop {

        private double globalMultipoles[][][] = null;
        private final BSplineRegion bSplines;

        public PermanentDensityLoop(SpatialDensityRegion region, BSplineRegion splines) {
            super(region, region.nSymm, region.actualCount);
            this.bSplines = splines;
        }

        public void setPermanent(double globalMultipoles[][][]) {
            this.globalMultipoles = globalMultipoles;
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            final double gm[] = globalMultipoles[iSymm][n];
            final double fm[] = fractionalMultipole[iSymm][n];
            // charge
            fm[0] = gm[0];
            // dipole
            for (int j = 1; j < 4; j++) {
                fm[j] = 0.0;
                for (int k = 1; k < 4; k++) {
                    fm[j] = fm[j] + tmm[j][k] * gm[k];
                }
            }
            // quadrupole
            for (int j = 4; j < 10; j++) {
                fm[j] = 0.0;
                for (int k = 4; k < 7; k++) {
                    fm[j] = fm[j] + tmm[j][k] * gm[k];
                }
                for (int k = 7; k < 10; k++) {
                    fm[j] = fm[j] + tmm[j][k] * 2.0 * gm[k];
                }
                /**
                 * Fractional quadrupole components are pre-multiplied by a
                 * factor of 1/3 that arises in their potential.
                 */
                fm[j] = fm[j] / 3.0;
            }
            final double[][] splx = bSplines.splineX[iSymm][n];
            final double[][] sply = bSplines.splineY[iSymm][n];
            final double[][] splz = bSplines.splineZ[iSymm][n];
            final int igrd0 = bSplines.initGrid[iSymm][n][0];
            final int jgrd0 = bSplines.initGrid[iSymm][n][1];
            int k0 = bSplines.initGrid[iSymm][n][2];
            final double c = fm[t000];
            final double dx = fm[t100];
            final double dy = fm[t010];
            final double dz = fm[t001];
            final double qxx = fm[t200];
            final double qyy = fm[t020];
            final double qzz = fm[t002];
            final double qxy = fm[t110];
            final double qxz = fm[t101];
            final double qyz = fm[t011];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double v2 = splzi[2];
                final double c0 = c * v0;
                final double dx0 = dx * v0;
                final double dy0 = dy * v0;
                final double dz1 = dz * v1;
                final double qxx0 = qxx * v0;
                final double qyy0 = qyy * v0;
                final double qzz2 = qzz * v2;
                final double qxy0 = qxy * v0;
                final double qxz1 = qxz * v1;
                final double qyz1 = qyz * v1;
                final int k = mod(++k0, fftZ);
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double u2 = splyi[2];
                    final double term0 = (c0 + dz1 + qzz2) * u0 + (dy0 + qyz1) * u1 + qyy0 * u2;
                    final double term1 = (dx0 + qxz1) * u0 + qxy0 * u1;
                    final double term2 = qxx0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double dq = splxi[0] * term0 + splxi[1] * term1 + splxi[2] * term2;
                        if (cudaFFT) {
                            floatGrid[ii] += dq;
                        } else {
                            densityGrid[ii] += dq;
                        }
                    }
                }
            }
        }
    }

    private class PolarizationDensityLoop extends SpatialDensityLoop {

        private double inducedDipole[][][] = null;
        private double inducedDipolep[][][] = null;
        private final BSplineRegion bSplines;

        public PolarizationDensityLoop(SpatialDensityRegion region, BSplineRegion splines) {
            super(region, region.nSymm, region.actualCount);
            this.bSplines = splines;
        }

        public void setPolarization(double inducedDipole[][][],
                                    double inducedDipolep[][][]) {
            this.inducedDipole = inducedDipole;
            this.inducedDipolep = inducedDipolep;
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            double ind[] = inducedDipole[iSymm][n];
            final double find[] = fractionalDipole[iSymm][n];



            find[0] = a[0][0] * ind[0] + a[0][1] * ind[1] + a[0][2] * ind[2];
            find[1] = a[1][0] * ind[0] + a[1][1] * ind[1] + a[1][2] * ind[2];
            find[2] = a[2][0] * ind[0] + a[2][1] * ind[1] + a[2][2] * ind[2];
            double inp[] = inducedDipolep[iSymm][n];
            final double finp[] = fractionalDipolep[iSymm][n];
            finp[0] = a[0][0] * inp[0] + a[0][1] * inp[1] + a[0][2] * inp[2];
            finp[1] = a[1][0] * inp[0] + a[1][1] * inp[1] + a[1][2] * inp[2];
            finp[2] = a[2][0] * inp[0] + a[2][1] * inp[1] + a[2][2] * inp[2];
            final double[][] splx = bSplines.splineX[iSymm][n];
            final double[][] sply = bSplines.splineY[iSymm][n];
            final double[][] splz = bSplines.splineZ[iSymm][n];
            final int igrd0 = bSplines.initGrid[iSymm][n][0];
            final int jgrd0 = bSplines.initGrid[iSymm][n][1];
            int k0 = bSplines.initGrid[iSymm][n][2];
            final double ux = find[0];
            final double uy = find[1];
            final double uz = find[2];
            final double px = finp[0];
            final double py = finp[1];
            final double pz = finp[2];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double dx0 = ux * v0;
                final double dy0 = uy * v0;
                final double dz1 = uz * v1;
                final double px0 = px * v0;
                final double py0 = py * v0;
                final double pz1 = pz * v1;
                final int k = mod(++k0, fftZ);
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double term0 = dz1 * u0 + dy0 * u1;
                    final double term1 = dx0 * u0;
                    final double termp0 = pz1 * u0 + py0 * u1;
                    final double termp1 = px0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double dq = splxi[0] * term0 + splxi[1] * term1;
                        final double pq = splxi[0] * termp0 + splxi[1] * termp1;
                        if (cudaFFT) {
                            floatGrid[ii] += dq;
                            floatGrid[ii + 1] += pq;
                        } else {
                            densityGrid[ii] += dq;
                            densityGrid[ii + 1] += pq;
                        }
                    }
                }
            }
        }
    }

    private class PermanentPhiRegion extends ParallelRegion {

        private final FractionalPhiLoop fractionalPhiLoop[];
        private final double splineX[][][][];
        private final double splineY[][][][];
        private final double splineZ[][][][];
        private final int initgrid[][][];

        public PermanentPhiRegion(BSplineRegion bSplineRegion) {
            this.initgrid = bSplineRegion.initGrid;
            this.splineX = bSplineRegion.splineX;
            this.splineY = bSplineRegion.splineY;
            this.splineZ = bSplineRegion.splineZ;
            fractionalPhiLoop = new FractionalPhiLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fractionalPhiLoop[i] = new FractionalPhiLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, fractionalPhiLoop[getThreadIndex()]);
            } catch (Exception e) {
                logger.severe(e.toString());
                e.printStackTrace();
            }
        }

        private class FractionalPhiLoop extends IntegerForLoop {

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
                for (int n = lb; n <= ub; n++) {
                    final double[][] splx = splineX[0][n];
                    final double[][] sply = splineY[0][n];
                    final double[][] splz = splineZ[0][n];
                    final int igrd[] = initgrid[0][n];
                    final int igrd0 = igrd[0];
                    final int jgrd0 = igrd[1];
                    int k0 = igrd[2];
                    double tuv000 = 0.0;
                    double tuv100 = 0.0;
                    double tuv010 = 0.0;
                    double tuv001 = 0.0;
                    double tuv200 = 0.0;
                    double tuv020 = 0.0;
                    double tuv002 = 0.0;
                    double tuv110 = 0.0;
                    double tuv101 = 0.0;
                    double tuv011 = 0.0;
                    double tuv300 = 0.0;
                    double tuv030 = 0.0;
                    double tuv003 = 0.0;
                    double tuv210 = 0.0;
                    double tuv201 = 0.0;
                    double tuv120 = 0.0;
                    double tuv021 = 0.0;
                    double tuv102 = 0.0;
                    double tuv012 = 0.0;
                    double tuv111 = 0.0;
                    for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                        final int k = mod(++k0, fftZ);
                        int j0 = jgrd0;
                        double tu00 = 0.0;
                        double tu10 = 0.0;
                        double tu01 = 0.0;
                        double tu20 = 0.0;
                        double tu11 = 0.0;
                        double tu02 = 0.0;
                        double tu30 = 0.0;
                        double tu21 = 0.0;
                        double tu12 = 0.0;
                        double tu03 = 0.0;
                        for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                            final int j = mod(++j0, fftY);
                            int i0 = igrd0;
                            double t0 = 0.0;
                            double t1 = 0.0;
                            double t2 = 0.0;
                            double t3 = 0.0;
                            for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                                final int i = mod(++i0, fftX);
                                final int ii = iComplex3D(i, j, k, fftX, fftY);
                                double tq;
                                if (cudaFFT) {
                                    tq = floatGrid[ii];
                                } else {
                                    tq = densityGrid[ii];
                                }
                                final double splxi[] = splx[ith1];
                                t0 += tq * splxi[0];
                                t1 += tq * splxi[1];
                                t2 += tq * splxi[2];
                                t3 += tq * splxi[3];
                            }
                            final double splyi[] = sply[ith2];
                            final double u0 = splyi[0];
                            final double u1 = splyi[1];
                            final double u2 = splyi[2];
                            final double u3 = splyi[3];
                            tu00 += t0 * u0;
                            tu10 += t1 * u0;
                            tu01 += t0 * u1;
                            tu20 += t2 * u0;
                            tu11 += t1 * u1;
                            tu02 += t0 * u2;
                            tu30 += t3 * u0;
                            tu21 += t2 * u1;
                            tu12 += t1 * u2;
                            tu03 += t0 * u3;
                        }
                        final double splzi[] = splz[ith3];
                        final double v0 = splzi[0];
                        final double v1 = splzi[1];
                        final double v2 = splzi[2];
                        final double v3 = splzi[3];
                        tuv000 += tu00 * v0;
                        tuv100 += tu10 * v0;
                        tuv010 += tu01 * v0;
                        tuv001 += tu00 * v1;
                        tuv200 += tu20 * v0;
                        tuv020 += tu02 * v0;
                        tuv002 += tu00 * v2;
                        tuv110 += tu11 * v0;
                        tuv101 += tu10 * v1;
                        tuv011 += tu01 * v1;
                        tuv300 += tu30 * v0;
                        tuv030 += tu03 * v0;
                        tuv003 += tu00 * v3;
                        tuv210 += tu21 * v0;
                        tuv201 += tu20 * v1;
                        tuv120 += tu12 * v0;
                        tuv021 += tu02 * v1;
                        tuv102 += tu10 * v2;
                        tuv012 += tu01 * v2;
                        tuv111 += tu11 * v1;
                    }
                    final double out[] = fractionalMultipolePhi[n];
                    out[t000] = tuv000;
                    out[t100] = tuv100;
                    out[t010] = tuv010;
                    out[t001] = tuv001;
                    out[t200] = tuv200;
                    out[t020] = tuv020;
                    out[t002] = tuv002;
                    out[t110] = tuv110;
                    out[t101] = tuv101;
                    out[t011] = tuv011;
                    out[t300] = tuv300;
                    out[t030] = tuv030;
                    out[t003] = tuv003;
                    out[t210] = tuv210;
                    out[t201] = tuv201;
                    out[t120] = tuv120;
                    out[t021] = tuv021;
                    out[t102] = tuv102;
                    out[t012] = tuv012;
                    out[t111] = tuv111;
                }
            }
        }
    }

    private class PolarizationPhiRegion extends ParallelRegion {

        private final PolarizationPhiInducedLoop polarizationPhiInducedLoop;
        private final double splineX[][][][];
        private final double splineY[][][][];
        private final double splineZ[][][][];
        private final int initgrid[][][];

        public PolarizationPhiRegion(BSplineRegion bSplineRegion) {
            this.initgrid = bSplineRegion.initGrid;
            this.splineX = bSplineRegion.splineX;
            this.splineY = bSplineRegion.splineY;
            this.splineZ = bSplineRegion.splineZ;
            polarizationPhiInducedLoop = new PolarizationPhiInducedLoop();
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, polarizationPhiInducedLoop);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class PolarizationPhiInducedLoop extends IntegerForLoop {

            private final IntegerSchedule schedule = IntegerSchedule.dynamic(1);
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(int lb, int ub) {
                for (int n = lb; n <= ub; n++) {
                    final double[][] splx = splineX[0][n];
                    final double[][] sply = splineY[0][n];
                    final double[][] splz = splineZ[0][n];
                    final int igrd[] = initgrid[0][n];
                    final int igrd0 = igrd[0];
                    final int jgrd0 = igrd[1];
                    int k0 = igrd[2];
                    double tuv000 = 0.0;
                    double tuv100 = 0.0;
                    double tuv010 = 0.0;
                    double tuv001 = 0.0;
                    double tuv200 = 0.0;
                    double tuv020 = 0.0;
                    double tuv002 = 0.0;
                    double tuv110 = 0.0;
                    double tuv101 = 0.0;
                    double tuv011 = 0.0;
                    double tuv300 = 0.0;
                    double tuv030 = 0.0;
                    double tuv003 = 0.0;
                    double tuv210 = 0.0;
                    double tuv201 = 0.0;
                    double tuv120 = 0.0;
                    double tuv021 = 0.0;
                    double tuv102 = 0.0;
                    double tuv012 = 0.0;
                    double tuv111 = 0.0;
                    double tuv000p = 0.0;
                    double tuv100p = 0.0;
                    double tuv010p = 0.0;
                    double tuv001p = 0.0;
                    double tuv200p = 0.0;
                    double tuv020p = 0.0;
                    double tuv002p = 0.0;
                    double tuv110p = 0.0;
                    double tuv101p = 0.0;
                    double tuv011p = 0.0;
                    double tuv300p = 0.0;
                    double tuv030p = 0.0;
                    double tuv003p = 0.0;
                    double tuv210p = 0.0;
                    double tuv201p = 0.0;
                    double tuv120p = 0.0;
                    double tuv021p = 0.0;
                    double tuv102p = 0.0;
                    double tuv012p = 0.0;
                    double tuv111p = 0.0;
                    for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                        final int k = mod(++k0, fftZ);
                        int j0 = jgrd0;
                        double tu00 = 0.0;
                        double tu10 = 0.0;
                        double tu01 = 0.0;
                        double tu20 = 0.0;
                        double tu11 = 0.0;
                        double tu02 = 0.0;
                        double tu30 = 0.0;
                        double tu21 = 0.0;
                        double tu12 = 0.0;
                        double tu03 = 0.0;
                        double tu00p = 0.0;
                        double tu10p = 0.0;
                        double tu01p = 0.0;
                        double tu20p = 0.0;
                        double tu11p = 0.0;
                        double tu02p = 0.0;
                        double tu30p = 0.0;
                        double tu21p = 0.0;
                        double tu12p = 0.0;
                        double tu03p = 0.0;
                        for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                            final int j = mod(++j0, fftY);
                            int i0 = igrd0;
                            double t0 = 0.0;
                            double t1 = 0.0;
                            double t2 = 0.0;
                            double t3 = 0.0;
                            double t0p = 0.0;
                            double t1p = 0.0;
                            double t2p = 0.0;
                            double t3p = 0.0;
                            for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                                final int i = mod(++i0, fftX);
                                final int ii = iComplex3D(i, j, k, fftX, fftY);
                                double tq, tp;
                                if (cudaFFT) {
                                    tq = floatGrid[ii];
                                    tp = floatGrid[ii + 1];
                                } else {
                                    tq = densityGrid[ii];
                                    tp = densityGrid[ii + 1];
                                }
                                final double splxi[] = splx[ith1];
                                t0 += tq * splxi[0];
                                t1 += tq * splxi[1];
                                t2 += tq * splxi[2];
                                t3 += tq * splxi[3];
                                t0p += tp * splxi[0];
                                t1p += tp * splxi[1];
                                t2p += tp * splxi[2];
                                t3p += tp * splxi[3];
                            }
                            final double splyi[] = sply[ith2];
                            final double u0 = splyi[0];
                            final double u1 = splyi[1];
                            final double u2 = splyi[2];
                            final double u3 = splyi[3];
                            tu00 += t0 * u0;
                            tu10 += t1 * u0;
                            tu01 += t0 * u1;
                            tu20 += t2 * u0;
                            tu11 += t1 * u1;
                            tu02 += t0 * u2;
                            tu30 += t3 * u0;
                            tu21 += t2 * u1;
                            tu12 += t1 * u2;
                            tu03 += t0 * u3;
                            tu00p += t0p * u0;
                            tu10p += t1p * u0;
                            tu01p += t0p * u1;
                            tu20p += t2p * u0;
                            tu11p += t1p * u1;
                            tu02p += t0p * u2;
                            tu30p += t3p * u0;
                            tu21p += t2p * u1;
                            tu12p += t1p * u2;
                            tu03p += t0p * u3;
                        }
                        final double splzi[] = splz[ith3];
                        final double v0 = splzi[0];
                        final double v1 = splzi[1];
                        final double v2 = splzi[2];
                        final double v3 = splzi[3];
                        tuv000 += tu00 * v0;
                        tuv100 += tu10 * v0;
                        tuv010 += tu01 * v0;
                        tuv001 += tu00 * v1;
                        tuv200 += tu20 * v0;
                        tuv020 += tu02 * v0;
                        tuv002 += tu00 * v2;
                        tuv110 += tu11 * v0;
                        tuv101 += tu10 * v1;
                        tuv011 += tu01 * v1;
                        tuv300 += tu30 * v0;
                        tuv030 += tu03 * v0;
                        tuv003 += tu00 * v3;
                        tuv210 += tu21 * v0;
                        tuv201 += tu20 * v1;
                        tuv120 += tu12 * v0;
                        tuv021 += tu02 * v1;
                        tuv102 += tu10 * v2;
                        tuv012 += tu01 * v2;
                        tuv111 += tu11 * v1;
                        tuv000p += tu00p * v0;
                        tuv100p += tu10p * v0;
                        tuv010p += tu01p * v0;
                        tuv001p += tu00p * v1;
                        tuv200p += tu20p * v0;
                        tuv020p += tu02p * v0;
                        tuv002p += tu00p * v2;
                        tuv110p += tu11p * v0;
                        tuv101p += tu10p * v1;
                        tuv011p += tu01p * v1;
                        tuv300p += tu30p * v0;
                        tuv030p += tu03p * v0;
                        tuv003p += tu00p * v3;
                        tuv210p += tu21p * v0;
                        tuv201p += tu20p * v1;
                        tuv120p += tu12p * v0;
                        tuv021p += tu02p * v1;
                        tuv102p += tu10p * v2;
                        tuv012p += tu01p * v2;
                        tuv111p += tu11p * v1;
                    }
                    double out[] = fractionalInducedDipolePhi[n];
                    out[t000] = tuv000;
                    out[t100] = tuv100;
                    out[t010] = tuv010;
                    out[t001] = tuv001;
                    out[t200] = tuv200;
                    out[t020] = tuv020;
                    out[t002] = tuv002;
                    out[t110] = tuv110;
                    out[t101] = tuv101;
                    out[t011] = tuv011;
                    out[t300] = tuv300;
                    out[t030] = tuv030;
                    out[t003] = tuv003;
                    out[t210] = tuv210;
                    out[t201] = tuv201;
                    out[t120] = tuv120;
                    out[t021] = tuv021;
                    out[t102] = tuv102;
                    out[t012] = tuv012;
                    out[t111] = tuv111;
                    out = fractionalInducedDipolepPhi[n];
                    out[t000] = tuv000p;
                    out[t100] = tuv100p;
                    out[t010] = tuv010p;
                    out[t001] = tuv001p;
                    out[t200] = tuv200p;
                    out[t020] = tuv020p;
                    out[t002] = tuv002p;
                    out[t110] = tuv110p;
                    out[t101] = tuv101p;
                    out[t011] = tuv011p;
                    out[t300] = tuv300p;
                    out[t030] = tuv030p;
                    out[t003] = tuv003p;
                    out[t210] = tuv210p;
                    out[t201] = tuv201p;
                    out[t120] = tuv120p;
                    out[t021] = tuv021p;
                    out[t102] = tuv102p;
                    out[t012] = tuv012p;
                    out[t111] = tuv111p;
                }
            }
        }
    }

    private double[] generalizedInfluenceFunction() {

        double influenceFunction[] = new double[complexFFT3DSpace / 2];


        double bsModX[] = new double[fftX];
        double bsModY[] = new double[fftY];
        double bsModZ[] = new double[fftZ];
        int maxfft = max(max(fftX, fftY), fftZ);
        double bsArray[] = new double[maxfft];
        double c[] = new double[bSplineOrder];

        bSpline(0.0, bSplineOrder, c);
        for (int i = 1; i < bSplineOrder + 1; i++) {
            bsArray[i] = c[i - 1];
        }
        discreteFTMod(bsModX, bsArray, fftX, bSplineOrder);
        discreteFTMod(bsModY, bsArray, fftY, bSplineOrder);
        discreteFTMod(bsModZ, bsArray, fftZ, bSplineOrder);

        double r00 = crystal.A[0][0];
        double r01 = crystal.A[0][1];
        double r02 = crystal.A[0][2];
        double r10 = crystal.A[1][0];
        double r11 = crystal.A[1][1];
        double r12 = crystal.A[1][2];
        double r20 = crystal.A[2][0];
        double r21 = crystal.A[2][1];
        double r22 = crystal.A[2][2];
        int ntot = fftX * fftY * fftZ;
        double piTerm = (PI / aewald) * (PI / aewald);
        double volTerm = PI * crystal.volume;
        int nfXY = fftX * fftY;
        int nX_2 = (fftX + 1) / 2;
        int nY_2 = (fftY + 1) / 2;
        int nZ_2 = (fftZ + 1) / 2;

        for (int i = 0; i < ntot - 1; i++) {
            int kZ = (i + 1) / nfXY;
            int j = i - kZ * nfXY + 1;
            int kY = j / fftX;
            int kX = j - kY * fftX;
            int h = kX;
            int k = kY;
            int l = kZ;
            if (kX >= nX_2) {
                h -= fftX;
            }
            if (kY >= nY_2) {
                k -= fftY;
            }
            if (kZ >= nZ_2) {
                l -= fftZ;
            }
            double sX = r00 * h + r01 * k + r02 * l;
            double sY = r10 * h + r11 * k + r12 * l;
            double sZ = r20 * h + r21 * k + r22 * l;
            double sSquared = sX * sX + sY * sY + sZ * sZ;
            double term = -piTerm * sSquared;
            double expterm = 0.0;
            if (term > -50.0) {
                double denom = sSquared * volTerm * bsModX[kX] * bsModY[kY] * bsModZ[kZ];
                expterm = exp(term) / denom;
                if (crystal.aperiodic()) {
                    expterm *= (1.0 - cos(PI * crystal.a * sqrt(sSquared)));
                }
            }
            int ii = iComplex3D(kX, kY, kZ, fftX, fftY) / 2;
            influenceFunction[ii] = expterm;
        }
        /**
         *  Account for the zeroth grid point for a periodic system.
         */
        influenceFunction[0] = 0.0;
        if (crystal.aperiodic()) {
            influenceFunction[0] = 0.5 * PI / crystal.a;
        }

        return influenceFunction;
    }

    private void fractionalToCartesianPhi(double frac[][], double cart[][]) {
        for (int i = 0; i < nAtoms; i++) {
            double in[] = frac[i];
            double out[] = cart[i];
            out[0] = tfm[0][0] * in[0];
            for (int j = 1; j < 4; j++) {
                out[j] = 0.0;
                for (int k = 1; k < 4; k++) {
                    out[j] += tfm[j][k] * in[k];
                }
            }
            for (int j = 4; j < 10; j++) {
                out[j] = 0.0;
                for (int k = 4; k < 10; k++) {
                    out[j] += tfm[j][k] * in[k];
                }
            }
        }
    }

    private void transformMultipoleMatrix(double mpole_xy[][]) {
        for (int i = 0; i < 3; i++) {
            a[0][i] = fftX * crystal.A[i][0];
            a[1][i] = fftY * crystal.A[i][1];
            a[2][i] = fftZ * crystal.A[i][2];
        }
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                mpole_xy[i][j] = 0.0;
            }
        }
        // charge
        mpole_xy[0][0] = 1.0;
        // dipole
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                mpole_xy[i][j] = a[i - 1][j - 1];
            }
        }
        // quadrupole
        for (int i1 = 0; i1 < 3; i1++) {
            int k = qi1[i1];
            for (int i2 = 0; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                mpole_xy[i1 + 4][i2 + 4] = a[k][i] * a[k][j];
            }
        }
        for (int i1 = 3; i1 < 6; i1++) {
            int k = qi1[i1];
            int l = qi2[i1];
            for (int i2 = 0; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                mpole_xy[i1 + 4][i2 + 4] = a[k][i] * a[l][j] + a[k][j] * a[l][i];
            }
        }
    }

    private void transformFieldMatrix(double field_xy[][]) {
        for (int i = 0; i < 3; i++) {
            a[i][0] = fftX * crystal.A[i][0];
            a[i][1] = fftY * crystal.A[i][1];
            a[i][2] = fftZ * crystal.A[i][2];
        }
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                field_xy[i][j] = 0.0;
            }
        }
        field_xy[0][0] = 1.0;
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                field_xy[i][j] = a[i - 1][j - 1];
            }
        }
        for (int i1 = 0; i1 < 3; i1++) {
            int k = qi1[i1];
            for (int i2 = 0; i2 < 3; i2++) {
                int i = qi1[i2];
                field_xy[i1 + 4][i2 + 4] = a[k][i] * a[k][i];
            }
            for (int i2 = 3; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                field_xy[i1 + 4][i2 + 4] = 2.0 * a[k][i] * a[k][j];
            }
        }
        for (int i1 = 3; i1 < 6; i1++) {
            int k = qi1[i1];
            int n = qi2[i1];
            for (int i2 = 0; i2 < 3; i2++) {
                int i = qi1[i2];
                field_xy[i1 + 4][i2 + 4] = a[k][i] * a[n][i];
            }
            for (int i2 = 3; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                field_xy[i1 + 4][i2 + 4] = a[k][i] * a[n][j] + a[n][i] * a[k][j];
            }
        }
    }

    /**
     * Computes the modulus of the discrete Fourier Transform of
     * "bsarray" and stores it in "bsmod".
     *
     * @param bsmod
     * @param bsarray
     * @param nfft
     * @param order
     */
    private static void discreteFTMod(double bsmod[], double bsarray[],
                                      int nfft, int order) {
        /**
         * Get the modulus of the discrete Fourier fft.
         */
        double factor = 2.0 * PI / nfft;
        for (int i = 0; i < nfft; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j = 0; j < nfft; j++) {
                double arg = factor * (i * j);
                sum1 = sum1 + bsarray[j] * cos(arg);
                sum2 = sum2 + bsarray[j] * sin(arg);
            }
            bsmod[i] = sum1 * sum1 + sum2 * sum2;
        }
        /**
         * Fix for exponential Euler spline interpolation failure.
         */
        double eps = 1.0e-7;
        if (bsmod[0] < eps) {
            bsmod[0] = 0.5 * bsmod[1];
        }
        for (int i = 1; i < nfft - 1; i++) {
            if (bsmod[i] < eps) {
                bsmod[i] = 0.5 * (bsmod[i - 1] + bsmod[i + 1]);
            }
        }
        if (bsmod[nfft - 1] < eps) {
            bsmod[nfft - 1] = 0.5 * bsmod[nfft - 2];
        }
        /**
         * Compute and apply the optimal zeta coefficient.
         */
        int jcut = 50;
        int order2 = 2 * order;
        for (int i = 0; i < nfft; i++) {
            int k = i;
            double zeta;
            if (i > nfft / 2) {
                k = k - nfft;
            }
            if (k == 0) {
                zeta = 1.0;
            } else {
                double sum1 = 1.0;
                double sum2 = 1.0;
                factor = PI * k / nfft;
                for (int j = 0; j < jcut; j++) {
                    double arg = factor / (factor + PI * (j + 1));
                    sum1 = sum1 + pow(arg, order);
                    sum2 = sum2 + pow(arg, order2);
                }
                for (int j = 0; j < jcut; j++) {
                    double arg = factor / (factor - PI * (j + 1));
                    sum1 = sum1 + pow(arg, order);
                    sum2 = sum2 + pow(arg, order2);
                }
                zeta = sum2 / sum1;
            }
            bsmod[i] = bsmod[i] * zeta * zeta;
        }
    }
    private static double toSeconds = 0.000000001;
    private final double a[][];
    private final double tfm[][] = new double[10][10];
    private final double tmm[][] = new double[10][10];
    /**
     * First lookup index to pack a 2D tensor into a 1D array.
     */
    private static final int qi1[] = {0, 1, 2, 0, 0, 1};
    /**
     * Second lookup index to pack a 2D tensor into a 1D array.
     */
    private static final int qi2[] = {0, 1, 2, 1, 2, 2};
    private static final int tensorCount = TensorRecursion.tensorCount(3);
}
