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
package ffx.numerics.fft;

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

/**
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs in parallel.
 * <p>
 * The location of the input point [i, j, k] within the input array must be:<br>
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ]<br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1]<br>
 * <br>
 * where<br>
 * int nextX = 2<br>
 * int nextY = 2*nX<br>
 * int nextZ = 2*nX*nY<br>
 * <p>
 * 
 * @author Michal J. Schnieders
 * @since 1.0
 */
public class RowMajorComplex3DParallel {

    private final int nX, nY, nZ;
    private final int nX2, nY2;
    private final int nextX, nextY, nextZ;
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final ParallelFFT parallelFFT;
    private final ParallelIFFT parallelIFFT;
    private final Convolution convolution;
    private final double[] recip;
    private final IntegerSchedule schedule;
    private static final Logger logger = Logger.getLogger(RowMajorComplex3DParallel.class.getName());

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX
     *            X-dimension.
     * @param nY
     *            Y-dimension.
     * @param nZ
     *            Z-dimension.
     * @param parallelTeam
     *            A ParallelTeam instance.
     *
     * @since 1.0
     */
    public RowMajorComplex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam) {
        this(nX, nY, nZ, parallelTeam, IntegerSchedule.fixed());
    }

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX
     *            X-dimension.
     * @param nY
     *            Y-dimension.
     * @param nZ
     *            Z-dimension.
     * @param parallelTeam
     *            A ParallelTeam instance.
     * @param integerSchedule
     *            The IntegerSchedule to use.
     *
     * @since 1.0
     */
    public RowMajorComplex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam,
                                     IntegerSchedule integerSchedule) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.parallelTeam = parallelTeam;
        recip = new double[nX * nY * nZ];
        nX2 = 2 * nX;
        nY2 = 2 * nY;
        nextZ = 2;
        nextY = 2 * nZ;
        nextX = nextY * nY;

        threadCount = parallelTeam.getThreadCount();
        if (integerSchedule != null) {
            schedule = integerSchedule;
        } else {
            schedule = IntegerSchedule.fixed();
        }

        parallelFFT = new ParallelFFT();
        parallelIFFT = new ParallelIFFT();
        convolution = new Convolution();
    }

    /**
     * Compute the 3D FFT in pararallel.
     *
     * @param input
     *            The input array must be of size 2 * nX * nY * nZ.
     *
     * @since 1.0
     */
    public void fft(final double input[]) {
        parallelFFT.input = input;
        try {
            parallelTeam.execute(parallelFFT);
        } catch (Exception e) {
            String message = "Fatal exception evaluating the FFT.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute the inverse 3D FFT in parallel.
     *
     * @param input
     *            The input array must be of size 2 * nX * nY * nZ.
     *
     * @since 1.0
     */
    public void ifft(final double input[]) {
        parallelIFFT.input = input;
        try {
            parallelTeam.execute(parallelIFFT);
        } catch (Exception e) {
            String message = "Fatal exception evaluating the inverse FFT.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute the 3D FFT, perfrom a multiplication in reciprocal space, and the
     * inverese 3D FFT all in parallel.
     *
     * @param input
     *            The input array must be of size 2 * nX * nY * nZ.
     *
     * @since 1.0
     */
    public void convolution(final double input[]) {
        convolution.input = input;
        try {
            parallelTeam.execute(convolution);
        } catch (Exception e) {
            String message = "Fatal exception evaluating a convolution.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void setRecip(double recip[]) {
        int offset, y, x, z, i;

        /**
         * Reorder the reciprocal space data into the order it is needed
         * by the convolution routine.
         */
        int index = 0;
        for (offset = 0, y = 0; y < nY; y++) {
            for (z = 0; z < nZ; z++, offset += 1) {
                for (i = 0, x = offset; i < nX2; i += 2, x += nZ * nY) {
                    this.recip[index++] = recip[x];
                }
            }
        }
    }

    private class ParallelFFT extends ParallelRegion {

        public double input[];
        private final int nXm1;
        private final int nZm1;
        private final FFTZYLoop fftZYLoop[];
        private final FFTXLoop fftXLoop[];

        public ParallelFFT() {
            nXm1 = nX - 1;
            nZm1 = nZ - 1;
            fftZYLoop = new FFTZYLoop[threadCount];
            fftXLoop = new FFTXLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftZYLoop[i] = new FFTZYLoop();
                fftXLoop[i] = new FFTXLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            fftZYLoop[threadIndex].input = input;
            fftXLoop[threadIndex].input = input;
            try {
                execute(0, nXm1, fftZYLoop[threadIndex]);
                execute(0, nZm1, fftXLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTZYLoop extends IntegerForLoop {

            public double input[];
            private int x, y, z, offset, stride;
            private final Complex fftZ = new Complex(nZ);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (x = lb; x <= ub; x++) {
                    for (y = 0, offset = x * nextX, stride = nextZ; y < nY; y++, offset += nextY) {
                        fftZ.fft(input, offset, stride);
                    }
                    for (z = 0, offset = x * nextX, stride = nextY; z < nZ; z++, offset += nextZ) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTXLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nX2];
            private int i, x, y, z, offset;
            private final Complex fft = new Complex(nX);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (z = lb, offset = lb * nY2; z <= ub; z++) {
                    for (y = 0; y < nY; y++, offset += 2) {
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            work[i] = input[x];
                            work[i + 1] = input[x + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            input[x] = work[i];
                            input[x + 1] = work[i + 1];
                        }
                    }
                }
            }
        }
    }

    private class ParallelIFFT extends ParallelRegion {

        public double input[];
        private final int nXm1;
        private final int nZm1;
        private final IFFTZYLoop ifftZYLoop[];
        private final IFFTXLoop ifftXLoop[];

        public ParallelIFFT() {
            nXm1 = nX - 1;
            nZm1 = nZ - 1;
            ifftZYLoop = new IFFTZYLoop[threadCount];
            ifftXLoop = new IFFTXLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                ifftZYLoop[i] = new IFFTZYLoop();
                ifftXLoop[i] = new IFFTXLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            ifftXLoop[threadIndex].input = input;
            ifftZYLoop[threadIndex].input = input;
            try {
                execute(0, nZm1, ifftXLoop[threadIndex]);
                execute(0, nXm1, ifftZYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class IFFTXLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nX2];
            private int i, x, y, z, offset;
            private final Complex fft = new Complex(nX);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (offset = lb * nY2, z = lb; z <= ub; z++) {
                    for (y = 0; y < nY; y++, offset += 2) {
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            work[i] = input[x];
                            work[i + 1] = input[x + 1];
                        }
                        fft.ifft(work, 0, 2);
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            input[x] = work[i];
                            input[x + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTZYLoop extends IntegerForLoop {

            public double input[];
            private int x, y, z, offset, stride;
            private final Complex fftZ = new Complex(nZ);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (x = lb; x <= ub; x++) {
                    for (offset = x * nextX, stride = nextY, z = 0; z < nZ; z++, offset += nextZ) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (offset = x * nextX, stride = nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftZ.ifft(input, offset, stride);
                    }
                }
            }
        }
    }

    private class Convolution extends ParallelRegion {

        public double input[];
        private final int nYm1;
        private final int nXm1;
        private final FFTZYLoop fftZYLoop[];
        private final FFTX_Multiply_IFFTXLoop fftX_Multiply_IFFTXLoop[];
        private final IFFTZYLoop ifftZYLoop[];

        public Convolution() {
            nYm1 = nY - 1;
            nXm1 = nX - 1;
            fftZYLoop = new FFTZYLoop[threadCount];
            fftX_Multiply_IFFTXLoop = new FFTX_Multiply_IFFTXLoop[threadCount];
            ifftZYLoop = new IFFTZYLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftZYLoop[i] = new FFTZYLoop();
                fftX_Multiply_IFFTXLoop[i] = new FFTX_Multiply_IFFTXLoop();
                ifftZYLoop[i] = new IFFTZYLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            fftZYLoop[threadIndex].input = input;
            fftX_Multiply_IFFTXLoop[threadIndex].input = input;
            ifftZYLoop[threadIndex].input = input;
            try {
                execute(0, nXm1, fftZYLoop[threadIndex]);
                execute(0, nYm1, fftX_Multiply_IFFTXLoop[threadIndex]);
                execute(0, nXm1, ifftZYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTZYLoop extends IntegerForLoop {

            public double input[];
            private final Complex fftZ = new Complex(nZ);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int x, y, z, offset, stride;
                for (x = lb; x <= ub; x++) {
                    for (offset = x * nextX, stride = nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftZ.fft(input, offset, stride);
                    }
                    for (offset = x * nextX, stride = nextY, z = 0; z < nZ; z++, offset += nextZ) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTX_Multiply_IFFTXLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nX2];
            private final Complex fft = new Complex(nZ);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int i, x, y, z, offset;
                int index = nX * nZ * lb;
                for (offset = lb * nextY, y = lb; y <= ub; y++) {
                    for (z = 0; z < nZ; z++, offset += 2) {
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            work[i] = input[x];
                            work[i + 1] = input[x + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (i = 0; i < nX2; i += 2) {
                            double r = recip[index++];
                            work[i] *= r;
                            work[i + 1] *= r;
                        }
                        fft.ifft(work, 0, 2);
                        for (i = 0, x = offset; i < nX2; i += 2, x += nextX) {
                            input[x] = work[i];
                            input[x + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTZYLoop extends IntegerForLoop {

            public double input[];
            private final Complex fftY = new Complex(nY);
            private final Complex fftZ = new Complex(nZ);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int x, y, z, offset, stride;
                for (x = lb; x <= ub; x++) {
                    for (offset = x * nextX, stride = nextY, z = 0; z < nZ; z++, offset += nextZ) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (offset = x * nextX, stride = nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftZ.ifft(input, offset, stride);
                    }
                }
            }
        }
    }

    /**
     * Test the Complex3DParallel FFT.
     *
     * @param args
     * @throws Exception
     *
     * @since 1.0
     */
    public static void main(String[] args) throws Exception {
        int dimNotFinal = 128;
        int ncpu = ParallelTeam.getDefaultThreadCount();
        int reps = 5;
        if (args != null) {
            try {
                dimNotFinal = Integer.parseInt(args[0]);
                if (dimNotFinal < 1) {
                    dimNotFinal = 100;
                }
                ncpu = Integer.parseInt(args[1]);
                if (ncpu < 1) {
                    ncpu = ParallelTeam.getDefaultThreadCount();
                }
                reps = Integer.parseInt(args[2]);
                if (reps < 1) {
                    reps = 5;
                }
            } catch (Exception e) {
            }
        }
        final int dim = dimNotFinal;
        System.out.println(String.format(
                "Initializing a %d cubed grid for %d CPUs.\n"
                + "The best timing out of %d repititions will be used.",
                dim, ncpu, reps));
        // One dimension of the serial array divided by the number of threads.
        Complex3D columnComplex3D = new Complex3D(dim, dim, dim);
        RowMajorComplex3D rowComplex3D = new RowMajorComplex3D(dim, dim, dim);

        ParallelTeam parallelTeam = new ParallelTeam(ncpu);
        Complex3DParallel columnComplex3DParallel = new Complex3DParallel(dim, dim, dim, parallelTeam);
        RowMajorComplex3DParallel rowComplex3DParallel = new RowMajorComplex3DParallel(
                dim, dim, dim, parallelTeam);

        final int dimCubed = dim * dim * dim;
        final double rowOrig[] = new double[dimCubed * 2];
        final double columnOrig[] = new double[dimCubed * 2];
        final double rowMajor[] = new double[dimCubed * 2];
        final double columnMajor[] = new double[dimCubed * 2];
        double rowRecip[] = new double[dimCubed];
        double columnRecip[] = new double[dimCubed];

        Random randomNumberGenerator = new Random();
        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    int columnMajorIndex = Complex3D.iComplex3D(x, y, z, dim, dim);
                    int rowMajorIndex = RowMajorComplex3D.iComplex3D(x, y, z, dim, dim);
                    double randomNumber = randomNumberGenerator.nextDouble();
                    rowOrig[rowMajorIndex] = randomNumber;
                    rowRecip[rowMajorIndex / 2] = randomNumber;
                    columnOrig[columnMajorIndex] = randomNumber;
                    columnRecip[columnMajorIndex / 2] = randomNumber;
                }
            }
        }

        rowComplex3D.setRecip(rowRecip);
        rowComplex3DParallel.setRecip(rowRecip);
        columnComplex3D.setRecip(columnRecip);
        columnComplex3DParallel.setRecip(columnRecip);
        rowRecip = null;
        columnRecip = null;

        System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
        System.arraycopy(columnOrig, 0, columnMajor, 0, dimCubed * 2);
        rowComplex3D.convolution(rowMajor);
        columnComplex3D.convolution(columnMajor);

        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    int columnMajorIndex = Complex3D.iComplex3D(x, y, z, dim, dim);
                    int rowMajorIndex = RowMajorComplex3D.iComplex3D(x, y, z, dim, dim);
                    double error = rowMajor[rowMajorIndex] - columnMajor[columnMajorIndex];
                    if (Math.abs(error) > 1.0e-6) {
                        logger.info(String.format("Error (%d,%d,%d): %12.6f", x, y, z, error));
                    }
                }
            }
        }
        logger.info("Single threaded column-major and row-major agree.");

        System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
        rowComplex3D.convolution(rowMajor);
        // note - using the columnMajor array for temporary storage
        System.arraycopy(rowOrig, 0, columnMajor, 0, dimCubed * 2);
        rowComplex3DParallel.convolution(columnMajor);

        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    int rowMajorIndex = RowMajorComplex3D.iComplex3D(x, y, z, dim, dim);
                    double error = rowMajor[rowMajorIndex] - columnMajor[rowMajorIndex];
                    if (Math.abs(error) > 1.0e-6) {
                        logger.info(String.format("Error (%d,%d,%d): %12.6f", x, y, z, error));
                    }
                }
            }
        }
        logger.info("Single threaded and SMP row-major agree.");

        System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
        System.arraycopy(columnOrig, 0, columnMajor, 0, dimCubed * 2);
        rowComplex3DParallel.convolution(rowMajor);
        columnComplex3DParallel.convolution(columnMajor);

        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    int columnMajorIndex = Complex3D.iComplex3D(x, y, z, dim, dim);
                    int rowMajorIndex = RowMajorComplex3D.iComplex3D(x, y, z, dim, dim);
                    double error = rowMajor[rowMajorIndex] - columnMajor[columnMajorIndex];
                    if (Math.abs(error) > 1.0e-6) {
                        logger.info(String.format("Error (%d,%d,%d): %12.6f", x, y, z, error));
                    }
                }
            }
        }
        logger.info("SMP column-major and row-major agree.");

        double toSeconds = 0.000000001;
        long parTime = Long.MAX_VALUE;
        long seqTime = Long.MAX_VALUE;
        for (int i = 0; i < reps; i++) {

            System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
            
            System.out.println(String.format("Iteration %d", i + 1));
            long time = System.nanoTime();
            rowComplex3D.fft(rowMajor);
            rowComplex3D.ifft(rowMajor);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f", toSeconds
                                                                  * time));
            if (time < seqTime) {
                seqTime = time;
            }

            System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
            time = System.nanoTime();
            rowComplex3D.convolution(rowMajor);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f (Convolution)",
                                             toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }

            System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
            time = System.nanoTime();
            rowComplex3DParallel.fft(rowMajor);
            rowComplex3DParallel.ifft(rowMajor);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Parallel:   %8.3f", toSeconds
                                                                  * time));
            if (time < parTime) {
                parTime = time;
            }

            System.arraycopy(rowOrig, 0, rowMajor, 0, dimCubed * 2);
            time = System.nanoTime();
            rowComplex3DParallel.convolution(rowMajor);
            time = (System.nanoTime() - time);
            System.out.println(String.format(
                    "Parallel:   %8.3f (Convolution)\n", toSeconds * time));
            if (time < parTime) {
                parTime = time;
            }
        }
        System.out.println(String.format("Best Sequential Time:  %8.3f",
                                         toSeconds * seqTime));
        System.out.println(String.format("Best Parallel Time:    %8.3f",
                                         toSeconds * parTime));
        System.out.println(String.format("Speedup: %15.5f", (double) seqTime
                                                            / parTime));
    }
}
