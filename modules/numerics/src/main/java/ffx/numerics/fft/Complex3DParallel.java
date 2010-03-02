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
public class Complex3DParallel {

    private final int nX, nY, nZ;
    private final int nY2, nZ2;
    private final int nextX, nextY, nextZ;
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final ParallelFFT parallelFFT;
    private final ParallelIFFT parallelIFFT;
    private final Convolution convolution;
    private final double[] recip;
    private final IntegerSchedule schedule;
    private static final Logger logger = Logger.getLogger(Complex3DParallel.class.getName());

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
    public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.parallelTeam = parallelTeam;
        recip = new double[nX * nY * nZ];
        nY2 = 2 * this.nY;
        nZ2 = 2 * this.nZ;
        nextX = 2;
        nextY = 2 * this.nX;
        nextZ = nextY * this.nY;
        threadCount = parallelTeam.getThreadCount();
        schedule = IntegerSchedule.fixed();
        parallelFFT = new ParallelFFT();
        parallelIFFT = new ParallelIFFT();
        convolution = new Convolution();
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
    public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam,
                             IntegerSchedule integerSchedule) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.parallelTeam = parallelTeam;
        recip = new double[nX * nY * nZ];
        nY2 = 2 * this.nY;
        nZ2 = 2 * this.nZ;
        nextX = 2;
        nextY = 2 * this.nX;
        nextZ = nextY * this.nY;
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
            System.exit(-1);
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
            System.exit(-1);
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
            for (x = 0; x < nX; x++, offset += 1) {
                for (i = 0, z = offset; i < nZ2; i += 2, z += nX * nY) {
                    this.recip[index++] = recip[z];
                }
            }
        }
    }

    private class ParallelFFT extends ParallelRegion {

        public double input[];
        private final int nXm1;
        private final int nZm1;
        private final FFTXYLoop fftXYLoop[];
        private final FFTZLoop fftZLoop[];

        public ParallelFFT() {
            nXm1 = nX - 1;
            nZm1 = nZ - 1;
            fftXYLoop = new FFTXYLoop[threadCount];
            fftZLoop = new FFTZLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftXYLoop[i] = new FFTXYLoop();
                fftZLoop[i] = new FFTZLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            fftXYLoop[threadIndex].input = input;
            fftZLoop[threadIndex].input = input;
            try {
                execute(0, nZm1, fftXYLoop[threadIndex]);
                execute(0, nXm1, fftZLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTXYLoop extends IntegerForLoop {

            public double input[];
            private int x, y, z, offset, stride;
            private final Complex fftX = new Complex(nX);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (z = lb; z <= ub; z++) {
                    for (y = 0, offset = z * nextZ, stride = nextX; y < nY; y++, offset += nextY) {
                        fftX.fft(input, offset, stride);
                    }
                    for (x = 0, offset = z * nextZ, stride = nextY; x < nX; x++, offset += nextX) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTZLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nZ2];
            private int i, x, y, z, offset;
            private final Complex fft = new Complex(nZ);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (x = lb, offset = lb * nY2; x <= ub; x++) {
                    for (y = 0; y < nY; y++, offset += 2) {
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
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
        private final IFFTXYLoop ifftXYLoop[];
        private final IFFTZLoop ifftZLoop[];

        public ParallelIFFT() {
            nXm1 = nX - 1;
            nZm1 = nZ - 1;
            ifftXYLoop = new IFFTXYLoop[threadCount];
            ifftZLoop = new IFFTZLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                ifftXYLoop[i] = new IFFTXYLoop();
                ifftZLoop[i] = new IFFTZLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            ifftZLoop[threadIndex].input = input;
            ifftXYLoop[threadIndex].input = input;
            try {
                execute(0, nXm1, ifftZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
                e.printStackTrace();
            }
        }

        private class IFFTZLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nZ2];
            private int i, x, y, z, offset;
            private final Complex fft = new Complex(nZ);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (offset = lb * nY2, x = lb; x <= ub; x++) {
                    for (y = 0; y < nY; y++, offset += 2) {
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.ifft(work, 0, 2);
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTXYLoop extends IntegerForLoop {

            public double input[];
            private int x, y, z, offset, stride;
            private final Complex fftX = new Complex(nX);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (z = lb; z <= ub; z++) {
                    for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                        fftX.ifft(input, offset, stride);
                    }
                }
            }
        }
    }

    private class Convolution extends ParallelRegion {

        public double input[];
        private final int nYm1;
        private final int nZm1;
        private final FFTXYLoop fftXYLoop[];
        private final FFTZ_Multiply_IFFTZLoop fftZ_Multiply_IFFTZLoop[];
        private final IFFTXYLoop ifftXYLoop[];

        public Convolution() {
            nYm1 = nY - 1;
            nZm1 = nZ - 1;
            fftXYLoop = new FFTXYLoop[threadCount];
            fftZ_Multiply_IFFTZLoop = new FFTZ_Multiply_IFFTZLoop[threadCount];
            ifftXYLoop = new IFFTXYLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftXYLoop[i] = new FFTXYLoop();
                fftZ_Multiply_IFFTZLoop[i] = new FFTZ_Multiply_IFFTZLoop();
                ifftXYLoop[i] = new IFFTXYLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            fftXYLoop[threadIndex].input = input;
            fftZ_Multiply_IFFTZLoop[threadIndex].input = input;
            ifftXYLoop[threadIndex].input = input;
            try {
                execute(0, nZm1, fftXYLoop[threadIndex]);
                execute(0, nYm1, fftZ_Multiply_IFFTZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Complex fftX = new Complex(nX);
            private final Complex fftY = new Complex(nY);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int x, y, z, offset, stride;
                for (z = lb; z <= ub; z++) {
                    for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                        fftX.fft(input, offset, stride);
                    }
                    for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTZ_Multiply_IFFTZLoop extends IntegerForLoop {

            public double input[];
            private double work[] = new double[nZ2];
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
                    for (x = 0; x < nX; x++, offset += 2) {
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (i = 0; i < nZ2; i += 2) {
                            double r = recip[index++];
                            work[i] *= r;
                            work[i + 1] *= r;
                        }
                        fft.ifft(work, 0, 2);
                        for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Complex fftY = new Complex(nY);
            private final Complex fftX = new Complex(nX);

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int x, y, z, offset, stride;
                for (z = lb; z <= ub; z++) {
                    for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                        fftX.ifft(input, offset, stride);
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
        Complex3D complexDoubleFFT3D = new Complex3D(dim, dim, dim);
        ParallelTeam parallelTeam = new ParallelTeam(ncpu);
        Complex3DParallel parallelComplexDoubleFFT3D = new Complex3DParallel(
                dim, dim, dim, parallelTeam);
        final int dimCubed = dim * dim * dim;
        final double data[] = new double[dimCubed * 2];
        final double work[] = new double[dimCubed * 2];
        // Parallel Array Initialization.
        try {
            parallelTeam.execute(new ParallelRegion() {

                @Override
                public void run() {
                    try {
                        execute(0, dim - 1, new IntegerForLoop() {

                            @Override
                            public void run(final int lb, final int ub) {
                                Random randomNumberGenerator = new Random(1);
                                int index = dim * dim * lb * 2;
                                for (int i = lb; i <= ub; i++) {
                                    for (int j = 0; j < dim; j++) {
                                        for (int k = 0; k < dim; k++) {
                                            double randomNumber = randomNumberGenerator.nextDouble();
                                            data[index] = randomNumber;
                                            index += 2;
                                        }
                                    }
                                }
                            }
                        });
                    } catch (Exception e) {
                        System.out.println(e.toString());
                        System.exit(-1);
                    }
                }
            });
        } catch (Exception e) {
            System.out.println(e.toString());
            System.exit(-1);
        }
        double toSeconds = 0.000000001;
        long parTime = Long.MAX_VALUE;
        long seqTime = Long.MAX_VALUE;
        complexDoubleFFT3D.setRecip(work);
        parallelComplexDoubleFFT3D.setRecip(work);
        for (int i = 0; i < reps; i++) {
            System.out.println(String.format("Iteration %d", i + 1));
            long time = System.nanoTime();
            complexDoubleFFT3D.fft(data);
            complexDoubleFFT3D.ifft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f", toSeconds
                                                                  * time));
            if (time < seqTime) {
                seqTime = time;
            }
            time = System.nanoTime();
            complexDoubleFFT3D.convolution(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f (Convolution)",
                                             toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }
            time = System.nanoTime();
            parallelComplexDoubleFFT3D.fft(data);
            parallelComplexDoubleFFT3D.ifft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Parallel:   %8.3f", toSeconds
                                                                  * time));
            if (time < parTime) {
                parTime = time;
            }
            time = System.nanoTime();
            parallelComplexDoubleFFT3D.convolution(data);
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
