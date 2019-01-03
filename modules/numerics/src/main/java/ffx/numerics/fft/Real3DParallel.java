/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
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
 * Compute the 3D FFT of real, double precision input of arbitrary dimensions in
 * parallel.
 * <p>
 *
 * @author Michal J. Schnieders
 *
 * @since 1.0
 *
 * @see Real
 */
public class Real3DParallel {

    private static final Logger logger = Logger.getLogger(Real3DParallel.class.getName());
    private final int nX, nY, nZ, nZ2, nX1;
    private final int n, nextX, nextY, nextZ;
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final ParallelIFFT parallelIFFT;
    private final ParallelFFT parallelFFT;
    private final ParallelConvolution parallelConvolution;
    private final double recip[];
    private final IntegerSchedule schedule;

    /**
     * Initialize the FFT for real input.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     * @since 1.0
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
     */
    public Real3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam) {
        this.nX = nX / 2;
        this.nY = nY;
        this.nZ = nZ;
        this.parallelTeam = parallelTeam;
        n = nX;
        nX1 = this.nX + 1;
        nZ2 = this.nZ * 2;
        nextX = 2;
        nextY = n + 2;
        nextZ = nextY * nY;
        recip = new double[nX1 * nY * nZ];
        threadCount = parallelTeam.getThreadCount();
        parallelFFT = new ParallelFFT();
        parallelIFFT = new ParallelIFFT();
        parallelConvolution = new ParallelConvolution();
        schedule = IntegerSchedule.fixed();
    }

    /**
     * Initialize the FFT for real input.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     * @param parallelTeam The ParallelTeam that will execute the transforms.
     * @param integerSchedule The IntegerSchedule to use.
     * @since 1.0
     */
    public Real3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam,
            IntegerSchedule integerSchedule) {
        this.nX = nX / 2;
        this.nY = nY;
        this.nZ = nZ;
        this.parallelTeam = parallelTeam;
        n = nX;
        nX1 = this.nX + 1;
        nZ2 = this.nZ * 2;
        nextX = 2;
        nextY = n + 2;
        nextZ = nextY * nY;
        recip = new double[nX1 * nY * nZ];
        threadCount = parallelTeam.getThreadCount();
        if (integerSchedule != null) {
            schedule = integerSchedule;
        } else {
            schedule = IntegerSchedule.fixed();
        }
        parallelFFT = new ParallelFFT();
        parallelIFFT = new ParallelIFFT();
        parallelConvolution = new ParallelConvolution();
    }

    /**
     * Compute the 3D FFT.
     *
     * @param input The input array must be of size (nX + 2) * nY * nZ.
     * @since 1.0
     */
    public void fft(final double input[]) {
        parallelFFT.input = input;
        try {
            parallelTeam.execute(parallelFFT);
        } catch (Exception e) {
            String message = "Fatal exception evaluating real 3D FFT in parallel.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    /**
     * Compute the inverese 3D FFT.
     *
     * @param input The input array must be of size (nX + 2) * nY * nZ.
     * @since 1.0
     */
    public void ifft(final double input[]) {
        parallelIFFT.input = input;
        try {
            parallelTeam.execute(parallelIFFT);
        } catch (Exception e) {
            String message = "Fatal exception evaluating real 3D inverse FFT in parallel.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    /**
     * Compute a convolution in parallel.
     *
     * @param input The input array must be of size (nX + 2) * nY * nZ.
     * @since 1.0
     */
    public void convolution(final double input[]) {
        parallelConvolution.input = input;
        try {
            parallelTeam.execute(parallelConvolution);
        } catch (Exception e) {
            String message = "Fatal exception evaluating a 3D convolution in parallel.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    /**
     * <p>
     * Setter for the field <code>recip</code>.</p>
     *
     * @param recip The recip array must be of size [(nX/2 + 1) * nY * nZ].
     */
    public void setRecip(double recip[]) {
        /**
         * Reorder the reciprocal space data into the order it is needed by the
         * convolution routine.
         */
        for (int index = 0, offset = 0, y = 0; y < nY; y++) {
            for (int x = 0; x < nX1; x++, offset += 1) {
                for (int i = 0, z = offset; i < nZ; i++, z += nX1 * nY) {
                    this.recip[index++] = recip[z];

                }
            }
        }
    }

    /**
     * Implement the 3D parallel FFT.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class ParallelFFT extends ParallelRegion {

        public double input[];
        private final int nZm1;
        private final FFTXYLoop fftXYLoop[];
        private final FFTZLoop fftZLoop[];

        public ParallelFFT() {
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
                // There are nX + 1 frequencies.
                execute(0, nX, fftZLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Real fftX;
            private final Complex fftY;

            private FFTXYLoop() {
                fftY = new Complex(nY);
                fftX = new Real(n);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int stride = nextY, z = lb; z <= ub; z++) {
                    for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftX.fft(input, offset);
                    }
                    for (int offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTZLoop extends IntegerForLoop {

            public double input[];
            private final double work[];
            private final Complex fft;

            private FFTZLoop() {
                work = new double[nZ2];
                fft = new Complex(nZ);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int x = lb; x <= ub; x++) {
                    for (int offset = x * 2, y = 0; y < nY; y++, offset += nextY) {
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
                        }
                    }
                }
            }
        }
    }

    /**
     * Implement the 3D parallel inverse FFT.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class ParallelIFFT extends ParallelRegion {

        public double input[];
        private final int nZm1;
        private final IFFTXYLoop ifftXYLoop[];
        private final IFFTZLoop ifftZLoop[];

        public ParallelIFFT() {
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
            ifftXYLoop[threadIndex].input = input;
            ifftZLoop[threadIndex].input = input;
            try {
                // There are xDim + 1 frequencies.
                execute(0, nX, ifftZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class IFFTZLoop extends IntegerForLoop {

            public double input[];
            private final double work[];
            private final Complex fft;

            private IFFTZLoop() {
                fft = new Complex(nZ);
                work = new double[nZ2];
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int x = lb; x <= ub; x++) {
                    for (int offset = x * 2, y = 0; y < nY; y++, offset += nextY) {
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.ifft(work, 0, 2);
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Real fftX;
            private final Complex fftY;

            private IFFTXYLoop() {
                fftX = new Real(n);
                fftY = new Complex(nY);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int stride = nextY, z = lb; z <= ub; z++) {
                    for (int offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftX.ifft(input, offset);
                    }
                }
            }
        }
    }

    /**
     * Implement the 3D parallel convolution.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class ParallelConvolution extends ParallelRegion {

        public double input[];
        private final int nZm1, nYm1, nX1nZ;
        private final FFTXYLoop fftXYLoop[];
        private final FFTZ_Multiply_IFFTZLoop fftZ_Multiply_ifftZLoop[];
        private final IFFTXYLoop ifftXYLoop[];

        public ParallelConvolution() {
            nZm1 = nZ - 1;
            nYm1 = nY - 1;
            nX1nZ = nX1 * nZ;
            fftXYLoop = new FFTXYLoop[threadCount];
            fftZ_Multiply_ifftZLoop = new FFTZ_Multiply_IFFTZLoop[threadCount];
            ifftXYLoop = new IFFTXYLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftXYLoop[i] = new FFTXYLoop();
                fftZ_Multiply_ifftZLoop[i] = new FFTZ_Multiply_IFFTZLoop();
                ifftXYLoop[i] = new IFFTXYLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            fftXYLoop[threadIndex].input = input;
            fftZ_Multiply_ifftZLoop[threadIndex].input = input;
            ifftXYLoop[threadIndex].input = input;
            try {
                execute(0, nZm1, fftXYLoop[threadIndex]);
                execute(0, nYm1, fftZ_Multiply_ifftZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Real fftX;
            private final Complex fftY;

            private FFTXYLoop() {
                fftY = new Complex(nY);
                fftX = new Real(n);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int stride = nextY, z = lb; z <= ub; z++) {
                    for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftX.fft(input, offset);
                    }
                    for (int offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                        fftY.fft(input, offset, stride);
                    }
                }
            }
        }

        private class FFTZ_Multiply_IFFTZLoop extends IntegerForLoop {

            public double input[];
            private final double work[];
            private final Complex fft;

            private FFTZ_Multiply_IFFTZLoop() {
                work = new double[nZ2];
                fft = new Complex(nZ);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                int index = lb * nX1nZ;
                for (int offset = lb * nextY, y = lb; y <= ub; y++) {
                    for (int x = 0; x < nX1; x++, offset += nextX) {
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            work[i] = input[z];
                            work[i + 1] = input[z + 1];
                        }
                        fft.fft(work, 0, 2);
                        for (int i = 0; i < nZ2; i += 2) {
                            double r = recip[index++];
                            work[i] *= r;
                            work[i + 1] *= r;
                        }
                        fft.ifft(work, 0, 2);
                        for (int z = offset, i = 0; i < nZ2; i += 2, z += nextZ) {
                            input[z] = work[i];
                            input[z + 1] = work[i + 1];
                        }
                    }
                }
            }
        }

        private class IFFTXYLoop extends IntegerForLoop {

            public double input[];
            private final Real fftX;
            private final Complex fftY;

            private IFFTXYLoop() {
                fftY = new Complex(nY);
                fftX = new Real(n);
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int stride = nextY, z = lb; z <= ub; z++) {
                    for (int offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                        fftY.ifft(input, offset, stride);
                    }
                    for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                        fftX.ifft(input, offset);
                    }
                }
            }
        }
    }

    /**
     * Test the real 3D FFT.
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
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
        if (dimNotFinal % 2 != 0) {
            dimNotFinal++;
        }
        final int dim = dimNotFinal;
        System.out.println(String.format(
                "Initializing a %d cubed grid for %d CPUs.\n"
                + "The best timing out of %d repititions will be used.",
                dim, ncpu, reps));

        Real3D real3D = new Real3D(dim, dim, dim);
        ParallelTeam parallelTeam = new ParallelTeam(ncpu);
        Real3DParallel real3DParallel = new Real3DParallel(dim, dim, dim,
                parallelTeam);

        final int dimCubed = (dim + 2) * dim * dim;
        final double data[] = new double[dimCubed];
        final double work[] = new double[dimCubed];

        // Parallel Array Initialization.
        try {
            parallelTeam.execute(new ParallelRegion() {
                @Override
                public void run() {
                    try {
                        execute(0, dim - 1, new IntegerForLoop() {
                            @Override
                            public void run(int lb, int ub) {
                                Random randomNumberGenerator = new Random(1);
                                int index = dim * dim * lb;
                                for (int z = lb; z <= ub; z++) {
                                    for (int y = 0; y < dim; y++) {
                                        for (int x = 0; x < dim; x++) {
                                            double randomNumber = randomNumberGenerator.nextDouble();
                                            data[index] = randomNumber;
                                            index++;
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
        real3D.setRecip(work);
        real3DParallel.setRecip(work);
        for (int i = 0; i < reps; i++) {
            System.out.println(String.format("Iteration %d", i + 1));
            long time = System.nanoTime();
            real3D.fft(data);
            real3D.ifft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f", toSeconds
                    * time));
            if (time < seqTime) {
                seqTime = time;
            }
            time = System.nanoTime();
            real3D.convolution(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Sequential: %8.3f (Convolution)",
                    toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }
            time = System.nanoTime();
            real3DParallel.fft(data);
            real3DParallel.ifft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("Parallel:   %8.3f", toSeconds
                    * time));
            if (time < parTime) {
                parTime = time;
            }
            time = System.nanoTime();
            real3DParallel.convolution(data);
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
