/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs in parallel.
 * <p>
 * The location of the input point [i, j, k] within the input array must be:
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ]
 * <br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1]
 * <br>
 * where
 * <br>
 * int nextX = 2
 * <br>
 * int nextY = 2*nX
 * <br>
 * int nextZ = 2*nX*nY
 * <br>
 *
 * @author Michal J. Schnieders
 *
 * @since 1.0
 *
 * @see Complex
 */
public class Complex3DParallel {

    private static final Logger logger = Logger.getLogger(Complex3DParallel.class.getName());
    private final int nX, nY, nZ;
    private final int nY2, nZ2;
    private final int strideX, strideY, strideZ;
    private final double[] recip;
    private final long convolutionTime[];
    private final int threadCount;
    private final ParallelTeam parallelTeam;
    private final Complex fftX[];
    private final Complex fftY[];
    private final Complex fftZ[];
    private final IntegerSchedule schedule;

    public double input[];
    public final int nXm1, nYm1, nZm1;
    public final FFTRegion fftRegion;
    public final IFFTRegion ifftRegion;
    public final ConvolutionRegion convRegion;

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     * @param parallelTeam A ParallelTeam instance.
     * @since 1.0
     */
    public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam) {
        this(nX, nY, nZ, parallelTeam, null);
    }

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     * @param parallelTeam A ParallelTeam instance.
     * @param integerSchedule The IntegerSchedule to use.
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
        strideX = 2;
        strideY = 2 * this.nX;
        strideZ = strideY * this.nY;
        nXm1 = this.nX - 1;
        nYm1 = this.nY - 1;
        nZm1 = this.nZ - 1;
        threadCount = parallelTeam.getThreadCount();
        if (integerSchedule != null) {
            schedule = integerSchedule;
        } else {
            schedule = IntegerSchedule.fixed();
        }
        fftX = new Complex[threadCount];
        fftY = new Complex[threadCount];
        fftZ = new Complex[threadCount];
        for (int i = 0; i < threadCount; i++) {
            fftX[i] = new Complex(nX);
            fftY[i] = new Complex(nY);
            fftZ[i] = new Complex(nZ);
        }
        fftRegion = new FFTRegion();
        ifftRegion = new IFFTRegion();
        convRegion = new ConvolutionRegion();
        convolutionTime = new long[threadCount];
    }

    public void initTiming() {
        for (int i = 0; i < threadCount; i++) {
            convolutionTime[i] = 0;
        }
    }

    public long[] getTimings() {
        return convolutionTime;
    }

    /**
     * Compute the 3D FFT in pararallel.
     *
     * @param input The input array must be of size 2 * nX * nY * nZ.
     * @since 1.0
     */
    public void fft(final double input[]) {
        this.input = input;
        try {
            parallelTeam.execute(fftRegion);
        } catch (Exception e) {
            String message = " Fatal exception evaluating the FFT.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute the inverse 3D FFT in parallel.
     *
     * @param input The input array must be of size 2 * nX * nY * nZ.
     * @since 1.0
     */
    public void ifft(final double input[]) {
        this.input = input;
        try {
            parallelTeam.execute(ifftRegion);
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
     * @param input The input array must be of size 2 * nX * nY * nZ.
     * @since 1.0
     */
    public void convolution(final double input[]) {
        this.input = input;
        try {
            parallelTeam.execute(convRegion);
        } catch (Exception e) {
            String message = "Fatal exception evaluating a convolution.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * <p>
     * Setter for the field <code>recip</code>.</p>
     *
     * @param recip an array of double.
     */
    public void setRecip(double recip[]) {
        int offset, y, x, z, i;

        /**
         * Reorder the reciprocal space data into the order it is needed by the
         * convolution routine.
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

    /**
     * An external ParallelRegion can be used as follows:
     *
     * <code>
     * start() {
     *  fftRegion.input = input;
     * }
     * run(){
     *  execute(0, nZm1, fftRegion.fftXYLoop[threadID]);
     *  execute(0, nXm1, fftRegion.fftZLoop[threadID]);
     * }
     * </code>
     */
    public class FFTRegion extends ParallelRegion {

        public final FFTXYLoop fftXYLoop[];
        public final FFTZLoop fftZLoop[];

        public FFTRegion() {
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
            try {
                execute(0, nZm1, fftXYLoop[threadIndex]);
                execute(0, nXm1, fftZLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }
    }

    /**
     * An external ParallelRegion can be used as follows:
     *
     * <code>
     * start() {
     *  ifftRegion.input = input;
     * }
     * run(){
     *  execute(0, nXm1, ifftRegion.ifftZLoop[threadID]);
     *  execute(0, nZm1, ifftRegion.ifftXYLoop[threadID]);
     * }
     * </code>
     */
    public class IFFTRegion extends ParallelRegion {

        public final IFFTXYLoop ifftXYLoop[];
        public final IFFTZLoop ifftZLoop[];

        public IFFTRegion() {
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
            try {
                execute(0, nXm1, ifftZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }
    }

    /**
     * An external ParallelRegion can be used as follows:
     *
     * <code>
     * start() {
     *  convRegion.input = input;
     * }
     * run(){
     *  execute(0, nZm1, convRegion.fftXYLoop[threadID]);
     *  execute(0, nYm1, convRegion.fftZIZLoop[threadID]);
     *  execute(0, nZm1, convRegion.ifftXYLoop[threadID]);
     * }
     * </code>
     */
    public class ConvolutionRegion extends ParallelRegion {

        private final FFTXYLoop fftXYLoop[];
        private final FFTZIZLoop fftZIZLoop[];
        private final IFFTXYLoop ifftXYLoop[];

        public ConvolutionRegion() {
            fftXYLoop = new FFTXYLoop[threadCount];
            fftZIZLoop = new FFTZIZLoop[threadCount];
            ifftXYLoop = new IFFTXYLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                fftXYLoop[i] = new FFTXYLoop();
                fftZIZLoop[i] = new FFTZIZLoop();
                ifftXYLoop[i] = new IFFTXYLoop();
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            convolutionTime[threadIndex] -= System.nanoTime();
            try {
                execute(0, nZm1, fftXYLoop[threadIndex]);
                execute(0, nYm1, fftZIZLoop[threadIndex]);
                execute(0, nZm1, ifftXYLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
            convolutionTime[threadIndex] += System.nanoTime();
        }
    }

    public class FFTXYLoop extends IntegerForLoop {

        private Complex localFFTX;
        private Complex localFFTY;

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            localFFTX = fftX[getThreadIndex()];
            localFFTY = fftY[getThreadIndex()];
        }

        @Override
        public void run(final int lb, final int ub) {
            for (int z = lb; z <= ub; z++) {
                for (int offset = z * strideZ, y = 0; y < nY; y++, offset += strideY) {
                    localFFTX.fft(input, offset, strideX);
                }
                for (int offset = z * strideZ, x = 0; x < nX; x++, offset += strideX) {
                    localFFTY.fft(input, offset, strideY);
                }
            }
        }
    }

    public class FFTZLoop extends IntegerForLoop {

        private final double work[];
        private Complex localFFTZ;

        private FFTZLoop() {
            work = new double[nZ2];
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            localFFTZ = fftZ[getThreadIndex()];
        }

        @Override
        public void run(final int lb, final int ub) {
            for (int x = lb, offset = lb * nY2; x <= ub; x++) {
                for (int y = 0; y < nY; y++, offset += 2) {
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        work[i] = input[z];
                        work[i + 1] = input[z + 1];
                    }
                    localFFTZ.fft(work, 0, 2);
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        input[z] = work[i];
                        input[z + 1] = work[i + 1];
                    }
                }
            }
        }
    }

    public class IFFTXYLoop extends IntegerForLoop {

        private Complex localFFTY;
        private Complex localFFTX;

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            localFFTX = fftX[getThreadIndex()];
            localFFTY = fftY[getThreadIndex()];
        }

        @Override
        public void run(final int lb, final int ub) {
            for (int z = lb; z <= ub; z++) {
                for (int offset = z * strideZ, x = 0; x < nX; x++, offset += strideX) {
                    localFFTY.ifft(input, offset, strideY);
                }
                for (int offset = z * strideZ, y = 0; y < nY; y++, offset += strideY) {
                    localFFTX.ifft(input, offset, strideX);
                }
            }
        }
    }

    public class IFFTZLoop extends IntegerForLoop {

        private final double work[];
        private Complex localFFTZ;

        private IFFTZLoop() {
            work = new double[nZ2];
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            localFFTZ = fftZ[getThreadIndex()];
        }

        @Override
        public void run(final int lb, final int ub) {
            for (int offset = lb * nY2, x = lb; x <= ub; x++) {
                for (int y = 0; y < nY; y++, offset += 2) {
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        work[i] = input[z];
                        work[i + 1] = input[z + 1];
                    }
                    localFFTZ.ifft(work, 0, 2);
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        input[z] = work[i];
                        input[z + 1] = work[i + 1];
                    }
                }
            }
        }
    }

    public class FFTZIZLoop extends IntegerForLoop {

        private final double work[];
        private Complex localFFTZ;

        private FFTZIZLoop() {
            work = new double[nZ2];
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            localFFTZ = fftZ[getThreadIndex()];
        }

        @Override
        public void run(final int lb, final int ub) {
            int index = nX * nZ * lb;
            for (int offset = lb * strideY, y = lb; y <= ub; y++) {
                for (int x = 0; x < nX; x++, offset += 2) {
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        work[i] = input[z];
                        work[i + 1] = input[z + 1];
                    }
                    localFFTZ.fft(work, 0, 2);
                    for (int i = 0; i < nZ2; i += 2) {
                        double r = recip[index++];
                        work[i] *= r;
                        work[i + 1] *= r;
                    }
                    localFFTZ.ifft(work, 0, 2);
                    for (int i = 0, z = offset; i < nZ2; i += 2, z += strideZ) {
                        input[z] = work[i];
                        input[z + 1] = work[i + 1];
                    }
                }
            }
        }
    }

    /**
     * Test the Complex3DParallel FFT.
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
        parallelTeam.shutdown();
    }
}
