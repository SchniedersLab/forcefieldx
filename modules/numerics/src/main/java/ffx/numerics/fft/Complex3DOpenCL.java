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

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;

/**
 * Compute a 3D FFT or Convolution using a Java Native Interface
 * bridge to OpenCL.
 *
 * @author Michal J. Schnieders
 * @since 1.0
 */
public class Complex3DOpenCL {

    private static final Logger logger = Logger.getLogger(
            ffx.numerics.fft.Complex3DOpenCL.class.getName());

    private native int init(int x, int y, int z,
            FloatBuffer real, FloatBuffer imag, FloatBuffer recip);

    public native int free();

    public native int convolution();

    /**
     * Load the native library.
     */
    static {
        try {
            System.loadLibrary("convolution");
        } catch (Error e) {
            logger.log(Level.INFO, System.getProperty("java.library.path"));
            String message = "Fatal error loading the OpenCL convolution library.";
            logger.log(Level.SEVERE, message, e);
        }
    }
    private final int nX, nY, nZ;

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     */
    public Complex3DOpenCL(int nX, int nY, int nZ, FloatBuffer real,
            FloatBuffer imag, FloatBuffer recip) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;

        int status = init(nX, nY, nZ, real, imag, recip);
        if (status != 1) {
            String message = " Java OpenCL init status: " + status;
            logger.severe(message);
            System.exit(-1);
        }
    }

    @Override
    protected void finalize() throws Throwable {
        try {
            int status = free();
            if (status != 1) {
                String message = " Java OpenCL finalize status: " + status;
                logger.severe(message);
            }
        } finally {
            super.finalize();

        }
    }

    public static void main(String[] args) throws Exception {
        int dimNotFinal = 64;
        int reps = 10;
        if (args != null) {
            try {
                dimNotFinal = Integer.parseInt(args[0]);
                if (dimNotFinal < 1) {
                    dimNotFinal = 100;
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
                "Initializing a %d cubed grid.\n"
                + "The best timing out of %d repititions will be used.",
                dim, reps));

        final int dimCubed = dim * dim * dim;
        int size = Float.SIZE / Byte.SIZE;
        ByteBuffer realTemp = ByteBuffer.allocateDirect(size * dimCubed);
        ByteBuffer imagTemp = ByteBuffer.allocateDirect(size * dimCubed);
        ByteBuffer recipTemp = ByteBuffer.allocateDirect(size * dimCubed);

        /**
         * Set the byte order to that of the native plaform.
         */
        float array[] = new float[1];
        FloatBuffer temp = FloatBuffer.wrap(array);
        realTemp.order(temp.order());
        imagTemp.order(temp.order());
        recipTemp.order(temp.order());

        /**
         * Finally create views of the ByteBuffers as FloatBuffers.
         */
        FloatBuffer real = realTemp.asFloatBuffer();
        FloatBuffer imag = imagTemp.asFloatBuffer();
        FloatBuffer recip = recipTemp.asFloatBuffer();

        /**
         * Create an array to save the initial input and result.
         */
        double orig[] = new double[dimCubed];
        double answer[] = new double[dimCubed];
        double dataArray[] = new double[dimCubed * 2];
        double recipArray[] = new double[dimCubed];

        Complex3D complex3D = new Complex3D(dim, dim, dim);
        Complex3DParallel complex3DParallel =
                new Complex3DParallel(dim, dim, dim, new ParallelTeam());
        Complex3DOpenCL complex3DOpenCL = new Complex3DOpenCL(dim, dim, dim,
                real, imag, recip);

        Random random = new Random(1);
        int index = 0;
        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < dim; j++) {
                for (int i = 0; i < dim; i++) {
                    orig[index] = random.nextFloat();
                    
                    real.put(index, (float) orig[index]);
                    dataArray[index * 2] = orig[index];

                    imag.put(index, 0.0f);
                    dataArray[index * 2 + 1] = 0.0;

                    recip.put(index, (float) orig[index]);
                    recipArray[index] = orig[index];

                    index++;
                }
            }
        }
        real.rewind();
        imag.rewind();
        recip.rewind();

        double toSeconds = 0.000000001;
        long parTime = Long.MAX_VALUE;
        long seqTime = Long.MAX_VALUE;
        long clTime = Long.MAX_VALUE;
        complex3D.setRecip(recipArray);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                dataArray[j * 2] = orig[j];
                dataArray[j * 2 + 1] = 0.0;
            }
            long time = System.nanoTime();
            complex3D.convolution(dataArray);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d Sequential: %8.3f", i + 1, toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }
        }

        for (int j = 0; j < dimCubed; j++) {
            answer[j] = dataArray[j * 2];
        }

        complex3DParallel.setRecip(recipArray);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                dataArray[j * 2] = orig[j];
                dataArray[j * 2 + 1] = 0.0;
            }
            long time = System.nanoTime();
            complex3DParallel.convolution(dataArray);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d Parallel:   %8.3f", i + 1, toSeconds * time));
            if (time < parTime) {
                parTime = time;
            }
        }

        double maxError = Double.MIN_VALUE;
        double rmse = 0.0;
        for (int i = 0; i < dimCubed; i++) {
            double error = Math.abs(answer[i] - dataArray[2 * i]);
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format("Parallel RMSE:   %12.10f, Max: %12.10f", rmse, maxError));

        float result[] = new float[dimCubed];
        for (int i = 0; i < reps; i++) {
            real.rewind();
            imag.rewind();
            for (int j = 0; j < dimCubed; j++) {
                real.put(j, (float) orig[j]);
                imag.put(j, 0.0f);
            }
            real.rewind();
            imag.rewind();
            long time = System.nanoTime();
            complex3DOpenCL.convolution();
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d OpenCL:     %8.3f", i + 1, toSeconds * time));
            if (time < clTime) {
                clTime = time;
            }
        }

        real.rewind();
        real.get(result, 0, dimCubed);
        maxError = Double.MIN_VALUE;
        double avg = 0.0;
        rmse = 0.0;
        for (int j = 0; j < result.length; j++) {
            double error = Math.abs((answer[j] - result[j]) / dimCubed);
            avg += error;
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        avg /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format("OpenCL RMSE:   %12.10f, Max: %12.10f, , Avg: %12.10f", rmse, maxError, avg));

        int status = complex3DOpenCL.free();
        logger.info("Java OpenCL finalize status: " + status);

        System.out.println(String.format("Best Sequential Time:  %8.3f",
                toSeconds * seqTime));
        System.out.println(String.format("Best Parallel Time:    %8.3f",
                toSeconds * parTime));
        System.out.println(String.format("Best OpenCL Time:      %8.3f",
                toSeconds * clTime));
        System.out.println(String.format("Parallel Speedup: %15.5f", (double) seqTime
                / parTime));
        System.out.println(String.format("OpenCL Speedup:   %15.5f", (double) seqTime
                / clTime));
    }
}
