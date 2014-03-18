/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.numerics.fft;

import java.io.File;
import java.net.URL;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.io.FileUtils;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;

import jcuda.LogLevel;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.*;
import jcuda.jcufft.*;

import static jcuda.driver.JCudaDriver.*;
import static jcuda.jcufft.JCufft.*;

/**
 * Compute a 3D Convolution using Java wrappers to the CUDA Driver API.
 *
 * @author Michal J. Schnieders
 * @since 1.0
 *
 */
public class Real3DCuda implements Runnable {

    private static final Logger logger = Logger.getLogger(Real3DCuda.class.getName());
    private final int nX, nY, nZ, len;
    private float data[], recip[];
    private boolean doConvolution = false;
    private boolean free = false;
    private boolean dead = false;
    CUfunction function;
    CUmodule module;
    cufftHandle planR2C, planC2R;
    Pointer dataPtr, recipPtr;
    CUdeviceptr dataDevice, recipDevice;
    Pointer dataDevicePtr, recipDevicePtr;

    /**
     * Blocking convolution method.
     *
     * @param data an array of float.
     * @return A status flag (0 for success, -1 for failure).
     */
    public int convolution(float data[]) {
        // This would be a programming error.
        if (dead || doConvolution) {
            return -1;
        }

        this.data = data;
        doConvolution = true;

        // Notify the CUDA thread and then block until it notifies us back.
        synchronized (this) {
            notify();
            while (doConvolution) {
                try {
                    wait();
                } catch (InterruptedException e) {
                    logger.severe(e.toString());
                }
            }
        }

        return 0;
    }

    /**
     * Blocking free method.
     *
     * @return A status flag (0 for success, -1 for failure).
     */
    public int free() {
        if (dead || doConvolution) {
            return -1;
        }

        free = true;

        // Notify the CUDA thread and then block until it notifies us back.
        synchronized (this) {
            notify();
            while (!dead) {
                try {
                    wait();
                } catch (InterruptedException e) {
                    logger.severe(e.toString());
                }
            }
        }
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        JCudaDriver.setExceptionsEnabled(true);
        JCudaDriver.setLogLevel(LogLevel.LOG_ERROR);
        JCufft.setExceptionsEnabled(true);
        JCufft.setLogLevel(LogLevel.LOG_ERROR);

        // Initialize the driver and create a context for the first device.
        cuInit(0);
        CUcontext pctx = new CUcontext();
        CUdevice dev = new CUdevice();
        CUdevprop prop = new CUdevprop();
        cuDeviceGetProperties(prop, dev);
        logger.info(" CUDA " + prop.toFormattedString());

        cuDeviceGet(dev, 0);
        cuCtxCreate(pctx, 0, dev);

        // Load the CUBIN file and obtain the "recipSummation" function.
        try {
            String bit = System.getProperty("sun.arch.data.model").trim();
            URL source = getClass().getClassLoader().getResource("ffx/numerics/fft/recipSummation-" + bit + ".cubin");
            File cubinFile = File.createTempFile("recipSummation", "cubin");
            FileUtils.copyURLToFile(source, cubinFile);
            module = new CUmodule();
            cuModuleLoad(module, cubinFile.getCanonicalPath());
            function = new CUfunction();
            cuModuleGetFunction(function, module, "recipSummation");
        } catch (Exception e) {
            String message = "Error loading the reciprocal summation kernel";
            logger.log(Level.SEVERE, message, e);
        }

        // Copy the data array to the device.
        dataDevice = new CUdeviceptr();
        cuMemAlloc(dataDevice, len * Sizeof.FLOAT);
        dataPtr = Pointer.to(data);
        cuMemcpyHtoD(dataDevice, dataPtr, len * Sizeof.FLOAT);

        // Copy the recip array to the device.
        recipDevice = new CUdeviceptr();
        cuMemAlloc(recipDevice, len * Sizeof.FLOAT);
        recipPtr = Pointer.to(recip);
        cuMemcpyHtoD(recipDevice, recipPtr, len * Sizeof.FLOAT);

        // Create a Real to Complex CUFFT plan
        planR2C = new cufftHandle();
        cufftPlan3d(planR2C, nX, nY, nZ, cufftType.CUFFT_R2C);
        cufftSetCompatibilityMode(planR2C, cufftCompatibility.CUFFT_COMPATIBILITY_FFTW_ALL);

        // Create a Complex to Real CUFFT plan
        planC2R = new cufftHandle();
        cufftPlan3d(planC2R, nX, nY, nZ, cufftType.CUFFT_C2R);
        cufftSetCompatibilityMode(planC2R, cufftCompatibility.CUFFT_COMPATIBILITY_FFTW_ALL);

        dataDevicePtr = Pointer.to(dataDevice);
        recipDevicePtr = Pointer.to(recipDevice);

        int threads = 512;
        int nBlocks = len / threads + (len % threads == 0 ? 0 : 1);
        int gridSize = (int) Math.floor(Math.sqrt(nBlocks)) + 1;

        logger.info(format(" CUDA thread initialized with %d threads per block", threads));
        logger.info(format(" Grid Size: (%d x %d x 1).", gridSize, gridSize));

        assert (gridSize * gridSize * threads >= len);

        synchronized (this) {
            while (!free) {
                if (doConvolution) {
                    cuMemcpyHtoD(dataDevice, dataPtr, len * Sizeof.FLOAT);
                    int ret = cufftExecR2C(planR2C, dataDevice, dataDevice);
                    if (ret != cufftResult.CUFFT_SUCCESS) {
                        logger.warning("R2C Result " + cufftResult.stringFor(ret));
                    }

                    // Set up the execution parameters for the kernel
                    cuFuncSetBlockShape(function, threads, 1, 1);
                    int offset = 0;
                    offset = align(offset, Sizeof.POINTER);
                    cuParamSetv(function, offset, dataDevicePtr, Sizeof.POINTER);
                    offset += Sizeof.POINTER;
                    offset = align(offset, Sizeof.POINTER);
                    cuParamSetv(function, offset, recipDevicePtr, Sizeof.POINTER);
                    offset += Sizeof.POINTER;
                    offset = align(offset, Sizeof.INT);
                    cuParamSeti(function, offset, len / 2);
                    offset += Sizeof.INT;
                    cuParamSetSize(function, offset);
                    // Call the kernel function.
                    cuLaunchGrid(function, gridSize, gridSize);

                    ret = cufftExecC2R(planC2R, dataDevice, dataDevice);
                    if (ret != cufftResult.CUFFT_SUCCESS) {
                        logger.warning("C2R Result " + cufftResult.stringFor(ret));
                    }
                    ret = cuMemcpyDtoH(dataPtr, dataDevice, len * Sizeof.FLOAT);
                    doConvolution = false;
                    notify();
                }
                try {
                    wait();
                } catch (InterruptedException e) {
                    logger.severe(e.toString());
                }
            }
            cufftDestroy(planR2C);
            cufftDestroy(planC2R);
            cuMemFree(dataDevice);
            cuMemFree(recipDevice);
            dead = true;
            notify();
        }
        logger.info(" CUDA Thread Done!");
    }

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     * @param data an array of float.
     * @param recip an array of float.
     */
    public Real3DCuda(int nX, int nY, int nZ, float data[], float recip[]) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.len = nX * nY * (nZ + 2);
        this.data = data;
        this.recip = recip;
        doConvolution = false;
        free = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void finalize() throws Throwable {
        try {
            free();
        } finally {
            super.finalize();
        }
    }

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
     */
    public static void main(String[] args) throws Exception {
        int dimNotFinal = 64;
        int reps = 10;
        if (args != null) {
            try {
                dimNotFinal = Integer.parseInt(args[0]);
                if (dimNotFinal < 1) {
                    dimNotFinal = 64;
                }
                reps = Integer.parseInt(args[1]);
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
                "Initializing a %d cubed grid.\n"
                + "The best timing out of %d repititions will be used.",
                dim, reps));

        final int dimCubed = dim * dim * dim;
        final int dimCubed2 = (dim + 2) * dim * dim;

        /**
         * Create an array to save the initial input and result.
         */
        final double orig[] = new double[dimCubed2];
        final double answer[] = new double[dimCubed2];
        final double data[] = new double[dimCubed2];
        final double recip[] = new double[dimCubed];

        final float origf[] = new float[dimCubed2];
        final float dataf[] = new float[dimCubed2];
        final float recipf[] = new float[dimCubed];

        Random randomNumberGenerator = new Random(1);
        int index = 0;
        int index2 = 0;

        /**
         * Row-major order.
         */
        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    float randomNumber = randomNumberGenerator.nextFloat();
                    orig[index] = randomNumber;
                    origf[index] = randomNumber;
                    index++;

                    recip[index2] = 1.0;
                    recipf[index2] = 1.0f;
                    index2++;
                }
                // Padding
                index += 2;
            }
        }

        Real3D real3D = new Real3D(dim, dim, dim);
        Real3DParallel real3DParallel =
                new Real3DParallel(dim, dim, dim, new ParallelTeam(),
                IntegerSchedule.fixed());
        Real3DCuda real3DCUDA = new Real3DCuda(dim, dim, dim, dataf, recipf);

        Thread cudaThread = new Thread(real3DCUDA);
        cudaThread.setPriority(Thread.MAX_PRIORITY);
        cudaThread.start();

        double toSeconds = 0.000000001;
        long parTime = Long.MAX_VALUE;
        long seqTime = Long.MAX_VALUE;
        long clTime = Long.MAX_VALUE;

        real3D.setRecip(recip);
        for (int i = 0; i < reps; i++) {
            System.arraycopy(orig, 0, data, 0, dimCubed2);
            long time = System.nanoTime();
            real3D.convolution(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d Sequential: %8.3f", i + 1, toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }
        }
        System.arraycopy(data, 0, answer, 0, dimCubed2);

        real3DParallel.setRecip(recip);
        for (int i = 0; i < reps; i++) {
            System.arraycopy(orig, 0, data, 0, dimCubed2);
            long time = System.nanoTime();
            real3DParallel.convolution(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d Parallel:   %8.3f", i + 1, toSeconds * time));
            if (time < parTime) {
                parTime = time;
            }
        }
        double maxError = Double.MIN_VALUE;
        double rmse = 0.0;
        index = 0;
        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    double error = Math.abs((orig[index] - data[index] / dimCubed));
                    if (error > maxError) {
                        maxError = error;
                    }
                    rmse += error * error;
                    index++;
                }
                index += 2;
            }
        }
        rmse /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format("Parallel RMSE:   %12.10f, Max: %12.10f", rmse, maxError));

        for (int i = 0; i < reps; i++) {
            System.arraycopy(origf, 0, dataf, 0, dimCubed2);
            long time = System.nanoTime();
            real3DCUDA.convolution(dataf);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d CUDA:     %8.3f", i + 1, toSeconds * time));
            if (time < clTime) {
                clTime = time;
            }
        }
        real3DCUDA.free();
        real3DCUDA = null;

        maxError = Double.MIN_VALUE;
        double avg = 0.0;
        rmse = 0.0;
        index = 0;
        for (int x = 0; x < dim; x++) {
            for (int y = 0; y < dim; y++) {
                for (int z = 0; z < dim; z++) {
                    if (Float.isNaN(dataf[index])) {
                        logger.info(String.format("Not a number %d %d %d", x, y, z));
                        System.exit(-1);
                    }
                    double error = Math.abs(origf[index] - dataf[index]);
                    avg += error;
                    if (error > maxError) {
                        maxError = error;
                    }
                    rmse += error * error;
                    index++;
                }
                index += 2;
            }
        }
        rmse /= dimCubed;
        avg /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format("CUDA RMSE:   %12.10f, Max: %12.10f, Avg: %12.10f", rmse, maxError, avg));

        System.out.println(String.format("Best Sequential Time:  %8.3f",
                toSeconds * seqTime));
        System.out.println(String.format("Best Parallel Time:    %8.3f",
                toSeconds * parTime));
        System.out.println(String.format("Best CUDA Time:        %8.3f",
                toSeconds * clTime));
        System.out.println(String.format("Parallel Speedup: %15.5f", (double) seqTime
                / parTime));
        System.out.println(String.format("CUDA Speedup:     %15.5f", (double) seqTime
                / clTime));
    }
}
