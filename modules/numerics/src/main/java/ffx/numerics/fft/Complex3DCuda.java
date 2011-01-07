/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import static java.lang.String.format;

import static jcuda.driver.JCudaDriver.*;
import static jcuda.jcufft.JCufft.*;

import java.io.File;
import java.net.URL;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FileUtils;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;

import jcuda.*;
import jcuda.driver.*;
import jcuda.jcufft.*;

/**
 * Compute a 3D Convolution using Java wrappers to the CUDA Driver API.
 *
 * @author Michal J. Schnieders
 * @since 1.0
 */
public class Complex3DCuda implements Runnable {

    private static final Logger logger = Logger.getLogger(Complex3DCuda.class.getName());
    private final int nX, nY, nZ, len;
    private float data[], recip[];
    private boolean doConvolution = false;
    private boolean free = false;
    private boolean dead = false;
    private int status;
    CUfunction function;
    CUmodule module;
    cufftHandle plan;
    Pointer dataPtr, recipPtr;
    CUdeviceptr dataDevice, recipDevice;
    Pointer dataDevicePtr, recipDevicePtr;

    /**
     * Blocking convolution method.
     * @param data
     * @return
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

        return status;
    }

    /**
     * Blocking free method.
     * @return
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

        return status;
    }

    @Override
    public void run() {
        JCudaDriver.setExceptionsEnabled(true);
        JCufft.setExceptionsEnabled(true);
        JCudaDriver.setLogLevel(LogLevel.LOG_ERROR);

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
        cuMemAlloc(dataDevice, len * 2 * Sizeof.FLOAT);
        dataPtr = Pointer.to(data);
        cuMemcpyHtoD(dataDevice, dataPtr, len * 2 * Sizeof.FLOAT);

        // Copy the recip array to the device.
        recipDevice = new CUdeviceptr();
        cuMemAlloc(recipDevice, len * Sizeof.FLOAT);
        recipPtr = Pointer.to(recip);
        cuMemcpyHtoD(recipDevice, recipPtr, len * Sizeof.FLOAT);

        // Create and execute a JCufft plan for the data
        plan = new cufftHandle();
        cufftPlan3d(plan, nZ, nY, nX, cufftType.CUFFT_C2C);

        dataDevicePtr = Pointer.to(dataDevice);
        recipDevicePtr = Pointer.to(recipDevice);
        
        int threads = 512;
        int nBlocks = len/threads + (len%threads == 0?0:1);
        int gridSize = (int) Math.floor(Math.sqrt(nBlocks)) + 1;

        logger.info(format(" CUDA thread initialized with %d threads per block", threads));
        logger.info(format(" Grid Size: (%d x %d x 1).", gridSize, gridSize));

        assert (gridSize * gridSize * threads >= len);

        synchronized (this) {
            while (!free) {
                if (doConvolution) {
                    cuMemcpyHtoD(dataDevice, dataPtr, len * 2 * Sizeof.FLOAT);
                    cufftExecC2C(plan, dataDevice, dataDevice, CUFFT_FORWARD);
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
                    cuParamSeti(function, offset, len);
                    offset += Sizeof.INT;
                    cuParamSetSize(function, offset);
                    // Call the kernel function.
                    cuLaunchGrid(function,gridSize,gridSize);
                    cufftExecC2C(plan, dataDevice, dataDevice, CUFFT_INVERSE);
                    cuMemcpyDtoH(dataPtr, dataDevice, len * 2 * Sizeof.FLOAT);
                    doConvolution = false;
                    notify();
                }
                try {
                    wait();
                } catch (InterruptedException e) {
                    logger.severe(e.toString());
                }
            }
            cufftDestroy(plan);
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
     */
    public Complex3DCuda(int nX, int nY, int nZ, float data[], float recip[]) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.len = nX * nY * nZ;
        this.data = data;
        this.recip = recip;
        doConvolution = false;
        free = false;
    }

    @Override
    protected void finalize() throws Throwable {
        try {
            free();
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
                    dimNotFinal = 64;
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


        /**
         * Create an array to save the initial input and result.
         */
        double orig[] = new double[dimCubed];
        double answer[] = new double[dimCubed];

        double dataArray[] = new double[dimCubed * 2];
        double recipArray[] = new double[dimCubed];

        float data[] = new float[dimCubed * 2];
        float recip[] = new float[dimCubed];

        Random random = new Random(1);
        int index = 0;
        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < dim; j++) {
                for (int i = 0; i < dim; i++) {
                    orig[index] = random.nextFloat();
                    dataArray[index * 2] = orig[index];
                    dataArray[index * 2 + 1] = 0.0;
                    recipArray[index] = orig[index];

                    data[index * 2] = (float) orig[index];
                    data[index * 2 + 1] = 0.0f;
                    recip[index] = (float) orig[index];

                    index++;
                }
            }
        }

        Complex3D complex3D = new Complex3D(dim, dim, dim);
        Complex3DParallel complex3DParallel =
                          new Complex3DParallel(dim, dim, dim, new ParallelTeam(), IntegerSchedule.fixed());
        Complex3DCuda complex3DCUDA = new Complex3DCuda(dim, dim, dim, data, recip);
        Thread cudaThread = new Thread(complex3DCUDA);
        cudaThread.setPriority(Thread.MAX_PRIORITY);
        cudaThread.start();

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
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                data[j * 2] = (float) orig[j];
                data[j * 2 + 1] = 0.0f;
            }
            long time = System.nanoTime();
            complex3DCUDA.convolution(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format("%2d CUDA:     %8.3f", i + 1, toSeconds * time));
            if (time < clTime) {
                clTime = time;
            }
        }

        maxError = Double.MIN_VALUE;
        double avg = 0.0;
        rmse = 0.0;
        for (int i = 0; i < dimCubed; i++) {
            double error = Math.abs((answer[i] - data[2 * i]) / dimCubed);
            avg += error;
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        avg /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format("CUDA RMSE:   %12.10f, Max: %12.10f, Avg: %12.10f", rmse, maxError, avg));

        complex3DCUDA.free();
        complex3DCUDA = null;

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
