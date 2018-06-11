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

import java.io.File;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
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
import jcuda.driver.CUcontext;
import jcuda.driver.CUctx_flags;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUdevprop;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.jcufft.JCufft;
import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import static jcuda.driver.JCudaDriver.CU_MEMHOSTALLOC_DEVICEMAP;
import static jcuda.driver.JCudaDriver.align;
import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuDeviceGetProperties;
import static jcuda.driver.JCudaDriver.cuFuncSetBlockShape;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLaunchGrid;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemFree;
import static jcuda.driver.JCudaDriver.cuMemFreeHost;
import static jcuda.driver.JCudaDriver.cuMemHostAlloc;
import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoad;
import static jcuda.driver.JCudaDriver.cuParamSetSize;
import static jcuda.driver.JCudaDriver.cuParamSeti;
import static jcuda.driver.JCudaDriver.cuParamSetv;
import static jcuda.jcufft.JCufft.CUFFT_FORWARD;
import static jcuda.jcufft.JCufft.CUFFT_INVERSE;
import static jcuda.jcufft.JCufft.cufftDestroy;
import static jcuda.jcufft.JCufft.cufftExecZ2Z;
import static jcuda.jcufft.JCufft.cufftPlan3d;

/**
 * Compute a 3D Convolution using Java wrappers to the CUDA Driver API.
 *
 * @author Michal J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("deprecation")
public class Complex3DCuda implements Runnable {

    private static final Logger logger = Logger.getLogger(Complex3DCuda.class.getName());
    private final int nX, nY, nZ, len;
    private MODE mode = null;
    private boolean free = false;
    private boolean dead = false;
    CUfunction function;
    CUmodule module;
    cufftHandle plan;
    Pointer pinnedMemory;
    DoubleBuffer pinnedMemoryBuffer;
    Pointer recipCPUPtr;
    Pointer dataGPUPtr, recipGPUPtr;
    CUdeviceptr dataDevice, recipDevice;
    private final boolean usePinnedMemory = false;

    private enum MODE {

        FFT, IFFT, CONVOLUTION, RECIP
    };

    public DoubleBuffer getDoubleBuffer() {
        return pinnedMemoryBuffer;
    }

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     */
    public Complex3DCuda(int nX, int nY, int nZ) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        this.len = nX * nY * nZ;
        mode = null;
        free = false;
    }

    /**
     * Compute the 3D FFT using CUDA.
     *
     * @param data The input data array must be of size 2 * nX * nY * nZ.
     * @since 1.0
     */
    public void fft(final double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        mode = MODE.FFT;
        execute();
    }

    /**
     * Compute the inverse 3D FFT using CUDA.
     *
     * @param data The input data array must be of size 2 * nX * nY * nZ.
     * @since 1.0
     */
    public void ifft(final double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        mode = MODE.IFFT;
        execute();
    }

    /**
     * Blocking convolution method.
     *
     * @param data Input/output data array.
     */
    public void convolution(double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        // Do the CUDA computation.
        mode = MODE.CONVOLUTION;
        execute();
    }

    /**
     * <p>
     * Setter for the field <code>recip</code>.</p>
     *
     * @param recip an array of double.
     */
    public void setRecip(double recip[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        recipCPUPtr = Pointer.to(recip);
        mode = MODE.RECIP;
        execute();
    }

    /**
     * Notify the CUDA thread and then wait until the requested operation
     * completes.
     *
     * @return Returns 0 upon success.
     */
    private int execute() {
        // Notify the CUDA thread and then wait until it notifies us back.
        synchronized (this) {
            notify();
            while (mode != null) {
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
        if (dead || mode != null) {
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
        JCufft.setExceptionsEnabled(true);
        JCudaDriver.setLogLevel(LogLevel.LOG_INFO);
        JCufft.setLogLevel(LogLevel.LOG_INFO);
        JCufft.initialize();

        // Initialize the driver and create a context for the first device.
        cuInit(0);
        CUcontext pctx = new CUcontext();
        CUdevice dev = new CUdevice();
        CUdevprop prop = new CUdevprop();
        cuDeviceGetProperties(prop, dev);
        logger.log(Level.INFO, "   CUDA {0}", prop.toFormattedString());
        cuDeviceGet(dev, 0);

        // Create a context that allows the GPU to map pinned host memory.
        if (usePinnedMemory) {
            cuCtxCreate(pctx, CUctx_flags.CU_CTX_MAP_HOST, dev);
        } else {
            // Create a context that does not allows the GPU to map pinned host memory.
            cuCtxCreate(pctx, 0, dev);
        }

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
            String message = " Error loading the reciprocal summation kernel";
            logger.log(Level.SEVERE, message, e);
        }

        pinnedMemory = new Pointer();

        if (usePinnedMemory) {
            // Allocate pinned memory mapped into the GPU address space.
            cuMemHostAlloc(pinnedMemory, len * 2 * Sizeof.DOUBLE, CU_MEMHOSTALLOC_DEVICEMAP);
        } else {
            // Allocate memory
            cuMemHostAlloc(pinnedMemory, len * 2 * Sizeof.DOUBLE, 0);
        }

        ByteBuffer byteBuffer = pinnedMemory.getByteBuffer(0, len * 2 * Sizeof.DOUBLE);
        byteBuffer.order(ByteOrder.nativeOrder());
        pinnedMemoryBuffer = byteBuffer.asDoubleBuffer();

        // Allocate a work array on the device.
        dataDevice = new CUdeviceptr();
        cuMemAlloc(dataDevice, len * 2 * Sizeof.DOUBLE);

        // Allocate memory on the device for the reciprocal space array.
        recipDevice = new CUdeviceptr();
        cuMemAlloc(recipDevice, len * Sizeof.DOUBLE);

        // Create and execute a JCufft plan for the data
        plan = new cufftHandle();

        cufftPlan3d(plan, nZ, nY, nX, cufftType.CUFFT_Z2Z);
        //cufftSetCompatibilityMode(plan, cufftCompatibility.CUFFT_COMPATIBILITY_FFTW_ALL);

        dataGPUPtr = Pointer.to(dataDevice);
        recipGPUPtr = Pointer.to(recipDevice);

        int threads = prop.maxThreadsPerBlock;
        int nBlocks = len / threads + (len % threads == 0 ? 0 : 1);
        int gridSize = (int) Math.floor(Math.sqrt(nBlocks)) + 1;

        logger.info(format("   CUDA thread initialized: %d threads per block", threads));
        logger.info(format("   Grid Size:                     (%3d,%3d,%3d)", gridSize, gridSize, 1));

        assert (gridSize * gridSize * threads >= len);

        synchronized (this) {
            while (!free) {
                if (mode != null) {
                    switch (mode) {

                        case RECIP:
                            cuMemcpyHtoD(recipDevice, recipCPUPtr, len * Sizeof.DOUBLE);
                            break;

                        case FFT:
                            // Zero Copy
                            if (usePinnedMemory) {
                                cufftExecZ2Z(plan, pinnedMemory, pinnedMemory, CUFFT_FORWARD);
                            } else {
                                cuMemcpyHtoD(dataDevice, pinnedMemory, 2 * len * Sizeof.DOUBLE);
                                cufftExecZ2Z(plan, dataDevice, dataDevice, CUFFT_FORWARD);
                                cuMemcpyDtoH(pinnedMemory, dataDevice, 2 * len * Sizeof.DOUBLE);
                            }
                            break;

                        case CONVOLUTION:

                            if (usePinnedMemory) {
                                // Zero Copy
                                cufftExecZ2Z(plan, pinnedMemory, dataDevice, CUFFT_FORWARD);
                            } else {
                                // Copy data to device and run forward FFT.
                                cuMemcpyHtoD(dataDevice, pinnedMemory, 2 * len * Sizeof.DOUBLE);
                                cufftExecZ2Z(plan, dataDevice, dataDevice, CUFFT_FORWARD);
                            }

                            // Set up the execution parameters for the kernel
                            cuFuncSetBlockShape(function, threads, 1, 1);
                            int offset = 0;
                            offset = align(offset, Sizeof.POINTER);
                            cuParamSetv(function, offset, dataGPUPtr, Sizeof.POINTER);
                            offset += Sizeof.POINTER;
                            offset = align(offset, Sizeof.POINTER);
                            cuParamSetv(function, offset, recipGPUPtr, Sizeof.POINTER);
                            offset += Sizeof.POINTER;
                            offset = align(offset, Sizeof.INT);
                            cuParamSeti(function, offset, len);
                            offset += Sizeof.INT;
                            cuParamSetSize(function, offset);
                            // Call the kernel function.
                            cuLaunchGrid(function, gridSize, gridSize);
                            if (usePinnedMemory) {
                                // Zero Copy
                                cufftExecZ2Z(plan, dataDevice, pinnedMemory, CUFFT_INVERSE);
                            } else {
                                // Perform inverse FFT and copy memory back to the CPU.
                                cufftExecZ2Z(plan, dataDevice, dataDevice, CUFFT_INVERSE);
                                cuMemcpyDtoH(pinnedMemory, dataDevice, 2 * len * Sizeof.DOUBLE);
                            }
                            break;

                        case IFFT:
                            // Zero Copy
                            if (usePinnedMemory) {
                                cufftExecZ2Z(plan, pinnedMemory, pinnedMemory, CUFFT_INVERSE);
                            } else {
                                cuMemcpyHtoD(dataDevice, pinnedMemory, 2 * len * Sizeof.DOUBLE);
                                cufftExecZ2Z(plan, dataDevice, dataDevice, CUFFT_INVERSE);
                                cuMemcpyDtoH(pinnedMemory, dataDevice, 2 * len * Sizeof.DOUBLE);
                            }
                            break;
                    }

                    // Block for the context's tasks to complete.
                    cuCtxSynchronize();

                    // Reset the mode to null and notify the calling thread.
                    mode = null;
                    notify();
                }
                // The CUDA thread will wait until it's notified again.
                try {
                    wait();
                } catch (InterruptedException e) {
                    logger.severe(e.toString());
                }
            }
            cufftDestroy(plan);
            cuMemFree(dataDevice);
            cuMemFree(recipDevice);
            cuMemFreeHost(pinnedMemory);
            dead = true;
            notify();
        }
        logger.info(" CUDA Thread Done!");
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
     * <p>
     * main</p>
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
                    reps = 10;
                }
            } catch (Exception e) {
            }
        }
        final int dim = dimNotFinal;

        System.out.println(String.format(
                " Initializing a %d cubed grid.\n"
                + " The best timing out of %d repititions will be used.",
                dim, reps));

        final int dimCubed = dim * dim * dim;

        /**
         * Create an array to save the initial input and result.
         */
        double orig[] = new double[dimCubed];
        double answer[] = new double[dimCubed];
        double data[] = new double[dimCubed * 2];
        double recip[] = new double[dimCubed];

        Random random = new Random(1);
        int index = 0;
        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < dim; j++) {
                for (int i = 0; i < dim; i++) {
                    orig[index] = random.nextDouble();
                    //recip[index] = orig[index];
                    recip[index] = 1.0;
                    index++;
                }
            }
        }

        Complex3D complex3D = new Complex3D(dim, dim, dim);
        Complex3DParallel complex3DParallel
                = new Complex3DParallel(dim, dim, dim, new ParallelTeam(), IntegerSchedule.fixed());
        complex3DParallel.setRecip(recip);

        Complex3DCuda complex3DCUDA = new Complex3DCuda(dim, dim, dim);
        Thread cudaThread = new Thread(complex3DCUDA);
        cudaThread.setPriority(Thread.MAX_PRIORITY);
        cudaThread.start();
        complex3DCUDA.setRecip(recip);

        double toSeconds = 0.000000001;
        long parTime = Long.MAX_VALUE;
        long seqTime = Long.MAX_VALUE;
        long clTime = Long.MAX_VALUE;

        complex3D.setRecip(recip);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                data[j * 2] = orig[j];
                data[j * 2 + 1] = 0.0;
            }
            long time = System.nanoTime();
            //complex3D.convolution(data);
            complex3D.fft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format(" %2d Sequential: %8.3f", i + 1, toSeconds * time));
            if (time < seqTime) {
                seqTime = time;
            }
        }

        for (int j = 0; j < dimCubed; j++) {
            answer[j] = data[j * 2];
        }

        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                data[j * 2] = orig[j];
                data[j * 2 + 1] = 0.0;
            }
            long time = System.nanoTime();
            //complex3DParallel.convolution(data);
            complex3DParallel.fft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format(" %2d Parallel:   %8.3f", i + 1, toSeconds * time));
            if (time < parTime) {
                parTime = time;
            }
        }

        double maxError = Double.MIN_VALUE;
        double rmse = 0.0;
        for (int i = 0; i < dimCubed; i++) {
            double error = Math.abs(answer[i] - data[2 * i]);
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format(" Parallel RMSE:   %12.10f, Max: %12.10f", rmse, maxError));

        DoubleBuffer cudaBuffer = complex3DCUDA.getDoubleBuffer();
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                // data[j * 2] = orig[j];
                // data[j * 2 + 1] = 0.0;
                cudaBuffer.put(j * 2, orig[j]);
                cudaBuffer.put(j * 2 + 1, 0.0);
            }
            long time = System.nanoTime();
            //complex3DCUDA.convolution(data);
            complex3DCUDA.fft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format(" %2d CUDA:     %8.3f", i + 1, toSeconds * time));
            if (time < clTime) {
                clTime = time;
            }
        }

        maxError = Double.MIN_VALUE;
        double avg = 0.0;
        rmse = 0.0;
        for (int i = 0; i < dimCubed; i++) {
            double error = Math.abs(answer[i] - cudaBuffer.get(2 * i));
            // double error = Math.abs(answer[i] / dimCubed -  data[2 * i]);
            avg += error;
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        avg /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format(" CUDA RMSE:   %12.10f, Max: %12.10f, Avg: %12.10f", rmse, maxError, avg));

        complex3DCUDA.free();
        complex3DCUDA = null;

        System.out.println(String.format(" Best Sequential Time:  %8.3f",
                toSeconds * seqTime));
        System.out.println(String.format(" Best Parallel Time:    %8.3f",
                toSeconds * parTime));
        System.out.println(String.format(" Best CUDA Time:        %8.3f",
                toSeconds * clTime));
        System.out.println(String.format(" Parallel Speedup: %15.5f", (double) seqTime
                / parTime));
        System.out.println(String.format(" CUDA Speedup:     %15.5f", (double) seqTime
                / clTime));
    }
}
