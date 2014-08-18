/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.numerics.fft;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLMemory;
import com.jogamp.opencl.CLPlatform;
import com.jogamp.opencl.CLProgram;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.DoubleBuffer;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;

/**
 * This class implements a Java wrapper for calling clMath functions.
 *
 * @author Michael J. Schnieders
 * @author Stephen LuCore
 */
public final class Complex3DOpenCL implements Runnable {

    static {
        System.loadLibrary("JclFFT");
    }

    private static final Logger logger = Logger.getLogger(Complex3DOpenCL.class.getName());
    private final int nX;
    private final int nY;
    private final int nZ;
    private final int len;
    private double data[];
    private double recip[];
    MODE mode;
    boolean free;
    boolean dead;

    public CLBuffer<DoubleBuffer> dataBuffer;
    public CLBuffer<DoubleBuffer> recipBuffer;
    private PlanHandle planHandle;

    /**
     * Constructor.
     *
     * @param x x-dimension.
     * @param y y-dimension.
     * @param z z-dimension.
     */
    public Complex3DOpenCL(int x, int y, int z) {
        nX = x;
        nY = y;
        nZ = z;
        this.len = nX * nY * nZ;
        mode = null;
        free = false;
        dead = false;
    }

    public void fft(final double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        this.data = data;
        mode = MODE.FFT;
        execute();
    }

    public void ifft(final double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        this.data = data;
        mode = MODE.IFFT;
        execute();
    }

    public void convolution(double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        this.data = data;
        mode = MODE.CONVOLUTION;
        execute();
    }

    public void setRecip(double recip[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        this.recip = recip;
        mode = MODE.RECIP;
        execute();
    }

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

    @Override
    public void run() {
        CLContext context = null;
        try {
            // Initialize the OpenCL Context
            context = CLContext.create();
            CLDevice device = context.getMaxFlopsDevice();
            logger.info(String.format(" Using device: %s\n", device));
            CLPlatform platform = device.getPlatform();
            CLCommandQueue queue = device.createCommandQueue();

            // Allocate memory on the device.
            int bufferSize = len * 2;
            int dims[] = {nX, nY, nZ};
            dataBuffer = context.createDoubleBuffer(bufferSize, CLMemory.Mem.READ_WRITE);
            DoubleBuffer doubleBuffer = dataBuffer.getBuffer();
            int MB = 1024 * 1024;
            logger.info(String.format(" FFT data buffer        [direct: %b, write: %b, size: %d MB]",
                    doubleBuffer.isDirect(), !doubleBuffer.isReadOnly(), dataBuffer.getCLSize() / MB));
            recipBuffer = context.createDoubleBuffer(len, CLMemory.Mem.READ_WRITE);
            doubleBuffer = recipBuffer.getBuffer();
            logger.info(String.format(" Reciprocal data buffer [direct: %b, write: %b, size: %d MB]",
                    doubleBuffer.isDirect(), !doubleBuffer.isReadOnly(), recipBuffer.getCLSize() / MB));

            // Initialize the OpenCL FFT library.
            setup();
            planHandle = createDefaultPlan(context, Complex3DOpenCL_DIMENSION.Complex3DOpenCL_3D, dims);
            setPlanPrecision(Complex3DOpenCL_PRECISION.DOUBLE);
            setLayout(Complex3DOpenCL_LAYOUT.Complex3DOpenCL_COMPLEX_INTERLEAVED,
                    Complex3DOpenCL_LAYOUT.Complex3DOpenCL_COMPLEX_INTERLEAVED);

            // Initialize the Reciprocal Space Multitply Kernal
            URL source = getClass().getClassLoader().getResource("ffx/numerics/fft/VectorMultiply.cl");
            InputStream input = source.openStream();
            CLProgram program = context.createProgram(input).build();

            // Get a reference to the kernel function with the name 'VectorMultiply'
            CLKernel kernel = program.createCLKernel("VectorMultiply");
            int localWorkSize = Math.min(device.getMaxWorkGroupSize(), 256);
            int globalWorkSize = roundUp(localWorkSize, len);

            synchronized (this) {
                while (!free) {
                    if (mode != null) {
                        switch (mode) {
                            case RECIP:
                                doubleBuffer = recipBuffer.getBuffer();
                                doubleBuffer.rewind();
                                doubleBuffer.put(recip);
                                doubleBuffer.rewind();
                                queue.putWriteBuffer(recipBuffer, true);
                                break;
                            case FFT:
                                doubleBuffer = dataBuffer.getBuffer();
                                doubleBuffer.rewind();
                                doubleBuffer.put(data);
                                doubleBuffer.rewind();
                                queue.putWriteBuffer(dataBuffer, false);
                                queue.putBarrier();
                                executeTransform(Complex3DOpenCL_DIRECTION.FORWARD, queue, dataBuffer, dataBuffer);
                                queue.finish();
                                queue.putReadBuffer(dataBuffer, true);
                                doubleBuffer.rewind();
                                doubleBuffer.get(data);
                                break;
                            case CONVOLUTION:
                                doubleBuffer = dataBuffer.getBuffer();
                                doubleBuffer.rewind();
                                doubleBuffer.put(data);
                                doubleBuffer.rewind();
                                queue.putWriteBuffer(dataBuffer, false);
                                queue.putBarrier();
                                // Forward FFT
                                executeTransform(Complex3DOpenCL_DIRECTION.FORWARD, queue, dataBuffer, dataBuffer);
                                queue.putBarrier();
                                // Reciprocal Space Multiply
                                kernel.rewind().putArgs(dataBuffer, recipBuffer).putArg(len);
                                queue.put1DRangeKernel(kernel, 0, globalWorkSize, localWorkSize);
                                queue.putBarrier();
                                // Backward FFT
                                executeTransform(Complex3DOpenCL_DIRECTION.BACKWARD, queue, dataBuffer, dataBuffer);
                                queue.finish();
                                queue.putReadBuffer(dataBuffer, true);
                                doubleBuffer.rewind();
                                doubleBuffer.get(data);
                                break;
                            case IFFT:
                                doubleBuffer = dataBuffer.getBuffer();
                                doubleBuffer.rewind();
                                doubleBuffer.put(data);
                                doubleBuffer.rewind();
                                queue.putWriteBuffer(dataBuffer, true);
                                executeTransform(Complex3DOpenCL_DIRECTION.BACKWARD, queue, dataBuffer, dataBuffer);
                                queue.finish();
                                queue.putReadBuffer(dataBuffer, true);
                                doubleBuffer.rewind();
                                doubleBuffer.get(data);
                                break;
                        }
                        // Reset the mode to null and notify the calling thread.
                        mode = null;
                        notify();
                    }
                    // The OpenCL thread will wait until it's notified again.
                    try {
                        wait();
                    } catch (InterruptedException e) {
                        logger.severe(e.toString());
                    }
                }
                destroyPlan();
                teardown();
                dead = true;
                notify();
            }
        } catch (IOException e) {
            logger.warning(e.toString());
        } finally {
            if (context != null) {
                context.release();
            }
        }
        logger.info(" OpenCL FFT/convolution thread is done.");
    }

    private enum MODE {

        FFT, IFFT, CONVOLUTION, RECIP
    };

    private class PlanHandle {

        public final long ID;

        public PlanHandle(long ID) {
            this.ID = ID;
        }

    }

    private enum Complex3DOpenCL_DIMENSION {

        Complex3DOpenCL_1D(1), Complex3DOpenCL_2D(2), Complex3DOpenCL_3D(3);

        public int ID;

        Complex3DOpenCL_DIMENSION(int ID) {
            this.ID = ID;
        }
    }

    private enum Complex3DOpenCL_LAYOUT {

        Complex3DOpenCL_COMPLEX_INTERLEAVED(0), Complex3DOpenCL_COMPLEX_PLANAR(1), Complex3DOpenCL_REAL(2);

        public int ID;

        Complex3DOpenCL_LAYOUT(int ID) {
            this.ID = ID;
        }
    }

    private enum Complex3DOpenCL_PRECISION {

        SINGLE(0), DOUBLE(1);

        public int ID;

        Complex3DOpenCL_PRECISION(int ID) {
            this.ID = ID;
        }
    }

    private enum Complex3DOpenCL_DIRECTION {

        FORWARD(1), BACKWARD(-1);

        public int ID;

        Complex3DOpenCL_DIRECTION(int ID) {
            this.ID = ID;
        }
    }

    private static void fillBuffer(double buffer[], int seed) {
        Random rnd = new Random(seed);
        for (int i = 0; i < buffer.length; i++) {
            buffer[i] = rnd.nextDouble() * 100;
        }
    }

    private static int roundUp(int groupSize, int globalSize) {
        int r = globalSize % groupSize;
        if (r == 0) {
            return globalSize;
        } else {
            return globalSize + groupSize - r;
        }
    }

    private long setup() {
        return (setupNative());
    }

    private PlanHandle createDefaultPlan(CLContext context, Complex3DOpenCL_DIMENSION dimension, int dimLengths[]) {
        int dimX = 0, dimY = 0, dimZ = 0;
        switch (dimLengths.length) {
            case 3:
                dimZ = dimLengths[2];
            case 2:
                dimY = dimLengths[1];
            case 1:
                dimX = dimLengths[0];
            default:
        }
        return new PlanHandle(createDefaultPlanNative(context.ID, dimension.ID, dimX, dimY, dimZ));
    }

    private int setPlanPrecision(Complex3DOpenCL_PRECISION precision) {
        return (setPlanPrecisionNative(planHandle.ID, precision.ID));
    }

    private int setLayout(Complex3DOpenCL_LAYOUT inLayout, Complex3DOpenCL_LAYOUT outLayout) {
        return (setLayoutNative(planHandle.ID, inLayout.ID, outLayout.ID));
    }

    private int executeTransform(Complex3DOpenCL_DIRECTION direction, CLCommandQueue queue,
            CLBuffer<DoubleBuffer> rBuffer, CLBuffer<DoubleBuffer> cBuffer) {
        return (executeTransformNative(planHandle.ID, direction.ID, queue.ID, rBuffer.ID, cBuffer.ID));
    }

    private int destroyPlan() {
        return (destroyPlanNative(planHandle.ID));
    }

    private int teardown() {
        return (teardownNative());
    }

    private static native long setupNative();

    private static native long createDefaultPlanNative(long context, int dim, int dimX, int dimY, int dimZ);

    private static native int setPlanPrecisionNative(long planHandle, int precision);

    private static native int setLayoutNative(long planHandle, int inLayout, int outLayout);

    private static native int executeTransformNative(long planHandle, int direction, long queue, long rBuffer, long cBuffer);

    private static native int destroyPlanNative(long planHandle);

    private static native int teardownNative();

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
                reps = Integer.parseInt(args[2]);
                if (reps < 1) {
                    reps = 5;
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
        Complex3DParallel complex3DParallel =
                new Complex3DParallel(dim, dim, dim, new ParallelTeam(), IntegerSchedule.fixed());
        complex3DParallel.setRecip(recip);

        Complex3DOpenCL complex3DOpenCL = new Complex3DOpenCL(dim, dim, dim);
        Thread openCLThread = new Thread(complex3DOpenCL);
        openCLThread.setPriority(Thread.MAX_PRIORITY);
        openCLThread.start();

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
            complex3D.convolution(data);
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
            complex3DParallel.convolution(data);
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

        complex3DOpenCL.setRecip(recip);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < dimCubed; j++) {
                data[j * 2] = orig[j];
                data[j * 2 + 1] = 0.0;
            }
            long time = System.nanoTime();
            complex3DOpenCL.convolution(data);
            //complex3DOpenCL.fft(data);
            time = (System.nanoTime() - time);
            System.out.println(String.format(" %2d OpenCL:     %8.3f", i + 1, toSeconds * time));
            if (time < clTime) {
                clTime = time;
            }
        }

        maxError = Double.MIN_VALUE;
        double avg = 0.0;
        rmse = 0.0;
        for (int i = 0; i < dimCubed; i++) {
            double error = Math.abs(answer[i] / dimCubed - data[2 * i]);
            avg += error;
            if (error > maxError) {
                maxError = error;
            }
            rmse += error * error;
        }
        rmse /= dimCubed;
        avg /= dimCubed;
        rmse = Math.sqrt(rmse);
        logger.info(String.format(" OpenCL RMSE:   %12.10f, Max: %12.10f, Avg: %12.10f", rmse, maxError, avg));

        complex3DOpenCL.free();
        complex3DOpenCL = null;

        System.out.println(String.format(" Best Sequential Time:  %8.3f",
                toSeconds * seqTime));
        System.out.println(String.format(" Best Parallel Time:    %8.3f",
                toSeconds * parTime));
        System.out.println(String.format(" Best OpenCL Time:      %8.3f",
                toSeconds * clTime));
        System.out.println(String.format(" Parallel Speedup: %15.5f", (double) seqTime
                / parTime));
        System.out.println(String.format(" OpenCL Speedup:   %15.5f", (double) seqTime
                / clTime));
    }

}
