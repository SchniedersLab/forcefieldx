/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.numerics.fft;

import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Random;
import java.util.logging.Logger;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLMemory;
import com.jogamp.opencl.CLPlatform;

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
    MODE mode;
    boolean free;
    boolean dead;

    //private CLCommandQueue queue;
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
        DoubleBuffer doubleBuffer = dataBuffer.getBuffer();
        doubleBuffer.put(data);
        mode = MODE.FFT;
        execute();
    }

    public void ifft(final double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        DoubleBuffer doubleBuffer = dataBuffer.getBuffer();
        doubleBuffer.put(data);
        mode = MODE.IFFT;
        execute();
    }

    public void convolution(double data[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        DoubleBuffer doubleBuffer = dataBuffer.getBuffer();
        doubleBuffer.put(data);
        mode = MODE.CONVOLUTION;
        execute();
    }

    public void setRecip(double recip[]) {
        // This would be a programming error.
        if (dead || mode != null) {
            return;
        }
        DoubleBuffer doubleBuffer = recipBuffer.getBuffer();
        doubleBuffer.put(recip);
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
        // Initialize the OpenCL Context
        CLContext context = CLContext.create();
        CLDevice device = context.getMaxFlopsDevice();
        CLPlatform platform = device.getPlatform();
        CLCommandQueue queue = device.createCommandQueue();
        logger.info(String.format(" Using device: %s\n", device));

        // Allocate memory on the device.
        int bufferSize = len * 2;
        int dims[] = {nX, nY, nZ};
        dataBuffer = context.createDoubleBuffer(bufferSize, CLMemory.Mem.READ_WRITE);
        recipBuffer = context.createDoubleBuffer(len, CLMemory.Mem.WRITE_ONLY);

        // Initialize the OpenCL FFT library.
        setup();
        planHandle = createDefaultPlan(context, CLFFT_DIMENSION.CLFFT_3D, dims);
        setPlanPrecision(CLFFT_PRECISION.DOUBLE);
        setLayout(CLFFT_LAYOUT.CLFFT_COMPLEX_INTERLEAVED, CLFFT_LAYOUT.CLFFT_COMPLEX_INTERLEAVED);

        synchronized (this) {
            while (!free) {
                if (mode != null) {
                    switch (mode) {
                        case RECIP:
                            queue.putWriteBuffer(recipBuffer, true);
                            break;
                        case FFT:
                            queue.putWriteBuffer(dataBuffer, true);
                            executeTransform(CLFFT_DIRECTION.FORWARD, queue, dataBuffer, dataBuffer);
                            queue.putReadBuffer(dataBuffer, true);
                            break;
                        case CONVOLUTION:
                            queue.putWriteBuffer(dataBuffer, true);
                            executeTransform(CLFFT_DIRECTION.FORWARD, queue, dataBuffer, dataBuffer);
                            // Reciprocal Space Multiply Needed...
                            executeTransform(CLFFT_DIRECTION.BACKWARD, queue, dataBuffer, dataBuffer);
                            queue.putReadBuffer(dataBuffer, true);
                            break;
                        case IFFT:
                            queue.putWriteBuffer(dataBuffer, true);
                            executeTransform(CLFFT_DIRECTION.BACKWARD, queue, dataBuffer, dataBuffer);
                            queue.putReadBuffer(dataBuffer, true);
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
            context.release();
            dead = true;
            notify();
        }
        logger.info(" OpenCL FFT/Convolution Thread Done!");
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

    private enum CLFFT_DIMENSION {

        CLFFT_1D(1), CLFFT_2D(2), CLFFT_3D(3);

        public int ID;

        CLFFT_DIMENSION(int ID) {
            this.ID = ID;
        }
    }

    private enum CLFFT_LAYOUT {

        CLFFT_COMPLEX_INTERLEAVED(0), CLFFT_COMPLEX_PLANAR(1), CLFFT_REAL(2);

        public int ID;

        CLFFT_LAYOUT(int ID) {
            this.ID = ID;
        }
    }

    private enum CLFFT_PRECISION {

        SINGLE(0), DOUBLE(1);

        public int ID;

        CLFFT_PRECISION(int ID) {
            this.ID = ID;
        }
    }

    private enum CLFFT_DIRECTION {

        FORWARD(1), BACKWARD(-1);

        public int ID;

        CLFFT_DIRECTION(int ID) {
            this.ID = ID;
        }
    }

    private static void fillBuffer(FloatBuffer buffer, int seed) {
        Random rnd = new Random(seed);
        while (buffer.remaining() != 0) {
            buffer.put(rnd.nextFloat() * 100);
        }
        buffer.rewind();
    }

    private static void fillBuffer(DoubleBuffer buffer, int seed) {

        Random rnd = new Random(seed);
        while (buffer.remaining() != 0) {
            buffer.put(rnd.nextFloat() * 100);
        }
        buffer.rewind();
    }

    private static int roundUp(int groupSize, int globalSize) {
        int r = globalSize % groupSize;
        if (r == 0) {
            return globalSize;
        } else {
            return globalSize + groupSize - r;
        }
    }

    private static void zeroBuffer(DoubleBuffer buffer) {
        while (buffer.remaining() != 0) {
            buffer.put(0);
        }
        buffer.rewind();
    }

    private long setup() {
        return (setupNative());
    }

    private PlanHandle createDefaultPlan(CLContext context, CLFFT_DIMENSION dimension, int dimLengths[]) {
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

    private int setPlanPrecision(CLFFT_PRECISION precision) {
        return (setPlanPrecisionNative(planHandle.ID, precision.ID));
    }

    private int setLayout(CLFFT_LAYOUT inLayout, CLFFT_LAYOUT outLayout) {
        return (setLayoutNative(planHandle.ID, inLayout.ID, outLayout.ID));
    }

    private int executeTransform(CLFFT_DIRECTION direction, CLCommandQueue queue,
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

    public static void main(String[] args) {

        int x = 160;
        int y = 160;
        int z = 160;
        int length = x * y * z * 2;
        double data[] = new double[length];

        // Create a Complex, 3D instance of FFT/Convolution operations, backed by an OpenCL implementation.
        Complex3DOpenCL complex3DOpenCL = new Complex3DOpenCL(x, y, z);

        // Run the operations in a thread.
        Thread thread = new Thread(complex3DOpenCL);
        thread.start();

        CLBuffer<DoubleBuffer> buffer = complex3DOpenCL.dataBuffer;
        fillBuffer(buffer.getBuffer(), 12345);
        buffer.getBuffer().get(data);

        System.out.format(" Original Data:------- \n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", data[i]);
        }

        complex3DOpenCL.fft(data);
        buffer.getBuffer().get(data);

        System.out.format(" Results from Forward Transform: \n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", data[i]);
        }

        complex3DOpenCL.ifft(data);
        buffer.getBuffer().get(data);

        System.out.format(" Results from Backward Transform: \n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", data[i]);
        }

        complex3DOpenCL.convolution(data);
        buffer.getBuffer().get(data);

        System.out.format(" Results from Convolution: \n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", data[i]);
        }

        complex3DOpenCL.free();
    }

}
