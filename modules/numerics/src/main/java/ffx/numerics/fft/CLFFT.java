/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.numerics.fft;

import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Random;

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
public final class CLFFT {

    static {
        System.loadLibrary("JclFFT");
    }

    public static native long clfftSetupNative();

    public static native long clfftCreateDefaultPlanNative(long context, int dim, int dimX, int dimY, int dimZ);

    public static native int clfftSetPlanPrecisionNative(long planHandle, int precision);

    public static native int clfftSetLayoutNative(long planHandle, int inLayout, int outLayout);

    public static native int clfftExecuteTransformNative(long planHandle, int direction, long queue, long rBuffer, long cBuffer);

    public static native int clfftDestroyPlanNative(long planHandle);

    public static native int clfftTeardownNative();

    private PlanHandle planHandle = null;

    /**
     * Constructor.
     *
     * @param context
     * @param dimension
     * @param dims
     */
    public CLFFT(CLContext context, CLFFT_DIMENSION dimension, int dims[]) {
        clfftSetup();
        clfftCreateDefaultPlan(context, CLFFT_DIMENSION.CLFFT_3D, dims);
    }

    private long clfftSetup() {
        return (clfftSetupNative());
    }

    private void clfftCreateDefaultPlan(CLContext context, CLFFT_DIMENSION dimension, int dimLengths[]) {
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
        planHandle = new PlanHandle(clfftCreateDefaultPlanNative(context.ID, dimension.ID, dimX, dimY, dimZ));
    }

    public int clfftSetPlanPrecision(CLFFT_PRECISION precision) {
        return (clfftSetPlanPrecisionNative(planHandle.ID, precision.ID));
    }

    public int clfftSetLayout(CLFFT_LAYOUT inLayout, CLFFT_LAYOUT outLayout) {
        return (clfftSetLayoutNative(planHandle.ID, inLayout.ID, outLayout.ID));
    }

    public int clfftExecuteTransform(CLFFT_DIRECTION direction, CLCommandQueue queue,
            CLBuffer<DoubleBuffer> rBuffer, CLBuffer<DoubleBuffer> cBuffer) {
        return (clfftExecuteTransformNative(planHandle.ID, direction.ID, queue.ID, rBuffer.ID, cBuffer.ID));
    }

    public int clfftDestroyPlan() {
        return (clfftDestroyPlanNative(planHandle.ID));
    }

    public int clfftTeardown() {
        return (clfftTeardownNative());
    }

    public class PlanHandle {

        public final long ID;

        public PlanHandle(long ID) {
            this.ID = ID;
        }

    }

    public enum CLFFT_DIMENSION {

        CLFFT_1D(1), CLFFT_2D(2), CLFFT_3D(3);

        public int ID;

        CLFFT_DIMENSION(int ID) {
            this.ID = ID;
        }
    }

    public enum CLFFT_LAYOUT {

        CLFFT_COMPLEX_INTERLEAVED(0), CLFFT_COMPLEX_PLANAR(1), CLFFT_REAL(2);

        public int ID;

        CLFFT_LAYOUT(int ID) {
            this.ID = ID;
        }
    }

    public enum CLFFT_PRECISION {

        SINGLE(0), DOUBLE(1);

        public int ID;

        CLFFT_PRECISION(int ID) {
            this.ID = ID;
        }
    }

    public enum CLFFT_DIRECTION {

        FORWARD(1), BACKWARD(-1);

        public int ID;

        CLFFT_DIRECTION(int ID) {
            this.ID = ID;
        }
    }

    public static void main(String[] args) {
        CLContext context = CLContext.create();
        CLDevice device = context.getMaxFlopsDevice();
        CLPlatform platform = device.getPlatform();
        CLCommandQueue queue = device.createCommandQueue();
        System.out.format(" Using device: %s\n", device);
        int elementCount = 160;     // must be power of 2,3,5
        //int localWorkSize = Math.min(device.getMaxWorkGroupSize(), 32);  // Local work size dimensions
        //int globalWorkSize = roundUp(localWorkSize, elementCount);   // rounded up to the nearest multiple of the localWorkSize
        int bufferSize = elementCount * elementCount * elementCount * 2; 
        CLBuffer<DoubleBuffer> rBuffer = context.createDoubleBuffer(bufferSize, CLMemory.Mem.READ_WRITE);
        CLBuffer<DoubleBuffer> cBuffer = context.createDoubleBuffer(bufferSize, CLMemory.Mem.READ_WRITE);
        fillBuffer(rBuffer.getBuffer(), 12345);
        zeroBuffer(cBuffer.getBuffer());
        queue.putWriteBuffer(rBuffer, true);
        queue.putWriteBuffer(cBuffer, true);

        System.out.format(" Original Data:------- \n");
        System.out.format(" Real:");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", rBuffer.getBuffer().get(i));
        }
        System.out.format(" Complex:\n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", cBuffer.getBuffer().get(i));
        }
        rBuffer.getBuffer().position(0);
        cBuffer.getBuffer().position(0);

        // forward 1D-FFT
        int dimX = elementCount;
        int dimY = elementCount;
        int dimZ = elementCount;
        int dims[] = {dimX, dimY, dimZ};

        CLFFT clFFT = new CLFFT(context, CLFFT_DIMENSION.CLFFT_3D, dims);
        clFFT.clfftSetPlanPrecision(CLFFT_PRECISION.DOUBLE);
        clFFT.clfftSetLayout(CLFFT_LAYOUT.CLFFT_COMPLEX_PLANAR, CLFFT_LAYOUT.CLFFT_COMPLEX_PLANAR);
        clFFT.clfftExecuteTransform(CLFFT_DIRECTION.FORWARD, queue, rBuffer, cBuffer);
        clFFT.clfftDestroyPlan();
        queue.putReadBuffer(rBuffer, true);
        queue.putReadBuffer(cBuffer, true);

        System.out.format(" Results in Java: \n");
        System.out.format(" Real:");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", rBuffer.getBuffer().get(i));
        }
        System.out.format(" Complex:\n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", cBuffer.getBuffer().get(i));
        }
        rBuffer.getBuffer().position(0);
        cBuffer.getBuffer().position(0);
        queue.putWriteBuffer(rBuffer, true);
        queue.putWriteBuffer(cBuffer, true);

        // reverse 1D-FFT
        clFFT.clfftCreateDefaultPlan(context, CLFFT_DIMENSION.CLFFT_3D, dims);
        clFFT.clfftSetPlanPrecision(CLFFT_PRECISION.DOUBLE);
        clFFT.clfftSetLayout(CLFFT_LAYOUT.CLFFT_COMPLEX_PLANAR, CLFFT_LAYOUT.CLFFT_COMPLEX_PLANAR);
        clFFT.clfftExecuteTransform(CLFFT_DIRECTION.BACKWARD, queue, rBuffer, cBuffer);
        clFFT.clfftDestroyPlan();
        clFFT.clfftTeardown();
        queue.putReadBuffer(rBuffer, true);
        queue.putReadBuffer(cBuffer, true);

        System.out.format(" Reverse Operation Results: \n");
        System.out.format(" Real:");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", rBuffer.getBuffer().get(i));
        }
        System.out.format(" Complex:\n");
        for (int i = 0; i < 10; i++) {
            System.out.format("\t%f\n", cBuffer.getBuffer().get(i));
        }
        rBuffer.getBuffer().position(0);
        cBuffer.getBuffer().position(0);
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

}
