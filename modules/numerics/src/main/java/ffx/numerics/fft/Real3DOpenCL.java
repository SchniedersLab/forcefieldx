/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.numerics.fft;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLMemory;
import com.jogamp.opencl.CLPlatform;
import java.nio.DoubleBuffer;
import jclfft.JclFFT;
import jclfft.JclFFT.*;
import jclfft.PlanHandle;

/**
 *
 * @author slucore
 */
public class Real3DOpenCL {
    
    private final int n, nextX, nextY, nextZ;
    private final int nX, nY, nZ;
    private final int nX1, nZ2;
    private final double work[];
    private final double recip[];
    private final Real fftX;
    private final Complex fftY, fftZ;

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     */
    public Real3DOpenCL(int nX, int nY, int nZ) {
        this.n = nX;
        this.nX = nX / 2;
        this.nY = nY;
        this.nZ = nZ;
        nX1 = this.nX + 1;
        nZ2 = 2 * nZ;
        nextX = 2;
        nextY = n + 2;
        nextZ = nextY * nY;
        work = new double[nZ2];
        recip = new double[nX1 * nY * nZ];
        fftX = new Real(n);
        fftY = new Complex(nY);
        fftZ = new Complex(nZ);
        
        JclFFT.clfftSetup();
    }
    
    /**
     * Compute the 3D FFT.
     *
     * @param input The input array must be of size (nX + 2) * nY * nZ.
     */
    public void fft(final double input[]) {
        CLContext context = CLContext.create();
        CLDevice device = context.getMaxFlopsDevice();
        CLPlatform platform = device.getPlatform();
        CLCommandQueue queue = device.createCommandQueue();
        System.out.format(" Using device: %s\n", device);
        int elementCount = 512;                                             // must be power of 2,3,5
        int localWorkSize = Math.min(device.getMaxWorkGroupSize(), 256);    // Local work size dimensions
        int globalWorkSize = roundUp(localWorkSize, elementCount);          // rounded up to the nearest multiple of the localWorkSize
        CLBuffer<DoubleBuffer> buffer = context.createDoubleBuffer(globalWorkSize, CLMemory.Mem.READ_WRITE);
        
        // write input data to device
        buffer.getBuffer().put(input);        
        queue.putWriteBuffer(buffer, true);

        // forward 1D-FFT
        int dimension = 1, direction = 1;
        int dimX = elementCount, dimY = 0, dimZ = 0;
        PlanHandle planHandle = new PlanHandle(JclFFT.clfftCreateDefaultPlan(context, CLFFT_DIMENSION.CLFFT_1D, dimX));
        planHandle.setPlanPrecision(CLFFT_PRECISION.DOUBLE);
        planHandle.setLayout(CLFFT_LAYOUT.CLFFT_COMPLEX_INTERLEAVED, CLFFT_LAYOUT.CLFFT_COMPLEX_INTERLEAVED);
        planHandle.executeTransform(CLFFT_DIRECTION.FORWARD, queue, buffer);
        planHandle.destroyPlan();
        
        queue.putReadBuffer(buffer, true);
        buffer.getBuffer().position(0);
    }

    private static int roundUp(int groupSize, int globalSize) {
        int r = globalSize % groupSize;
        if (r == 0) {
            return globalSize;
        } else {
            return globalSize + groupSize - r;
        }
    }
    
}
