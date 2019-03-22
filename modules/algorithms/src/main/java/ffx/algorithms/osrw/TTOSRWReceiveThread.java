package ffx.algorithms.osrw;

import java.io.IOException;
import java.io.InterruptedIOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.stream;

import edu.rit.mp.DoubleBuf;

/**
 * The TTOSRWReceiveThread accumulates TT-OSRW statistics from multiple asynchronous walkers.
 */
class TTOSRWReceiveThread extends Thread {

    private static final Logger logger = Logger.getLogger(TTOSRWReceiveThread.class.getName());


    /**
     * Private reference to the TTOSRW instance.
     */
    private TransitionTemperedOSRW transitionTemperedOSRW;

    /**
     * Storage to recieve a recursion count.
     */
    private final double[] recursionCount;

    /**
     * DoubleBuf to wrap the recursion count.
     */
    private final DoubleBuf recursionCountBuf;

    /**
     * Constructor for a thread to asynchronously receive recursion counts.
     *
     * @param transitionTemperedOSRW Private TransitionTemperedOSRW instance of.
     */
    TTOSRWReceiveThread(TransitionTemperedOSRW transitionTemperedOSRW) {
        this.transitionTemperedOSRW = transitionTemperedOSRW;
        recursionCount = new double[3];
        recursionCountBuf = DoubleBuf.buffer(recursionCount);
    }

    /**
     * Run the TTOSRWReceiveThread receive thread.
     */
    @Override
    public void run() {
        while (true) {
            try {
                transitionTemperedOSRW.world.receive(null, recursionCountBuf);
            } catch (InterruptedIOException ioe) {
                String message = " TTOSRWReceiveThread was interrupted at world.receive; " +
                        "future message passing may be in an error state.";
                logger.log(Level.WARNING, message, ioe);
                break;
            } catch (IOException e) {
                String message = e.getMessage();
                logger.log(Level.WARNING, message, e);
            }

            // 3x NaN is a message (usually sent by the same process) indicating that it is time to shut down.
            boolean terminateSignal = stream(recursionCount).allMatch(Double::isNaN);
            if (terminateSignal) {
                logger.fine(" Termination signal (3x NaN) received; TTOSRWReceiveThread shutting down.");
                break;
            }

            // Check that the FLambda range of the Recursion kernel includes both the minimum and maximum FLambda value.
            transitionTemperedOSRW.checkRecursionKernelSize(recursionCount[1]);

            // Increment the Recursion Kernel based on the input of current walker.
            int walkerLambda = transitionTemperedOSRW.binForLambda(recursionCount[0]);
            int walkerFLambda = transitionTemperedOSRW.binForFLambda(recursionCount[1]);
            double weight = recursionCount[2];

            // If the weight is less than 1.0, then a walker has activated tempering.
            if (!transitionTemperedOSRW.isTempering() && weight < 1.0) {
                transitionTemperedOSRW.setTempering(true);
                logger.info(format(" Tempering activated due to received weight of (%8.6f)", weight));
            }

            if (transitionTemperedOSRW.resetStatistics && recursionCount[0] > transitionTemperedOSRW.lambdaResetValue) {
                transitionTemperedOSRW.allocateRecursionKernel();
                transitionTemperedOSRW.resetStatistics = false;
                logger.info(format(" Cleared OSRW histogram (Lambda = %6.4f).", recursionCount[0]));
            }

            // Increase the Recursion Kernel based on the input of current walker.
            transitionTemperedOSRW.addToRecursionKernelValue(walkerLambda, walkerFLambda, weight);
            if (isInterrupted()) {
                logger.log(Level.FINE, " TTOSRWReceiveThread was interrupted; ceasing execution.");
                // No pending message receipt, so no warning.
                break;
            }
        }
    }

}
