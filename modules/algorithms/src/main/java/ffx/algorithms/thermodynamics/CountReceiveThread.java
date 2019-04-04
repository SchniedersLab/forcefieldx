//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.thermodynamics;

import java.io.IOException;
import java.io.InterruptedIOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.stream;

import edu.rit.mp.DoubleBuf;

/**
 * The CountReceiveThread accumulates TT-OSRW statistics from multiple asynchronous walkers.
 */
class CountReceiveThread extends Thread {

    private static final Logger logger = Logger.getLogger(CountReceiveThread.class.getName());


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
    CountReceiveThread(TransitionTemperedOSRW transitionTemperedOSRW) {
        this.transitionTemperedOSRW = transitionTemperedOSRW;
        recursionCount = new double[3];
        recursionCountBuf = DoubleBuf.buffer(recursionCount);
    }

    /**
     * Run the CountReceiveThread receive thread.
     */
    @Override
    public void run() {
        while (true) {
            try {
                transitionTemperedOSRW.world.receive(null, recursionCountBuf);
            } catch (InterruptedIOException ioe) {
                String message = " CountReceiveThread was interrupted at world.receive; " +
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
                logger.fine(" Termination signal (3x NaN) received; CountReceiveThread shutting down.");
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
                logger.log(Level.FINE, " CountReceiveThread was interrupted; ceasing execution.");
                // No pending message receipt, so no warning.
                break;
            }
        }
    }

}
