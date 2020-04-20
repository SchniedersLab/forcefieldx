// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
// ******************************************************************************
package ffx.algorithms.thermodynamics;

import static java.lang.String.format;
import static java.util.Arrays.stream;
import static org.apache.commons.math3.util.FastMath.round;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import java.io.IOException;
import java.io.InterruptedIOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The CountReceiveThread accumulates OST statistics from multiple asynchronous walkers.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
class AsynchronousSend extends Thread {

  private static final Logger logger = Logger.getLogger(AsynchronousSend.class.getName());
  /** Storage to send a recursion count. */
  private final double[] myRecursionWeight;
  /** DoubleBuf to wrap the recursion count for sending. */
  private final DoubleBuf myRecursionWeightBuf;
  /** Storage to receive a recursion count. */
  private final double[] recursionCount;
  /** DoubleBuf to wrap the recursion count. */
  private final DoubleBuf recursionCountBuf;
  /** World communicator. */
  private final Comm world = Comm.world();
  /** Rank. */
  private final int rank = world.rank();
  /** Number of processes. */
  private final int numProc = world.size();
  /** Private reference to the Histogram instance. */
  private final Histogram histogram;

  /**
   * Constructor for a thread to asynchronously receive recursion counts.
   *
   * @param histogram Histogram instance.
   */
  AsynchronousSend(Histogram histogram) {
    this.histogram = histogram;
    // Send.
    myRecursionWeight = new double[4];
    myRecursionWeightBuf = DoubleBuf.buffer(myRecursionWeight);
    // Receive.
    recursionCount = new double[4];
    recursionCountBuf = DoubleBuf.buffer(recursionCount);
  }

  /** Run the AsynchronousSend receive thread. */
  @Override
  public void run() {
    while (true) {
      try {
        histogram.world.receive(null, recursionCountBuf);
      } catch (InterruptedIOException ioe) {
        String message =
            " CountReceiveThread was interrupted at world.receive; "
                + "future message passing may be in an error state.";
        logger.log(Level.WARNING, message, ioe);
        break;
      } catch (IOException e) {
        String message = e.getMessage();
        logger.log(Level.WARNING, message, e);
      }

      // 4x NaN is a message (usually sent by the same process) indicating that it is time to shut
      // down.
      boolean terminateSignal = stream(recursionCount).allMatch(Double::isNaN);
      if (terminateSignal) {
        logger.fine(" Termination signal received; CountReceiveThread shutting down.");
        break;
      }

      int rank = (int) round(recursionCount[0]);
      double lambda = recursionCount[1];
      double fLambda = recursionCount[2];

      // If independent, only add bias values from this walker
      if (histogram.getIndependentWalkers() && histogram.getRank() != rank) {
        continue;
      }

      // Check that the FLambda range of the Recursion kernel includes both the minimum and maximum
      // FLambda value.
      histogram.checkRecursionKernelSize(fLambda);

      // Increment the Recursion Kernel based on the input of current walker.
      int walkerLambda = histogram.indexForLambda(lambda);
      int walkerFLambda = histogram.binForFLambda(fLambda);
      double weight = recursionCount[3];

      if (histogram.getResetStatistics() && lambda > histogram.getLambdaResetValue()) {
        histogram.allocateRecursionKernel();
        histogram.disableResetStatistics();
        logger.info(format(" Cleared OST histogram (Lambda = %6.4f).", lambda));
      }

      // Increase the Recursion Kernel based on the input of current walker.
      // Guaranteed to be from a different process.
      histogram.addToRecursionKernelValue(walkerLambda, walkerFLambda, weight, true);
      if (isInterrupted()) {
        logger.log(Level.FINE, " CountReceiveThread was interrupted; ceasing execution.");
        // No pending message receipt, so no warning.
        break;
      }
    }
  }

  /**
   * Send an OST count to all other processes.
   *
   * @param lambda Current value of lambda.
   * @param dUdL Current value of dU/dL.
   */
  public void send(double lambda, double dUdL, double temperingWeight) {
    myRecursionWeight[0] = rank;
    myRecursionWeight[1] = lambda;
    myRecursionWeight[2] = dUdL;
    myRecursionWeight[3] = temperingWeight;

    histogram.setLastReceivedLambda(lambda);
    histogram.setLastReceiveddUdL(dUdL);

    for (int i = 0; i < numProc; i++) {
      try {
        world.send(i, myRecursionWeightBuf);
      } catch (Exception ex) {
        String message = " Asynchronous Multiwalker OST send failed.";
        logger.log(Level.SEVERE, message, ex);
      }
    }
  }
}
