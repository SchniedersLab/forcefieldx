package ffx.algorithms.thermodynamics;

import static edu.rit.mp.DoubleBuf.buffer;
import static java.lang.String.format;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Synchronous (blocking) communication of OST counts.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class SendSynchronous {

  private static final Logger logger = Logger.getLogger(SendSynchronous.class.getName());

  /** Parallel Java world communicator. */
  protected final Comm world = Comm.world();
  /** Rank of this process. */
  protected final int rank = world.rank();
  /** Number of processes. */
  private final int numProc = world.size();
  /**
   * The counts array stores [Lambda, dU/dL, temper] for each process. Therefore the array is of size
   * [numProc][3].
   *
   * <p>Each 3 entry array must be wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final double[][] counts;
  /**
   * myCounts is a convenience pointer for this process to counts[rank].
   */
  private final double[] myCounts;
  /**
   * countsBuf wraps the counts arrays for each process.
   */
  private final DoubleBuf[] countsBuf;
  /**
   * myCountsBuf is a convenience pointer for this process to countsBuf[rank].
   */
  private final DoubleBuf myCountsBuf;
  /**
   * The histograms to update.
   */
  private Histogram[] histograms;
  /**
   * Map from ranks to histograms.
   */
  private int[] rankToHistogramMap;

  /**
   * Synchronous Communication Constructor.
   *
   * @param histograms An array of Bias Histograms.
   * @param rankToHistogramMap A map from process rank to Histogram.
   */
  public SendSynchronous(Histogram[] histograms, int[] rankToHistogramMap) {
    counts = new double[numProc][3];
    countsBuf = new DoubleBuf[numProc];
    for (int i = 0; i < numProc; i++) {
      countsBuf[i] = buffer(counts[i]);
    }
    myCounts = counts[rank];
    myCountsBuf = countsBuf[rank];

    this.histograms = histograms;
    this.rankToHistogramMap = rankToHistogramMap;
  }

  public int getHistogramIndex() {
    return rankToHistogramMap[rank];
  }

  /**
   * Send an OST count to all other processes while also receiving an OST count from all other
   * processes.
   *
   * @param lambda Current value of lambda.
   * @param dUdL Current value of dU/dL.
   * @param temperingWeight Current value of the temperingWeight.
   */
  public void send(double lambda, double dUdL, double temperingWeight) {
    // All-Gather counts from each walker.
    myCounts[0] = lambda;
    myCounts[1] = dUdL;
    myCounts[2] = temperingWeight;
    try {
      world.allGather(myCountsBuf, countsBuf);
    } catch (IOException ex) {
      String message = " Multi-walker OST allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }

    // Increment the Recursion Kernel(s) based on the input of each walker.
    for (int i = 0; i < numProc; i++) {

      int his = rankToHistogramMap[i];
      Histogram currentHistogram = histograms[his];

      // Only include this walkers bias.
      if (currentHistogram.getIndependentWalkers() && i != rank) {
        continue;
      }

      double walkerLambda = counts[i][0];
      double walkerdUdL = counts[i][1];
      double weight = counts[i][2];

      currentHistogram.setLastReceivedLambda(walkerLambda);
      currentHistogram.setLastReceiveddUdL(walkerdUdL);

      boolean resetStatistics = currentHistogram.getResetStatistics();
      double lambdaResetValue = currentHistogram.getLambdaResetValue();
      if (resetStatistics && walkerLambda > lambdaResetValue) {
        currentHistogram.allocateRecursionKernel();
        logger.info(format(" Cleared OST histogram (Lambda = %6.4f).", walkerLambda));
      }

      // For i == rank, the addBias method will handle updating FLambda (and optionally printing).
      currentHistogram.addToRecursionKernelValue(walkerLambda, walkerdUdL, weight, i != rank);
    }
  }

  /**
   * Update the synchronous communication histograms.
   *
   * @param histograms Histograms in use.
   * @param rankToHistogramMap Map from rank to histogram.
   */
  public void setHistograms(Histogram[] histograms, int[] rankToHistogramMap) {
    this.histograms = histograms;
    this.rankToHistogramMap = rankToHistogramMap;
  }

  /**
   * Update the map of rank-to-histogram.
   *
   * @param updatedRankToHisto Updated rank-to-histogram mappings.
   */
  public void updateRanks(int[] updatedRankToHisto) {
    assert updatedRankToHisto.length == rankToHistogramMap.length;
    System.arraycopy(updatedRankToHisto, 0, rankToHistogramMap, 0, rankToHistogramMap.length);
  }
}
