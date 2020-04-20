package ffx.algorithms.thermodynamics;

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
public class SynchronousSend {

  private static final Logger logger = Logger.getLogger(SynchronousSend.class.getName());

  /** Parallel Java world communicator. */
  protected final Comm world;
  /** Rank of this process. */
  protected final int rank;
  /** Number of processes. */
  private final int numProc;
  /**
   * The recursionWeights stores the [Lambda, FLambda] weight for each process. Therefore the array
   * is of size [number of Processes][2].
   *
   * <p>Each 2 entry array must be wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final double[][] recursionWeights;

  private final double[] myRecursionWeight;
  /** These DoubleBufs wrap the recursionWeight arrays. */
  private final DoubleBuf[] recursionWeightsBuf;

  private final DoubleBuf myRecursionWeightBuf;
  private boolean independentWalkers;
  /** The histograms to update. */
  private Histogram[] histograms;
  /** Map from ranks to histograms */
  private int[] rankToHistogramMap;

  /**
   * Constructor.
   *
   * @param histograms An array of Bias Histograms.
   * @param rankToHistogramMap A map from process rank to Histogram.
   * @param independentWalkers If true, each walker has its own Histogram.
   */
  public SynchronousSend(
      Histogram[] histograms, int[] rankToHistogramMap, boolean independentWalkers) {
    world = Comm.world();
    numProc = world.size();
    rank = world.rank();
    // Use synchronous communication.
    recursionWeights = new double[numProc][3];
    recursionWeightsBuf = new DoubleBuf[numProc];
    for (int i = 0; i < numProc; i++) {
      recursionWeightsBuf[i] = DoubleBuf.buffer(recursionWeights[i]);
    }
    myRecursionWeight = recursionWeights[rank];
    myRecursionWeightBuf = recursionWeightsBuf[rank];

    this.histograms = histograms;
    this.rankToHistogramMap = rankToHistogramMap;
    this.independentWalkers = independentWalkers;
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
    myRecursionWeight[0] = lambda;
    myRecursionWeight[1] = dUdL;
    myRecursionWeight[2] = temperingWeight;

    try {
      world.allGather(myRecursionWeightBuf, recursionWeightsBuf);
    } catch (IOException ex) {
      String message = " Multi-walker OST allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }

    // Find the minimum and maximum FLambda bin for the gathered counts.
    for (int i = 0; i < numProc; i++) {
      // Only include this walkers bias.
      if (independentWalkers && i != rank) {
        continue;
      }
      int his = rankToHistogramMap[i];
      histograms[his].checkRecursionKernelSize(recursionWeights[i][1]);
    }

    // Increment the Recursion Kernel(s) based on the input of each walker.
    for (int i = 0; i < numProc; i++) {

      // Only include this walkers bias.
      if (independentWalkers && i != rank) {
        continue;
      }

      int his = rankToHistogramMap[i];
      Histogram currentHistogram = histograms[his];

      currentHistogram.setLastReceivedLambda(recursionWeights[i][0]);
      currentHistogram.setLastReceiveddUdL(recursionWeights[i][1]);

      int walkerLambda = currentHistogram.indexForLambda(recursionWeights[i][0]);
      int walkerFLambda = currentHistogram.binForFLambda(recursionWeights[i][1]);
      double weight = recursionWeights[i][2];

      boolean resetStatistics = currentHistogram.getResetStatistics();
      double lambdaResetValue = currentHistogram.getLambdaResetValue();
      if (resetStatistics && recursionWeights[i][0] > lambdaResetValue) {
        currentHistogram.allocateRecursionKernel();
        logger.info(format(" Cleared OST histogram (Lambda = %6.4f).", recursionWeights[i][0]));
      }

      // For i == rank, the addBias method will handle updating FLambda (and optionally printing).
      currentHistogram.addToRecursionKernelValue(walkerLambda, walkerFLambda, weight, i != rank);
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

  public void setIndependentWalkers(boolean independentWalkers) {
    this.independentWalkers = independentWalkers;
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
