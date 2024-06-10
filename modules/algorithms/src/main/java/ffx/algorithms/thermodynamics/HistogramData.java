// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import ffx.utilities.FileUtils;
import ffx.utilities.HistogramXmlAdapter;
import jakarta.xml.bind.JAXBContext;
import jakarta.xml.bind.Marshaller;
import jakarta.xml.bind.Unmarshaller;
import jakarta.xml.bind.annotation.XmlAccessOrder;
import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorOrder;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlElement;
import jakarta.xml.bind.annotation.XmlRootElement;
import jakarta.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.io.File;
import java.util.logging.Logger;

import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

@XmlRootElement(name = "HistogramData")
@XmlAccessorOrder(XmlAccessOrder.ALPHABETICAL)
@XmlAccessorType(XmlAccessType.NONE)
public class HistogramData {

  private static final Logger logger = Logger.getLogger(HistogramData.class.getName());

  /**
   * The default Gaussian bias magnitude is 0.05 kcal/mol.
   */
  private static final double DEFAULT_BIAS_MAG = 0.05;

  /**
   * The default dU/dL bias cutoff is 5 bins.
   */
  private static final int DEFAULT_DUDL_BIAS_CUTOFF = 5;

  /**
   * The default lambda bias cutoff is 5 bins for continuous lambda and 0 bins for discrete lambda.
   */
  private static final int DEFAULT_LAMBDA_BIAS_CUTOFF = 5;

  /**
   * The default count interval is 10 MD steps.
   */
  private static final int DEFAULT_COUNT_INTERVAL = 10;

  /**
   * The default lambda bin width is 0.005.
   * There are (1.0 / width - 1) full size bins (e.g., 199 for width = 0.005).
   * The first and last bins are half size (e.g., 0.0025 for width = 0.005).
   */
  private static final double DEFAULT_LAMBDA_BIN_WIDTH = 0.005;

  /**
   * The default number of lambda bins is 201 based on DEFAULT_LAMBDA_BIN_WIDTH = 0.005.
   * The first and last bin are of size DEFAULT_LAMBDA_BIN_WIDTH / 2.
   * Bin 0:   0.0 to DEFAULT_LAMBDA_BIN_WIDTH / 2.
   * Bin 200: 1.0 - DEFAULT_LAMBDA_BIN_WIDTH / 2 to 1.0.
   * Bins 1 to 199 have size DEFAULT_LAMBDA_BIN_WIDTH.
   */
  private static final int DEFAULT_LAMBDA_BINS = 201;

  /**
   * The default dU/dL bin width is 2.0 kcal/mol.
   */
  private static final double DEFAULT_DUDL_BIN_WIDTH = 2.0;

  /**
   * Default number of dU/dL bins is 101.
   */
  private static final int DEFAULT_DUDL_BINS = 101;

  /**
   * The default dU/dL range is from -DEFAULT_DUDL_BINS to DEFAULT_DUDL_BINS.
   * This depends on the 1) default bin width of 2.0 kcal/mol and 2) the central bin range of -1.0 to 1.0.
   * The default dU/dL minimum is -101 kcal/mol and the maximum is 101 kcal/mol.
   */
  private static final double DEFAULT_DUDL_MINIMUM = -DEFAULT_DUDL_BINS;

  /**
   * Default is that meta-dynamics is not used.
   */
  private static final boolean DEFAULT_METADYNAMICS = false;

  /**
   * The default tempering factor is 2.0 (in units of kBT).
   */
  private static final double DEFAULT_TEMPERING_FACTOR = 2.0;

  /**
   * The default lambda reset value is at lambda = 0.99.
   */
  private static final double DEFAULT_RESET_HISTOGRAM_AT_LAMBDA = 0.99;

  /**
   * The default reset statistics flag is false.
   */
  private static final boolean DEFAULT_RESET_HISTOGRAM = false;

  /**
   * The default tempering offset is DEFAULT_BIAS_TO_TEMPERING_OFFSET * DEFAULT_BIAS_MAG = 1.0 kcal/mol.
   */
  private final double DEFAULT_BIAS_TO_TEMPERING_OFFSET = 20.0;

  /**
   * By default, each walker contributes to the same histogram.
   */
  private final boolean DEFAULT_INDEPENDENT_WALKERS = false;

  /**
   * By default, only the Rank 0 process writes the histogram.
   * <p>
   * Note that if independent walkers are used, each walker will write its own histogram.
   */
  private final boolean DEFAULT_INDEPENDENT_WRITE = false;

  /**
   * By default, the lambda state variable is continuous. If discreteLambda is true, the lambda state is discrete.
   */
  private static final boolean DEFAULT_DISCRETE_LAMBDA = false;

  /**
   * Magnitude of each hill (not including tempering). The default biasMag = 0.05 (kcal/mol).
   */
  @XmlElement(name = "BiasMag", defaultValue = "" + DEFAULT_BIAS_MAG)
  double biasMag = DEFAULT_BIAS_MAG;

  /**
   * When evaluating the biasing potential, contributions from a Gaussian centered on a bin more
   * than lambdaBiasCutoff away will be neglected.
   * The continuous lambda simulations, the default lambdaBiasCutoff = 5.
   * For discrete lambda simulations, the default lambdaBiasCutoff = 0.
   */
  @XmlElement(name = "LambdaBiasCutoff", defaultValue = "" + DEFAULT_LAMBDA_BIAS_CUTOFF)
  int lambdaBiasCutoff = DEFAULT_LAMBDA_BIAS_CUTOFF;

  /**
   * When evaluating the biasing potential, contributions from a Gaussian centered on a bin more
   * than dUdLBiasCutoff" away will be neglected. The default dUdLBiasCutoff = 5.
   */
  @XmlElement(name = "dUdLBiasCutoff", defaultValue = "" + DEFAULT_DUDL_BIAS_CUTOFF)
  int dUdLBiasCutoff = DEFAULT_DUDL_BIAS_CUTOFF;

  /**
   * Interval between adding a count to the Recursion kernel in MD steps. The default countInterval = 10.
   */
  @XmlElement(name = "CountInterval", defaultValue = "" + DEFAULT_COUNT_INTERVAL)
  int countInterval = DEFAULT_COUNT_INTERVAL;

  /**
   * Width of a lambda bin, or the distance between discrete lambda values.
   * The default is 1.0 / (lambdaBins - 1).
   */
  @XmlElement(name = "LambdaBinWidth", defaultValue = "" + DEFAULT_LAMBDA_BIN_WIDTH)
  double lambdaBinWidth = DEFAULT_LAMBDA_BIN_WIDTH;

  /**
   * If true, use discrete lambda values instead of continuous lambda values.
   */
  @XmlElement(name = "DiscreteLambda", defaultValue = "" + DEFAULT_DISCRETE_LAMBDA)
  boolean discreteLambda = DEFAULT_DISCRETE_LAMBDA;

  /**
   * The width of the FLambda bin. The default is 2.0 (kcal/mol).
   */
  @XmlElement(name = "dUdLBinWidth", defaultValue = "" + DEFAULT_DUDL_BIN_WIDTH)
  double dUdLBinWidth = DEFAULT_DUDL_BIN_WIDTH;

  /**
   * The minimum value of the first dU/dL bin.
   * Initially this is set to: mindUdL = -(dUdLBinWidth * dUdLBins) / 2.0.
   * However, mindUdL may be updated as the count matrix is dynamically resized.
   */
  @XmlElement(name = "dUdLMinimum", defaultValue = "" + DEFAULT_DUDL_MINIMUM)
  double dUdLMinimum = DEFAULT_DUDL_MINIMUM;

  /**
   * Flag to indicate use of 1D meta-dynamics instead of OST.
   */
  @XmlElement(name = "MetaDynamics", defaultValue = "" + DEFAULT_METADYNAMICS)
  boolean metaDynamics = DEFAULT_METADYNAMICS;

  /**
   * The Dama et al. transition-tempering rate parameter. A reasonable value is about 2 to 8 kT,
   * with larger values being resulting in slower decay. The default temperingFactor = 2.0.
   */
  @XmlElement(name = "TemperingFactor", defaultValue = "" + DEFAULT_TEMPERING_FACTOR)
  double temperingFactor = DEFAULT_TEMPERING_FACTOR;

  /**
   * An offset applied before calculating tempering weight.
   * <p>
   * First, for every Lambda, we calculate the maximum bias at that lambda by searching all
   * populated dU/dL bins: maxdUdL(L) = max[ G(L,F_L) ] where the max operator searches over all F_L
   * bins.
   * <p>
   * Then, the minimum bias coverage is determined by searching the maxdUdL(L) array over Lambda.
   * minBias = min[ maxdUdL(L) ] where the min operator searches all entries in the array.
   * <p>
   * Then the temperOffset is applied to the minBias: biasHeight = max[minBias - temperingOffset, 0]
   * <p>
   * The default temperingOffset = 20x the Gaussian bias magnitude.
   */
  @XmlElement(name = "TemperingOffset", defaultValue = "-1.0")
  private double temperingOffset = -1.0;

  /**
   * Once the lambda reset value is reached, OST statistics are reset.
   */
  @XmlElement(name = "ResetHistogramAtLambda", defaultValue = "" + DEFAULT_RESET_HISTOGRAM_AT_LAMBDA)
  double resetHistogramAtLambda = DEFAULT_RESET_HISTOGRAM_AT_LAMBDA;

  /**
   * Flag set to false once OST statistics are reset at lambdaResetValue.
   */
  @XmlElement(name = "ResetHistogram", defaultValue = "" + DEFAULT_RESET_HISTOGRAM)
  boolean resetHistogram = DEFAULT_RESET_HISTOGRAM;

  @XmlElement(name = "IndependentWalkers", defaultValue = "" + DEFAULT_INDEPENDENT_WALKERS)
  boolean independentWalkers = DEFAULT_INDEPENDENT_WALKERS;

  @XmlElement(name = "IndependentWrite", defaultValue = "" + DEFAULT_INDEPENDENT_WRITE)
  boolean independentWrite = DEFAULT_INDEPENDENT_WRITE;

  /**
   * Flag to indicate if OST should send and receive counts between processes synchronously or
   * asynchronously. The latter can be faster by ~40% because simulation with Lambda &gt; 0.75 must
   * compute two condensed phase self-consistent fields to interpolate polarization.
   */
  @XmlElement(name = "Asynchronous", defaultValue = "" + true)
  boolean asynchronous = true;

  /**
   * It is useful to have an odd number of bins, so that there is a bin from FL=-dFL/2 to dFL/2 so
   * that as FL approaches zero its contribution to thermodynamic integration goes to zero.
   * <p>
   * Otherwise a contribution of zero from a L bin can only result from equal sampling of the
   * ranges -dFL to 0 and 0 to dFL.
   * <p>
   * The default FLambdaBins = 101.
   */
  @XmlElement(name = "dUdLBins", defaultValue = "" + DEFAULT_DUDL_BINS)
  public int dUdLBins = DEFAULT_DUDL_BINS;

  /**
   * The number of hills added to the recursion kernel.
   */
  @XmlElement(name = "Counts", defaultValue = "0")
  public int counts = 0;

  /**
   * The recursion kernel stores the weight of each [lambda][dU/dL] bin.
   * Note that the variable name begins with Z so that its last in the histogram XML file.
   */
  @XmlElement(name = "Data", required = true)
  @XmlJavaTypeAdapter(HistogramXmlAdapter.class)
  double[][] zHistogram = new double[DEFAULT_LAMBDA_BINS][DEFAULT_DUDL_BINS];

  /**
   * The histogram file to write to / read from.
   */
  private File histogramFile = null;

  /**
   * The name of the histogram file.
   */
  private String histogramFileName = null;

  /**
   * Flag indicating if this histogram data came from reading a file.
   */
  private boolean histogramRead = false;

  /**
   * For continuous lambda: The first Lambda bin is centered on 0.0 (-0.0025 to 0.0025). The final
   * Lambda bin is centered on 1.0 ( 0.9975 to 1.0025).
   * <p>
   * With this scheme, the maximum of biasing Gaussian hills is at exactly 0.0 and 1.0.
   * <p>
   * For discrete lambda: The first value of lambda is 0.0 and last value is 1.0.
   * <p>
   * The default lambdaBins = 201.
   */
  public int lambdaBins = DEFAULT_LAMBDA_BINS;

  /**
   * Half the width of a lambda bin, or zero for discrete lambda values.
   */
  public double lambdaBinWidth_2 = DEFAULT_LAMBDA_BIN_WIDTH / 2.0;

  /**
   * The minimum value of the first lambda bin.
   * <p>
   * minLambda = -dL_2 for continuous lambda.
   * <p>
   * minLambda = 0 for discrete lambda.
   */
  public double minLambda = -DEFAULT_LAMBDA_BIN_WIDTH / 2.0;

  /**
   * Either the discrete lambda values used, or null (continuous lambda).
   */
  public double[] lambdaLadder = null;

  /**
   * The variance for the Gaussian bias in the lambda dimension.
   * lambdaVariance = 2.0 * dL * 2.0 * dL;
   */
  public double lambdaVariance = pow(2.0 * DEFAULT_LAMBDA_BIN_WIDTH, 2);


  /**
   * The maximum value of the last dUdL bin.
   * <p>
   * maxdUdL = mindUdL + dUdLBins * dUdLBinWidth.
   * <p>
   * default = -101 + 101 * 2.0 = 101 = DEFAULT_DUDL_BINS.
   */
  public double dUdLMaximum = DEFAULT_DUDL_BINS;

  /**
   * Half the width of the F_lambda bin.
   */
  public double dUdLBinWidth_2 = DEFAULT_DUDL_BIN_WIDTH / 2.0;

  /**
   * The variance for the Gaussian bias in the dU/dL dimension.
   * dUdLVariance = 2.0 * dFL * 2.0 * dFL;
   */
  public double dUdLVariance = pow(2.0 * dUdLBinWidth, 2);

  /**
   * Marshall the histogram data to a file.
   *
   * @param histogramData The Histogram data.
   * @param file          The file to write to.
   */
  public static void writeHistogram(HistogramData histogramData, File file) {
    try {
      JAXBContext jaxbContext = JAXBContext.newInstance(HistogramData.class);
      Marshaller jaxbMarshaller = jaxbContext.createMarshaller();
      jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
      jaxbMarshaller.marshal(histogramData, file);
    } catch (Exception e) {
      logger.warning(" Exception writing histogram:\n " + e);
    }
  }

  /**
   * Unmarshall the histogram data.
   *
   * @param file The file to read from.
   * @return The Histogram data.
   */
  public static HistogramData readHistogram(@Nullable File file) {
    if (file != null && file.exists() && file.canRead()) {
      try {
        JAXBContext jaxbContext = JAXBContext.newInstance(HistogramData.class);
        Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
        HistogramData histogramData = (HistogramData) unmarshaller.unmarshal(file);
        histogramData.setHistogramFile(file);
        histogramData.setHistogramRead(true);
        histogramData.rectify();
        return histogramData;
      } catch (Exception e) {
        logger.warning(" Exception reading histogram:\n " + e);
      }
    }
    HistogramData histogramData = new HistogramData();
    histogramData.setHistogramFile(file);
    histogramData.rectify();
    return histogramData;
  }

  private void rectify() {
    // Check that histogram size matches the number of lambda and dU/dL bins.
    int hisLambdaBins = zHistogram.length;
    if (lambdaBins != hisLambdaBins) {
      logger.fine(" Current LambdaBins: " + lambdaBins);
      logger.fine(" Updated to:         " + hisLambdaBins);
      lambdaBins = hisLambdaBins;
      lambdaBinWidth = 1.0 / (lambdaBins - 1);
    }
    int hisDUDLBins = zHistogram[0].length;
    if (dUdLBins != hisDUDLBins) {
      logger.fine(" Current dUdLBins:   " + dUdLBins);
      logger.fine(" Updated to:         " + hisDUDLBins);
      dUdLBins = hisDUDLBins;
    }
    rectifyLambdaVariables();
    rectifyDUDLVariables();
  }

  /**
   * All transient variables that depend on dUdLBinWidth are set here.
   */
  private void rectifyDUDLVariables() {
    // The center of the central bin is at 0.
    dUdLBinWidth_2 = dUdLBinWidth / 2.0;
    dUdLMaximum = dUdLMinimum + (dUdLBins * dUdLBinWidth);
    dUdLVariance = pow(2.0 * dUdLBinWidth, 2);
  }

  /**
   * All transient variables that depend on lambdaBinWidth or discreteLambda are set here.
   */
  private void rectifyLambdaVariables() {
    // Get the closest number of full width bins based on the current lambdaBinWidth.
    int n = (int) round(1.0 / lambdaBinWidth);
    // Now set the lambdaBinWidth exactly.
    lambdaBinWidth = 1.0 / n;
    // Divide a bin into two half size bins (the first and last bins are centered on 0 and 1, respectively).
    lambdaBins = n + 1;
    if (discreteLambda) {
      lambdaLadder = new double[lambdaBins];
      lambdaLadder[0] = 0.0;
      lambdaLadder[lambdaBins - 1] = 1.0;
      for (int i = 1; i < lambdaBins - 1; i++) {
        lambdaLadder[i] = lambdaBinWidth * i;
      }
      lambdaBinWidth_2 = 0.0;
      minLambda = 0.0;
      lambdaBiasCutoff = 0;
    } else {
      lambdaLadder = null;
      lambdaBinWidth_2 = lambdaBinWidth / 2.0;
      minLambda = -lambdaBinWidth_2;
      // Leave the lambda bias cutoff at it current value, unless its zero.
      if (lambdaBiasCutoff == 0) {
        // In this case, set it back to the default value.
        lambdaBiasCutoff = DEFAULT_LAMBDA_BIAS_CUTOFF;
      }
    }
    lambdaVariance = pow(2.0 * lambdaBinWidth, 2);
  }

  public void setCountInterval(int countInterval) {
    this.countInterval = countInterval;
  }

  /**
   * Marshall this histogram data to a file.
   */
  public void writeHistogram() {
    writeHistogram(this, histogramFile);
  }

  public String getHistogramFileName() {
    return histogramFileName;
  }

  public File getHistogramFile() {
    return histogramFile;
  }

  public void setHistogramFile(@Nullable File histogramFile) {
    this.histogramFile = histogramFile;
    if (histogramFile != null) {
      histogramFileName = FileUtils.relativePathTo(histogramFile).toString();
    }
  }

  public boolean wasHistogramRead() {
    return histogramRead;
  }

  public void setHistogramRead(boolean histogramRead) {
    this.histogramRead = histogramRead;
  }

  /**
   * Sets the value of independentWalkers; if true, it also sets writeIndependent to true.
   *
   * @param independentWalkers Value to set independentWalkers to.
   */
  public void setIndependentWalkers(boolean independentWalkers) {
    this.independentWalkers = independentWalkers;
    if (independentWalkers) {
      setWriteIndependent(true);
    }
  }

  /**
   * Sets the value of writeIndependent.
   *
   * @param independentWrite Value to set writeIndependent to.
   * @throws IllegalArgumentException If independentWrite is false and independentWalkers is true.
   */
  public void setWriteIndependent(boolean independentWrite) throws IllegalArgumentException {
    if (!independentWrite && independentWalkers) {
      throw new IllegalArgumentException(" Independent walkers implies independent writing.");
    }
    this.independentWrite = independentWrite;
  }

  /**
   * Returns the value of independentWrite.
   *
   * @return independentWrite.
   */
  public boolean independentWrite() {
    return independentWrite;
  }

  /**
   * Gets the Gaussian bias magnitude in kcal/mol.
   *
   * @return Gaussian bias magnitude in kcal/mol.
   */
  public double getBiasMag() {
    return biasMag;
  }

  /**
   * Set the bias magnitude. If temperingOffset has not been explicitly set,
   * the temperOffset is set to bias magnitude * 20.
   *
   * @param biasMag Gaussian bias magnitude in kcal/mol
   */
  public void setBiasMag(double biasMag) {
    this.biasMag = biasMag;
  }

  public void setTemperingFactor(double temperingFactor) {
    this.temperingFactor = temperingFactor;
  }

  public void setMetaDynamics(boolean metaDynamics) {
    this.metaDynamics = metaDynamics;
  }

  public double getLambdaBinWidth() {
    return lambdaBinWidth;
  }

  public int getLambdaBins() {
    return lambdaBins;
  }

  public int getDUDLBins() {
    return dUdLBins;
  }

  /**
   * Sets lambdaBinWidth; if an invalid value is provided (not 0-1), resets it to default 0.005.
   *
   * @param lambdaBinWidth Lambda bin width.
   */
  public void setLambdaBinWidth(double lambdaBinWidth) {
    if (lambdaBinWidth <= 0.0 || lambdaBinWidth > 1.0) {
      this.lambdaBinWidth = DEFAULT_LAMBDA_BIN_WIDTH;
    } else {
      this.lambdaBinWidth = lambdaBinWidth;
    }
    rectifyLambdaVariables();
  }

  public void setdUdLBinWidth(double dUdLBinWidth) {
    // The scale the minimum value of the first dU/dL bin.
    double ratio = dUdLBinWidth / this.dUdLBinWidth;
    this.dUdLMinimum *= ratio;
    this.dUdLBinWidth = dUdLBinWidth;
    rectifyDUDLVariables();
  }

  /**
   * Gets the tempering offset.
   *
   * @return Gets the tempering offset in kcal/mol.
   */
  public double getTemperingOffset() {
    if (temperingOffset >= 0.0) {
      return temperingOffset;
    }
    return biasMag * DEFAULT_BIAS_TO_TEMPERING_OFFSET;
  }

  /**
   * Sets an explicit tempering offset in kcal/mol.
   *
   * @param temperingOffset Tempering offset in kcal/mol.
   */
  public void setTemperingOffset(double temperingOffset) {
    this.temperingOffset = temperingOffset;
  }

  public void setAsynchronous(boolean asynchronous) {
    this.asynchronous = asynchronous;
  }

  public double resetHistogramAtLambda() {
    return resetHistogramAtLambda;
  }

  /**
   * Set the lambda value at which to reset the histogram.
   *
   * @param resetHistogramAtLambda The lambda value at which to reset the histogram.
   */
  public void setResetHistogramAtLambda(double resetHistogramAtLambda) {
    if (resetHistogramAtLambda < 0.0 || resetHistogramAtLambda > 1.0) {
      this.resetHistogramAtLambda = DEFAULT_RESET_HISTOGRAM_AT_LAMBDA;
      this.resetHistogram = false;
    } else {
      this.resetHistogramAtLambda = resetHistogramAtLambda;
      this.resetHistogram = true;
    }
  }

  public void setDiscreteLambda(boolean discreteLambda) {
    this.discreteLambda = discreteLambda;
    rectifyLambdaVariables();
  }

  public void applyProperties(CompositeConfiguration properties) {
    if (properties.containsKey("lambda-bias-cutoff")) {
      lambdaBiasCutoff = properties.getInt("lambda-bias-cutoff");
    }
    if (properties.containsKey("dudl-bias-cutoff")) {
      dUdLBiasCutoff = properties.getInt("dudl-bias-cutoff");
    }
    if (properties.containsKey("ost-bias-mag")) {
      double bias = properties.getDouble("ost-bias-mag");
      setBiasMag(bias);
    }
    if (properties.containsKey("ost-temperOffset")) {
      double offset = properties.getDouble("ost-temperOffset");
      setTemperingOffset(offset);
    }
    if (properties.containsKey("lambda-bin-width")) {
      double width = properties.getDouble("lambda-bin-width");
      setLambdaBinWidth(width);
    }
    if (properties.containsKey("flambda-bin-width")) {
      double dUdLBinWidth = properties.getDouble("flambda-bin-width");
      setdUdLBinWidth(dUdLBinWidth);
    }
    if (properties.containsKey("discrete-lambda")) {
      boolean discreteLambda = properties.getBoolean("discrete-lambda");
      setDiscreteLambda(discreteLambda);
    }
  }

  public String toString() {
    return format("  Lambda bins:      %6d\n", lambdaBins)
        + format("  Lambda bin width: %6.3f\n", lambdaBinWidth)
        + format("  dU/dL bins:       %6d\n", dUdLBins)
        + format("  dU/dL bin width:  %6.3f (kcal/mol)\n", dUdLBinWidth)
        + format("  Bias magnitude:   %6.3f (kcal/mol)\n", biasMag)
        + format("  Tempering offset: %6.3f (kcal/mol)\n", getTemperingOffset())
        + format("  Tempering rate:   %6.3f\n", temperingFactor)
        + format("  Discrete lambda:  %6B\n", discreteLambda)
        + format("  Meta-dynamics:    %6B\n\n", metaDynamics);
  }

}
