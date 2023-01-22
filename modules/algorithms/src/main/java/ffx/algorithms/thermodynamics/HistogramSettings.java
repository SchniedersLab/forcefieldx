// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import ffx.numerics.integrate.Integrate1DNumeric;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * HistogramSettings is a mutable settings class for OST histograms. Many fields in Histogram are
 * (or should be) final at construction, but there are too many to reasonably put in a constructor.
 * As such, HistogramSettings stores these values and their defaults, so a Histogram can have many
 * final fields with a small-signature constructor.
 *
 * <p>While most fields are package-private and can be set easily (technically breaking
 * encapsulation), some fields are fully private because they are involved in side-effects. For
 * example, dL may be adjusted slightly downwards if continuous lambda bins are in use and we need
 * to arrive at an odd number of lambda bins.
 */
public class HistogramSettings {
  private static final boolean DEFAULT_DISCRETE_LAMBDA = false;
  private static final double DEFAULT_DL = 0.005;
  private static final double DEFAULT_DFL = 2.0;
  private static final int DEFAULT_BIAS_CUTOFF = 5;
  private static final double DEFAULT_BIAS_MAG = 0.05;
  private static final double DEFAULT_TEMPERATURE = 298.15; // TODO: Consider setting this in Constants.
  private static final double DEFAULT_TEMPERING_FACTOR = 2.0;
  private static final int DEFAULT_FLAMDA_PRINT_INTERVAL = 25;
  private static final Integrate1DNumeric.IntegrationType DEFAULT_INTEGRATION =
      Integrate1DNumeric.IntegrationType.SIMPSONS;
  private static final double DEFAULT_TIMESTEP = 1.0;
  private static final double DEFAULT_THETA_FRICTION = 1.0e-19;
  private static final double DEFAULT_THETA_MASS = 1.0e-18;
  private static final int DEFAULT_COUNT_INTERVAL = 10;
  private static final double DEFAULT_LAMBDA_RESET = 0.99;
  private static final boolean DEFAULT_RESET_STATISTICS = false;
  private static final boolean DEFAULT_TEMPERING = true;
  /**
   * If true, use discrete lambda values instead of continuous lambda values. Is a final field to
   * avoid strange interactions with dL and number of lambda bins.
   */
  public final boolean discreteLambda;
  /**
   * Each walker reads the same histogram restart file. Only the walker of rank 0 writes the
   * histogram restart file.
   */
  private final File histogramFile;
  /**
   * Relative path to the histogram restart file. Assumption: a Histogram object will never change
   * its histogram or lambda files.
   */
  private final String hisFileName;
  /**
   * Relative path to the lambda restart file. Assumption: a Histogram object will never change its
   * histogram or lambda files.
   */
  private final String lambdaFileName;
  /**
   * Temperature in Kelvin.
   *
   * <p>The default is 298.15.
   */
  public double temperature = DEFAULT_TEMPERATURE;
  /**
   * The Dama et al. transition-tempering rate parameter. A reasonable value is about 2 to 8 kT,
   * with larger values being resulting in slower decay.
   *
   * <p>The default temperingFactor = 2.0.
   */
  public double temperingFactor = DEFAULT_TEMPERING_FACTOR;
  /** Time step in picoseconds. */
  public double dt = DEFAULT_TIMESTEP;
  /** Reasonable thetaFriction is ~60 per picosecond (1.0e-12). */
  public double thetaFriction = DEFAULT_THETA_FRICTION;
  /** Reasonable thetaMass is ~100 a.m.u. (100 a.m.u is 1.6605e-22 grams). */
  public double thetaMass = DEFAULT_THETA_MASS;
  /**
   * Interval between adding a count to the Recursion kernel in MD steps.
   *
   * <p>The default countInterval = 10.
   */
  public int countInterval = DEFAULT_COUNT_INTERVAL;
  /**
   * Flag to indicate if OST should send and receive counts between processes synchronously or
   * asynchronously. The latter can be faster by ~40% because simulation with Lambda &gt; 0.75 must
   * compute two condensed phase self-consistent fields to interpolate polarization.
   */
  public boolean asynchronous = true;

  double dFL;
  /**
   * When evaluating the biasing potential, contributions from Gaussians centered on bins more the
   * "biasCutoff" away will be neglected.
   *
   * <p>The default biasCutoff = 5.
   */
  int biasCutoff;
  /**
   * When evaluating the biasing potential, contributions from Gaussians centered on bins more than
   * "lambdaBiasCutoff" in the lambda dimension away will be neglected.
   *
   * <p>The default biasCutoff = 5 (continuous) or 0 (discrete).
   */
  int lambdaBiasCutoff;
  /**
   * Interval between how often the 1D histogram is printed to screen versus silently updated in
   * background.
   *
   * <p>The fLambdaPrintInterval is 25.
   */
  final int fLambdaPrintInterval = DEFAULT_FLAMDA_PRINT_INTERVAL;
  /** The integration algorithm used for thermodynamic integration. */
  Integrate1DNumeric.IntegrationType integrationType = DEFAULT_INTEGRATION;

  boolean tempering = DEFAULT_TEMPERING;
  /**
   * Width of standard lambda bins. Must be accessed by method because of an interaction with
   * discreteLambda; if that field is false, the Histogram needs an odd number of bins, which may
   * affect dL.
   */
  private double dL;
  /**
   * Magnitude of each hill (not including tempering). Must be set by method because the default
   * behavior for temperingThreshold is to be 20x this value.
   *
   * <p>The default biasMag = 0.05 (kcal/mol).
   */
  private double biasMag = DEFAULT_BIAS_MAG;

  private final double DEFAULT_BIAS_TO_OFFSET = 20.0;
  /**
   * An offset applied before calculating tempering weight.
   *
   * <p>First, for every Lambda, we calculate the maximum bias at that lambda by searching all
   * populated dU/dL bins: maxdUdL(L) = max[ G(L,F_L) ] where the max operator searches over all F_L
   * bins.
   *
   * <p>Then, the minimum bias coverage is determined by searching the maxdUdL(L) array over Lambda.
   * minBias = min[ maxdUdL(L) ] where the min operator searches all entries in the array.
   *
   * <p>Then the temperOffset is applied to the minBias:
   *
   * <p>biasHeight = max[minBias - temperOffset, 0]
   *
   * <p>The default temperOffset = 1.0 kcal/mol.
   *
   * <p>Must be set/accessed programmatically due to the default being 20x the Gaussian bias
   * magnitude.
   */
  private double temperOffset = biasMag * DEFAULT_BIAS_TO_OFFSET;
  // If set true (either by method or property setting temperOffset), no longer calculate based on
  // bias magnitude.
  private boolean temperOffsetSet = false;
  /** Once the lambda reset value is reached, OST statistics are reset. */
  private double lambdaResetValue = DEFAULT_LAMBDA_RESET;
  /** Flag set to false once OST statistics are reset at lambdaResetValue. */
  private boolean resetStatistics = DEFAULT_RESET_STATISTICS;

  private boolean metaDynamics = false;
  private boolean writeIndependent = false;
  private boolean independentWalkers = false;
  /** Flag indicating if a histogram file was read in. */
  final boolean histogramRead;

  public HistogramSettings(File hisFile, String lamFileName, CompositeConfiguration properties)
      throws IOException {
    this(
        hisFile,
        lamFileName,
        properties,
        properties.getBoolean("discrete-lambda", DEFAULT_DISCRETE_LAMBDA));
  }

  public HistogramSettings(
      File hisFile, String lamFileName, CompositeConfiguration properties, boolean discreteLambda)
      throws IOException {
    histogramFile = hisFile;
    hisFileName = hisFile.toString();
    this.lambdaFileName = lamFileName;

    biasCutoff = properties.getInt("lambda-bias-cutoff", DEFAULT_BIAS_CUTOFF);
    biasMag = properties.getDouble("bias-gaussian-mag", DEFAULT_BIAS_MAG);
    setDL(properties.getDouble("lambda-bin-width", DEFAULT_DL));
    dFL = properties.getDouble("flambda-bin-width", DEFAULT_DFL);
    this.discreteLambda = discreteLambda;
    // TODO: Strongly consider just eliminating the tempering flag, a relic of our earlier tempering scheme.
    tempering = properties.getBoolean("ost-alwaysTemper", DEFAULT_TEMPERING);
    if (properties.containsKey("ost-temperOffset")) {
      temperOffsetSet = true;
      temperOffset = properties.getDouble("ost-temperOffset");
    }

    if (histogramFile.exists()) {
      try (HistogramReader hr =
          new HistogramReader(null, new BufferedReader(new FileReader(histogramFile)))) {
        hr.readHistogramFile();
        temperature = hr.getTemperature();
        thetaMass = hr.getThetaMass();
        thetaFriction = hr.getThetaFriction();
        biasMag = hr.getBiasMag();
        biasCutoff = hr.getBiasCutoff();
        countInterval = hr.getCountInterval();
        setDL(hr.getLambdaBins());
        dFL = hr.getdUdLBinWidth();
        histogramRead = true;
      }
    } else {
      histogramRead = false;
    }
    if (properties.containsKey("lambda-bias-cutoff")) {
      lambdaBiasCutoff = properties.getInt("lambda-bias-cutoff");
    } else if (this.discreteLambda) {
      lambdaBiasCutoff = 0;
    } else {
      lambdaBiasCutoff = biasCutoff;
    }
  }

  // If setting a value can have side effects, or can be the side effect of another variable,
  // set it by method rather than raw access.

  /**
   * Gets the Gaussian bias magnitude in kcal/mol.
   *
   * @return Gaussian bias magnitude in kcal/mol.
   */
  public double getBiasMag() {
    return biasMag;
  }

  /**
   * Sets the bias magnitude, and if temperOffset has not previously been set by method or property,
   * sets temperOffset to 20x that value.
   *
   * @param biasMag Gaussian bias magnitude in kcal/mol
   */
  public void setBiasMag(double biasMag) {
    this.biasMag = biasMag;
    if (!temperOffsetSet) {
      // TODO: Logger.
      temperOffset = DEFAULT_BIAS_TO_OFFSET * biasMag;
    }
  }

  public double getDL() {
    continuousLambdaBins();
    return dL;
  }

  /**
   * Sets dLambda; if an invalid value is provided (not 0-1), resets it to default 0.005.
   *
   * @param dLambda Lambda bin width.
   */
  public void setDL(double dLambda) {
    if (dLambda < 0.0 || dLambda > 1.0) {
      dL = DEFAULT_DL;
    } else {
      dL = dLambda;
    }
  }

  /**
   * Sets dL based on a provided number of lambda bins.
   *
   * @param lambdaBins Number of lambda bins.
   */
  public void setDL(int lambdaBins) {
    assert lambdaBins > 2;
    dL = 1.0 / (lambdaBins - 1);
  }

  public String getHisFileName() {
    return hisFileName;
  }

  public File getHistogramFile() {
    return histogramFile;
  }

  public String getLambdaFileName() {
    return lambdaFileName;
  }

  public double getLambdaResetValue() {
    return lambdaResetValue;
  }

  public void setLambdaResetValue(double rsVal) {
    lambdaResetValue = rsVal;
    // TODO: Should this be a (0-1] check instead of [0-1]?
    resetStatistics = rsVal >= 0.0 && rsVal <= 1.0;
  }

  /**
   * Gets the tempering offset.
   *
   * @return Gets the tempering offset in kcal/mol.
   */
  public double getTemperOffset() {
    return temperOffset;
  }

  /**
   * Sets the tempering offset in kcal/mol.
   *
   * @param temperingOffset Tempering offset in kcal/mol.
   */
  public void setTemperOffset(double temperingOffset) {
    temperOffset = temperingOffset;
    temperOffsetSet = true;
  }

  /**
   * Returns the value of metaDynamics.
   *
   * @return metaDynamics.
   */
  public boolean isMetaDynamics() {
    return metaDynamics;
  }

  public void setMetaDynamics(boolean metaDynamics) {
    this.metaDynamics = metaDynamics;
  }

  /**
   * Returns the value of independentWalkers.
   *
   * @return independentWalkers.
   */
  public boolean independentWalkers() {
    return independentWalkers;
  }

  public boolean resetStatistics() {
    return resetStatistics;
  }

  /**
   * Sets the value of independentWalkers; if true, it also sets writeIndependent to true.
   *
   * @param iw Value to set independentWalkers to.
   */
  public void setIndependentWalkers(boolean iw) {
    independentWalkers = iw;
    if (iw) {
      setWriteIndependent(true);
    }
  }

  /**
   * Sets the value of writeIndependent.
   *
   * @param wi Value to set writeIndependent to.
   * @throws IllegalArgumentException If wi is false and independentWalkers is true.
   */
  public void setWriteIndependent(boolean wi) throws IllegalArgumentException {
    if (!wi && independentWalkers) {
      throw new IllegalArgumentException(" Independent walkers implies independent writing!");
    }
    writeIndependent = wi;
  }

  /**
   * Returns the value of writeIndependent.
   *
   * @return writeIndependent.
   */
  public boolean writeIndependent() {
    return writeIndependent;
  }

  /**
   * If continuous lambda bins are in use, rectify dL to have an odd number of bins. If (1.0 / dL)
   * would produce an odd number, dL is increased to make (1.0 / dL) the next even number.
   */
  private void continuousLambdaBins() {
    if (!discreteLambda) {
      // This is not really "lambda bins" so much as it's just 1.0 / dL.
      int fullWidthBins = (int) Math.round(1.0 / dL);
      if (fullWidthBins % 2 == 0) {
        ++fullWidthBins;
      }
      dL = 1.0 / (fullWidthBins - 1);
    }
  }
}
