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

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.cli.LambdaParticleOptions;
import ffx.algorithms.dynamics.Barostat;
import ffx.algorithms.dynamics.integrators.Stochastic;
import ffx.algorithms.optimize.Minimize;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.numerics.integrate.DataSet;
import ffx.numerics.integrate.DoublesDataSet;
import ffx.numerics.integrate.Integrate1DNumeric;
import ffx.numerics.integrate.Integrate1DNumeric.IntegrationType;
import ffx.potential.MolecularAssembly;
import ffx.potential.SystemState;
import ffx.potential.Utilities;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import javax.annotation.Nullable;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.integrate.Integrate1DNumeric.IntegrationType.SIMPSONS;
import static ffx.utilities.Constants.R;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * An implementation of the Orthogonal Space Tempering algorithm.
 *
 * <p>This only partially implements the LambdaInterface, since it does not return 2nd lambda
 * derivatives. The 2nd derivatives of the bias require 3rd derivatives of the underlying Hamiltonian
 * (not these derivatives are not needed for OST MD).
 *
 * @author Michael J. Schnieders, James Dama, Wei Yang and Pengyu Ren
 * @since 1.0
 */
public class OrthogonalSpaceTempering implements CrystalPotential, LambdaInterface {

  private static final Logger logger = Logger.getLogger(OrthogonalSpaceTempering.class.getName());
  /**
   * The potential energy of the system.
   */
  protected final CrystalPotential potential;
  /**
   * Reference to the Barostat in use; if present this must be turned off during optimization.
   */
  protected final Barostat barostat;
  /**
   * The AlgorithmListener is called each time a count is added.
   */
  protected final AlgorithmListener algorithmListener;
  /**
   * Print detailed energy information.
   */
  protected final boolean print = false;
  /**
   * Number of variables.
   */
  protected final int nVariables;
  /**
   * A potential energy that implements the LambdaInterface.
   */
  private final LambdaInterface lambdaInterface;
  /**
   * List of additional Histograms this OST can switch to.
   */
  private final List<Histogram> allHistograms = new ArrayList<>();
  /**
   * Parameters to control saving local optimizations.
   */
  private final OptimizationParameters optimizationParameters;
  /**
   * Properties.
   */
  private final CompositeConfiguration properties;
  /**
   * The MolecularAssembly being simulated.
   */
  protected MolecularAssembly molecularAssembly;
  /**
   * Are FAST varying energy terms being computed, SLOW varying energy terms, or BOTH. OST is not
   * active when only FAST varying energy terms are being propagated.
   */
  protected Potential.STATE state = Potential.STATE.BOTH;
  /**
   * Force Field Potential Energy (i.e. with no bias terms added).
   */
  protected double forceFieldEnergy;
  /**
   * Contains counts for the OST bias.
   */
  private Histogram histogram;
  /**
   * Index of the current Histogram.
   */
  private int histogramIndex;
  /**
   * Flag to indicate that the Lambda particle should be propagated.
   */
  private boolean propagateLambda = true;
  /**
   * Mixed second partial derivative with respect to coordinates and lambda.
   */
  private final double[] dUdXdL;

  /**
   * Gradient array needed when the OST Energy method is called.
   */
  private final double[] tempGradient;

  /**
   * Partial derivative of the force field energy with respect to lambda.
   */
  private double dForceFieldEnergydL;
  /**
   * Magnitude of the 2D orthogonal space bias G(L,dE/dL).
   */
  private double gLdEdL = 0.0;
  /**
   * OST Bias energy.
   */
  private double biasEnergy;
  /**
   * Total system energy.
   */
  private double totalEnergy;
  /**
   * Total partial derivative of the potential (U) being sampled with respect to lambda.
   */
  private double dUdLambda;
  /**
   * Second partial derivative of the potential being sampled with respect to lambda.
   */
  private double d2UdL2;
  /**
   * If true, values of (lambda, dU/dL) that have not been observed are rejected.
   */
  private boolean hardWallConstraint = false;

  private final DynamicsOptions dynamicsOptions;

  private final LambdaParticleOptions lambdaParticleOptions;

  /**
   * OST Constructor.
   *
   * @param lambdaInterface       defines Lambda and dU/dL.
   * @param potential             defines the Potential energy.
   * @param histogramData         contains histogram restart data.
   * @param lambdaData            contains lambda restart data.
   * @param properties            defines System properties.
   * @param dynamicsOptions       defines molecular dynamics parameters.
   * @param lambdaParticleOptions defines lambda particle parameters.
   * @param algorithmListener     the AlgorithmListener to be notified of progress.
   */
  public OrthogonalSpaceTempering(LambdaInterface lambdaInterface, CrystalPotential potential,
                                  HistogramData histogramData, LambdaData lambdaData, CompositeConfiguration properties,
                                  DynamicsOptions dynamicsOptions, LambdaParticleOptions lambdaParticleOptions,
                                  AlgorithmListener algorithmListener) {
    this.lambdaInterface = lambdaInterface;
    this.potential = potential;
    this.properties = properties;
    this.dynamicsOptions = dynamicsOptions;
    this.lambdaParticleOptions = lambdaParticleOptions;
    this.algorithmListener = algorithmListener;
    nVariables = potential.getNumberOfVariables();

    if (potential instanceof Barostat) {
      barostat = (Barostat) potential;
    } else {
      barostat = null;
    }

    dUdXdL = new double[nVariables];
    tempGradient = new double[nVariables];

    // Init the Histogram.
    histogram = new Histogram(properties, histogramData, lambdaData);
    histogramIndex = 0;
    allHistograms.add(histogram);

    // Configure optimization parameters.
    optimizationParameters = new OptimizationParameters(properties);
  }

  /**
   * Add an alternate Histogram this OST can use.
   *
   * @param histogramData Settings to use for the new Histogram.
   */
  public void addHistogram(HistogramData histogramData, LambdaData lambdaData) {
    Histogram newHisto = new Histogram(properties, histogramData, lambdaData);
    histogramData.asynchronous = this.histogram.hd.asynchronous;
    allHistograms.add(newHisto);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean dEdLZeroAtEnds() {
    return false;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    // Shut down the CountReceiveThread.
    histogram.destroy();
    return potential.destroy();
  }

  /**
   * Compute the force field + bias energy.
   */
  public double energy(double[] x) {

    // OST is propagated with the slowly varying terms.
    if (state == Potential.STATE.FAST) {
      forceFieldEnergy = potential.energy(x);
      return forceFieldEnergy;
    }

    // Stochastic dynamics pre-force propagation.
    if (propagateLambda) {
      histogram.stochasticPreForce();
    }

    // We have to compute the energy and gradient, since we need dU/dL to be computed.
    fill(tempGradient, 0.0);
    forceFieldEnergy = potential.energyAndGradient(x, tempGradient);

    dForceFieldEnergydL = lambdaInterface.getdEdL();
    d2UdL2 = lambdaInterface.getd2EdL2();
    int lambdaBin = histogram.indexForLambda();
    dUdLambda = dForceFieldEnergydL;

    gLdEdL = 0.0;
    double bias1D;
    // Update the free energy difference.
    histogram.updateFreeEnergyDifference(false, false);
    if (histogram.hd.metaDynamics) {
      bias1D = histogram.energyAndGradientMeta(true);
    } else {
      // Calculate recursion kernel G(L, F_L) and its derivatives with respect to L and F_L.
      if (histogram.hd.biasMag > 0.0) {
        double[] chainRule = new double[2];
        gLdEdL = histogram.energyAndGradient2D(dUdLambda, chainRule);
        double dGdLambda = chainRule[0];
        double dGdFLambda = chainRule[1];
        dUdLambda += dGdLambda + dGdFLambda * d2UdL2;
      }

      // Compute the energy and gradient for the recursion worker at F(L) using interpolation.
      bias1D = histogram.energyAndGradient1D(true);
    }

    // The total bias energy is the sum of the 1D and 2D terms.
    biasEnergy = bias1D + gLdEdL;

    if (print) {
      logger.info(format(" Bias Energy        %16.8f", biasEnergy));
      logger.info(format(" %s %16.8f  (Kcal/mole)", "OST Potential    ", forceFieldEnergy + biasEnergy));
    }

    if (propagateLambda) {
      histogram.stochasticPostForce();

      histogram.ld.stepsTaken++;
      long energyCount = histogram.ld.stepsTaken;

      // Log the current Lambda state.
      int printFrequency = dynamicsOptions.getReportFrequency(100);
      if (energyCount % printFrequency == 0) {
        double dBdL = dUdLambda - dForceFieldEnergydL;
        double lambda = histogram.ld.lambda;
        int lambdaBins = histogram.hd.getLambdaBins();
        if (lambdaBins < 1000) {
          logger.info(format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f", lambda,
              lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, histogram.ld.thetaVelocity));
        } else {
          logger.info(format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f", lambda,
              lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, histogram.ld.thetaVelocity));
        }
      }

      // Metadynamics grid counts (every 'countInterval' steps).
      if (energyCount % histogram.hd.countInterval == 0) {
        histogram.addBias(dForceFieldEnergydL);

        // Locally optimize the current state.
        if (optimizationParameters.doOptimization) {
          optimizationParameters.optimize(forceFieldEnergy, x, tempGradient);
        }

        if (algorithmListener != null) {
          algorithmListener.algorithmUpdate(molecularAssembly);
        }
      }
    }

    totalEnergy = forceFieldEnergy + biasEnergy;

    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] gradient) {

    // Stochastic dynamics pre-force propagation.
    if (state != STATE.FAST && propagateLambda) {
      histogram.stochasticPreForce();
    }

    // Compute the force field energy and gradient.
    forceFieldEnergy = potential.energyAndGradient(x, gradient);

    // OST is propagated with the slowly varying terms.
    if (state == STATE.FAST) {
      return forceFieldEnergy;
    }

    // dU/dL is the partial derivative of the force field energy with respect to lambda.
    dForceFieldEnergydL = lambdaInterface.getdEdL();
    dUdLambda = dForceFieldEnergydL;
    d2UdL2 = lambdaInterface.getd2EdL2();

    gLdEdL = 0.0;
    double bias1D;
    // Update the free energy difference.
    histogram.updateFreeEnergyDifference(false, false);
    if (histogram.hd.metaDynamics) {
      bias1D = histogram.energyAndGradientMeta(true);
    } else {
      if (histogram.hd.biasMag > 0) {
        double[] chainRule = new double[2];
        gLdEdL = histogram.energyAndGradient2D(dUdLambda, chainRule);
        double dGdLambda = chainRule[0];
        double dGdFLambda = chainRule[1];

        // Lambda gradient due to recursion kernel G(L, F_L).
        dUdLambda += dGdLambda + dGdFLambda * d2UdL2;

        // Cartesian coordinate gradient due to recursion kernel G(L, F_L).
        fill(dUdXdL, 0.0);
        lambdaInterface.getdEdXdL(dUdXdL);
        for (int i = 0; i < nVariables; i++) {
          gradient[i] += dGdFLambda * dUdXdL[i];
        }
      }

      // Compute the energy and gradient for the recursion worker at F(L) using interpolation.
      bias1D = histogram.energyAndGradient1D(true);
    }

    // The total bias is the sum of 1D and 2D terms.
    biasEnergy = bias1D + gLdEdL;

    if (print) {
      logger.info(format(" %s %16.8f", "Bias Energy       ", biasEnergy));
      logger.info(format(" %s %16.8f  %s", "OST Potential    ", forceFieldEnergy + biasEnergy, "(Kcal/mole)"));
    }

    if (propagateLambda) {
      histogram.stochasticPostForce();

      histogram.ld.stepsTaken++;
      long energyCount = histogram.ld.stepsTaken;

      // Log the current Lambda state.
      int printFrequency = dynamicsOptions.getReportFrequency(100);
      if (energyCount % printFrequency == 0) {
        double dBdL = dUdLambda - dForceFieldEnergydL;
        int lambdaBin = histogram.indexForLambda();
        double lambda = histogram.ld.lambda;
        int lambdaBins = histogram.hd.getLambdaBins();
        if (lambdaBins < 1000) {
          logger.info(format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f", lambda,
              lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, histogram.ld.thetaVelocity));
        } else {
          logger.info(format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f", lambda,
              lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, histogram.ld.thetaVelocity));
        }
      }

      // Metadynamics grid counts (every 'countInterval' steps).
      if (energyCount % histogram.hd.countInterval == 0) {
        histogram.addBias(dForceFieldEnergydL);

        // Locally optimize the current state.
        if (optimizationParameters.doOptimization) {
          optimizationParameters.optimize(forceFieldEnergy, x, gradient);
        }

        if (algorithmListener != null) {
          algorithmListener.algorithmUpdate(molecularAssembly);
        }
      }
    }

    totalEnergy = forceFieldEnergy + biasEnergy;

    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    return potential.getAcceleration(acceleration);
  }

  public Histogram[] getAllHistograms() {
    int nHisto = allHistograms.size();
    Histogram[] ret = new Histogram[nHisto];
    ret = allHistograms.toArray(ret);
    return ret;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] doubles) {
    return potential.getCoordinates(doubles);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Crystal getCrystal() {
    return potential.getCrystal();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    potential.setCrystal(crystal);
  }

  /**
   * Returns the number of energy evaluations performed by this OST, including those picked up in the
   * lambda file.
   *
   * @return Number of energy steps taken by this walker.
   */
  public long getEnergyCount() {
    return histogram.ld.stepsTaken;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Potential.STATE getEnergyTermState() {
    return state;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setEnergyTermState(Potential.STATE state) {
    this.state = state;
    potential.setEnergyTermState(state);
  }

  /**
   * Getter for the field <code>forceFieldEnergy</code>.
   *
   * @return a double.
   */
  public double getForceFieldEnergy() {
    return forceFieldEnergy;
  }

  /**
   * Return the current 2D Histogram of counts.
   *
   * @return the Histogram.
   */
  public Histogram getHistogram() {
    return histogram;
  }

  /**
   * Getter for the field <code>lambda</code>.
   *
   * @return a double.
   */
  public double getLambda() {
    return histogram.ld.lambda;
  }

  /**
   * Setter for the field <code>lambda</code>.
   *
   * @param lambda a double.
   */
  public void setLambda(double lambda) {
    if (histogram != null) {
      lambda = histogram.mapLambda(lambda);
    } else {
      logger.warning(" OrthogonalSpaceTempering.setLambda was called before histogram constructed!");
      logger.info(Utilities.stackTraceToString(new RuntimeException()));
    }
    lambdaInterface.setLambda(lambda);
    histogram.ld.lambda = lambda;
    histogram.ld.theta = asin(sqrt(lambda));
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    return potential.getMass();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return potential.getNumberOfVariables();
  }

  /**
   * Return the OST optimization information.
   *
   * @return The OST optimization parameters.
   */
  public OptimizationParameters getOptimizationParameters() {
    return optimizationParameters;
  }

  /**
   * getPotentialEnergy.
   *
   * @return a {@link ffx.numerics.Potential} object.
   */
  public Potential getPotentialEnergy() {
    return potential;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    return potential.getPreviousAcceleration(previousAcceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return potential.getScaling();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(double[] scaling) {
    potential.setScaling(scaling);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * getTotaldEdLambda.
   *
   * @return a double.
   */
  public double getTotaldEdLambda() {
    return dUdLambda;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reference to each variables type.
   */
  @Override
  public Potential.VARIABLE_TYPE[] getVariableTypes() {
    return potential.getVariableTypes();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    return potential.getVelocity(velocity);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    throw new UnsupportedOperationException(
        " Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    return getTotaldEdLambda();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] gradient) {
    throw new UnsupportedOperationException(
        " Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
  }

  // TODO: Delete method when debugging of RepexOST is done.
  public void logOutputFiles(int index) {
    logger.info(format(" OST: Lambda file %s, histogram %s", histogram.ld.getLambdaFile(),
        allHistograms.get(index).hd.getHistogramFile()));
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    potential.setAcceleration(acceleration);
  }

  /**
   * If this flag is true, (lambda, dU/dL) Monte Carlo samples that have no weight in the Histogram
   * are rejected.
   *
   * @param hardWallConstraint If true, MC samples outside the current range are rejected.
   */
  public void setHardWallConstraint(boolean hardWallConstraint) {
    this.hardWallConstraint = hardWallConstraint;
  }

  public void setMolecularAssembly(MolecularAssembly molecularAssembly) {
    this.molecularAssembly = molecularAssembly;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    potential.setPreviousAcceleration(previousAcceleration);
  }

  /**
   * Indicate if the Lambda extended system particle should be propagated using Langevin dynamics.
   *
   * @param propagateLambda If true, Lambda will be propagated using Langevin dynamics.
   */
  public void setPropagateLambda(boolean propagateLambda) {
    this.propagateLambda = propagateLambda;
  }

  /**
   * If true, the Lambda extended system particle is propagated using Langevin dynamics.
   */
  public boolean getPropagateLambda() {
    return propagateLambda;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    potential.setVelocity(velocity);
  }

  /**
   * Switch to an alternate Histogram.
   *
   * @param index Index of the Histogram to use.
   */
  public void switchHistogram(int index) {
    histogramIndex = index;
    histogram = allHistograms.get(histogramIndex);
    logger.info(" OST switching to histogram " + histogramIndex);
  }

  @Override
  public void writeAdditionalRestartInfo(boolean recursive) {
    histogram.writeRestart();
    if (recursive) {
      potential.writeAdditionalRestartInfo(recursive);
    }
  }

  double getLambdaWriteOut() {
    return histogram.hd.resetHistogramAtLambda;
  }

  /**
   * getForceFielddEdL.
   *
   * @return a double.
   */
  double getForceFielddEdL() {
    return dForceFieldEnergydL;
  }

  /**
   * Getter for the field <code>biasEnergy</code>.
   *
   * @return a double.
   */
  double getBiasEnergy() {
    return biasEnergy;
  }

  /**
   * If the dUdLHardWall flag is set to true, this method will return false if the (lambda, dU/dL)
   * sample is has not been seen.
   *
   * @param lambda The proposed lambda value.
   * @param dUdL   The proposed dU/dL value.
   * @return Returns false only if the dUdLHardWall flag is true, and the (lambda, dU/dL) sample has
   * not been seen.
   */
  boolean insideHardWallConstraint(double lambda, double dUdL) {
    if (hardWallConstraint) {
      double weight = histogram.getRecursionKernelValue(lambda, dUdL);
      return weight > 0.0;
    }
    return true;
  }

  /**
   * Parameters for running local optimizations during OST sampling.
   */
  public class OptimizationParameters {

    /**
     * Flag to turn on OST optimization.
     *
     * <p>The default doOptimization = false.
     */
    private boolean doOptimization = false;
    /**
     * Reset unit cell parameters, molecular orientation and translation.
     */
    private final boolean doUnitCellReset;
    /**
     * OST optimization only runs if Lambda is greater than the lambdaCutoff.
     *
     * <p>The default lambdaCutoff = 0.8.
     */
    private final double lambdaCutoff;
    /**
     * The lowest energy found via optimizations.
     *
     * <p>The optimumEnergy is initially set to Double.MAX_VALUE.
     */
    private double optimumEnergy = Double.MAX_VALUE;
    /**
     * The OST optimization frequency
     *
     * <p>The default is once every 10,000 steps.
     */
    private final int frequency;
    /**
     * The OST optimization convergence criteria.
     *
     * <p>The default eps = 0.1.
     */
    private final double eps;
    /**
     * The OST tolerance when checking for equal energy after coordinate reversion.
     *
     * <p>The default is 1.0e-8 kcal/mol.
     */
    private final double tolerance;
    /**
     * The OST optimization energy window.
     *
     * <p>The default is 4.0 kcal/mol, which is convenient for small organic crystals.
     */
    private final double energyWindow;
    /**
     * Holds the lowest potential energy coordinates.
     */
    private double[] optimumCoords;
    /**
     * File instance used for saving optimized structures.
     */
    private File optimizationFile;
    /**
     * SystemFilter used to save optimized structures.
     */
    private SystemFilter optimizationFilter;

    /**
     * Empty constructor.
     */
    OptimizationParameters(CompositeConfiguration properties) {
      energyWindow = properties.getDouble("ost-opt-energy-window", 10.0);
      eps = properties.getDouble("ost-opt-eps", 0.1);
      tolerance = properties.getDouble("ost-opt-tolerance", 1.0e-8);
      frequency = properties.getInt("ost-opt-frequency", 10000);
      lambdaCutoff = properties.getDouble("ost-opt-lambda-cutoff", 0.8);
      doUnitCellReset = properties.getBoolean("ost-opt-unitcell-reset", false);
    }

    private void log() {
      logger.info("\n Optimization Parameters");
      logger.info(format("  Energy Window:                  %6.3f (kcal/mol)", energyWindow));
      logger.info(format("  EPS:                            %6.4f RMS (kcal/mol/Å)", eps));
      logger.info(format("  Tolerance:                      %6.4f (kcal/mol)", tolerance));
      logger.info(format("  Frequency:                      %6d (steps)", frequency));
      logger.info(format("  Lambda Cutoff:                  %6.4f", lambdaCutoff));
      logger.info(format("  Unit Cell Reset:                %6B", doUnitCellReset));
      logger.info(format("  File:                           %s", optimizationFile.getName()));
    }

    /**
     * getOptimumCoordinates.
     *
     * @return an array of {@link double} objects.
     */
    public double[] getOptimumCoordinates() {
      if (optimumEnergy < Double.MAX_VALUE) {
        return optimumCoords;
      } else {
        logger.info(
            "Lambda optimization cutoff was not reached. Try increasing the number of timesteps.");
        return null;
      }
    }

    /**
     * getOptimumEnergy.
     *
     * @return a double.
     */
    public double getOptimumEnergy() {
      if (optimumEnergy == Double.MAX_VALUE) {
        logger.info(
            "Lambda optimization cutoff was not reached. Try increasing the number of timesteps.");
      }
      return optimumEnergy;
    }

    /**
     * Run a local optimization.
     *
     * @param e        Current energy.
     * @param x        Current atomic coordinates.
     * @param gradient Work array for collecting the gradient.
     */
    public void optimize(double e, double[] x, @Nullable double[] gradient) {

      // Return if the optimization flag is not set, or if lambda is not beyond the cutoff.
      double lambda = histogram.ld.lambda;
      long stepsTaken = histogram.ld.stepsTaken;
      if (doOptimization && lambda > lambdaCutoff && stepsTaken % frequency == 0) {
        if (gradient == null) {
          gradient = new double[x.length];
        }
      } else {
        return;
      }

      logger.info(format("\n OST Minimization (Step %d)", stepsTaken));

      // Set the underlying Potential's Lambda value to 1.0.
      lambdaInterface.setLambda(1.0);

      // Use all energy terms.
      potential.setEnergyTermState(Potential.STATE.BOTH);

      // Turn off the Barostat.
      boolean origBaroActive = true;
      if (barostat != null) {
        origBaroActive = barostat.isActive();
        barostat.setActive(false);
      }

      // Optimize the system.
      try {
        double startingEnergy = potential.energy(x);
        Minimize minimize = new Minimize(molecularAssembly, potential, algorithmListener);
        minimize.minimize(eps);
        // Collect the minimum energy.
        double minEnergy = potential.getTotalEnergy();
        // Check for a new minimum within an energy window of the lowest energy structure found.
        if (minEnergy < optimumEnergy + energyWindow) {
          if (minEnergy < optimumEnergy) {
            optimumEnergy = minEnergy;
          }
          int n = potential.getNumberOfVariables();
          optimumCoords = new double[n];
          optimumCoords = potential.getCoordinates(optimumCoords);
          double mass = molecularAssembly.getMass();
          Crystal crystal = molecularAssembly.getCrystal();
          double density = crystal.getDensity(mass);
          optimizationFilter.writeFile(optimizationFile, true);
          Crystal uc = crystal.getUnitCell();
          logger.info(format(" Minimum: %12.6f %s (%12.6f g/cc) optimized from %12.6f at step %d.",
              minEnergy, uc.toShortString(), density, startingEnergy, stepsTaken));
        }
      } catch (Exception ex) {
        String message = ex.getMessage();
        logger.info(
            format(" Exception minimizing coordinates at lambda=%8.6f\n %s.", lambda, message));
        logger.info(" Sampling will continue.");
      }

      // Set the underlying Potential's Lambda value back to current lambda value.
      lambdaInterface.setLambda(lambda);

      // Reset the Potential State
      potential.setEnergyTermState(state);

      // Reset the Barostat
      if (barostat != null) {
        barostat.setActive(origBaroActive);
      }

      if (doUnitCellReset) {
        logger.info("\n Resetting Unit Cell");
        double mass = molecularAssembly.getMass();
        double density = molecularAssembly.getCrystal().getDensity(mass);
        molecularAssembly.applyRandomDensity(density);
        molecularAssembly.applyRandomSymOp(0.0);
        lambda = 0.0;
        lambdaInterface.setLambda(lambda);
      } else {
        // Revert to the coordinates and gradient prior to optimization.
        double eCheck = potential.energyAndGradient(x, gradient);

        if (abs(eCheck - e) > tolerance) {
          logger.warning(
              format(" Optimization could not revert coordinates %16.8f vs. %16.8f.", e, eCheck));
        }
      }
    }

    /**
     * setOptimization.
     *
     * @param doOptimization    a boolean.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     */
    public void setOptimization(boolean doOptimization, MolecularAssembly molecularAssembly) {
      this.doOptimization = doOptimization;
      OrthogonalSpaceTempering.this.molecularAssembly = molecularAssembly;
      File file = molecularAssembly.getFile();
      String fileName = FilenameUtils.removeExtension(file.getAbsolutePath());
      String ext = FilenameUtils.getExtension(file.getAbsolutePath());
      if (optimizationFilter == null) {
        if (ext.toUpperCase().contains("XYZ")) {
          optimizationFile = new File(fileName + "_opt.arc");
          optimizationFilter = new XYZFilter(optimizationFile, molecularAssembly,
              molecularAssembly.getForceField(), molecularAssembly.getProperties());
        } else {
          optimizationFile = new File(fileName + "_opt.pdb");
          PDBFilter pdbFilter = new PDBFilter(optimizationFile, molecularAssembly,
              molecularAssembly.getForceField(), molecularAssembly.getProperties());
          int models = pdbFilter.countNumModels();
          pdbFilter.setModelNumbering(models);
          optimizationFilter = pdbFilter;
        }
      }

      if (this.doOptimization) {
        log();
      }
    }
  }

  /**
   * Store and operate on a 2D Histogram of (Lambda, dU/dL) observations to produce an OST bias.
   *
   * @author Michael J. Schnieders
   * @since 1.0
   */
  public class Histogram {

    /**
     * If a 2D bias is in use (i.e. lambda-bias-cutoff > 0), do not normalize bias height.
     */
    private static final double TWO_D_NORMALIZATION = 1.0;
    /**
     * If a 1D bias is in use (i.e. lambda-bias-cutoff == 0), normalize bias height to deposit the same
     * volume of bias.
     */
    private static final double ONE_D_NORMALIZATION = Math.sqrt(2.0 * Math.PI);
    /**
     * If the real bias magnitude is 0 (for no 2D bias), temporarily set biasMag to this value to calculate the ensemble average dU/dL.
     *
     * <p>Any value that does not overflow/underflow double precision summations of the count matrix
     * will give identical results. For example, values of 1.0e-20, 1.0 and 1.0e20 were tested and
     * found to give identical results.
     *
     * <p>Thus, a value of 1.0 is good choice.
     */
    private static final double PSEUDO_BIAS_MAGNITUDE = 1.0;
    /**
     * Parallel Java world communicator.
     */
    protected final Comm world;
    /**
     * Rank of this process.
     */
    protected final int rank;
    /**
     * 1D PMF with respect to lambda F(L).
     */
    final double[] ensembleAveragedUdL;
    /**
     * This deltaT is used to determine the tempering weight as described below for the
     * temperingWeight variable.
     *
     * <p>deltaT = temperingFactor * kB * T.
     */
    private final double deltaT;
    /**
     * The integration algorithm used for thermodynamic integration.
     */
    private final IntegrationType integrationType;
    /**
     * Send OST counts asynchronously.
     */
    private final SendAsynchronous sendAsynchronous;
    /**
     * Send OST counts synchronously.
     */
    private final SendSynchronous sendSynchronous;
    /**
     * The Dama et al. transition-tempering weight. The initial temperingWeight = 1.0,
     * and more generally temperingWeight = exp(-biasHeight/deltaT)
     */
    private double temperingWeight = 1.0;
    /**
     * If the recursion kernel becomes too large or too small for some combinations of (Lambda,
     * dU/dL), then its statistical weight = exp(kernel * beta) cannot be represented by a double
     * value.
     */
    private double[] kernelValues;
    /**
     * Most recent lambda values for each Walker.
     */
    private double lastReceivedLambda;
    /**
     * Most recent dU/dL value for each walker.
     */
    private double lastReceiveddUdL;
    private final boolean spreadBias;
    /**
     * The HistogramData object contains the parameters for the Histogram that are saved to disk.
     */
    final HistogramData hd;
    /**
     * The LambdaData object contains the parameters for the Lambda particle that are saved to disk.
     */
    final LambdaData ld;
    /**
     * The Stochastic integrator used to propagate the lambda particle.
     */
    private final Stochastic stochastic;
    /**
     * The state of the lambda particle.
     */
    private final SystemState lambdaState;
    /**
     * Number of hills when the last free energy estimate was made.
     */
    private int currentNumberOfHills = 0;
    /**
     * The current free energy difference.
     */
    private double currentFreeEnergyDifference = 0.0;

    /**
     * Histogram constructor.
     *
     * @param properties    a CompositeConfiguration used to configure the Histogram.
     * @param histogramData An object containing the values this Histogram will be set to.
     */
    Histogram(CompositeConfiguration properties, HistogramData histogramData, LambdaData lambdaData) {
      hd = histogramData;
      ld = lambdaData;

      /*
      double gaussNormalization;
      if (hd.lambdaBiasCutoff == 0 && !hd.wasHistogramRead()) {
        gaussNormalization = ONE_D_NORMALIZATION;
        logger.info(format(" Bias magnitude multiplied by a factor of %.4f "
            + "sqrt(2*pi) to match 1D Gaussian volume to 2D Gaussian volume.", gaussNormalization));
      } else {
        gaussNormalization = TWO_D_NORMALIZATION;
      }
      biasMag = settings.getBiasMag() * gaussNormalization;
       */

      deltaT = hd.temperingFactor * R * dynamicsOptions.getTemperature();

      // Allocate space to compute the <dU/dL>
      ensembleAveragedUdL = new double[hd.getLambdaBins()];

      // Allocate space to regularize kernel values.
      kernelValues = new double[hd.getDUDLBins()];

      String propString = properties.getString("ost-integrationType", "SIMPSONS");
      IntegrationType testType;
      try {
        testType = IntegrationType.valueOf(propString.toUpperCase());
      } catch (Exception ex) {
        logger.warning(format(" Invalid argument %s to ost-integrationType; resetting to SIMPSONS", propString));
        testType = SIMPSONS;
      }
      integrationType = testType;

      spreadBias = properties.getBoolean("ost-spread-bias", false);

      /*
       Set up the multi-walker communication variables for Parallel Java
       communication between nodes.
      */
      world = Comm.world();
      int numProc = world.size();
      rank = world.rank();
      if (hd.asynchronous) {
        // Use asynchronous communication.
        sendAsynchronous = new SendAsynchronous(this);
        sendAsynchronous.start();
        sendSynchronous = null;
      } else {
        Histogram[] histograms = new Histogram[numProc];
        int[] rankToHistogramMap = new int[numProc];
        for (int i = 0; i < numProc; i++) {
          histograms[i] = this;
          rankToHistogramMap[i] = 0;
        }
        sendSynchronous = new SendSynchronous(histograms, rankToHistogramMap);
        sendAsynchronous = null;
      }

      lastReceivedLambda = ld.lambda;
      if (hd.discreteLambda) {
        lastReceivedLambda = mapLambda(lastReceivedLambda);
        logger.info(format(" Discrete lambda: initializing lambda to nearest bin %.5f", lastReceivedLambda));
        ld.lambda = lastReceivedLambda;
        ld.theta = asin(sqrt(lastReceivedLambda));
        lambdaInterface.setLambda(lastReceivedLambda);
        lambdaState = null;
        stochastic = null;
      } else {
        // Configure the Stochastic integrator.
        lambdaState = new SystemState(1);
        stochastic = new Stochastic(lambdaParticleOptions.getLambdaFriction(), lambdaState);
        stochastic.setTemperature(dynamicsOptions.getTemperature());
        stochastic.setTimeStep(dynamicsOptions.getDtPsec());
        double[] mass = lambdaState.getMass();
        double[] thetaPosition = lambdaState.x();
        double[] thetaVelocity = lambdaState.v();
        double[] thetaAccel = lambdaState.a();
        mass[0] = lambdaParticleOptions.getLambdaMass();
        thetaPosition[0] = ld.theta;
        thetaVelocity[0] = ld.thetaVelocity;
        thetaAccel[0] = ld.thetaAcceleration;
      }

      lastReceiveddUdL = getdEdL();

      if (hd.wasHistogramRead()) {
        updateFreeEnergyDifference(true, false);
      }
    }

    public void disableResetStatistics() {
      hd.resetHistogram = false;
    }

    /**
     * evaluateTotalBias.
     *
     * @param bias If false, return the negative of the Total OST bias.
     * @return A StringBuffer The total OST Bias.
     */
    public StringBuffer evaluateTotalOSTBias(boolean bias) {
      StringBuffer sb = new StringBuffer();
      double[] chainRule = new double[2];
      for (int dUdLBin = 0; dUdLBin < hd.getDUDLBins(); dUdLBin++) {
        double currentdUdL = dUdLforBin(dUdLBin);
        sb.append(format(" %16.8f", currentdUdL));
        for (int lambdaBin = 0; lambdaBin < hd.getLambdaBins(); lambdaBin++) {
          double currentLambda = lambdaForIndex(lambdaBin);
          double bias1D = energyAndGradient1D(currentLambda, false);
          double bias2D = energyAndGradient2D(currentLambda, currentdUdL, chainRule, hd.biasMag);
          double totalBias = bias1D + bias2D;
          // If the bias flag is false, turn the total bias into the PMF.
          if (!bias) {
            totalBias = -totalBias;
          }
          sb.append(format(" %16.8f", totalBias));
        }
        sb.append("\n");
      }
      return sb;
    }

    /**
     * evaluate2DOSTBias.
     *
     * @param bias If false, return the negative of the 2D OST bias.
     * @return A StringBuffer with the 2D OST bias.
     */
    public StringBuffer evaluate2DOSTBias(boolean bias) {
      StringBuffer sb = new StringBuffer();
      double[] chainRule = new double[2];
      for (int dUdLBin = 0; dUdLBin < hd.getDUDLBins(); dUdLBin++) {
        double currentdUdL = dUdLforBin(dUdLBin);
        sb.append(format(" %16.8f", currentdUdL));
        for (int lambdaBin = 0; lambdaBin < hd.getLambdaBins(); lambdaBin++) {
          double currentLambda = lambdaForIndex(lambdaBin);
          double bias2D = energyAndGradient2D(currentLambda, currentdUdL, chainRule, hd.biasMag);
          if (!bias) {
            bias2D = -bias2D;
          }
          sb.append(format(" %16.8f", bias2D));
        }
        sb.append("\n");
      }
      return sb;
    }

    /**
     * For MPI parallel jobs, returns true if the walkers are independent (i.e. contribute to only
     * their own histogram).
     *
     * @return True if the walkers are independent.
     */
    public boolean getIndependentWalkers() {
      return hd.independentWalkers;
    }

    public double getLambdaResetValue() {
      return hd.resetHistogramAtLambda;
    }

    /**
     * For MPI parallel jobs, return the rank of this process.
     *
     * @return The rank of this process.
     */
    public int getRank() {
      return rank;
    }

    public boolean getResetStatistics() {
      return hd.resetHistogram;
    }

    /**
     * Return the SynchronousSend associated with this Histogram, if any.
     *
     * @return The SynchronousSend, if any.
     */
    public Optional<SendSynchronous> getSynchronousSend() {
      return Optional.ofNullable(sendSynchronous);
    }

    public void setLastReceiveddUdL(double lastReceiveddUdL) {
      this.lastReceiveddUdL = lastReceiveddUdL;
    }

    public double updateFreeEnergyDifference(boolean print, boolean save) {
      if (hd.metaDynamics) {
        return updateMetaDynamicsFreeEnergyDifference(print, save);
      } else {
        return updateOSTFreeEnergyDifference(print, save);
      }
    }

    /**
     * Update the free energy estimate for Meta Dynamics.
     *
     * @param print Whether to write the histogram to screen.
     * @param save  Whether to write the histogram to disk.
     * @return Free energy (via integration of ensemble-average dU/dL)
     */
    public double updateMetaDynamicsFreeEnergyDifference(boolean print, boolean save) {
      synchronized (this) {
        if (!print && currentNumberOfHills == hd.counts) {
          return currentFreeEnergyDifference;
        }
        double freeEnergyDifferenceMeta = 0.0;
        double freeEnergyDifferenceTI = 0.0;

        double minFL = Double.MAX_VALUE;

        int lambdaBins = hd.getLambdaBins();
        double[] metaFreeEnergy = new double[lambdaBins];

        // Total histogram weight.
        double totalWeight = 0;
        StringBuilder stringBuilder = new StringBuilder();

        // Loop over lambda bins, computing <dU/dL> for each bin.
        for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
          int firstdUdLBin = firstdUdLBin(lambdaBin);
          int lastdUdLBin = lastdUdLBin(lambdaBin);
          double lambdaCount = 0;
          // The dUdL range sampled for lambda.
          double mindUdLForLambda = 0.0;
          double maxdUdLforLambda = 0.0;
          double maxBias = 0;
          if (firstdUdLBin == -1 || lastdUdLBin == -1) {
            ensembleAveragedUdL[lambdaBin] = 0.0;
            minFL = 0.0;
          } else {
            double ensembleAverageFLambda = 0.0;
            double partitionFunction = 0.0;

            // Evaluate and regularize all kernel values for this value of lambda.
            double offset = evaluateKernelForLambda(lambdaBin, firstdUdLBin, lastdUdLBin);

            for (int dUdLBin = firstdUdLBin; dUdLBin <= lastdUdLBin; dUdLBin++) {
              double kernel = kernelValues[dUdLBin];

              if (kernel - offset > maxBias) {
                maxBias = kernel - offset;
              }

              // The weight is just the kernel value for no 2D bias.
              partitionFunction += kernel;

              double currentdUdL = dUdLforBin(dUdLBin);
              ensembleAverageFLambda += currentdUdL * kernel;
              lambdaCount += getRecursionKernelValue(lambdaBin, dUdLBin);
            }
            if (minFL > maxBias) {
              minFL = maxBias;
            }
            ensembleAveragedUdL[lambdaBin] = (partitionFunction == 0) ? 0 : ensembleAverageFLambda / partitionFunction;
            mindUdLForLambda = hd.dUdLMinimum + firstdUdLBin * hd.dUdLBinWidth;
            maxdUdLforLambda = hd.dUdLMinimum + (lastdUdLBin + 1) * hd.dUdLBinWidth;
          }

          double deltaFreeEnergy = ensembleAveragedUdL[lambdaBin] * deltaForLambdaBin(lambdaBin);
          freeEnergyDifferenceTI += deltaFreeEnergy;
          totalWeight += lambdaCount;


          double lambdaBinWidth = hd.lambdaBinWidth;
          double llL = lambdaBin * lambdaBinWidth - hd.lambdaBinWidth_2;
          double ulL = llL + lambdaBinWidth;
          if (llL < 0.0) {
            llL = 0.0;
          }
          if (ulL > 1.0) {
            ulL = 1.0;
          }

          double midLambda = llL;
          if (!hd.discreteLambda) {
            midLambda = (llL + ulL) / 2.0;
          }
          double bias1D = energyAndGradient1D(midLambda, false);
          metaFreeEnergy[lambdaBin] = energyAndGradientMeta(midLambda, false);
          freeEnergyDifferenceMeta = -(metaFreeEnergy[lambdaBin] - metaFreeEnergy[0]);

          if (print || save) {
            stringBuilder.append(format(" %6.2e %7.5f %7.1f %7.1f %9.2f %9.2f %9.2f %9.2f %9.2f\n", lambdaCount,
                midLambda, mindUdLForLambda, maxdUdLforLambda, ensembleAveragedUdL[lambdaBin],
                bias1D, freeEnergyDifferenceTI, metaFreeEnergy[lambdaBin], freeEnergyDifferenceMeta));
          }
        }

        double temperingOffset = hd.getTemperingOffset();
        double temperEnergy = (minFL > temperingOffset) ? temperingOffset - minFL : 0;
        temperingWeight = exp(temperEnergy / deltaT);

        if (print) {
          logger.info("  Weight   Lambda      dU/dL Bins   <dU/dL>      g(L) dG_OST(L)  Meta(L) dG_Meta(L)");
          logger.info(stringBuilder.toString());
          logger.info(" Histogram Evaluation");
          logger.info(format("  Free Energy Difference: %12.4f kcal/mol", freeEnergyDifferenceMeta));
          logger.info(format("  Number of Hills:        %12d", hd.counts));
          logger.info(format("  Total Bias Added:       %12.4f kcal/mol", totalWeight * hd.biasMag));
          logger.info(format("  Minimum Bias:           %12.4f kcal/mol", minFL));
          logger.info(format("  Tempering Percentage:   %12.4f %%\n", temperingWeight * 100));
        }

        if (save) {
          String modelFilename = molecularAssembly.getFile().getAbsolutePath();
          File saveDir = new File(FilenameUtils.getFullPath(modelFilename));
          String dirName = saveDir + File.separator;
          String fileName = dirName + "histogram.txt";
          try {
            logger.info(" Writing " + fileName);
            PrintWriter printWriter = new PrintWriter(fileName);
            printWriter.write(stringBuilder.toString());
            printWriter.close();
          } catch (Exception e) {
            logger.info(format(" Failed to write %s.", fileName));
          }
        }

        currentNumberOfHills = hd.counts;
        currentFreeEnergyDifference = freeEnergyDifferenceMeta;
        return currentFreeEnergyDifference;
      }
    }

    /**
     * Eqs. 7 and 8 from the 2012 Crystal Thermodynamics paper.
     *
     * @param print Whether to write the histogram to screen.
     * @param save  Whether to write the histogram to disk.
     * @return Free energy (via integration of ensemble-average dU/dL)
     */
    public double updateOSTFreeEnergyDifference(boolean print, boolean save) {
      synchronized (this) {
        if (!print && currentNumberOfHills == hd.counts) {
          return currentFreeEnergyDifference;
        }

        double freeEnergyDifferenceOST = 0.0;
        double minFL = Double.MAX_VALUE;

        // If the bias magnitude is zero, computing <dU/dL> from
        // counts will not be correct. Assign a temporary non-zero bias magnitude.
        boolean biasMagZero = hd.biasMag <= 0.0;

        // Total histogram weight.
        double totalWeight = 0;
        double beta = 1.0 / (R * dynamicsOptions.getTemperature());
        StringBuilder stringBuilder = new StringBuilder();

        // Loop over lambda bins, computing <dU/dL> for each bin.
        for (int lambdaBin = 0; lambdaBin < hd.lambdaBins; lambdaBin++) {
          int firstdUdLBin = firstdUdLBin(lambdaBin);
          int lastdUdLBin = lastdUdLBin(lambdaBin);
          double lambdaCount = 0;
          // The FL range sampled for lambda bin [iL*dL .. (iL+1)*dL]
          double mindUdLforLambda = 0.0;
          double maxdUdLforLambda = 0.0;
          double maxBias = 0;
          if (firstdUdLBin == -1 || lastdUdLBin == -1) {
            ensembleAveragedUdL[lambdaBin] = 0.0;
            minFL = 0.0;
          } else {
            double ensembleAverage = 0.0;
            double partitionFunction = 0.0;

            // Evaluate and regularize all kernel values for this value of lambda.
            double offset = evaluateKernelForLambda(lambdaBin, firstdUdLBin, lastdUdLBin);

            for (int dUdLBin = firstdUdLBin; dUdLBin <= lastdUdLBin; dUdLBin++) {
              double kernel = kernelValues[dUdLBin];

              if (kernel - offset > maxBias) {
                maxBias = kernel - offset;
              }

              double weight;
              if (biasMagZero) {
                // The weight is just the kernel value for no 2D bias.
                weight = kernel;
              } else {
                // The Boltzmann weight based on temperature and the 2D bias height.
                weight = exp(kernel * beta);
              }
              partitionFunction += weight;

              double currentdUdL = dUdLforBin(dUdLBin);
              ensembleAverage += currentdUdL * weight;
              lambdaCount += getRecursionKernelValue(lambdaBin, dUdLBin);
            }
            if (minFL > maxBias) {
              minFL = maxBias;
            }
            ensembleAveragedUdL[lambdaBin] = (partitionFunction == 0) ? 0 : ensembleAverage / partitionFunction;
            mindUdLforLambda = hd.dUdLMinimum + firstdUdLBin * hd.dUdLBinWidth;
            maxdUdLforLambda = hd.dUdLMinimum + (lastdUdLBin + 1) * hd.dUdLBinWidth;
          }

          double deltaFreeEnergy = ensembleAveragedUdL[lambdaBin] * deltaForLambdaBin(lambdaBin);
          freeEnergyDifferenceOST += deltaFreeEnergy;
          totalWeight += lambdaCount;

          if (print || save) {
            double llL = lambdaBin * hd.lambdaBinWidth - hd.lambdaBinWidth_2;
            double ulL = llL + hd.lambdaBinWidth;
            if (llL < 0.0) {
              llL = 0.0;
            }
            if (ulL > 1.0) {
              ulL = 1.0;
            }

            double midLambda = llL;
            if (!hd.discreteLambda) {
              midLambda = (llL + ulL) / 2.0;
            }
            double bias1D = energyAndGradient1D(midLambda, false);

            double bias2D = 0.0;
            if (!biasMagZero) {
              bias2D = computeBiasEnergy(midLambda, ensembleAveragedUdL[lambdaBin]) - bias1D;
            }

            stringBuilder.append(
                format(" %6.2e %7.5f %7.1f %7.1f %8.2f %8.2f %8.2f %8.2f %8.2f   %8.2f\n", lambdaCount,
                    midLambda, mindUdLforLambda, maxdUdLforLambda, ensembleAveragedUdL[lambdaBin],
                    bias1D, bias2D, bias1D + bias2D, freeEnergyDifferenceOST, bias1D + bias2D + freeEnergyDifferenceOST));
          }
        }

        if (!biasMagZero) {
          double temperingOffset = hd.getTemperingOffset();
          double temperEnergy = (minFL > temperingOffset) ? temperingOffset - minFL : 0;
          temperingWeight = exp(temperEnergy / deltaT);
        }

        freeEnergyDifferenceOST = integrateNumeric(ensembleAveragedUdL, integrationType);

        if (print) {
          logger.info("  Weight   Lambda      dU/dL Bins  <dU/dL>    g(L)  f(L,<dU/dL>) Bias    dG(L) Bias+dG(L)");
          logger.info(stringBuilder.toString());
          logger.info(" Histogram Evaluation");
          logger.info(format("  Free Energy Difference: %12.4f kcal/mol", freeEnergyDifferenceOST));
          logger.info(format("  Number of Hills:        %12d", hd.counts));
          logger.info(format("  Total Bias Added:       %12.4f kcal/mol", totalWeight * hd.biasMag));
          logger.info(format("  Minimum Bias:           %12.4f kcal/mol", minFL));
          logger.info(format("  Tempering Percentage:   %12.4f %%\n", temperingWeight * 100));
        }

        if (save) {
          String modelFilename = molecularAssembly.getFile().getAbsolutePath();
          File saveDir = new File(FilenameUtils.getFullPath(modelFilename));
          String dirName = saveDir + File.separator;
          String fileName = dirName + "histogram.txt";
          try {
            logger.info("\n Writing " + fileName);
            PrintWriter printWriter = new PrintWriter(fileName);
            printWriter.write(stringBuilder.toString());
            printWriter.close();
          } catch (Exception e) {
            logger.info(format(" Failed to write %s.", fileName));
          }
        }

        currentNumberOfHills = hd.counts;
        currentFreeEnergyDifference = freeEnergyDifferenceOST;
        return currentFreeEnergyDifference;
      }
    }

    /**
     * For thermodynamic integration, return the integration width for the given Lambda lambdaBin.
     *
     * @param lambdaBin The lambda lambdaBin.
     * @return The integration width.
     */
    private double deltaForLambdaBin(int lambdaBin) {
      if (!hd.discreteLambda && (lambdaBin == 0 || lambdaBin == hd.lambdaBins - 1)) {
        // The first and last lambda bins are half size for continuous lambda.
        return hd.lambdaBinWidth_2;
      } else if (hd.discreteLambda && lambdaBin == 0) {
        // The free energy change to move from L=0 to L=0 is zero.
        return 0.0;
      }
      // All other cases.
      return hd.lambdaBinWidth;
    }

    /**
     * Returns the index of the first dU/dL bin with counts, or -1 if there are no counts for the
     * given lambda bin.
     *
     * @param lambdaBin Lambda bin to consider.
     * @return Index of the first dUdL bin with counts.
     */
    private int firstdUdLBin(int lambdaBin) {
      // Synchronize use of the recursion kernel.
      synchronized (this) {
        // Find the smallest FL bin that has counts.
        for (int jFL = 0; jFL < hd.dUdLBins; jFL++) {
          double count = hd.zHistogram[lambdaBin][jFL];
          if (count > 0) {
            return jFL;
          }
        }
        return -1;
      }
    }

    /**
     * Returns the index of the last dU/dL bin with counts, or -1 if there are no counts for the
     * given lambda bin.
     *
     * @param lambdaBin Lambda bin to consider.
     * @return Index of the last dUdL bin with counts.
     */
    private int lastdUdLBin(int lambdaBin) {
      // Synchronize use of the recursion kernel.
      synchronized (this) {
        // Find the largest FL bin that has counts.
        for (int jFL = hd.dUdLBins - 1; jFL >= 0; jFL--) {
          double count = hd.zHistogram[lambdaBin][jFL];
          if (count > 0) {
            return jFL;
          }
        }
        return -1;
      }
    }

    /**
     * Write histogram and lambda restart files.
     * <p>
     * For MPI parallel jobs, only rank 0 writes a histogram restart file (unless the "independentWrite" flag is set).
     */
    void writeRestart() {
      if (rank == 0 || hd.independentWrite) {
        updateFreeEnergyDifference(true, false);
        try {
          hd.writeHistogram();
          logger.info(format(" Wrote histogram restart to: %s.", hd.getHistogramFileName()));
        } catch (Exception ex) {
          String message = format(" Exception writing histogram restart file %s.", hd.getHistogramFileName());
          logger.log(Level.INFO, Utilities.stackTraceToString(ex));
          logger.log(Level.SEVERE, message, ex);
        }
      }

      // All ranks write a lambda restart file.
      try {
        ld.writeLambdaData();
        logger.info(format(" Wrote lambda restart to:    %s.", ld.getLambdaFileName()));
      } catch (Exception ex) {
        String message = format(" Exception writing lambda restart file %s.", ld.getLambdaFileName());
        logger.log(Level.INFO, Utilities.stackTraceToString(ex));
        logger.log(Level.SEVERE, message, ex);
      }
    }

    private double mapLambda(double lambda) {
      if (hd.discreteLambda) {
        return hd.lambdaLadder[indexForDiscreteLambda(lambda)];
      } else {
        return lambda;
      }
    }

    /**
     * For continuous lambda, the returned value is the lambda bin. For discrete lambda, the returned
     * value is the discrete lambda index.
     *
     * @param lambda a double.
     * @return a int.
     */
    int indexForLambda(double lambda) {
      if (hd.discreteLambda) {
        return indexForDiscreteLambda(lambda);
      } else {
        return indexForContinuousLambda(lambda);
      }
    }

    int indexForLambda() {
      return indexForLambda(ld.lambda);
    }

    private double lambdaForIndex(int bin) {
      if (hd.discreteLambda) {
        return hd.lambdaLadder[bin];
      } else {
        return bin * hd.lambdaBinWidth;
      }
    }

    private int indexForContinuousLambda(double lambda) {
      int lambdaBin = (int) floor((lambda - hd.minLambda) / hd.lambdaBinWidth);
      if (lambdaBin < 0) {
        lambdaBin = 0;
      }
      if (lambdaBin >= hd.lambdaBins) {
        lambdaBin = hd.lambdaBins - 1;
      }
      return lambdaBin;
    }

    private int indexForDiscreteLambda(double lambda) {
      assert hd.discreteLambda && hd.lambdaLadder != null && hd.lambdaLadder.length > 0;

      int initialGuess = indexForContinuousLambda(lambda);
      double minErr = Double.MAX_VALUE;
      int minErrBin = -1;
      for (int i = -1; i < 2; i++) {
        int guessBin = i + initialGuess;
        if (guessBin < 0 || guessBin >= hd.lambdaBins) {
          continue;
        }
        double guessLam = hd.lambdaLadder[guessBin];
        double guessErr = Math.abs(guessLam - lambda);
        if (guessErr < minErr) {
          minErr = guessErr;
          minErrBin = guessBin;
        }
      }

      assert minErr < 1.0E-6;
      return minErrBin;
    }

    /**
     * Find the bin for the supplied dEdLambda.
     * <p>
     * If the supplied dEdL is outside the range of the count matrix, then -1 is returned.
     *
     * @param dUdL a double.
     * @return The dUdL bin.
     */
    int binFordUdL(double dUdL) {

      // No counts below mindUdL.
      if (dUdL < hd.dUdLMinimum) {
        return -1;
      }

      // No counts above the maxdUdL.
      if (dUdL > hd.dUdLMaximum) {
        return -1;
      }

      int bin = (int) floor((dUdL - hd.dUdLMinimum) / hd.dUdLBinWidth);

      if (bin == hd.dUdLBins) {
        bin = hd.dUdLBins - 1;
      }

      if (bin < 0) {
        bin = 0;
      }

      return bin;
    }

    /**
     * Return the value of a recursion kernel bin. Mirror conditions are applied to the lambda bin.
     * If the dUdLBin is outside the range of the count matrix, zero is returned.
     *
     * @param lambdaBin The lambda bin.
     * @param dUdLBin   The dU/dL bin.
     * @return The value of the bin.
     */
    double getRecursionKernelValue(int lambdaBin, int dUdLBin) {
      // Synchronize use of the recursion kernel.
      synchronized (this) {
        // Apply lambda mirror condition.
        lambdaBin = lambdaMirror(lambdaBin);

        // For dUdL outside the count matrix the weight is 0.
        if (dUdLBin < 0 || dUdLBin >= hd.dUdLBins) {
          return 0.0;
        }

        return hd.zHistogram[lambdaBin][dUdLBin];
      }
    }

    /**
     * Return the value of a recursion kernel bin. Mirror conditions are applied to the lambda bin.
     * If the dUdL is outside the range of the count matrix, zero is returned.
     *
     * @param lambda the lambda value.
     * @param dUdL   the dU/dL value.
     * @return The value of the Histogram.
     */
    double getRecursionKernelValue(double lambda, double dUdL) {
      return getRecursionKernelValue(indexForLambda(lambda), binFordUdL(dUdL));
    }

    /**
     * Add to the value of a recursion kernel bin.
     *
     * @param currentLambda The lambda bin.
     * @param currentdUdL   The dU/dL bin.
     * @param value         The value of the bin.
     */
    void addToRecursionKernelValue(double currentLambda, double currentdUdL, double value) {
      synchronized (this) {
        if (spreadBias) {
          // Expand the recursion kernel if necessary.
          int dUdLbiasCutoff = hd.dUdLBiasCutoff;
          checkRecursionKernelSize(dUdLforBin(binFordUdL(currentdUdL) - dUdLbiasCutoff));
          checkRecursionKernelSize(dUdLforBin(binFordUdL(currentdUdL) + dUdLbiasCutoff));

          int currentLambdaBin = indexForLambda(currentLambda);
          int currentdUdLBin = binFordUdL(currentdUdL);

          // Variances are only used when dividing by twice their value.
          double invLs2 = 0.5 / hd.lambdaVariance;
          double invFLs2 = 0.5 / hd.dUdLVariance;

          // Compute the normalization.
          double normalize = 0.0;
          int lambdaBiasCutoff = hd.lambdaBiasCutoff;
          for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
            int lambdaBin = currentLambdaBin + iL;
            double deltaL = currentLambda - lambdaBin * hd.lambdaBinWidth;
            double deltaL2 = deltaL * deltaL;
            // Pre-compute the lambda bias.
            double L2exp = exp(-deltaL2 * invLs2);
            for (int jFL = -dUdLbiasCutoff; jFL <= dUdLbiasCutoff; jFL++) {
              int dUdLBin = currentdUdLBin + jFL;
              double deltaFL = currentdUdL - dUdLforBin(dUdLBin);
              double deltaFL2 = deltaFL * deltaFL;
              normalize += L2exp * exp(-deltaFL2 * invFLs2);
            }
          }

          for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
            int lambdaBin = currentLambdaBin + iL;
            double deltaL = currentLambda - lambdaBin * hd.lambdaBinWidth;
            double deltaL2 = deltaL * deltaL;
            // Pre-compute the lambda bias.
            double L2exp = exp(-deltaL2 * invLs2);
            lambdaBin = lambdaMirror(lambdaBin);
            for (int jFL = -dUdLbiasCutoff; jFL <= dUdLbiasCutoff; jFL++) {
              int dUdLBin = currentdUdLBin + jFL;
              double deltaFL = currentdUdL - dUdLforBin(dUdLBin);
              double deltaFL2 = deltaFL * deltaFL;
              double weight = value / normalize * L2exp * exp(-deltaFL2 * invFLs2);
              hd.zHistogram[lambdaBin][dUdLBin] += weight;
            }
          }
          hd.counts++;
        } else {
          // Check the recursion kernel size.
          checkRecursionKernelSize(currentdUdL);
          int lambdaBin = indexForLambda(currentLambda);
          int dUdLBin = binFordUdL(currentdUdL);
          try {
            hd.zHistogram[lambdaBin][dUdLBin] += value;
            hd.counts++;
          } catch (IndexOutOfBoundsException e) {
            logger.warning(format(
                " Count skipped in addToRecursionKernelValue due to an index out of bounds exception.\n L=%10.8f (%d), dU/dL=%10.8f (%d) and count=%10.8f",
                currentLambda, lambdaBin, currentdUdL, dUdLBin, value));
          }
        }
      }
    }

    /**
     * Allocate memory for the recursion kernel.
     */
    void allocateRecursionKernel() {
      // Synchronize updates of the recursion kernel.
      synchronized (this) {
        hd.zHistogram = new double[hd.lambdaBins][hd.dUdLBins];
        kernelValues = new double[hd.dUdLBins];
      }
    }

    /**
     * Integrates dUdL over lambda using more sophisticated techniques than midpoint rectangular
     * integration.
     *
     * <p>The ends (from 0 to dL and 1-dL to 1) are integrated with trapezoids for continuous
     * lambda.
     *
     * @param dUdLs dUdL at the midpoint of each bin.
     * @param type  Integration type to use.
     * @return Current delta-G estimate.
     */
    private double integrateNumeric(double[] dUdLs, IntegrationType type) {
      double val;
      if (hd.discreteLambda) {
        // Integrate between the first bin and the last bin.
        double[] lams = Integrate1DNumeric.generateXPoints(0.0, 1.0, hd.lambdaBins, false);
        DataSet dSet = new DoublesDataSet(lams, dUdLs, false);
        val = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, type);
      } else {
        // Integrate between the second bin midpoint and the second-to-last bin midpoint.
        double[] midLams = Integrate1DNumeric.generateXPoints(hd.lambdaBinWidth, 1.0 - hd.lambdaBinWidth,
            (hd.lambdaBins - 2), false);
        double[] midVals = Arrays.copyOfRange(dUdLs, 1, (hd.lambdaBins - 1));
        DataSet dSet = new DoublesDataSet(midLams, midVals, false);

        val = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, type);

        // Everything after this is just adding in the endpoint contributions.

        double dL_4 = hd.lambdaBinWidth_2 * 0.5;

        // Initially, assume dU/dL is exactly 0 at the endpoints. This is sometimes a true
        // assumption.
        double val0 = 0;
        double val1 = 0;

        // If we cannot guarantee that dUdL is exactly 0 at the endpoints, interpolate.
        if (!lambdaInterface.dEdLZeroAtEnds()) {
          double recipSlopeLen = 1.0 / (hd.lambdaBinWidth * 0.75);

          double slope = dUdLs[0] - dUdLs[1];
          slope *= recipSlopeLen;
          val0 = dUdLs[0] + (slope * dL_4);

          slope = dUdLs[hd.lambdaBins - 1] - dUdLs[hd.lambdaBins - 2];
          slope *= recipSlopeLen;
          val1 = dUdLs[hd.lambdaBins - 1] + (slope * dL_4);
          logger.fine(format(" Inferred dU/dL values at 0 and 1: %10.5g , %10.5g", val0, val1));
        }

        // Integrate trapezoids from 0 to the second bin midpoint, and from second-to-last bin
        // midpoint to 1.
        val += trapezoid(0, dL_4, val0, dUdLs[0]);
        val += trapezoid(dL_4, hd.lambdaBinWidth, dUdLs[0], dUdLs[1]);
        val += trapezoid(1.0 - hd.lambdaBinWidth, 1.0 - dL_4, dUdLs[hd.lambdaBins - 2],
            dUdLs[hd.lambdaBins - 1]);
        val += trapezoid(1.0 - dL_4, 1.0, dUdLs[hd.lambdaBins - 1], val1);
      }

      return val;
    }

    /**
     * Integrates a trapezoid.
     *
     * @param x0  First x point
     * @param x1  Second x point
     * @param fx0 First f(x) point
     * @param fx1 Second f(x) point
     * @return The area under a trapezoid.
     */
    private double trapezoid(double x0, double x1, double fx0, double fx1) {
      double val = 0.5 * (fx0 + fx1);
      val *= (x1 - x0);
      return val;
    }

    /**
     * Evaluate the kernel across dU/dL values at given value of lambda. The largest kernel value V
     * is used to define an offset (-V), which is added to all to kernel values. Then, the largest
     * kernel value is zero, and will result in a statistical weight of 1.
     *
     * <p>The offset avoids the recursion kernel becoming too large for some combinations of
     * (Lambda, dU/dL). This can result in its statistical weight = exp(kernel * beta) exceeding the
     * maximum representable double value.
     *
     * @param lambda Value of Lambda to evaluate the kernel for.
     * @param llFL   Lower FLambda bin.
     * @param ulFL   Upper FLambda bin.
     * @return the applied offset.
     */
    private double evaluateKernelForLambda(int lambda, int llFL, int ulFL) {
      double maxKernel = Double.MIN_VALUE;

      double gaussianBiasMagnitude = hd.biasMag;
      if (hd.biasMag <= 0.0) {
        gaussianBiasMagnitude = PSEUDO_BIAS_MAGNITUDE;
      }

      for (int jFL = llFL; jFL <= ulFL; jFL++) {
        double value = evaluateKernel(lambda, jFL, gaussianBiasMagnitude);
        kernelValues[jFL] = value;
        if (value > maxKernel) {
          maxKernel = value;
        }
      }

      // Only offset the kernel values for use with the 2D OST bias.
      double offset = 0.0;
      if (hd.biasMag > 0.0 && !hd.metaDynamics) {
        offset = -maxKernel;
        for (int jFL = llFL; jFL <= ulFL; jFL++) {
          kernelValues[jFL] += offset;
        }
      }
      return offset;
    }

    /**
     * Mirror the lambda bin if its less < 0 or greater than the last bin.
     *
     * @param bin Lambda bin to mirror.
     * @return The mirrored lambda bin.
     */
    private int lambdaMirror(int bin) {
      if (bin < 0) {
        return -bin;
      }
      int lastBin = hd.lambdaBins - 1;
      if (bin > lastBin) {
        // Number of bins past the last bin.
        bin -= lastBin;
        // Return Mirror bin
        return lastBin - bin;
      }
      // No mirror condition.
      return bin;
    }

    /**
     * For continuous lambda, the width of the first and last bins is dLambda_2, so the mirror
     * condition is to double their counts.
     *
     * @param bin Current lambda bin.
     * @return The mirror factor (either 1.0 or 2.0).
     */
    private double mirrorFactor(int bin) {
      if (hd.discreteLambda) {
        return 1.0;
      }
      if (bin == 0 || bin == hd.lambdaBins - 1) {
        return 2.0;
      }
      return 1.0;
    }

    /**
     * Compute the value of dU/dL for the given Histogram bin.
     *
     * @param dUdLBin The bin index in the dU/dL dimension.
     * @return The value of dU/dL at the center of the bin.
     */
    private double dUdLforBin(int dUdLBin) {
      return hd.dUdLMinimum + dUdLBin * hd.dUdLBinWidth + hd.dUdLBinWidth_2;
    }

    /**
     * Evaluate the bias at [currentLambdaBin, cF_lambda]
     */
    private double evaluateKernel(int currentLambdaBin, int currentdUdLBin, double gaussianBiasMagnitude) {
      // Compute the value of L and FL for the center of the current bin.
      double currentLambda = currentLambdaBin * hd.lambdaBinWidth;
      double currentdUdL = dUdLforBin(currentdUdLBin);

      // Variances are only used when dividing by twice their value.
      double invLs2 = 0.5 / hd.lambdaVariance;
      double invFLs2 = 0.5 / hd.dUdLVariance;

      double sum = 0.0;
      int lambdaBiasCutoff = hd.lambdaBiasCutoff;
      for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
        int lambdaBin = currentLambdaBin + iL;
        double deltaL = currentLambda - lambdaBin * hd.lambdaBinWidth;
        double deltaL2 = deltaL * deltaL;

        // Pre-compute the lambda bias.
        double L2exp = exp(-deltaL2 * invLs2);

        // Mirror condition for Lambda counts.
        lambdaBin = lambdaMirror(lambdaBin);
        double mirrorFactor = mirrorFactor(lambdaBin);

        int dUdLbiasCutoff = hd.dUdLBiasCutoff;
        for (int jFL = -dUdLbiasCutoff; jFL <= dUdLbiasCutoff; jFL++) {
          int dUdLBin = currentdUdLBin + jFL;
          double weight = mirrorFactor * getRecursionKernelValue(lambdaBin, dUdLBin);
          if (weight <= 0.0) {
            continue;
          }
          double deltaFL = currentdUdL - dUdLforBin(dUdLBin);
          double deltaFL2 = deltaFL * deltaFL;
          double e = weight * L2exp * exp(-deltaFL2 * invFLs2);
          sum += e;
        }
      }
      return gaussianBiasMagnitude * sum;
    }

    /**
     * Compute the total Bias energy at (currentLambda, currentdUdL).
     *
     * @param currentLambda The value of lambda.
     * @param currentdUdL   The value of dU/dL.
     * @return The bias energy.
     */
    double computeBiasEnergy(double currentLambda, double currentdUdL) {
      synchronized (this) {
        int currentLambdaBin = indexForLambda(currentLambda);
        int currentdUdLBin = binFordUdL(currentdUdL);

        double bias1D;
        double bias2D = 0.0;

        if (!hd.metaDynamics) {
          if (hd.biasMag > 0.0) {
            int lambdaBiasCutoff = hd.lambdaBiasCutoff;
            for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {

              int lambdaBin = currentLambdaBin + iL;
              double deltaL = currentLambda - (lambdaBin * hd.lambdaBinWidth);
              double deltaL2 = deltaL * deltaL;
              double expL2 = exp(-deltaL2 / (2.0 * hd.lambdaVariance));

              // Mirror conditions for recursion kernel counts.
              lambdaBin = lambdaMirror(lambdaBin);
              double mirrorFactor = mirrorFactor(lambdaBin);

              int dUdLbiasCutoff = hd.dUdLBiasCutoff;
              for (int iFL = -dUdLbiasCutoff; iFL <= dUdLbiasCutoff; iFL++) {
                int dUdLBin = currentdUdLBin + iFL;

                double weight = mirrorFactor * getRecursionKernelValue(lambdaBin, dUdLBin);
                if (weight <= 0.0) {
                  continue;
                }

                double deltaFL = currentdUdL - dUdLforBin(dUdLBin);
                double deltaFL2 = deltaFL * deltaFL;
                double bias = weight * hd.biasMag * expL2 * exp(-deltaFL2 / (2.0 * hd.dUdLVariance));
                bias2D += bias;
              }
            }
          }
          // Compute the energy for the recursion worker at F(L) using interpolation.
          bias1D = energyAndGradient1D(currentLambda, false);
        } else {
          bias1D = energyAndGradientMeta(currentLambda, false);
        }

        // Return the total bias.
        return bias1D + bias2D;
      }
    }

    /**
     * Compute the total Bias energy at (currentLambda, currentdUdL) and its gradient.
     *
     * @param currentdUdLambda The value of dU/dL.
     * @param chainRule        The chain rule terms for the bias.
     * @return The bias energy.
     */
    double energyAndGradient2D(double currentdUdLambda, double[] chainRule) {
      return energyAndGradient2D(ld.lambda, currentdUdLambda, chainRule, hd.biasMag);
    }

    /**
     * Compute the total Bias energy at (currentLambda, currentdUdL) and its gradient.
     *
     * @param currentLambda         The value of lambda.
     * @param currentdUdLambda      The value of dU/dL.
     * @param chainRule             The chain rule terms for the bias.
     * @param gaussianBiasMagnitude The magnitude of the Gaussian bias.
     * @return The bias energy.
     */
    double energyAndGradient2D(double currentLambda, double currentdUdLambda, double[] chainRule,
                               double gaussianBiasMagnitude) {

      if (gaussianBiasMagnitude <= 0.0) {
        chainRule[0] = 0.0;
        chainRule[1] = 0.0;
        return 0;
      }

      // Calculate recursion kernel G(L, F_L) and its derivatives with respect to L and F_L.
      double gLdEdL = 0.0;
      double dGdLambda = 0.0;
      double dGdFLambda = 0.0;
      int currentLambdaBin = indexForLambda(currentLambda);
      int currentdUdLBin = binFordUdL(currentdUdLambda);

      int lambdaBiasCutoff = hd.lambdaBiasCutoff;
      for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
        int lambdaBin = currentLambdaBin + iL;

        // Pass in the bin offset to the weight instance.
        double deltaL = currentLambda - (lambdaBin * hd.lambdaBinWidth);
        double deltaL2 = deltaL * deltaL;
        double expL2 = exp(-deltaL2 / (2.0 * hd.lambdaVariance));

        // Mirror conditions for recursion kernel counts.
        double mirrorFactor = mirrorFactor(lambdaBin);
        int dUdLbiasCutoff = hd.dUdLBiasCutoff;
        for (int iFL = -dUdLbiasCutoff; iFL <= dUdLbiasCutoff; iFL++) {
          int dUdLBin = currentdUdLBin + iFL;

          double weight = mirrorFactor * getRecursionKernelValue(lambdaBin, dUdLBin);
          if (weight <= 0.0) {
            continue;
          }

          double deltaFL = currentdUdLambda - dUdLforBin(dUdLBin);
          double deltaFL2 = deltaFL * deltaFL;
          double bias = weight * gaussianBiasMagnitude * expL2 * exp(-deltaFL2 / (2.0 * hd.dUdLVariance));
          gLdEdL += bias;
          dGdLambda -= deltaL / hd.lambdaVariance * bias;
          dGdFLambda -= deltaFL / hd.dUdLVariance * bias;
        }
      }

      chainRule[0] = dGdLambda;
      chainRule[1] = dGdFLambda;
      return gLdEdL;
    }


    /**
     * This calculates the 1D OST bias and its derivative with respect to Lambda.
     * <p>
     * See Equation 8 in <a href="http://doi.org/10.1021/ct300035u">
     * The Structure, Thermodynamics, and Solubility of Organic Crystals from Simulation with a
     * Polarizable Force Field</a>.
     * <p>
     *
     * @param gradient If true, compute the gradient.
     * @return a double.
     */
    private double energyAndGradient1D(boolean gradient) {
      return energyAndGradient1D(ld.lambda, gradient);
    }

    /**
     * This calculates the 1D OST bias and its derivative with respect to Lambda.
     * <p>
     * See Equation 8 in <a href="http://doi.org/10.1021/ct300035u">
     * The Structure, Thermodynamics, and Solubility of Organic Crystals from Simulation with a
     * Polarizable Force Field</a>.
     * <p>
     *
     * @param currentLambda The current value of lambda.
     * @param gradient      If true, compute the gradient.
     * @return a double.
     */
    private double energyAndGradient1D(double currentLambda, boolean gradient) {
      synchronized (this) {
        double biasEnergy = 0.0;
        for (int iL0 = 0; iL0 < hd.lambdaBins - 1; iL0++) {
          int iL1 = iL0 + 1;

          // Find bin centers and values for interpolation / extrapolation points.
          double L0 = iL0 * hd.lambdaBinWidth;
          double L1 = L0 + hd.lambdaBinWidth;
          double FL0 = ensembleAveragedUdL[iL0];
          double FL1 = ensembleAveragedUdL[iL1];
          double deltaFL = FL1 - FL0;
          /*
           If the lambda is less than or equal to the upper limit, this is
           the final interval. Set the upper limit to L, compute the partial
           derivative and break.
          */
          boolean done = false;
          if (currentLambda <= L1) {
            done = true;
            L1 = currentLambda;
          }

          // Upper limit - lower limit of the integral of the extrapolation / interpolation.
          biasEnergy += (FL0 * L1 + deltaFL * L1 * (0.5 * L1 - L0) / hd.lambdaBinWidth);
          biasEnergy -= (FL0 * L0 + deltaFL * L0 * (-0.5 * L0) / hd.lambdaBinWidth);
          if (done) {
            // Compute the gradient d F(L) / dL at L.
            if (gradient) {
              dUdLambda -= FL0 + (L1 - L0) * deltaFL / hd.lambdaBinWidth;
            }
            break;
          }
        }
        return -biasEnergy;
      }
    }

    /**
     * This calculates the 1D Metadynamics bias and its derivative with respect to Lambda.
     *
     * @param gradient If true, compute the gradient.
     * @return a double.
     */
    private double energyAndGradientMeta(boolean gradient) {
      return energyAndGradientMeta(ld.lambda, gradient);
    }

    /**
     * This calculates the 1D Metadynamics bias and its derivative with respect to Lambda.
     *
     * @param currentLambda The current value of lambda.
     * @param gradient      If true, compute the gradient.
     * @return a double.
     */
    private double energyAndGradientMeta(double currentLambda, boolean gradient) {
      // Zero out the metadynamics bias energy.
      double biasEnergy = 0.0;

      // Current lambda bin.
      int currentBin = indexForLambda(currentLambda);

      // Loop over bins within the cutoff.
      int lambdaBiasCutoff = hd.lambdaBiasCutoff;
      for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
        int lambdaBin = currentBin + iL;
        double deltaL = currentLambda - (lambdaBin * hd.lambdaBinWidth);
        double deltaL2 = deltaL * deltaL;
        // Mirror conditions for the lambda bin and count magnitude.
        lambdaBin = lambdaMirror(lambdaBin);
        double mirrorFactor = mirrorFactor(lambdaBin);
        double weight = mirrorFactor * hd.biasMag * countsForLambda(lambdaBin);
        if (weight > 0) {
          double e = weight * exp(-deltaL2 / (2.0 * hd.lambdaVariance));
          biasEnergy += e;
          // Add the dMeta/dL contribution.
          if (gradient) {
            dUdLambda -= deltaL / hd.lambdaVariance * e;
          }
        }
      }

      return biasEnergy;
    }

    /**
     * Marginalize over dU/dL counts for the given lambda bin.
     *
     * @param lambdaBin Lambda bin to marginalize for.
     * @return Total number of counts.
     */
    private double countsForLambda(int lambdaBin) {
      synchronized (this) {
        double count = 0.0;
        for (int i = 0; i < hd.dUdLBins; i++) {
          count += hd.zHistogram[lambdaBin][i];
        }
        return count;
      }
    }

    /**
     * If necessary, allocate more space.
     */
    void checkRecursionKernelSize(double currentdUdL) {
      // Synchronize updates of the recursion kernel.
      synchronized (this) {
        if (currentdUdL > hd.dUdLMaximum) {
          logger.info(format(" Current F_lambda %8.2f > maximum histogram size %8.2f.", currentdUdL, hd.dUdLMaximum));

          double origDeltaG = updateFreeEnergyDifference(false, false);

          int newdUdLBins = hd.dUdLBins;
          while (hd.dUdLMinimum + newdUdLBins * hd.dUdLBinWidth < currentdUdL) {
            newdUdLBins += 100;
          }
          double[][] newRecursionKernel = new double[hd.lambdaBins][newdUdLBins];

          // We have added bins above the indices of the current counts just copy them into the new array.
          for (int i = 0; i < hd.lambdaBins; i++) {
            arraycopy(hd.zHistogram[i], 0, newRecursionKernel[i], 0, hd.dUdLBins);
          }

          hd.zHistogram = newRecursionKernel;
          hd.dUdLBins = newdUdLBins;
          kernelValues = new double[hd.dUdLBins];
          hd.dUdLMaximum = hd.dUdLMinimum + hd.dUdLBinWidth * hd.dUdLBins;
          logger.info(format(" New histogram %8.2f to %8.2f with %d bins.\n", hd.dUdLMinimum, hd.dUdLMaximum, hd.dUdLBins));

          double newFreeEnergy = updateFreeEnergyDifference(true, false);
          assert (origDeltaG == newFreeEnergy);
        }
        if (currentdUdL < hd.dUdLMinimum) {
          logger.info(format(" Current F_lambda %8.2f < minimum histogram size %8.2f.", currentdUdL, hd.dUdLMinimum));

          double origDeltaG = updateFreeEnergyDifference(false, false);

          int offset = 100;
          while (currentdUdL < hd.dUdLMinimum - offset * hd.dUdLBinWidth) {
            offset += 100;
          }
          int newdUdLBins = hd.dUdLBins + offset;
          double[][] newRecursionKernel = new double[hd.lambdaBins][newdUdLBins];

          // We have added bins below the current counts,
          // so their indices must be increased by: offset = newFLBins - FLBins
          for (int i = 0; i < hd.lambdaBins; i++) {
            arraycopy(hd.zHistogram[i], 0, newRecursionKernel[i], offset, hd.dUdLBins);
          }

          hd.zHistogram = newRecursionKernel;
          hd.dUdLMinimum = hd.dUdLMinimum - offset * hd.dUdLBinWidth;
          hd.dUdLBins = newdUdLBins;
          kernelValues = new double[hd.dUdLBins];
          logger.info(format(" New histogram %8.2f to %8.2f with %d bins.\n", hd.dUdLMinimum, hd.dUdLMaximum, hd.dUdLBins));

          double newFreeEnergy = updateFreeEnergyDifference(true, false);
          assert (origDeltaG == newFreeEnergy);
        }
      }
    }

    /**
     * Add a Gaussian hill to the Histogram at (lambda, dUdL).
     *
     * @param dUdL The value of dU/dL.
     */
    void addBias(double dUdL) {
      // Communicate adding the bias to all walkers.
      if (hd.asynchronous) {
        sendAsynchronous.send(ld.lambda, dUdL, temperingWeight);
      } else {
        sendSynchronous.send(ld.lambda, dUdL, temperingWeight);
      }
    }

    /**
     * Pre-force portion of the stochastic integrator.
     */
    private void stochasticPreForce() {
      // Update theta.
      double[] thetaPosition = lambdaState.x();
      thetaPosition[0] = ld.theta;
      stochastic.preForce(OrthogonalSpaceTempering.this);
      ld.theta = thetaPosition[0];

      // Maintain theta in the interval -PI to PI.
      if (ld.theta > PI) {
        ld.theta -= 2.0 * PI;
      } else if (ld.theta <= -PI) {
        ld.theta += 2.0 * PI;
      }

      // Update lambda as sin(theta)^2.
      double sinTheta = sin(ld.theta);
      ld.lambda = sinTheta * sinTheta;
      lambdaInterface.setLambda(ld.lambda);
    }

    /**
     * Post-force portion of the stochastic integrator.
     */
    private void stochasticPostForce() {
      double[] thetaPosition = lambdaState.x();
      double[] gradient = lambdaState.gradient();
      // Update theta and dU/dtheta.
      gradient[0] = dUdLambda * sin(2 * ld.theta);
      thetaPosition[0] = ld.theta;
      stochastic.postForce(gradient);

      // Load the new theta, velocity and acceleration.
      double[] thetaVelocity = lambdaState.v();
      double[] thetaAcceleration = lambdaState.a();
      ld.theta = thetaPosition[0];
      ld.thetaVelocity = thetaVelocity[0];
      ld.thetaAcceleration = thetaAcceleration[0];

      // Maintain theta in the interval -PI to PI.
      if (ld.theta > PI) {
        ld.theta -= 2.0 * PI;
      } else if (ld.theta <= -PI) {
        ld.theta += 2.0 * PI;
      }

      // Update lambda as sin(theta)^2.
      double sinTheta = sin(ld.theta);
      ld.lambda = sinTheta * sinTheta;
      lambdaInterface.setLambda(ld.lambda);
    }

    /**
     * Gets the last lambda value received by this Histogram. This can be out-of-date w.r.t. the
     * OST's current lambda!
     *
     * @return Lambda value of the last bias added to this Histogram.
     */
    double getLastReceivedLambda() {
      return lastReceivedLambda;
    }

    public void setLastReceivedLambda(double lastReceivedLambda) {
      this.lastReceivedLambda = lastReceivedLambda;
    }

    /**
     * Gets the last dU/dL value received by this Histogram. This can be out-of-date w.r.t. the OST's
     * current dU/dL!
     *
     * @return dU/dL value of the last bias added to this Histogram.
     */
    double getLastReceivedDUDL() {
      return lastReceiveddUdL;
    }

    void destroy() {
      if (sendAsynchronous != null && sendAsynchronous.isAlive()) {
        double[] killMessage = new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        DoubleBuf killBuf = DoubleBuf.buffer(killMessage);
        try {
          logger.fine(" Sending the termination message.");
          world.send(rank, killBuf);
          logger.fine(" Termination message was sent successfully.");
          logger.fine(format(" Receive thread alive %b status %s", sendAsynchronous.isAlive(),
              sendAsynchronous.getState()));
        } catch (Exception ex) {
          String message = format(" Asynchronous Multiwalker OST termination signal "
              + "failed to be sent for process %d.", rank);
          logger.log(Level.SEVERE, message, ex);
        }
      } else {
        logger.fine(
            " CountReceiveThread was either not initialized, or is not alive. This is the case for the Histogram script.");
      }
    }

  }

}