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
import static org.apache.commons.math3.util.FastMath.round;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.Barostat;
import ffx.algorithms.optimize.Minimize;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.numerics.integrate.DataSet;
import ffx.numerics.integrate.DoublesDataSet;
import ffx.numerics.integrate.Integrate1DNumeric;
import ffx.numerics.integrate.Integrate1DNumeric.IntegrationType;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.utils.EnergyException;
import ffx.utilities.Constants;
import ffx.utilities.FileUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

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
  /** The potential energy of the system. */
  protected final CrystalPotential potential;
  /** Reference to the Barostat in use; if present this must be turned off during optimization. */
  protected final Barostat barostat;
  /** The AlgorithmListener is called each time a count is added. */
  protected final AlgorithmListener algorithmListener;
  /** Print detailed energy information. */
  protected final boolean print = false;
  /** Number of variables. */
  protected final int nVariables;
  /** A potential energy that implements the LambdaInterface. */
  private final LambdaInterface lambdaInterface;
  /** List of additional Histograms this OST can switch to. */
  private final List<Histogram> allHistograms = new ArrayList<>();
  /** Parameters to control saving local optimizations. */
  private final OptimizationParameters optimizationParameters;
  /** Each walker has a unique lambda restart file. */
  private final File lambdaFile;
  /** Write out structures only for lambda values greater than or equal to this threshold. */
  private final double lambdaWriteOut;
  /** Properties. */
  private final CompositeConfiguration properties;
  /** Temperature in Kelvins. */
  private final double temperature;
  /** Timestep. */
  private final double dt;
  /** Whether to use asynchronous communications. */
  private final boolean asynchronous;
  /** The MolecularAssembly being simulated. */
  protected MolecularAssembly molecularAssembly;
  /** State variable lambda ranges from 0.0 .. 1.0. */
  protected double lambda;
  /**
   * Are FAST varying energy terms being computed, SLOW varying energy terms, or BOTH. OST is not
   * active when only FAST varying energy terms are being propagated.
   */
  protected Potential.STATE state = Potential.STATE.BOTH;
  /** Force Field Potential Energy (i.e. with no bias terms added). */
  protected double forceFieldEnergy;
  /**
   * Interval between writing an OST restart file in steps.
   *
   * <p>The default saveFrequency = 1000.
   */
  int saveFrequency;
  /**
   * Number of times the OST biasing potential has been evaluated with the "propagateLambda" flag
   * true.
   */
  long energyCount;
  /** Contains counts for the OST bias. */
  private Histogram histogram;
  /** Index of the current Histogram. */
  private int histogramIndex;
  /** Flag to indicate that the Lambda particle should be propagated. */
  private boolean propagateLambda = true;
  /**
   * Interval between printing information on the lambda particle in steps.
   *
   * <p>The default printFrequency = 100.
   */
  private int printFrequency;
  /** Mixed second partial derivative with respect to coordinates and lambda. */
  private final double[] dUdXdL;
  /** Partial derivative of the force field energy with respect to lambda. */
  private double dForceFieldEnergydL;
  /** Magnitude of the 2D orthogonal space bias G(L,dE/dL). */
  private double gLdEdL = 0.0;
  /** OST Bias energy. */
  private double biasEnergy;
  /** Total system energy. */
  private double totalEnergy;
  /** Total partial derivative of the potential (U) being sampled with respect to lambda. */
  private double dUdLambda;
  /** Second partial derivative of the potential being sampled with respect to lambda. */
  private double d2UdL2;
  /**
   * Save the previous free energy, in order to limit logging to time points where the free energy
   * has changed.
   */
  private double previousFreeEnergy = 0.0;
  /** If true, values of (lambda, dU/dL) that have not been observed are rejected. */
  private boolean hardWallConstraint = false;
  /**
   * If a 2D bias is in use (i.e. lambda-bias-cutoff > 0), do not normalize bias height.
   */
  private static final double TWO_D_NORMALIZATION = 1.0;
  /**
   * If a 1D bias is in use (i.e. lambda-bias-cutoff == 0), normalize bias height to deposit the same volume of bias.
   */
  private static final double ONE_D_NORMALIZATION = Math.sqrt(2.0 * Math.PI);

  /**
   * OST Constructor.
   *
   * @param lambdaInterface defines Lambda and dU/dL.
   * @param potential defines the Potential energy.
   * @param lambdaFile contains the current Lambda particle position and velocity.
   * @param histoSettings contains histogram-centric options.
   * @param properties defines System properties.
   * @param temperature the simulation temperature.
   * @param dt the time step in femtoseconds.
   * @param printInterval number of steps between logging updates.
   * @param saveInterval number of steps between restart file updates.
   * @param asynchronous set to true if walkers run asynchronously.
   * @param resetNumSteps whether to reset energy counts to 0
   * @param algorithmListener the AlgorithmListener to be notified of progress.
   */
  public OrthogonalSpaceTempering(
      LambdaInterface lambdaInterface,
      CrystalPotential potential,
      File lambdaFile,
      HistogramSettings histoSettings,
      CompositeConfiguration properties,
      double temperature,
      double dt,
      double printInterval,
      double saveInterval,
      boolean asynchronous,
      boolean resetNumSteps,
      AlgorithmListener algorithmListener) {
    this(
        lambdaInterface,
        potential,
        lambdaFile,
        histoSettings,
        properties,
        temperature,
        dt,
        printInterval,
        saveInterval,
        asynchronous,
        resetNumSteps,
        algorithmListener,
        0.0);
  }

  /**
   * OST Constructor.
   *
   * @param lambdaInterface defines Lambda and dU/dL.
   * @param potential defines the Potential energy.
   * @param lambdaFile contains the current Lambda particle position and velocity.
   * @param histoSettings contains histogram-centric options.
   * @param properties defines System properties.
   * @param temperature the simulation temperature.
   * @param dt the time step in femtoseconds.
   * @param printInterval number of steps between logging updates.
   * @param saveInterval number of steps between restart file updates.
   * @param asynchronous set to true if walkers run asynchronously.
   * @param resetNumSteps whether to reset energy counts to 0
   * @param algorithmListener the AlgorithmListener to be notified of progress.
   * @param lambdaWriteOut Minimum lambda value to print out snapshots.
   */
  public OrthogonalSpaceTempering(
      LambdaInterface lambdaInterface,
      CrystalPotential potential,
      File lambdaFile,
      HistogramSettings histoSettings,
      CompositeConfiguration properties,
      double temperature,
      double dt,
      double printInterval,
      double saveInterval,
      boolean asynchronous,
      boolean resetNumSteps,
      AlgorithmListener algorithmListener,
      double lambdaWriteOut) {
    this.lambdaInterface = lambdaInterface;
    this.potential = potential;
    this.lambdaFile = lambdaFile;
    this.algorithmListener = algorithmListener;
    nVariables = potential.getNumberOfVariables();

    if (potential instanceof Barostat) {
      barostat = (Barostat) potential;
    } else {
      barostat = null;
    }

    this.lambdaWriteOut = lambdaWriteOut;

    // Convert the time step to picoseconds.
    dt *= Constants.FSEC_TO_PSEC;

    // Convert the print interval to a print frequency.
    printFrequency = 100;
    if (printInterval >= dt) {
      printFrequency = (int) (printInterval / dt);
    }

    // Convert the save interval to a save frequency.
    saveFrequency = 1000;
    if (saveInterval >= dt) {
      saveFrequency = (int) (saveInterval / dt);
    }

    energyCount = -1;
    dUdXdL = new double[nVariables];

    // Init the Histogram and read a restart file if it exists.
    this.properties = properties;
    this.temperature = temperature;
    this.dt = dt;
    this.asynchronous = asynchronous;
    histoSettings.temperature = temperature;
    histoSettings.dt = dt;
    histoSettings.asynchronous = asynchronous;
    histogram = new Histogram(properties, histoSettings);
    histogramIndex = 0;
    allHistograms.add(histogram);

    // Load the OST lambda restart file if it exists.
    if (lambdaFile != null && lambdaFile.exists()) {
      try {
        LambdaReader lambdaReader = new LambdaReader(new FileReader(lambdaFile));
        lambdaReader.readLambdaFile(resetNumSteps);
        lambdaReader.setVariables(this);
        lambdaReader.close();
        logger.info(format("\n Continuing OST lambda from %s.", histogram.lambdaFileName));
      } catch (FileNotFoundException ex) {
        logger.info(" Lambda restart file could not be found and will be ignored.");
      } catch (IOException ioe) {
        logger.warning(" Could not close lambda restart file reader!");
      }
    }

    // Configure optimization parameters.
    optimizationParameters = new OptimizationParameters(properties);

    // Log parameters.
    logger.info("\n Orthogonal Space Random Walk Parameters");
    logger.info(format("  Gaussian Bias Magnitude:        %6.4f (kcal/mol)", histogram.biasMag));
    logger.info(format("  Gaussian dU/dL Bias Cutoff:     %6d bins", histogram.biasCutoff));
    logger.info(format("  Gaussian Lambda Bias Cutoff:    %6d bins", histogram.lambdaBiasCutoff));
    logger.info(format("  Print Interval:                 %6.3f psec", printInterval));
    logger.info(format("  Save Interval:                  %6.3f psec", saveInterval));
  }

  /**
   * Add an alternate Histogram this OST can use.
   *
   * @param settings Settings to use for the new Histogram.
   */
  public void addHistogram(HistogramSettings settings) {
    Histogram newHisto = new Histogram(properties, settings);
    settings.temperature = temperature;
    settings.dt = dt;
    settings.asynchronous = asynchronous;
    allHistograms.add(newHisto);
  }

  /** {@inheritDoc} */
  @Override
  public boolean dEdLZeroAtEnds() {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public boolean destroy() {
    // Shut down the CountReceiveThread.
    histogram.destroy();
    return potential.destroy();
  }

  /** Compute the force field + bias energy. */
  public double energy(double[] x) {

    forceFieldEnergy = potential.energy(x);

    // OST is propagated with the slowly varying terms.
    if (state == Potential.STATE.FAST) {
      return forceFieldEnergy;
    }

    dForceFieldEnergydL = lambdaInterface.getdEdL();
    d2UdL2 = lambdaInterface.getd2EdL2();
    int lambdaBin = histogram.indexForLambda(lambda);
    dUdLambda = dForceFieldEnergydL;
    gLdEdL = 0.0;

    // Calculate recursion kernel G(L, F_L) and its derivatives with respect to L and F_L.
    if (histogram.biasMag > 0.0) {
      double[] chainRule = new double[2];
      gLdEdL = histogram.energyAndGradient2D(lambda, dUdLambda, chainRule, histogram.biasMag);
      double dGdLambda = chainRule[0];
      double dGdFLambda = chainRule[1];
      dUdLambda += dGdLambda + dGdFLambda * d2UdL2;
    }

    // Compute the energy and gradient for the recursion worker at F(L) using interpolation.
    double bias1D = histogram.energyAndGradient1D(lambda, true);

    // The total bias energy is the sum of the 1D and 2D terms.
    biasEnergy = bias1D + gLdEdL;

    if (print) {
      logger.info(format(" Bias Energy        %16.8f", biasEnergy));
      logger.info(
          format(" %s %16.8f  (Kcal/mole)", "OST Potential    ", forceFieldEnergy + biasEnergy));
    }

    if (propagateLambda) {
      energyCount++;

      // Log the current Lambda state.
      if (energyCount % printFrequency == 0) {
        double dBdL = dUdLambda - dForceFieldEnergydL;
        if (histogram.lambdaBins < 1000) {
          logger.info(
              format(
                  " L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                  lambda,
                  lambdaBin,
                  dForceFieldEnergydL,
                  dBdL,
                  dUdLambda,
                  histogram.halfThetaVelocity));
        } else {
          logger.info(
              format(
                  " L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                  lambda,
                  lambdaBin,
                  dForceFieldEnergydL,
                  dBdL,
                  dUdLambda,
                  histogram.halfThetaVelocity));
        }
      }

      // Metadynamics grid counts (every 'countInterval' steps).
      if (energyCount % histogram.countInterval == 0) {
        histogram.addBias(dForceFieldEnergydL, x, null);
      }

      histogram.langevin();
    }

    totalEnergy = forceFieldEnergy + biasEnergy;

    return totalEnergy;
  }

  /** {@inheritDoc} */
  @Override
  public double energyAndGradient(double[] x, double[] gradient) {

    forceFieldEnergy = potential.energyAndGradient(x, gradient);

    // OST is propagated with the slowly varying terms.
    if (state == STATE.FAST) {
      return forceFieldEnergy;
    }

    dForceFieldEnergydL = lambdaInterface.getdEdL();
    dUdLambda = dForceFieldEnergydL;
    d2UdL2 = lambdaInterface.getd2EdL2();
    gLdEdL = 0.0;

    if (histogram.biasMag > 0) {
      double[] chainRule = new double[2];
      gLdEdL = histogram.energyAndGradient2D(lambda, dUdLambda, chainRule, histogram.biasMag);
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
    double bias1D = histogram.energyAndGradient1D(lambda, true);

    // The total bias is the sum of 1D and 2D terms.
    biasEnergy = bias1D + gLdEdL;

    if (print) {
      logger.info(format(" %s %16.8f", "Bias Energy       ", biasEnergy));
      logger.info(
          format(
              " %s %16.8f  %s", "OST Potential    ", forceFieldEnergy + biasEnergy, "(Kcal/mole)"));
    }

    if (propagateLambda) {
      energyCount++;

      // Log the current Lambda state.
      if (energyCount % printFrequency == 0) {
        double dBdL = dUdLambda - dForceFieldEnergydL;
        int lambdaBin = histogram.indexForLambda(lambda);
        if (histogram.lambdaBins < 1000) {
          logger.info(
              format(
                  " L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                  lambda,
                  lambdaBin,
                  dForceFieldEnergydL,
                  dBdL,
                  dUdLambda,
                  histogram.halfThetaVelocity));
        } else {
          logger.info(
              format(
                  " L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                  lambda,
                  lambdaBin,
                  dForceFieldEnergydL,
                  dBdL,
                  dUdLambda,
                  histogram.halfThetaVelocity));
        }
      }

      // Metadynamics grid counts (every 'countInterval' steps).
      if (energyCount % histogram.countInterval == 0) {
        histogram.addBias(dForceFieldEnergydL, x, gradient);
      }

      histogram.langevin();
    }

    totalEnergy = forceFieldEnergy + biasEnergy;

    return totalEnergy;
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
  @Override
  public double[] getCoordinates(double[] doubles) {
    return potential.getCoordinates(doubles);
  }

  /** {@inheritDoc} */
  @Override
  public Crystal getCrystal() {
    return potential.getCrystal();
  }

  /** {@inheritDoc} */
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
    return energyCount;
  }

  /**
   * Set the number of counts.
   *
   * @param counts Number of counts.
   */
  public void setEnergyCount(long counts) {
    this.energyCount = counts;
  }

  /** {@inheritDoc} */
  @Override
  public Potential.STATE getEnergyTermState() {
    return state;
  }

  /** {@inheritDoc} */
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
   * Return a 2D Histogram of counts.
   *
   * @param index Index of the Histogram to fetch.
   * @return the index'th Histogram.
   */
  public Histogram getHistogram(int index) {
    return allHistograms.get(index);
  }

  /**
   * Getter for the field <code>lambda</code>.
   *
   * @return a double.
   */
  public double getLambda() {
    return lambda;
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
    this.lambda = lambda;
    histogram.theta = asin(sqrt(lambda));
  }

  /**
   * Getter for the field <code>lambdaInterface</code>.
   *
   * @return a {@link ffx.potential.bonded.LambdaInterface} object.
   */
  public LambdaInterface getLambdaInterface() {
    return lambdaInterface;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getMass() {
    return potential.getMass();
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    return potential.getPreviousAcceleration(previousAcceleration);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getScaling() {
    return potential.getScaling();
  }

  /** {@inheritDoc} */
  @Override
  public void setScaling(double[] scaling) {
    potential.setScaling(scaling);
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
  @Override
  public double[] getVelocity(double[] velocity) {
    return potential.getVelocity(velocity);
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    throw new UnsupportedOperationException(
        " Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {
    return getTotaldEdLambda();
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] gradient) {
    throw new UnsupportedOperationException(
        " Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
  }

  // TODO: Delete method when debugging of RepexOST is done.
  public void logOutputFiles(int index) {
    logger.info(
        format(
            " OST: Lambda file %s, histogram %s",
            histogram.lambdaFileName, allHistograms.get(index).hisFileName));
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
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

  /**
   * Switch to an alternate Histogram.
   *
   * @param histo Histogram to switch to.
   */
  public void switchHistogram(Histogram histo) {
    if (allHistograms.contains(histo)) {
      histogram = histo;
      histogramIndex = allHistograms.indexOf(histo);
    } else {
      logger.warning(
          " Likely unintended behavior: switching to a Histogram this OST does not yet know about!");
      histogram = histo;
      allHistograms.add(histogram);
      histogramIndex = allHistograms.size() - 1;
    }
  }

  @Override
  public void writeAdditionalRestartInfo(boolean recursive) {
    writeRestart();
    if (recursive) {
      potential.writeAdditionalRestartInfo(true);
    }
  }

  /** Write histogram and lambda restart files. */
  private void writeRestart() {
    if (algorithmListener != null) {
      algorithmListener.algorithmUpdate(molecularAssembly);
    }

    // Only the rank 0 process writes a histogram restart file.
    histogram.writeRestart();

    // All ranks write a lambda restart file.
    try {
      LambdaWriter lambdaWriter =
          new LambdaWriter(this, new BufferedWriter(new FileWriter(lambdaFile)));
      lambdaWriter.writeLambdaFile();
      lambdaWriter.flush();
      lambdaWriter.close();
      logger.info(format(" Wrote lambda restart file to %s.", histogram.lambdaFileName));
    } catch (IOException ex) {
      String message = format(" Exception writing lambda restart file %s.", lambdaFile);
      logger.log(Level.INFO, Utilities.stackTraceToString(ex));
      logger.log(Level.SEVERE, message, ex);
    }
  }

  double getLambdaWriteOut() {
    return lambdaWriteOut;
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
   * Returns basic information about the current histogram.
   *
   * @return Tempering parameter, threshold, and bias magnitude.
   */
  String histoInfo() {
    return histoInfo(histogramIndex);
  }

  /**
   * Returns basic information about the histogram.
   *
   * @param index Index of the histogram to print information about.
   * @return Tempering parameter, threshold, and bias magnitude.
   */
  String histoInfo(int index) {
    Histogram histo = allHistograms.get(index);
    return format(
        " Tempering parameter: %.5f kBT. Threshold: %.5f kcal/mol. Bias magnitude: %.5f kcal/mol",
        histo.temperingFactor, histo.temperOffset, histo.biasMag);
  }

  /**
   * If the dUdLHardWall flag is set to true, this method will return false if the (lambda, dU/dL)
   * sample is has not been seen.
   *
   * @param lambda The proposed lambda value.
   * @param dUdL The proposed dU/dL value.
   * @return Returns false only if the dUdLHardWall flag is true, and the (lambda, dU/dL) sample has
   *     not been seen.
   */
  boolean insideHardWallConstraint(double lambda, double dUdL) {
    if (hardWallConstraint) {
      double weight = histogram.evaluateHistogram(lambda, dUdL);
      return weight > 0.0;
    }
    return true;
  }

  /** Parameters for running local optimizations during OST sampling. */
  public class OptimizationParameters {

    /**
     * Flag to turn on OST optimization.
     *
     * <p>The default doOptimization = false.
     */
    private boolean doOptimization = false;
    /** Reset unit cell parameters, molecular orientation and translation. */
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
    /** Holds the lowest potential energy coordinates. */
    private double[] optimumCoords;
    /** File instance used for saving optimized structures. */
    private File optimizationFile;
    /** SystemFilter used to save optimized structures. */
    private SystemFilter optimizationFilter;

    /** Empty constructor. */
    OptimizationParameters(CompositeConfiguration properties) {
      energyWindow = properties.getDouble("ost-opt-energy-window", 10.0);
      eps = properties.getDouble("ost-opt-eps", 0.1);
      tolerance = properties.getDouble("ost-opt-tolerance", 1.0e-8);
      frequency = properties.getInt("ost-opt-frequency", 10000);
      lambdaCutoff = properties.getDouble("ost-opt-lambda-cutoff", 0.8);
      doUnitCellReset = properties.getBoolean("ost-opt-unitcell-reset", false);

      logger.info("\n Optimization Parameters");
      logger.info(format("  Energy Window:                  %6.4f (kcal/mol)", energyWindow));
      logger.info(format("  EPS:                            %6.5f ", eps));
      logger.info(format("  Tolerance:                      %6.4f (kcal/mol)", tolerance));
      logger.info(format("  Frequency:                      %6d (steps)", frequency));
      logger.info(format("  Lambda Cutoff:                  %6.4f ", lambdaCutoff));
      logger.info(format("  Unit Cell Reset:                %B ", doUnitCellReset));
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
     * @param e Current energy.
     * @param x Current atomic coordinates.
     * @param gradient Work array for collecting the gradient.
     */
    public void optimize(double e, double[] x, double[] gradient) {

      // Return if the optimization flag is not set, or if lambda is not beyond the cutoff.
      if (doOptimization && lambda > lambdaCutoff && energyCount % frequency == 0) {
        if (gradient == null) {
          gradient = new double[x.length];
        }
      } else {
        return;
      }

      logger.info(format("\n OST Minimization (Step %d)", energyCount));

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
          optimizationFilter.writeFile(optimizationFile, false);
          Crystal uc = crystal.getUnitCell();
          logger.info(
              format(
                  " Minimum: %12.6f %s (%12.6f g/cc) optimized from %12.6f at step %d.",
                  minEnergy, uc.toShortString(), density, startingEnergy, energyCount));
        }
      } catch (EnergyException ex) {
        String message = ex.getMessage();
        logger.info(
            format(
                " Energy exception minimizing coordinates at lambda=%8.6f\n %s.", lambda, message));
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
        lambdaInterface.setLambda(0.0);
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
     * @param doOptimization a boolean.
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
          optimizationFile = new File(fileName + "_opt.xyz");
          optimizationFilter = new XYZFilter(optimizationFile, molecularAssembly, null, null);
        } else {
          optimizationFile = new File(fileName + "_opt.pdb");
          optimizationFilter = new PDBFilter(optimizationFile, molecularAssembly, null, null);
        }
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
     * If "realBiasMagnitude" is 0, temporarily set biasMag to this value to calculate the the
     * ensemble average dU/dL.
     *
     * <p>Any value that does not overflow/underflow double precision summations of the count matrix
     * will give identical results. For example, values of 1.0e-20, 1.0 and 1.0e20 were tested and
     * found to give identical results.
     *
     * <p>Thus, a value of 1.0 is good choice.
     */
    private static final double PSEUDO_BIAS_MAGNITUDE = 1.0;
    /**
     * Temperature in Kelvin.
     *
     * <p>The default is 298.15.
     */
    protected final double temperature;
    /** Time step in picoseconds. */
    protected final double dt;
    /** Parallel Java world communicator. */
    protected final Comm world;
    /** Rank of this process. */
    protected final int rank;
    /** If true, use discrete lambda values instead of continuous lambda values. */
    final boolean discreteLambda;
    /**
     * Width of a lambda bin, or the distance between discrete lambda values.
     *
     * <p>The default dL = (1.0 / (lambdaBins - 1).
     */
    final double dL;
    /** Half the width of a lambda bin, or zero for discrete lambda values. */
    final double dL_2;
    /**
     * The width of the FLambda bin.
     *
     * <p>The default dFL = 2.0 (kcal/mol).
     */
    final double dFL;
    /** Half the width of the F_lambda bin. */
    final double dFL_2;
    /**
     * For continuous lambda: The first Lambda bin is centered on 0.0 (-0.005 .. 0.005). The final
     * Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     *
     * <p>With this scheme, the maximum of biasing Gaussians is at the edges.
     *
     * <p>For discrete lambda: The first value of lambda is 0.0 and last value is 1.0.
     *
     * <p>The default lambdaBins = 201.
     */
    final int lambdaBins;
    /**
     * When evaluating the biasing potential, contributions from Gaussians centered on bins more the
     * "biasCutoff" away will be neglected.
     *
     * <p>The default biasCutoff = 5.
     */
    final int biasCutoff;
    /**
     * When evaluating the biasing potential, contributions from Gaussians centered on bins more the
     * "biasCutoff" away will be neglected.
     *
     * <p>The continuous lambda simulations, the default lambdaBiasCutoff = biasCutoff.
     * <p>For discrete lambda simulations, the default lambdaBiasCutoff = 0.
     */
    final int lambdaBiasCutoff;
    /** 1D PMF with respect to lambda F(L). */
    final double[] FLambda;
    /** Reasonable thetaFriction is ~60 per picosecond (1.0e-12). */
    final double thetaFriction;
    /** Reasonable thetaMass is ~100 a.m.u. (100 a.m.u is 1.6605e-22 grams). */
    final double thetaMass;
    /**
     * Interval between adding a count to the Recursion kernel in MD steps.
     *
     * <p>The default countInterval = 10.
     */
    final int countInterval;
    /**
     * Each walker reads the same histogram restart file. Only the walker of rank 0 writes the
     * histogram restart file.
     */
    private final File histogramFile;
    /**
     * The minimum value of the first lambda bin.
     *
     * <p>minLambda = -dL_2 for continuous lambda.
     *
     * <p>minLambda = 0 for discrete lambda.
     */
    private final double minLambda;
    /**
     * The Dama et al. transition-tempering rate parameter. A reasonable value is about 2 to 8 kT,
     * with larger values being resulting in slower decay.
     *
     * <p>The default temperingFactor = 2.0.
     */
    private final double temperingFactor;
    /**
     * This deltaT is used to determine the tempering weight as described below for the
     * temperingWeight variable.
     *
     * <p>deltaT = temperingFactor * kB * T.
     */
    private final double deltaT;
    /**
     * An offset applied before calculating tempering weight.
     *
     * <p>First, for every Lambda, we calculate the maximum bias at that lambda by searching all
     * populated dU/dL bins: maxdUdL(L) = max[ G(L,F_L) ] where the max operator searches over all
     * F_L bins.
     *
     * <p>Then, the minimum bias coverage is determined by searching the maxdUdL(L) array over
     * Lambda. minBias = min[ maxdUdL(L) ] where the min operator searches all entries in the array.
     *
     * <p>Then the temperOffset is applied to the minBias:
     *
     * <p>biasHeight = max[minBias - temperOffset, 0]
     *
     * <p>The default temperOffset = 1.0 kcal/mol.
     */
    private final double temperOffset;
    /**
     * Interval between how often the 1D histogram is printed to screen versus silently updated in
     * background.
     *
     * <p>The fLambdaPrintInterval is 25.
     */
    private final int fLambdaPrintInterval;
    /** The integration algorithm used for thermodynamic integration. */
    private final IntegrationType integrationType;
    /**
     * Random force conversion to kcal/mol/A; Units: Sqrt (4.184 Joule per calorie) / (nanometers per
     * meter)
     */
    private final double randomConvert = sqrt(4.184) / 10e9;
    /** randomConvert squared. Units: Joule per calorie / (nanometer per meter)^2 */
    private final double randomConvert2 = randomConvert * randomConvert;
    /** Random variable for stochastic lambda particle generation. */
    private final Random stochasticRandom;
    /** Once the lambda reset value is reached, OST statistics are reset. */
    private final double lambdaResetValue;

    private final boolean writeIndependent;
    private final boolean independentWalkers;
    /**
     * Flag to indicate if OST should send and receive counts between processes synchronously or
     * asynchronously. The latter can be faster by ~40% because simulation with Lambda &gt; 0.75 must
     * compute two condensed phase self-consistent fields to interpolate polarization.
     */
    private final boolean asynchronous;
    /** The CountReceiveThread accumulates OST statistics from multiple asynchronous walkers. */
    private final AsynchronousSend asynchronousSend;

    private final SynchronousSend synchronousSend;
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
     * The minimum value of the first F_lambda bin.
     *
     * <p>minFLambda = -(dFL * FLambdaBins) / 2.0.
     */
    double minFLambda;
    /**
     * It is useful to have an odd number of bins, so that there is a bin from FL=-dFL/2 to dFL/2 so
     * that as FL approaches zero its contribution to thermodynamic integration goes to zero.
     *
     * <p>Otherwise a contribution of zero from a L bin can only result from equal sampling of the
     * ranges -dFL to 0 and 0 to dFL.
     *
     * <p>The default FLambdaBins = 401.
     */
    int FLambdaBins;

    double halfThetaVelocity = 0.0;
    /** The recursion kernel stores the weight of each [lambda][Flambda] bin. */
    private double[][] recursionKernel;
    /**
     * The maximum value of the last F_lambda bin.
     *
     * <p>maxFLambda = minFLambda + FLambdaBins * dFL.
     */
    private double maxFLambda;
    /**
     * Magnitude of each hill (not including tempering).
     *
     * <p>The default biasMag = 0.05 (kcal/mol).
     */
    private final double biasMag;
    /**
     * The Dama et al. transition-tempering weight.
     *
     * <p>The initial temperingWeight = 1.0, and more generally temperingWeight =
     * exp(-biasHeight/deltaT)
     */
    private double temperingWeight = 1.0;
    /** A count of FLambdaUpdates. */
    private int fLambdaUpdates = 0;
    /**
     * If the recursion kernel becomes too large or too small for some combinations of (Lambda,
     * dU/dL), then its statistical weight = exp(kernel * beta) will cannot be represented by a
     * double value.
     */
    private double[] kernelValues;
    /** Number of times a Gaussian has been added. */
    private int biasCount = 0;
    /**
     * Map lambda to a periodic variable theta. <code>theta = asin(sqrt(lambda))</code> <code> lambda
     * = sin^2 (theta).</code>
     */
    private double theta;
    /** Flag set to false once OST statistics are reset at lambdaResetValue. */
    private boolean resetStatistics;
    /** Most recent lambda values for each Walker. */
    private double lastReceivedLambda;
    /** Most recent dU/dL value for each walker. */
    private double lastReceiveddUdL;
    /** Either the discrete lambda values used, or null (continuous lambda). */
    private final double[] lambdaLadder;

    /**
     * Histogram constructor.
     *
     * @param properties a CompositeConfiguration used to configure the Histogram.
     * @param settings An object containing the values this Histogram will be set to.
     */
    Histogram(CompositeConfiguration properties, HistogramSettings settings) {
      this.temperature = settings.temperature;
      this.dt = settings.dt;
      this.histogramFile = settings.getHistogramFile();
      this.asynchronous = settings.asynchronous;

      discreteLambda = settings.discreteLambda;
      biasCutoff = settings.biasCutoff;
      /*if (discreteLambda) {
        lambdaBiasCutoff = 0;
      } else {
        lambdaBiasCutoff = biasCutoff;
      }*/
      lambdaBiasCutoff = settings.lambdaBiasCutoff;
      double gaussNormalization;
      if (lambdaBiasCutoff == 0 && !settings.histogramRead) {
        gaussNormalization = ONE_D_NORMALIZATION;
        logger.info(format(" Bias magnitude multiplied by a factor of %.4f " +
                "sqrt(2*pi) to match 1D Gaussian volume to 2D Gaussian volume.",
                gaussNormalization));
      } else {
        gaussNormalization = TWO_D_NORMALIZATION;
      }
      biasMag = settings.getBiasMag() * gaussNormalization;
      dFL = settings.dFL;
      temperingFactor = settings.temperingFactor;
      temperOffset = settings.getTemperOffset();

      deltaT = temperingFactor * R * temperature;

      dL = settings.getDL();
      lambdaBins = 1 + (int) round(1.0 / dL);
      if (discreteLambda) {
        lambdaLadder = new double[lambdaBins];
        lambdaLadder[0] = 0.0;
        lambdaLadder[lambdaBins - 1] = 1.0;
        for (int i = 1; i < lambdaBins - 1; i++) {
          lambdaLadder[i] = dL * i;
        }
        dL_2 = 0.0;
        minLambda = 0.0;
      } else {
        lambdaLadder = null;
        dL_2 = dL * 0.5;
        minLambda = -dL_2;
      }

      // The center of the central bin is at 0.
      FLambdaBins = 101;
      minFLambda = -(dFL * FLambdaBins) / 2.0;
      dFL_2 = dFL / 2.0;
      maxFLambda = minFLambda + (FLambdaBins * dFL);
      FLambda = new double[lambdaBins];

      writeIndependent = settings.writeIndependent();
      independentWalkers = settings.independentWalkers();
      assert !(independentWalkers && !writeIndependent);

      lambdaResetValue = settings.getLambdaResetValue();
      resetStatistics = settings.resetStatistics();
      countInterval = settings.countInterval;
      thetaFriction = settings.thetaFriction;
      thetaMass = settings.thetaMass;
      fLambdaPrintInterval = settings.fLambdaPrintInterval;

      // Allocate space for the recursion kernel that stores weights.
      recursionKernel = new double[lambdaBins][FLambdaBins];
      // Allocate space to regularize kernel values.
      kernelValues = new double[FLambdaBins];

      // Random numbers for MD-OST.
      stochasticRandom = new Random();

      String propString = properties.getString("ost-integrationType", "SIMPSONS");
      IntegrationType testType;
      try {
        testType = IntegrationType.valueOf(propString.toUpperCase());
      } catch (Exception ex) {
        logger.warning(
            format(
                " Invalid argument %s to ost-integrationType; resetting to SIMPSONS", propString));
        testType = SIMPSONS;
      }
      integrationType = testType;

      /*
       Set up the multi-walker communication variables for Parallel Java
       communication between nodes.
      */
      world = Comm.world();
      int numProc = world.size();
      rank = world.rank();
      if (asynchronous) {
        // Use asynchronous communication.
        asynchronousSend = new AsynchronousSend(this);
        asynchronousSend.start();
        synchronousSend = null;
      } else {
        Histogram[] histograms = new Histogram[numProc];
        int[] rankToHistogramMap = new int[numProc];
        for (int i = 0; i < numProc; i++) {
          histograms[i] = this;
          rankToHistogramMap[i] = 0;
        }
        synchronousSend = new SynchronousSend(histograms, rankToHistogramMap, independentWalkers);
        asynchronousSend = null;
      }
      lastReceivedLambda = getLambda();
      if (discreteLambda) {
        lastReceivedLambda = discretizedLambda(lastReceivedLambda);
        logger.info(
            format(
                " Discrete lambda: initializing lambda to nearest bin %.5f", lastReceivedLambda));
        lambda = mapLambda(lastReceivedLambda);
        theta = asin(sqrt(lambda));
        lambdaInterface.setLambda(lastReceivedLambda);
      }
      lastReceiveddUdL = getdEdL();

      // Attempt to load a restart file if one exists.
      readRestart();

      hisFileName = FileUtils.relativePathTo(histogramFile).toString();
      if (lambdaFile != null) {
        lambdaFileName = FileUtils.relativePathTo(lambdaFile).toString();
      } else {
        lambdaFileName = FilenameUtils.removeExtension(hisFileName) + ".lam";
      }

      logger.info("\n" + toString());
    }

    /** checkRecursionKernelSize. */
    public void checkRecursionKernelSize() {
      double[] x = new double[nVariables];
      x = potential.getCoordinates(x);
      double[] g = new double[nVariables];
      potential.energyAndGradient(x, g, false);
      double dudl = lambdaInterface.getdEdL();
      checkRecursionKernelSize(dudl);
    }

    public void disableResetStatistics() {
      resetStatistics = false;
    }

    /**
     * evaluate2DPMF.
     *
     * @return A StringBuffer with 2D Bias PMF.
     */
    public StringBuffer evaluate2DPMF() {
      StringBuffer sb = new StringBuffer();
      for (int fLambdaBin = 0; fLambdaBin < FLambdaBins; fLambdaBin++) {
        double currentFL = minFLambda + fLambdaBin * dFL + dFL_2;
        sb.append(format(" %16.8f", currentFL));
        for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
          double bias = -evaluateKernel(lambdaBin, fLambdaBin, biasMag);
          sb.append(format(" %16.8f", bias));
        }
        sb.append("\n");
      }
      return sb;
    }

    /**
     * evaluateTotalPMF.
     *
     * @return A StringBuffer the total 2D PMF.
     */
    public StringBuffer evaluateTotalPMF() {
      StringBuffer sb = new StringBuffer();
      for (int fLambdaBin = 0; fLambdaBin < FLambdaBins; fLambdaBin++) {

        double currentFL = minFLambda + fLambdaBin * dFL + dFL_2;
        sb.append(format(" %16.8f", currentFL));

        for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
          setLambda(lambdaForIndex(lambdaBin));
          double bias1D = -energyAndGradient1D(lambda, false);
          double totalBias = bias1D - evaluateKernel(lambdaBin, fLambdaBin, biasMag);
          sb.append(format(" %16.8f", totalBias));
        }
        sb.append("\n");
      }
      return sb;
    }

    public double getBiasMagnitude() {
      return biasMag;
    }

    /**
     * For MPI parallel jobs, returns true if the walkers are independent (i.e. contribute to only
     * their own histogram).
     *
     * @return True if the walkers are independent.
     */
    public boolean getIndependentWalkers() {
      return independentWalkers;
    }

    public double getLambdaResetValue() {
      return lambdaResetValue;
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
      return resetStatistics;
    }

    /**
     * Return the SynchronousSend associated with this Histogram, if any.
     *
     * @return The SynchronousSend, if any.
     */
    public Optional<SynchronousSend> getSynchronousSend() {
      return Optional.ofNullable(synchronousSend);
    }

    public void setHalfThetaVelocity(double halfThetaV) {
      halfThetaVelocity = halfThetaV;
    }

    public void setLastReceiveddUdL(double lastReceiveddUdL) {
      this.lastReceiveddUdL = lastReceiveddUdL;
    }

    public String toString() {
      StringBuilder sb = new StringBuilder(" OST Histogram");
      sb.append(format("\n  Gaussian bias magnitude: %6.3f (kcal/mol)", biasMag));
      sb.append(format("\n  Tempering offset:        %6.3f (kcal/mol)", temperOffset));
      sb.append(format("\n  Tempering rate:          %6.3f", temperingFactor));
      sb.append(format("\n  Discrete lambda:         %6B", discreteLambda));
      sb.append(format("\n  Number of Lambda bins:   %6d", lambdaBins));
      sb.append(format("\n  Lambda bin width:        %6.3f", dL));
      sb.append(format("\n  Number of dU/dL bins:    %6d", FLambdaBins));
      sb.append(format("\n  dU/dL bin width:         %6.3f (kcal/mol)", dFL));
      sb.append(format("\n  Histogram restart:       %s",
          FileUtils.relativePathTo(histogramFile).toString()));
      return sb.toString();
    }

    /**
     * Eqs. 7 and 8 from the 2012 Crystal Thermodynamics paper.
     *
     * @param print Whether to write the histogram to screen.
     * @param save Whether to write the histogram to disk.
     * @return Free energy (via integration of ensemble-average dU/dL)
     */
    public synchronized double updateFLambda(boolean print, boolean save) {
      double freeEnergy = 0.0;
      double minFL = Double.MAX_VALUE;

      // If the bias magnitude is zero, computing <dU/dL> from
      // counts will not be correct. Assign a temporary non-zero bias magnitude.
      boolean biasMagZero = false;
      if (biasMag <= 0.0) {
        biasMagZero = true;
      }

      // Total histogram weight.
      double totalWeight = 0;
      double beta = 1.0 / (R * temperature);
      StringBuilder stringBuilder = new StringBuilder();

      // Loop over lambda bins, computing <dU/dL> for each bin.
      for (int iL = 0; iL < lambdaBins; iL++) {
        int ulFL = -1;
        int llFL = -1;

        // Find the smallest FL bin that has counts.
        for (int jFL = 0; jFL < FLambdaBins; jFL++) {
          double count = recursionKernel[iL][jFL];
          if (count > 0) {
            llFL = jFL;
            break;
          }
        }

        // Find the largest FL bin that has counts.
        for (int jFL = FLambdaBins - 1; jFL >= 0; jFL--) {
          double count = recursionKernel[iL][jFL];
          if (count > 0) {
            ulFL = jFL;
            break;
          }
        }

        double lambdaCount = 0;
        // The FL range sampled for lambda bin [iL*dL .. (iL+1)*dL]
        double lla = 0.0;
        double ula = 0.0;
        double maxBias = 0;
        if (llFL == -1 || ulFL == -1) {
          FLambda[iL] = 0.0;
          minFL = 0.0;
        } else {
          double ensembleAverageFLambda = 0.0;
          double partitionFunction = 0.0;

          // Evaluate and regularize all kernel values for this value of lambda.
          double offset = evaluateKernelforLambda(iL, llFL, ulFL);

          for (int jFL = llFL; jFL <= ulFL; jFL++) {
            double kernel = kernelValues[jFL];

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

            double currentFLambda = minFLambda + jFL * dFL + dFL_2;
            ensembleAverageFLambda += currentFLambda * weight;
            lambdaCount += recursionKernel[iL][jFL];
          }
          if (minFL > maxBias) {
            minFL = maxBias;
          }
          FLambda[iL] = (partitionFunction == 0) ? 0 : ensembleAverageFLambda / partitionFunction;
          lla = minFLambda + llFL * dFL;
          ula = minFLambda + (ulFL + 1) * dFL;
        }

        double delta = dL;
        if (!discreteLambda && (iL == 0 || iL == lambdaBins - 1)) {
          // The first and last lambda bins are half size.
          delta = dL_2;
        } else if (discreteLambda && iL == 0) {
          // The free energy change to move from L=0 to L=0 is zero.
          delta = 0.0;
        }
        double deltaFreeEnergy = FLambda[iL] * delta;
        freeEnergy += deltaFreeEnergy;
        totalWeight += lambdaCount;

        if (print || save) {
          double llL = iL * dL - dL_2;
          double ulL = llL + dL;
          if (llL < 0.0) {
            llL = 0.0;
          }
          if (ulL > 1.0) {
            ulL = 1.0;
          }

          double midLambda = llL;
          if (!discreteLambda) {
            midLambda = (llL + ulL) / 2.0;
          }
          double bias1D = energyAndGradient1D(midLambda, false);

          double bias2D = 0.0;
          if (!biasMagZero) {
            bias2D = computeBiasEnergy(midLambda, FLambda[iL]) - bias1D;
          }

          stringBuilder.append(
              format(
                  " %6.2e %7.5f %7.1f %7.1f %8.2f %8.2f %8.2f %8.2f %8.2f   %8.2f\n",
                  lambdaCount,
                  midLambda,
                  lla,
                  ula,
                  FLambda[iL],
                  bias1D,
                  bias2D,
                  bias1D + bias2D,
                  freeEnergy,
                  bias1D + bias2D + freeEnergy));
        }
      }

      if (!biasMagZero) {
        double temperEnergy = (minFL > temperOffset) ? temperOffset - minFL : 0;
        temperingWeight = exp(temperEnergy / deltaT);
      }

      freeEnergy = integrateNumeric(FLambda, integrationType);

      if (print && abs(freeEnergy - previousFreeEnergy) > 0.001) {
        logger.info(
            "  Weight   Lambda      dU/dL Bins  <dU/dL>    g(L)  f(L,<dU/dL>) Bias    dG(L) Bias+dG(L)");
        logger.info(stringBuilder.toString());
        logger.info(
            format(
                " The free energy is %12.4f kcal/mol (Total Weight: %6.2e, Tempering: %6.4f, Counts: %12d).",
                freeEnergy, totalWeight, temperingWeight, biasCount));
        logger.info(format(" Minimum Bias %8.3f", minFL));
        previousFreeEnergy = freeEnergy;
      } else if (!save && (print || biasCount % printFrequency == 0)) {
        logger.info(
            format(
                " The free energy is %12.4f kcal/mol (Total Weight: %6.2e, Tempering: %6.4f, Counts: %12d).",
                freeEnergy, totalWeight, temperingWeight, biasCount));
      }

      if (save) {
        String modelFilename = molecularAssembly.getFile().getAbsolutePath();
        File saveDir = new File(FilenameUtils.getFullPath(modelFilename));
        String dirName = saveDir.toString() + File.separator;
        String fileName = dirName + "histogram.txt";
        try {
          logger.info(" Writing " + fileName);
          PrintWriter printWriter = new PrintWriter(new File(fileName));
          printWriter.write(stringBuilder.toString());
          printWriter.close();
        } catch (Exception e) {
          logger.info(format(" Failed to write %s.", fileName));
        }
      }

      return freeEnergy;
    }

    /**
     * Converts a continuous lambda value into its nearest discretized value.
     *
     * @param lambda Lambda to discretize.
     * @return Discretized lambda.
     */
    private double discretizedLambda(double lambda) {
      assert discreteLambda;
      int bin = indexForLambda(lambda);
      return bin * dL;
    }

    /** Write a Histogram restart file (sometimes skipped for rank > 0). */
    void writeRestart() {
      if (rank == 0 || writeIndependent) {
        try {
          HistogramWriter histogramWriter =
              new HistogramWriter(this, new BufferedWriter(new FileWriter(histogramFile)));
          histogramWriter.writeHistogramFile();
          histogramWriter.flush();
          histogramWriter.close();
          logger.info(format(" Wrote histogram restart file to %s.", hisFileName));
        } catch (IOException ex) {
          String message = " Exception writing OST histogram restart file.";
          logger.log(Level.INFO, Utilities.stackTraceToString(ex));
          logger.log(Level.SEVERE, message, ex);
        }
      }
    }

    /** Read a Histogram restart file (all ranks). */
    void readRestart() {
      File histogramToRead = histogramFile;
      // Independent walkers will write Histogram re-starts into the subdirectories,
      // but we might be restarting from a Histogram in the parent directory.
      if (histogramFile != null && !histogramFile.exists()) {
        // Try to find a restart in the parent directory:
        // Update /path/to/filename.his --> /path/to/../filename.his
        String histogramString = histogramFile.getAbsolutePath();
        String path = FilenameUtils.getPath(histogramString);
        String name = FilenameUtils.getName(histogramString);
        String parentHistogram = File.separator + path + ".." + File.separator + name;
        histogramToRead = new File(parentHistogram);
        if (histogramToRead.exists()) {
          logger.info(" Reading parent histogram: " + parentHistogram);
        }
      }
      if (histogramToRead != null && histogramToRead.exists()) {
        try {
          HistogramReader histogramReader =
              new HistogramReader(this, new FileReader(histogramToRead));
          histogramReader.readHistogramFile();
          // Currently, the restart format does not include info on using discrete lambda values.
          updateFLambda(true, false);
          logger.info(format("\n Read OST histogram from %s.", hisFileName));
        } catch (FileNotFoundException ex) {
          logger.info(" Histogram restart file could not be found and will be ignored.");
        }
      }
    }

    private double mapLambda(double lambda) {
      if (discreteLambda) {
        return lambdaLadder[indexForDiscreteLambda(lambda)];
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
      // I give up on making this a nice pretty method reference.
      if (discreteLambda) {
        return indexForDiscreteLambda(lambda);
      } else {
        return indexForContinuousLambda(lambda);
      }
    }

    private double lambdaForIndex(int bin) {
      if (discreteLambda) {
        return lambdaLadder[bin];
      } else {
        return (bin * dL) - dL_2;
      }
    }

    private int indexForContinuousLambda(double lambda) {
      int lambdaBin = (int) floor((lambda - minLambda) / dL);
      if (lambdaBin < 0) {
        lambdaBin = 0;
      }
      if (lambdaBin >= lambdaBins) {
        lambdaBin = lambdaBins - 1;
      }
      return lambdaBin;
    }

    private int indexForDiscreteLambda(double lambda) {
      assert discreteLambda && lambdaLadder != null && lambdaLadder.length > 0;

      int initialGuess = indexForContinuousLambda(lambda);
      double minErr = Double.MAX_VALUE;
      int minErrBin = -1;
      for (int i = -1; i < 2; i++) {
        int guessBin = i + initialGuess;
        if (guessBin < 0 || guessBin >= lambdaBins) {
          continue;
        }
        double guessLam = lambdaLadder[guessBin];
        double guessErr = Math.abs(guessLam - lambda);
        if (guessErr < minErr) {
          minErr = guessErr;
          minErrBin = guessBin;
        }
      }

      assert minErr < 1.0E-6 && minErrBin > -1;
      return minErrBin;
    }

    /**
     * binForFLambda.
     *
     * @param dEdLambda a double.
     * @return a int.
     */
    int binForFLambda(double dEdLambda) {
      int FLambdaBin = (int) floor((dEdLambda - minFLambda) / dFL);
      if (FLambdaBin == FLambdaBins) {
        FLambdaBin = FLambdaBins - 1;
      }
      assert (FLambdaBin < FLambdaBins);
      assert (FLambdaBin >= 0);
      return FLambdaBin;
    }

    /**
     * Return the value of a recursion kernel bin.
     *
     * @param lambdaBin The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @return The value of the bin.
     */
    double getRecursionKernelValue(int lambdaBin, int fLambdaBin) {
      return recursionKernel[lambdaBin][fLambdaBin];
    }

    /**
     * Set the value of a recursion kernel bin.
     *
     * @param lambdaBin The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @param value The value of the bin.
     */
    void setRecursionKernelValue(int lambdaBin, int fLambdaBin, double value) {
      recursionKernel[lambdaBin][fLambdaBin] = value;
    }

    /**
     * Add to the value of a recursion kernel bin.
     *
     * @param lambdaBin The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @param value The value of the bin.
     * @param updateFLambda Whether to update the 1D bias (typically true for biases received
     *     from other processes)
     */
    synchronized void addToRecursionKernelValue(
        int lambdaBin, int fLambdaBin, double value, boolean updateFLambda) {
      recursionKernel[lambdaBin][fLambdaBin] += value;
      if (updateFLambda) {
        updateFLambda(false, false);
      }
      ++biasCount;
    }

    /** Allocate memory for the recursion kernel. */
    void allocateRecursionKernel() {
      recursionKernel = new double[lambdaBins][FLambdaBins];
      kernelValues = new double[FLambdaBins];
    }

    /**
     * Integrates dUdL over lambda using more sophisticated techniques than midpoint rectangular
     * integration.
     *
     * <p>The ends (from 0 to dL and 1-dL to 1) are integrated with trapezoids for continuous
     * lambda.
     *
     * @param dUdLs dUdL at the midpoint of each bin.
     * @param type Integration type to use.
     * @return Current delta-G estimate.
     */
    private double integrateNumeric(double[] dUdLs, IntegrationType type) {
      double val;
      if (discreteLambda) {
        // Integrate between the first bin and the last bin.
        double[] lams = Integrate1DNumeric.generateXPoints(0.0, 1.0, lambdaBins, false);
        DataSet dSet = new DoublesDataSet(lams, dUdLs, false);
        val = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, type);
      } else {
        // Integrate between the second bin midpoint and the second-to-last bin midpoint.
        double[] midLams =
            Integrate1DNumeric.generateXPoints(dL, 1.0 - dL, (lambdaBins - 2), false);
        double[] midVals = Arrays.copyOfRange(dUdLs, 1, (lambdaBins - 1));
        DataSet dSet = new DoublesDataSet(midLams, midVals, false);

        val = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, type);

        // Everything after this is just adding in the endpoint contributions.

        double dL_4 = dL_2 * 0.5;

        // Initially, assume dU/dL is exactly 0 at the endpoints. This is sometimes a true
        // assumption.
        double val0 = 0;
        double val1 = 0;

        // If we cannot guarantee that dUdL is exactly 0 at the endpoints, interpolate.
        if (!lambdaInterface.dEdLZeroAtEnds()) {
          double recipSlopeLen = 1.0 / (dL * 0.75);

          double slope = dUdLs[0] - dUdLs[1];
          slope *= recipSlopeLen;
          val0 = dUdLs[0] + (slope * dL_4);

          slope = dUdLs[lambdaBins - 1] - dUdLs[lambdaBins - 2];
          slope *= recipSlopeLen;
          val1 = dUdLs[lambdaBins - 1] + (slope * dL_4);
          logger.fine(format(" Inferred dU/dL values at 0 and 1: %10.5g , %10.5g", val0, val1));
        }

        // Integrate trapezoids from 0 to the second bin midpoint, and from second-to-last bin
        // midpoint to 1.
        val += trapezoid(0, dL_4, val0, dUdLs[0]);
        val += trapezoid(dL_4, dL, dUdLs[0], dUdLs[1]);
        val += trapezoid(1.0 - dL, 1.0 - dL_4, dUdLs[lambdaBins - 2], dUdLs[lambdaBins - 1]);
        val += trapezoid(1.0 - dL_4, 1.0, dUdLs[lambdaBins - 1], val1);
      }

      return val;
    }

    /**
     * Integrates a trapezoid.
     *
     * @param x0 First x point
     * @param x1 Second x point
     * @param fX1 First f(x) point
     * @param fX2 Second f(x) point
     * @return The area under a trapezoid.
     */
    private double trapezoid(double x0, double x1, double fX1, double fX2) {
      double val = 0.5 * (fX1 + fX2);
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
     * @param lambda Value of Lambda to evaluate the kernal for.
     * @param llFL Lower FLambda bin.
     * @param ulFL Upper FLambda bin.
     * @return the applied offset.
     */
    private double evaluateKernelforLambda(int lambda, int llFL, int ulFL) {
      double maxKernel = Double.MIN_VALUE;

      double gaussianBiasMagnitude = biasMag;
      if (biasMag <= 0.0) {
        gaussianBiasMagnitude = PSEUDO_BIAS_MAGNITUDE;
      }

      for (int jFL = llFL; jFL <= ulFL; jFL++) {
        double value = evaluateKernel(lambda, jFL, gaussianBiasMagnitude);
        kernelValues[jFL] = value;
        if (value > maxKernel) {
          maxKernel = value;
        }
      }

      // Only offset the kernel values for non-zero bias magnitude.
      double offset = 0.0;
      if (biasMag > 0.0) {
        offset = -maxKernel;
        for (int jFL = llFL; jFL <= ulFL; jFL++) {
          kernelValues[jFL] += offset;
        }
      }
      return offset;
    }

    /** Evaluate the bias at [cLambda, cF_lambda] */
    private double evaluateKernel(int cLambda, int cF_Lambda, double gaussianBiasMagnitude) {
      // Compute the value of L and FL for the center of the current bin.
      double vL = cLambda * dL;
      double vFL = minFLambda + cF_Lambda * dFL + dFL_2;

      // Set the variances for the Gaussian bias.
      double Ls2 = 2.0 * dL * 2.0 * dL;
      double FLs2 = 2.0 * dFL * 2.0 * dFL;

      // Variances are only used when dividing by twice their value.
      double invLs2 = 0.5 / Ls2;
      double invFLs2 = 0.5 / FLs2;

      double sum = 0.0;
      for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
        int Lcenter = cLambda + iL;
        double deltaL = vL - Lcenter * dL;
        double deltaL2 = deltaL * deltaL;

        // Pre-compute the lambda-width exponent.
        double L2exp = exp(-deltaL2 * invLs2);

        // Mirror condition for Lambda counts.
        int lcount = Lcenter;
        double mirrorFactor = 1.0;
        if (lcount == 0 || lcount == lambdaBins - 1) {
          // For continuous lambda, the width of the first and last bins is dLambda_2,
          // so the mirror condition is to double their counts.
          if (!discreteLambda) {
            mirrorFactor = 2.0;
          }
        } else if (lcount < 0) {
          lcount = -lcount;
        } else if (lcount > lambdaBins - 1) {
          // number of bins past the last bin.
          lcount -= (lambdaBins - 1);
          // mirror bin
          lcount = lambdaBins - 1 - lcount;
        }

        for (int jFL = -biasCutoff; jFL <= biasCutoff; jFL++) {
          int FLcenter = cF_Lambda + jFL;
          // For FLambda outside the count matrix the weight is 0 so we continue.
          if (FLcenter < 0 || FLcenter >= FLambdaBins) {
            continue;
          }
          double rc = recursionKernel[lcount][FLcenter];
          if (rc <= 0.0) {
            continue;
          }
          double deltaFL = vFL - (minFLambda + FLcenter * dFL + dFL_2);
          double deltaFL2 = deltaFL * deltaFL;
          double weight = mirrorFactor * rc;
          if (weight > 0) {
            double e = weight * gaussianBiasMagnitude * L2exp * exp(-deltaFL2 * invFLs2);
            sum += e;
          }
        }
      }
      return sum;
    }

    /**
     * Compute the total Bias energy at (currentLambda, currentdUdL).
     *
     * @param currentLambda The value of lambda.
     * @param currentdUdL The value of dU/dL.
     * @return The bias energy.
     */
    double computeBiasEnergy(double currentLambda, double currentdUdL) {

      int lambdaBin = indexForLambda(currentLambda);
      int FLambdaBin = binForFLambda(currentdUdL);

      double bias2D = 0.0;
      if (biasMag > 0.0) {
        double ls2 = (2.0 * dL) * (2.0 * dL);
        double FLs2 = (2.0 * dFL) * (2.0 * dFL);
        for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
          int lcenter = lambdaBin + iL;
          double deltaL = currentLambda - (lcenter * dL);
          double deltaL2 = deltaL * deltaL;
          // Mirror conditions for recursion kernel counts.
          int lcount = lcenter;
          double mirrorFactor = 1.0;
          if (lcount == 0 || lcount == lambdaBins - 1) {
            // For continuous lambda, the width of the first and last bins is dLambda_2,
            // so the mirror condition is to double their counts.
            if (!discreteLambda) {
              mirrorFactor = 2.0;
            }
          } else if (lcount < 0) {
            lcount = -lcount;
          } else if (lcount > lambdaBins - 1) {
            // Number of bins past the last bin
            lcount -= (lambdaBins - 1);
            // Mirror bin
            lcount = lambdaBins - 1 - lcount;
          }
          for (int iFL = -biasCutoff; iFL <= biasCutoff; iFL++) {
            int FLcenter = FLambdaBin + iFL;

            // If either of the following FL edge conditions are true,
            // then there are no counts and we continue.
            if (FLcenter < 0 || FLcenter >= FLambdaBins) {
              continue;
            }
            double rc = recursionKernel[lcount][FLcenter];
            if (rc <= 0.0) {
              continue;
            }
            double currentFL = minFLambda + FLcenter * dFL + dFL_2;
            double deltaFL = currentdUdL - currentFL;
            double deltaFL2 = deltaFL * deltaFL;
            double weight = mirrorFactor * rc;
            double bias =
                weight * biasMag * exp(-deltaL2 / (2.0 * ls2)) * exp(
                    -deltaFL2 / (2.0 * FLs2));
            bias2D += bias;
          }
        }
      }

      // Compute the energy for the recursion worker at F(L) using interpolation.
      double bias1D = energyAndGradient1D(currentLambda, false);

      // Return the total bias.
      return bias1D + bias2D;
    }

    double energyAndGradient2D(
        double currentLambda, double currentdUdLambda, double[] chainRule,
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
      int lambdaBin = indexForLambda(currentLambda);
      int FLambdaBin = binForFLambda(currentdUdLambda);

      double ls2 = (2.0 * dL) * (2.0 * dL);
      double FLs2 = (2.0 * dFL) * (2.0 * dFL);
      for (int iL = -lambdaBiasCutoff; iL <= lambdaBiasCutoff; iL++) {
        int lcenter = lambdaBin + iL;
        // double deltaL = lambda - (lcenter * dL);
        double deltaL = currentLambda - (lcenter * dL);
        double deltaL2 = deltaL * deltaL;
        // Mirror conditions for recursion kernel counts.
        int lcount = lcenter;
        double mirrorFactor = 1.0;
        if (lcount == 0 || lcount == lambdaBins - 1) {
          // For continuous lambda, the width of the first and last bins is dLambda_2,
          // so the mirror condition is to double their counts.
          if (!discreteLambda) {
            mirrorFactor = 2.0;
          }
        } else if (lcount < 0) {
          lcount = -lcount;
        } else if (lcount > lambdaBins - 1) {
          // Number of bins past the last bin
          lcount -= (lambdaBins - 1);
          // Mirror bin
          lcount = lambdaBins - 1 - lcount;
        }
        for (int iFL = -biasCutoff; iFL <= biasCutoff; iFL++) {
          int FLcenter = FLambdaBin + iFL;

          // If either of the following FL edge conditions are true, then there are no counts and we
          // continue.
          if (FLcenter < 0 || FLcenter >= FLambdaBins) {
            continue;
          }
          double rc = recursionKernel[lcount][FLcenter];
          if (rc <= 0.0) {
            continue;
          }
          // double deltaFL = dUdLambda - (minFLambda + FLcenter * dFL + dFL_2);
          double deltaFL = currentdUdLambda - (minFLambda + FLcenter * dFL + dFL_2);
          double deltaFL2 = deltaFL * deltaFL;
          double weight = mirrorFactor * rc;
          double bias =
              weight * gaussianBiasMagnitude * exp(-deltaL2 / (2.0 * ls2)) * exp(
                  -deltaFL2 / (2.0 * FLs2));
          gLdEdL += bias;
          dGdLambda -= deltaL / ls2 * bias;
          dGdFLambda -= deltaFL / FLs2 * bias;
        }
      }

      chainRule[0] = dGdLambda;
      chainRule[1] = dGdFLambda;
      return gLdEdL;
    }

    /**
     * This calculates the 1D OST bias and its derivative with respect to Lambda.
     *
     * <p>See Equation 8 in http://doi.org/10.1021/ct300035u.
     *
     * @return a double.
     */
    private double energyAndGradient1D(double currentLambda, boolean gradient) {
      double biasEnergy = 0.0;
      for (int iL0 = 0; iL0 < lambdaBins - 1; iL0++) {
        int iL1 = iL0 + 1;

        // Find bin centers and values for interpolation / extrapolation points.
        double L0 = iL0 * dL;
        double L1 = L0 + dL;
        double FL0 = FLambda[iL0];
        double FL1 = FLambda[iL1];
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
        biasEnergy += (FL0 * L1 + deltaFL * L1 * (0.5 * L1 - L0) / dL);
        biasEnergy -= (FL0 * L0 + deltaFL * L0 * (-0.5 * L0) / dL);
        if (done) {
          // Compute the gradient d F(L) / dL at L.
          if (gradient) {
            dUdLambda -= FL0 + (L1 - L0) * deltaFL / dL;
          }
          break;
        }
      }
      return -biasEnergy;
    }

    /**
     * evaluateHistogram.
     *
     * @param lambda the lambda value.
     * @param dUdL the dU/dL value.
     * @return The value of the Histogram.
     */
    double evaluateHistogram(double lambda, double dUdL) {
      int lambdaBin = indexForLambda(lambda);
      int dUdLBin = binForFLambda(dUdL);
      try {
        return recursionKernel[lambdaBin][dUdLBin];
      } catch (Exception e) {
        // Catch an index out of bounds exception.
        return 0.0;
      }
    }

    /** If necessary, allocate more space. */
    void checkRecursionKernelSize(double dEdLambda) {
      if (dEdLambda > maxFLambda) {
        logger.info(
            format(
                " Current F_lambda %8.2f > maximum histogram size %8.2f.", dEdLambda, maxFLambda));

        double origDeltaG = updateFLambda(false, false);

        int newFLambdaBins = FLambdaBins;
        while (minFLambda + newFLambdaBins * dFL < dEdLambda) {
          newFLambdaBins += 100;
        }
        double[][] newRecursionKernel = new double[lambdaBins][newFLambdaBins];

        // We have added bins above the indeces of the current counts just copy them into the new
        // array.
        for (int i = 0; i < lambdaBins; i++) {
          arraycopy(recursionKernel[i], 0, newRecursionKernel[i], 0, FLambdaBins);
        }
        recursionKernel = newRecursionKernel;
        FLambdaBins = newFLambdaBins;
        kernelValues = new double[FLambdaBins];
        maxFLambda = minFLambda + dFL * FLambdaBins;
        logger.info(
            format(
                " New histogram %8.2f to %8.2f with %d bins.\n",
                minFLambda, maxFLambda, FLambdaBins));

        double newFreeEnergy = updateFLambda(false, false);
        assert (origDeltaG == newFreeEnergy);
      }
      if (dEdLambda < minFLambda) {
        logger.info(
            format(
                " Current F_lambda %8.2f < minimum histogram size %8.2f.", dEdLambda, minFLambda));

        double origDeltaG = updateFLambda(false, false);

        int offset = 100;
        while (dEdLambda < minFLambda - offset * dFL) {
          offset += 100;
        }
        int newFLambdaBins = FLambdaBins + offset;
        double[][] newRecursionKernel = new double[lambdaBins][newFLambdaBins];

        // We have added bins below the current counts,
        // so their indeces must be increased by: offset = newFLBins - FLBins
        for (int i = 0; i < lambdaBins; i++) {
          arraycopy(recursionKernel[i], 0, newRecursionKernel[i], offset, FLambdaBins);
        }
        recursionKernel = newRecursionKernel;
        minFLambda = minFLambda - offset * dFL;
        FLambdaBins = newFLambdaBins;
        kernelValues = new double[FLambdaBins];

        logger.info(
            format(
                " New histogram %8.2f to %8.2f with %d bins.\n",
                minFLambda, maxFLambda, FLambdaBins));

        double newFreeEnergy = updateFLambda(false, false);
        assert (origDeltaG == newFreeEnergy);
      }
    }

    /** Add a Gaussian hill to the Histogram at (lambda, dEdU). */
    void addBias(double dEdU, double[] x, double[] gradient) {
      // Communicate adding the bias to all walkers.
      if (asynchronous) {
        asynchronousSend.send(lambda, dEdU, temperingWeight);
      } else {
        synchronousSend.send(lambda, dEdU, temperingWeight);
      }

      // Update F(L)
      fLambdaUpdates++;
      boolean printFLambda = fLambdaUpdates % fLambdaPrintInterval == 0;
      updateFLambda(printFLambda, false);

      // Locally optimize the current state.
      optimizationParameters.optimize(forceFieldEnergy, x, gradient);

      // Write out restart files.
      if (energyCount > 0 && energyCount % saveFrequency == 0) {
        OrthogonalSpaceTempering.this.writeRestart();
      }
    }

    /** Propagate Lambda using Langevin dynamics. */
    private void langevin() {
      // Compute the random force pre-factor (kcal/mol * psec^-2).
      double rt2 = 2.0 * Constants.R * histogram.temperature * thetaFriction / dt;

      // Compute the random force.
      double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / randomConvert;

      // Compute dEdL (kcal/mol).
      double dEdL = -dUdLambda * sin(2.0 * theta);

      // Update halfThetaVelocity (ps-1).
      halfThetaVelocity =
          (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
              + randomConvert2 * 2.0 * dt * (dEdL + randomForce))
              / (2.0 * thetaMass + thetaFriction * dt);

      // Update theta.
      theta = theta + dt * halfThetaVelocity;

      // Maintain theta in the interval PI to -PI.
      if (theta > PI) {
        theta -= 2.0 * PI;
      } else if (theta <= -PI) {
        theta += 2.0 * PI;
      }

      // Compute the sin(theta).
      double sinTheta = sin(theta);

      // Compute lambda as sin(theta)^2.
      setLambda(sinTheta * sinTheta);
      lambdaInterface.setLambda(lambda);
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
      if (asynchronousSend != null && asynchronousSend.isAlive()) {
        double[] killMessage = new double[] {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        DoubleBuf killBuf = DoubleBuf.buffer(killMessage);
        try {
          logger.fine(" Sending the termination message.");
          world.send(rank, killBuf);
          logger.fine(" Termination message was sent successfully.");
          logger.fine(
              format(
                  " Receive thread alive %b status %s",
                  asynchronousSend.isAlive(), asynchronousSend.getState()));
        } catch (Exception ex) {
          String message =
              format(
                  " Asynchronous Multiwalker OST termination signal "
                      + "failed to be sent for process %d.",
                  rank);
          logger.log(Level.SEVERE, message, ex);
        }
      } else {
        logger.fine(
            " CountReceiveThread was either not initialized, or is not alive. This is the case for the Histogram script.");
      }
    }

    int getHistogramIndex() {
      if (asynchronous) {
        return writeIndependent ? rank : 0;
      } else {
        return synchronousSend.getHistogramIndex();
      }
    }
  }
}
