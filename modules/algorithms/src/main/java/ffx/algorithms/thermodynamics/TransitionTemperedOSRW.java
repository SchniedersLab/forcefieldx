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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
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
import static ffx.numerics.integrate.Integrate1DNumeric.IntegrationType.SIMPSONS;
import static ffx.utilities.Constants.R;

/**
 * An implementation of Transition-Tempered Orthogonal Space Random Walk
 * algorithm.
 * <p>
 * This only partially implements the LambdaInterface, since it does not return 2nd lambda
 * derivatives. The 2nd derivatives of the bias require 3rd derivatives of
 * the underlying Hamiltonian (not these derivatives are not needed for OSRW MD).
 *
 * @author Michael J. Schnieders, James Dama, Wei Yang and Pengyu Ren
 * @since 1.0
 */
public class TransitionTemperedOSRW implements CrystalPotential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(TransitionTemperedOSRW.class.getName());

    /**
     * The MolecularAssembly being simulated.
     */
    protected MolecularAssembly molecularAssembly;
    /**
     * A potential energy that implements the LambdaInterface.
     */
    private final LambdaInterface lambdaInterface;
    /**
     * The potential energy of the system.
     */
    protected final CrystalPotential potential;
    /**
     * Reference to the Barostat in use; if present this must be turned off
     * during optimization.
     */
    protected final Barostat barostat;
    /**
     * The AlgorithmListener is called each time a count is added.
     */
    protected final AlgorithmListener algorithmListener;
    /**
     * Number of variables.
     */
    protected final int nVariables;
    /**
     * Each walker has a unique lambda restart file.
     */
    private final File lambdaFile;
    /**
     * Each walker reads the same histogram restart file. Only the walker of
     * rank 0 writes the histogram restart file.
     */
    private final File histogramFile;
    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    protected double lambda;
    /**
     * Flag to indicate that the Lambda particle should be propagated.
     */
    private boolean propagateLambda = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "propagateLambda" flag true.
     */
    int energyCount;
    /**
     * Number of times a Gaussian has been added.
     */
    private int biasCount = 0;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005). The final
     * Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     * <p>
     * With this scheme, the maximum of biasing Gaussians is at the edges.
     * <p>
     * The default lambdaBins = 201.
     */
    int lambdaBins;
    /**
     * It is useful to have an odd number of bins, so that there is a bin from
     * FL=-dFL/2 to dFL/2 so that as FL approaches zero its contribution to
     * thermodynamic integration goes to zero.
     * <p>
     * Otherwise a contribution of zero from a L bin can only result from equal
     * sampling of the ranges -dFL to 0 and 0 to dFL.
     * <p>
     * The default FLambdaBins = 401.
     */
    int FLambdaBins;
    /**
     * Parallel Java world communicator.
     */
    protected final Comm world;
    /**
     * Number of processes.
     */
    private final int numProc;
    /**
     * Rank of this process.
     */
    protected final int rank;
    /**
     * When evaluating the biasing potential, contributions from Gaussians
     * centered on bins more the "biasCutoff" away will be neglected.
     * <p>
     * The default biasCutoff = 5.
     */
    int biasCutoff;
    /**
     * Width of the lambda bin.
     * <p>
     * The default dL = (1.0 / (lambdaBins - 1).
     */
    protected double dL;
    /**
     * Half the width of a lambda bin.
     */
    double dL_2;
    /**
     * The width of the FLambda bin.
     * <p>
     * The default dFL = 2.0 (kcal/mol).
     */
    double dFL;
    /**
     * Half the width of the F_lambda bin.
     */
    double dFL_2;
    /**
     * The minimum value of the first lambda bin.
     * <p>
     * minLambda = -dL_2.
     */
    private double minLambda;
    /**
     * The minimum value of the first F_lambda bin.
     * <p>
     * minFLambda = -(dFL * FLambdaBins) / 2.0.
     */
    double minFLambda;
    /**
     * The maximum value of the last F_lambda bin.
     * <p>
     * maxFLambda = minFLambda + FLambdaBins * dFL.
     */
    private double maxFLambda;
    /**
     * Force Field Potential Energy (i.e. with no bias terms added).
     */
    protected double forceFieldEnergy;
    /**
     * Partial derivative of the force field energy with respect to lambda.
     */
    private double dForceFieldEnergydL;
    /**
     * OSRW Bias energy.
     */
    private double biasEnergy;
    /**
     * Total partial derivative of the potential (U) being sampled
     * with respect to lambda.
     */
    private double dUdLambda;
    /**
     * Second partial derivative of the potential being sampled with respect to lambda.
     */
    private double d2UdL2;
    /**
     * Mixed second partial derivative with respect to coordinates and lambda.
     */
    private double[] dUdXdL;
    /**
     * Magnitude of each hill (not including tempering).
     * <p>
     * The default biasMag = 0.05 (kcal/mol).
     */
    double biasMag;
    /**
     * 1D PMF with respect to lambda F(L).
     */
    double[] FLambda;
    /**
     * Magnitude of the 2D orthogonal space bias G(L,dE/dL).
     */
    private double gLdEdL = 0.0;
    /**
     * Reasonable thetaFriction is ~60 per picosecond (1.0e-12).
     */
    double thetaFriction = 1.0e-19;
    /**
     * Reasonable thetaMass is ~100 a.m.u. (100 a.m.u is 1.6605e-22 grams).
     */
    double thetaMass = 1.0e-18;
    double halfThetaVelocity = 0.0;
    /**
     * Time step in picoseconds.
     */
    protected final double dt;
    /**
     * Temperature in Kelvin.
     * <p>
     * The default is 298.15.
     */
    protected double temperature;
    /**
     * Interval between adding a count to the Recursion kernel in steps.
     * <p>
     * The default countInterval = 10.
     */
    int countInterval = 10;
    /**
     * Interval between printing information on the lambda particle in steps.
     * <p>
     * The default printFrequency = 100.
     */
    private int printFrequency;
    /**
     * Once the lambda reset value is reached, Transition-Tempered OSRW
     * statistics are reset.
     */
    final double lambdaResetValue = 0.99;
    /**
     * Flag set to false once Transition-Tempered OSRW statistics are reset at
     * lambdaResetValue.
     */
    boolean resetStatistics = false;
    /**
     * Flag to turn on OSRW optimization.
     * <p>
     * The default osrwOptimization = false.
     */
    private boolean osrwOptimization = false;
    /**
     * Holds the lowest potential-energy parameters for loopBuilder runs from
     * all visits to lambda &gt; osrwOptimizationLambdaCutoff.
     */
    private double[] osrwOptimumCoords;
    /**
     * The lowest energy found via optimizations.
     * <p>
     * The osrwOptimum is initially set to Double.MAX_VALUE.
     */
    private double osrwOptimum = Double.MAX_VALUE;
    /**
     * OSRW optimization only runs if Lambda is greater than the osrwOptimizationLambdaCutoff.
     * <p>
     * The default osrwOptimizationLambdaCutoff = 0.8.
     */
    private double osrwOptimizationLambdaCutoff = 0.8;
    /**
     * The OSRW optimization frequency.
     */
    private int osrwOptimizationFrequency = 10000;
    /**
     * The OSRW optimization convergence criteria.
     * <p>
     * The default osrwOptimizationEps = 0.1.
     */
    private double osrwOptimizationEps = 0.1;
    /**
     * The OSRW optimization tolerance.
     */
    private double osrwOptimizationTolerance = 1.0e-8;
    /**
     * The OSRW optimization energy window.
     */
    private double osrwOptimizationEnergyWindow = 2.0;
    /**
     * File instance used for saving optimized structures.
     */
    private File osrwOptimizationFile;
    /**
     * SystemFilter used to save optimized structures.
     */
    private SystemFilter osrwOptimizationFilter;
    /**
     * Interval between how often the 1D histogram is printed to screen versus
     * silently updated in background.
     * <p>
     * The fLambdaPrintInterval is 25.
     */
    private final int fLambdaPrintInterval = 25;
    /**
     * A count of FLambdaUpdates.
     */
    private int fLambdaUpdates = 0;
    /**
     * Interval between writing an OSRW restart file in steps.
     * <p>
     * The default saveFrequency = 1000.
     */
    int saveFrequency;
    /**
     * Print detailed energy information.
     */
    protected final boolean print = false;
    /**
     * Total system energy.
     */
    private double totalEnergy;
    /**
     * Save the previous free energy, in order to limit logging to time points
     * where the free energy has changed.
     */
    private double previousFreeEnergy = 0.0;
    /**
     * Are FAST varying energy terms being computed, SLOW varying energy terms,
     * or BOTH. OSRW is not active when only FAST varying energy terms are being
     * propagated.
     */
    protected Potential.STATE state = Potential.STATE.BOTH;
    /**
     * Flag to indicate if OSRW should send and receive counts between processes
     * synchronously or asynchronously. The latter can be faster by ~40% because
     * simulation with Lambda &gt; 0.75 must compute two condensed phase
     * self-consistent fields to interpolate polarization.
     */
    private final boolean asynchronous;
    /**
     * Write out structures only for lambda values greater than or equal to this threshold.
     */
    double lambdaWriteOut = 0.0;
    /**
     * Should the 1D OSRW bias be included in the target function.
     */
    private boolean include1DBias;
    /**
     * Map lambda to a periodic variable theta.
     *
     * <code>theta = asin(sqrt(lambda))</code>
     *
     * <code>lambda = sin^2 (theta).</code>
     */
    private double theta;
    private final Random stochasticRandom;
    /**
     * Random force conversion to kcal/mol/A;
     * Units: Sqrt (4.184 Joule per calorie) / (nanometers per meter)
     */
    private static final double randomConvert = sqrt(4.184) / 10e9;
    /**
     * randomConvert squared.
     * Units: Joule per calorie / (nanometer per meter)^2
     */
    private static final double randomConvert2 = randomConvert * randomConvert;
    /**
     * The recursion kernel stores the weight of each [lambda][Flambda] bin.
     */
    private double[][] recursionKernel;
    /**
     * The integration algorithm used for thermodynamic integration.
     */
    private final IntegrationType integrationType;
    /**
     * A flag to indicate if the transition has been crossed and Dama et al.
     * transition-tempering should begin.
     * <p>
     * Currently "tempering" is always true, but a temperOffset (default 1.0 kcal/mol) prevents tempering
     * prior to there being a bias everywhere along the lambda path of at least temperOffset.
     */
    private boolean tempering = true;
    /**
     * The Dama et al. transition-tempering rate parameter. A reasonable value
     * is about 2 to 8 kT, with larger values being resulting in slower decay.
     * <p>
     * The default temperingFactor = 2.0.
     */
    private double temperingFactor = 2.0;
    /**
     * This deltaT is used to determine the tempering weight as described below for the temperingWeight variable.
     * <p>
     * deltaT = temperingFactor * kB * T.
     */
    private double deltaT;
    /**
     * The Dama et al. transition-tempering weight.
     * <p>
     * The initial temperingWeight = 1.0, and more generally temperingWeight = exp(-biasHeight/deltaT)
     */
    private double temperingWeight = 1.0;
    /**
     * An offset applied before calculating tempering weight.
     * <p>
     * First, for every Lambda, we calculate the maximum bias at that lambda by searching all populated dU/dL bins:
     * maxdUdL(L) = max[ G(L,F_L) ] where the max operator searches over all F_L bins.
     * <p>
     * Then, the minimum bias coverage is determined by searching the maxdUdL(L) array over Lambda.
     * minBias = min[ maxdUdL(L) ] where the min operator searches all entries in the array.
     * <p>
     * Then the temperOffset is applied to the minBias:
     * <p>
     * biasHeight = max[minBias - temperOffset, 0]
     * <p>
     * The default temperOffset = 1.0 kcal/mol.
     */
    private double temperOffset;
    /**
     * If this flag is true, (lambda, dU/dL) samples that have no weight in the Histogram are rejected.
     */
    private boolean hardWallConstraint = false;
    /**
     * If the recursion kernel becomes too large or too small for some combinations of (Lambda, dU/dL),
     * then its statistical weight = exp(kernel * beta) will cannot be represented by a double value.
     */
    private double[] kernelValues;
    /**
     * The CountReceiveThread accumulates OSRW statistics from multiple asynchronous walkers.
     */
    private final CountReceiveThread receiveThread;
    /**
     * The recursionWeights stores the [Lambda, FLambda] weight for each
     * process. Therefore the array is of size [number of Processes][2].
     * <p>
     * Each 2 entry array must be wrapped inside a Parallel Java DoubleBuf for the
     * All-Gather communication calls.
     */
    private final double[][] recursionWeights;
    private final double[] myRecursionWeight;
    /**
     * These DoubleBufs wrap the recursionWeight arrays.
     */
    private final DoubleBuf[] recursionWeightsBuf;
    private final DoubleBuf myRecursionWeightBuf;
    /**
     * Most recent lambda values for each Walker.
     */
    private final double[] currentLambdaValues;

    /**
     * Transition Tempered OSRW Constructor.
     *
     * @param lambdaInterface   defines Lambda and dU/dL.
     * @param potential         defines the Potential energy.
     * @param lambdaFile        contains the current Lambda particle position and velocity.
     * @param histogramFile     contains the Lambda and dU/dL histogram.
     * @param properties        defines System properties.
     * @param temperature       the simulation temperature.
     * @param dt                the time step.
     * @param printInterval     number of steps between logging updates.
     * @param saveInterval      number of steps between restart file updates.
     * @param asynchronous      set to true if walkers run asynchronously.
     * @param algorithmListener the AlgorithmListener to be notified of
     *                          progress.
     */
    public TransitionTemperedOSRW(LambdaInterface lambdaInterface, CrystalPotential potential,
                                  File lambdaFile, File histogramFile, CompositeConfiguration properties,
                                  double temperature, double dt, double printInterval,
                                  double saveInterval, boolean asynchronous,
                                  AlgorithmListener algorithmListener) {
        this(lambdaInterface, potential, lambdaFile, histogramFile, properties,
                temperature, dt, printInterval, saveInterval, asynchronous,
                true, algorithmListener);
    }

    /**
     * Transition Tempered OSRW Constructor.
     *
     * @param lambdaInterface   defines Lambda and dU/dL.
     * @param potential         defines the Potential energy.
     * @param lambdaFile        contains the current Lambda particle position and velocity.
     * @param histogramFile     contains the Lambda and dU/dL histogram.
     * @param properties        defines System properties.
     * @param temperature       the simulation temperature.
     * @param dt                the time step.
     * @param printInterval     number of steps between logging updates.
     * @param saveInterval      number of steps between restart file updates.
     * @param asynchronous      set to true if walkers run asynchronously.
     * @param resetNumSteps     whether to reset energy counts to 0
     * @param algorithmListener the AlgorithmListener to be notified of
     *                          progress.
     */
    public TransitionTemperedOSRW(LambdaInterface lambdaInterface, CrystalPotential potential,
                                  File lambdaFile, File histogramFile, CompositeConfiguration properties,
                                  double temperature, double dt, double printInterval,
                                  double saveInterval, boolean asynchronous, boolean resetNumSteps,
                                  AlgorithmListener algorithmListener) {

        this.lambdaInterface = lambdaInterface;
        this.potential = potential;
        this.lambdaFile = lambdaFile;
        this.histogramFile = histogramFile;
        this.temperature = temperature;
        this.asynchronous = asynchronous;
        this.algorithmListener = algorithmListener;
        nVariables = potential.getNumberOfVariables();

        if (potential instanceof Barostat) {
            barostat = (Barostat) potential;
        } else {
            barostat = null;
        }

        // Convert the time step to picoseconds.
        this.dt = dt * 0.001;

        // Convert the print interval to a print frequency.
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        // Convert the save interval to a save frequency.
        saveFrequency = 1000;
        if (saveInterval >= this.dt) {
            saveFrequency = (int) (saveInterval / this.dt);
        }

        biasCutoff = properties.getInt("lambda-bias-cutoff", 5);
        biasMag = properties.getDouble("bias-gaussian-mag", 0.05);
        dL = properties.getDouble("lambda-bin-width", 0.005);
        dFL = properties.getDouble("flambda-bin-width", 2.0);

        // Require modest sampling of the lambda path.
        if (dL > 0.1) {
            dL = 0.1;
        }

        /*
          Many lambda bin widths do not evenly divide into 1.0; here we correct
          for this by computing an integer number of bins, then re-setting the
          lambda variable appropriately. Note that we also choose to have an
          odd number of lambda bins, so that the centers of the first and last
          bin are at 0 and 1.
         */
        lambdaBins = (int) (1.0 / dL);
        if (lambdaBins % 2 == 0) {
            lambdaBins++;
        }

        /*
          The initial number of FLambda bins does not really matter, since a
          larger number is automatically allocated as needed. The center of the
          central bin is at 0.
         */
        FLambdaBins = 101;
        minFLambda = -(dFL * FLambdaBins) / 2.0;

        energyCount = -1;

        dL = 1.0 / (lambdaBins - 1);
        dL_2 = dL / 2.0;
        minLambda = -dL_2;
        dFL_2 = dFL / 2.0;
        maxFLambda = minFLambda + FLambdaBins * dFL;
        FLambda = new double[lambdaBins];
        dUdXdL = new double[nVariables];
        stochasticRandom = new Random();

        /*
          Set up the multi-walker communication variables for Parallel Java
          communication between nodes.
         */
        world = Comm.world();
        numProc = world.size();
        rank = world.rank();

        include1DBias = properties.getBoolean("osrw-1D-bias", true);

        deltaT = temperingFactor * R * this.temperature;

        // Allocate space for the recursion kernel that stores weights.
        recursionKernel = new double[lambdaBins][FLambdaBins];

        // Allocate space to regularize kernel values.
        kernelValues = new double[FLambdaBins];

        // Load the OSRW histogram restart file if it exists.
        boolean readHistogramRestart = false;
        if (histogramFile != null && histogramFile.exists()) {
            try {
                HistogramReader osrwHistogramReader = new HistogramReader(this, new FileReader(histogramFile));
                osrwHistogramReader.readHistogramFile();
                logger.info(format("\n Continuing OSRW histogram from %s.", histogramFile.getName()));
                readHistogramRestart = true;
            } catch (FileNotFoundException ex) {
                logger.info(" Histogram restart file could not be found and will be ignored.");
            }
        }

        // Load the OSRW lambda restart file if it exists.
        if (lambdaFile != null && lambdaFile.exists()) {
            try {
                LambdaReader osrwLambdaReader = new LambdaReader(this, new FileReader(lambdaFile));
                osrwLambdaReader.readLambdaFile(resetNumSteps);
                logger.info(format("\n Continuing OSRW lambda from %s.", lambdaFile.getName()));
            } catch (FileNotFoundException ex) {
                logger.info(" Lambda restart file could not be found and will be ignored.");
            }
        }

        currentLambdaValues = new double[world.size()];
        if (asynchronous) {
            // Use asynchronous communication.
            myRecursionWeight = new double[4];
            myRecursionWeightBuf = DoubleBuf.buffer(myRecursionWeight);
            receiveThread = new CountReceiveThread(this);
            receiveThread.start();
            recursionWeights = null;
            recursionWeightsBuf = null;
        } else {
            // Use synchronous communication.
            recursionWeights = new double[numProc][3];
            recursionWeightsBuf = new DoubleBuf[numProc];
            for (int i = 0; i < numProc; i++) {
                recursionWeightsBuf[i] = DoubleBuf.buffer(recursionWeights[i]);
            }
            myRecursionWeight = recursionWeights[rank];
            myRecursionWeightBuf = recursionWeightsBuf[rank];
            receiveThread = null;
        }

        double defaultOffset = 20.0 * biasMag;
        temperOffset = properties.getDouble("ttosrw-temperOffset", defaultOffset);
        if (temperOffset < 0.0) {
            temperOffset = 0.0;
        }
        logger.info(format("  Coverage before tempering:     %7.4g kcal/mol", temperOffset));

        String propString = properties.getString("ttosrw-integrationType", "SIMPSONS");
        IntegrationType testType;
        try {
            testType = IntegrationType.valueOf(propString.toUpperCase());
        } catch (Exception ex) {
            logger.warning(format(" Invalid argument %s to ttosrw-integrationType; resetting to SIMPSONS", propString));
            testType = SIMPSONS;
        }
        integrationType = testType;

        // Update and print out the recursion slave.
        if (readHistogramRestart) {
            updateFLambda(true, false);
        }

        // Log parameters.
        logger.info("\n Orthogonal Space Random Walk Parameters");
        logger.info(format("  Gaussian Bias Magnitude:       %6.4f (kcal/mol)", biasMag));
        logger.info(format("  Gaussian Bias Cutoff:           %6d bins", biasCutoff));
        logger.info(format("  Print Interval:                 %6.3f psec", printInterval));
        logger.info(format("  Save Interval:                  %6.3f psec", saveInterval));
        logger.info(format("  Include 1D OSRW Bias:           %6b", include1DBias));
    }

    /**
     * Sets the Dama et al tempering parameter, as a multiple of kBT.
     *
     * @param temper a double.
     */
    public void setTemperingParameter(double temper) {
        temperingFactor = temper;
        if (temperingFactor > 0.0) {
            deltaT = temperingFactor * R * temperature;
        } else {
            deltaT = Double.MAX_VALUE;
        }
    }

    /**
     * Check if tempering has started.
     *
     * @return True if tempering has begun.
     */
    boolean isTempering() {
        return tempering;
    }

    /**
     * Set the tempering flag.
     *
     * @param tempering True to indicate tempering.
     */
    void setTempering(boolean tempering) {
        this.tempering = tempering;
    }

    /**
     * Return the value of a recursion kernel bin.
     *
     * @param lambdaBin  The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @return The value of the bin.
     */
    double getRecursionKernelValue(int lambdaBin, int fLambdaBin) {
        return recursionKernel[lambdaBin][fLambdaBin];
    }

    /**
     * Set the value of a recursion kernel bin.
     *
     * @param lambdaBin  The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @param value      The value of the bin.
     */
    void setRecursionKernelValue(int lambdaBin, int fLambdaBin, double value) {
        recursionKernel[lambdaBin][fLambdaBin] = value;
    }

    /**
     * Add to the value of a recursion kernel bin.
     *
     * @param lambdaBin  The lambda bin.
     * @param fLambdaBin The dU/dL bin.
     * @param value      The value of the bin.
     */
    void addToRecursionKernelValue(int lambdaBin, int fLambdaBin, double value) {
        recursionKernel[lambdaBin][fLambdaBin] += value;
    }

    /**
     * Allocate memory for the recursion kernel.
     */
    void allocateRecursionKernel() {
        recursionKernel = new double[lambdaBins][FLambdaBins];
        kernelValues = new double[FLambdaBins];
    }

    /**
     * Set the OSRW Gaussian biasing potential magnitude (kcal/mol).
     *
     * @param biasMag Gaussian biasing potential magnitude (kcal/mol)
     */
    public void setBiasMagnitude(double biasMag) {
        this.biasMag = biasMag;
        logger.info(format("\n  Gaussian Bias Magnitude:       %6.4f (kcal/mol)", biasMag));

        double defaultOffset = 20.0 * biasMag;
        String propString = System.getProperty("ttosrw-temperOffset", Double.toString(defaultOffset));
        temperOffset = defaultOffset;
        try {
            temperOffset = Double.parseDouble(propString);
        } catch (NumberFormatException ex) {
            logger.info(format(" Exception in parsing ttosrw-temperOffset, resetting to 1.0 kcal/mol: %s", ex.toString()));
            temperOffset = defaultOffset;
        }
        if (temperOffset < 0.0) {
            temperOffset = 0.0;
        }
        logger.info(format("  Coverage before tempering:     %6.4f (kcal/mol)", temperOffset));
    }

    /**
     * Compute the force field + bias energy.
     */
    public double energy(double[] x) {

        forceFieldEnergy = potential.energy(x);

        // OSRW is propagated with the slowly varying terms.
        if (state == Potential.STATE.FAST) {
            return forceFieldEnergy;
        }

        gLdEdL = 0.0;
        dUdLambda = lambdaInterface.getdEdL();
        d2UdL2 = lambdaInterface.getd2EdL2();
        int lambdaBin = binForLambda(lambda);
        int FLambdaBin = binForFLambda(dUdLambda);
        dForceFieldEnergydL = dUdLambda;

        // Calculate recursion kernel G(L, F_L) and its derivatives with respect to L and F_L.
        double dGdLambda = 0.0;
        double dGdFLambda = 0.0;
        double ls2 = (2.0 * dL) * (2.0 * dL);
        double FLs2 = (2.0 * dFL) * (2.0 * dFL);
        int kernelCount = 0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int lcenter = lambdaBin + iL;
            double deltaL = lambda - (lcenter * dL);
            double deltaL2 = deltaL * deltaL;
            // Mirror conditions for recursion kernel counts.
            int lcount = lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                mirrorFactor = 2.0;
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

                // If either of the following FL edge conditions are true, then there are no counts and we continue.
                if (FLcenter < 0 || FLcenter >= FLambdaBins) {
                    continue;
                }
                double deltaFL = dUdLambda - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double rc = recursionKernel[lcount][FLcenter];
                if (rc > 0.0) {
                    kernelCount++;
                } else {
                    continue;
                }
                double weight = mirrorFactor * rc;
                double bias = weight * biasMag
                        * exp(-deltaL2 / (2.0 * ls2))
                        * exp(-deltaFL2 / (2.0 * FLs2));
                gLdEdL += bias;
                dGdLambda -= deltaL / ls2 * bias;
                dGdFLambda -= deltaFL / FLs2 * bias;
            }
        }

        // Up-weight the 2D bias to account for bins outside the Wall who have no weight.
//        if (hardWallConstraint && kernelCount > 0) {
//            if (hardWallConstraint && propagateLambda) {
//                logger.severe(" The Hard Wall constraint is not supported for MD-OST.");
//            }
//            int maxCount = (2 * biasCutoff + 1);
//            maxCount *= maxCount;
//            double upweightFactor = (double) maxCount / (double) kernelCount;
//            dGdLambda *= upweightFactor;
//            dGdFLambda *= upweightFactor;
//        }

        // Lambda gradient due to recursion kernel G(L, F_L).
        dUdLambda += dGdLambda + dGdFLambda * d2UdL2;

        // Compute the energy and gradient for the recursion slave at F(L) using interpolation.
        double bias1D = 0.0;
        if (include1DBias) {
            bias1D = current1DBiasEnergy(lambda, true);
        }
        biasEnergy = bias1D + gLdEdL;

        if (print) {
            logger.info(format(" Bias Energy        %16.8f", biasEnergy));
            logger.info(format(" %s %16.8f  (Kcal/mole)", "OSRW Potential    ", forceFieldEnergy + biasEnergy));
        }

        if (propagateLambda) {
            energyCount++;

            // Log the current Lambda state.
            if (energyCount % printFrequency == 0) {
                double dBdL = dUdLambda - dForceFieldEnergydL;
                if (lambdaBins < 1000) {
                    logger.info(format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                            lambda, lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, halfThetaVelocity));
                } else {
                    logger.info(format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                            lambda, lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, halfThetaVelocity));
                }
            }

            // Metadynamics grid counts (every 'countInterval' steps).
            if (energyCount % countInterval == 0) {
                addBias(dForceFieldEnergydL, x, null);
            }

            langevin();
        }

        totalEnergy = forceFieldEnergy + biasEnergy;

        return totalEnergy;
    }

    /**
     * Compute the Bias energy at (currentLambda, currentdUdL).
     *
     * @param currentLambda The value of lambda.
     * @param currentdUdL   The value of dU/dL.
     * @return The bias energy.
     */
    double computeBiasEnergy(double currentLambda, double currentdUdL) {

        int lambdaBin = binForLambda(currentLambda);
        int FLambdaBin = binForFLambda(currentdUdL);

        double bias2D = 0.0;
        double ls2 = (2.0 * dL) * (2.0 * dL);
        double FLs2 = (2.0 * dFL) * (2.0 * dFL);
        int kernelCount = 0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int lcenter = lambdaBin + iL;
            double deltaL = currentLambda - (lcenter * dL);
            double deltaL2 = deltaL * deltaL;
            // Mirror conditions for recursion kernel counts.
            int lcount = lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                mirrorFactor = 2.0;
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
                double currentFL = minFLambda + FLcenter * dFL + dFL_2;
                double deltaFL = currentdUdL - currentFL;
                double deltaFL2 = deltaFL * deltaFL;
                double rc = recursionKernel[lcount][FLcenter];
                if (rc > 0.0) {
                    kernelCount++;
                } else {
                    continue;
                }
                double weight = mirrorFactor * rc;
                double bias = weight * biasMag
                        * exp(-deltaL2 / (2.0 * ls2))
                        * exp(-deltaFL2 / (2.0 * FLs2));

                // logger.info(format("(L=%6.4f FL=%8.2f) L=%6.4f Bin=%3d; FL=%8.3f Bin=%6d; Bias: %8.6f",
                //        currentLambda, currentdUdL, lcenter * dL, lcount, currentFL, FLcenter, bias));

                bias2D += bias;
            }
        }

        // Up-weight the 2D bias to account for bins outside the Wall who have no weight.
//        if (hardWallConstraint && kernelCount > 0) {
//            int maxCount = (2 * biasCutoff + 1);
//            maxCount *= maxCount;
//            double upweightFactor = (double) maxCount / (double) kernelCount;
//            bias2D *= upweightFactor;
//        }

        // Compute the energy for the recursion slave at F(L) using interpolation.
        double bias1D = 0.0;
        if (include1DBias) {
            bias1D = current1DBiasEnergy(currentLambda, false);
        }

        return bias1D + bias2D;
    }

    /**
     * Run a local optimization beginning at the given coordinates.
     *
     * @param e        Current energy before optimization.
     * @param x        Current coordinates.
     * @param gradient Work array to store the gradient.
     */
    private void optimization(double e, double[] x, double[] gradient) {
        if (energyCount % osrwOptimizationFrequency == 0) {
            logger.info(format(" OSRW Minimization (Step %d)", energyCount));

            // Set the underlying Potential's Lambda value to 1.0.
            lambdaInterface.setLambda(1.0);

            potential.setEnergyTermState(Potential.STATE.BOTH);

            boolean origBaroActive = true;
            if (barostat != null) {
                origBaroActive = barostat.isActive();
                barostat.setActive(false);
            }

            try {
                // Optimize the system.
                double startingEnergy = potential.energy(x);

                Minimize minimize = new Minimize(molecularAssembly, potential, algorithmListener);
                minimize.minimize(osrwOptimizationEps);

                // Collect the minimum energy.
                double minEnergy = potential.getTotalEnergy();
                // Check for a new minimum within an energy window of the lowest energy structure found.
                if (minEnergy < osrwOptimum + osrwOptimizationEnergyWindow) {
                    if (minEnergy < osrwOptimum) {
                        osrwOptimum = minEnergy;
                    }
                    int n = potential.getNumberOfVariables();
                    osrwOptimumCoords = new double[n];
                    osrwOptimumCoords = potential.getCoordinates(osrwOptimumCoords);
                    double mass = molecularAssembly.getMass();
                    Crystal crystal = molecularAssembly.getCrystal();
                    double density = crystal.getDensity(mass);
                    osrwOptimizationFilter.writeFile(osrwOptimizationFile, false);
                    Crystal uc = crystal.getUnitCell();
                    logger.info(format(" Minimum: %12.6f %s (%12.6f g/cc) optimized from %12.6f at step %d.",
                            minEnergy, uc.toShortString(), density, startingEnergy, energyCount));
                }
            } catch (EnergyException ex) {
                String message = ex.getMessage();
                logger.info(format(" Energy exception minimizing coordinates at lambda=%8.6f\n %s.", lambda, message));
                logger.info(" TT-OSRW sampling will continue.");
            }

            // Set the underlying Potential's Lambda value back to current lambda value.
            lambdaInterface.setLambda(lambda);

            // Reset the Potential State
            potential.setEnergyTermState(state);

            // Reset the Barostat
            if (barostat != null) {
                barostat.setActive(origBaroActive);
            }

            // Revert to the coordinates and gradient prior to optimization.
            double eCheck = potential.energyAndGradient(x, gradient);

            if (abs(eCheck - e) > osrwOptimizationTolerance) {
                logger.warning(format(
                        " TT-OSRW optimization could not revert coordinates %16.8f vs. %16.8f.", e, eCheck));
            }
        }
    }

    /**
     * Evaluate the kernel across dU/dL values at given value of lambda. The largest kernel value V
     * is used to define an offset (-V), which is added to all to kernel values. Then,
     * the largest kernel value is zero, and will result in a statistical weight of 1.
     * <p>
     * The offset avoids the recursion kernel becoming too large for some combinations of (Lambda, dU/dL).
     * This can result in its statistical weight = exp(kernel * beta)
     * exceeding the maximum representable double value.
     *
     * @param lambda Value of Lambda to evaluate the kernal for.
     * @param llFL   Lower FLambda bin.
     * @param ulFL   Upper FLambda bin.
     * @return the applied offset.
     */
    private double evaluateKernelforLambda(int lambda, int llFL, int ulFL) {
        double maxKernel = Double.MIN_VALUE;
        for (int jFL = llFL; jFL <= ulFL; jFL++) {
            double value = evaluateKernel(lambda, jFL);
            kernelValues[jFL] = value;
            if (value > maxKernel) {
                maxKernel = value;
            }
        }
        double offset = -maxKernel;
        for (int jFL = llFL; jFL <= ulFL; jFL++) {
            kernelValues[jFL] += offset;
        }
        return offset;
    }

    /**
     * Integrates dUdL over lambda using more sophisticated techniques than midpoint rectangular integration.
     * <p>
     * The ends (from 0 to dL and 1-dL to 1) are integrated with trapezoids
     *
     * @param dUdLs dUdL at the midpoint of each bin.
     * @param type  Integration type to use.
     * @return Current delta-G estimate.
     */
    private double integrateNumeric(double[] dUdLs, IntegrationType type) {
        // Integrate between the second bin midpoint and the second-to-last bin midpoint.
        double[] midLams = Integrate1DNumeric.generateXPoints(dL, 1.0 - dL, (lambdaBins - 2), false);
        double[] midVals = Arrays.copyOfRange(dUdLs, 1, (lambdaBins - 1));
        DataSet dSet = new DoublesDataSet(midLams, midVals, false);

        double val = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, type);

        double dL_4 = dL_2 * 0.5;

        // Initially, assume dU/dL is exactly 0 at the endpoints. This is sometimes a true assumption.
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

        // Integrate trapezoids from 0 to the second bin midpoint, and from second-to-last bin midpoint to 1.
        val += trapezoid(0, dL_4, val0, dUdLs[0]);
        val += trapezoid(dL_4, dL, dUdLs[0], dUdLs[1]);
        val += trapezoid(1.0 - dL, 1.0 - dL_4, dUdLs[lambdaBins - 2], dUdLs[lambdaBins - 1]);
        val += trapezoid(1.0 - dL_4, 1.0, dUdLs[lambdaBins - 1], val1);

        return val;
    }

    /**
     * Integrates a trapezoid.
     *
     * @param x0  First x point
     * @param x1  Second x point
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
     * Update a local array of current lambda values for each walker.
     *
     * @param rank   Walker rank.
     * @param lambda Walker's current lambda value.
     */
    void setCurrentLambdaforRank(int rank, double lambda) {
        currentLambdaValues[rank] = lambda;
    }

    /**
     * Send an OSRW count to all other processes while also receiving an OSRW
     * count from all other processes.
     *
     * @param lambda Current value of lambda.
     * @param dUdL   Current value of dU/dL.
     */
    private void synchronousSend(double lambda, double dUdL) {

        // All-Gather counts from each walker.
        myRecursionWeight[0] = lambda;
        myRecursionWeight[1] = dUdL;
        myRecursionWeight[2] = temperingWeight;
        try {
            world.allGather(myRecursionWeightBuf, recursionWeightsBuf);
        } catch (IOException ex) {
            String message = " Multi-walker OSRW allGather failed.";
            logger.log(Level.SEVERE, message, ex);
        }

        // Find the minimum and maximum FLambda bin for the gathered counts.
        double minRequired = Double.MAX_VALUE;
        double maxRequired = Double.MIN_VALUE;
        for (int i = 0; i < numProc; i++) {
            minRequired = min(minRequired, recursionWeights[i][1]);
            maxRequired = max(maxRequired, recursionWeights[i][1]);
        }

        // Check that the FLambda range of the Recursion kernel includes both the minimum and maximum FLambda value.
        checkRecursionKernelSize(minRequired);
        checkRecursionKernelSize(maxRequired);

        // Increment the Recursion Kernel based on the input of each walker.
        for (int i = 0; i < numProc; i++) {
            currentLambdaValues[i] = recursionWeights[i][0];
            int walkerLambda = binForLambda(recursionWeights[i][0]);
            int walkerFLambda = binForFLambda(recursionWeights[i][1]);
            double weight = recursionWeights[i][2];

            // If the weight is less than 1.0, then a walker has activated tempering.
            if (!tempering && weight < 1.0) {
                tempering = true;
                logger.info(format(" Tempering activated due to received weight of (%8.6f)", weight));
            }

            if (resetStatistics && recursionWeights[i][0] > lambdaResetValue) {
                recursionKernel = new double[lambdaBins][FLambdaBins];
                resetStatistics = false;
                logger.info(format(" Cleared OSRW histogram (Lambda = %6.4f).", recursionWeights[i][0]));
            }

            recursionKernel[walkerLambda][walkerFLambda] += weight;
        }
    }

    /**
     * Send an OSRW count to all other processes.
     *
     * @param lambda Current value of lambda.
     * @param dUdL   Current value of dU/dL.
     */
    private void asynchronousSend(double lambda, double dUdL) {
        myRecursionWeight[0] = world.rank();
        myRecursionWeight[1] = lambda;
        myRecursionWeight[2] = dUdL;
        myRecursionWeight[3] = temperingWeight;

        for (int i = 0; i < numProc; i++) {
            try {
                world.send(i, myRecursionWeightBuf);
            } catch (Exception ex) {
                String message = " Asynchronous Multiwalker OSRW send failed.";
                logger.log(Level.SEVERE, message, ex);
            }
        }
    }

    /**
     * Evaluate the bias at [cLambda, cF_lambda]
     */
    private double evaluateKernel(int cLambda, int cF_Lambda) {
        // Compute the value of L and FL for the center of the current bin.
        double vL = cLambda * dL;
        double vFL = minFLambda + cF_Lambda * dFL + dFL_2;

        // Set the variances for the Gaussian bias.
        double Ls2 = 2.0 * dL * 2.0 * dL;
        double FLs2 = 2.0 * dFL * 2.0 * dFL;

        // Variances are only used when dividing by twice their value, so pre-compute!
        double invLs2 = 0.5 / Ls2;
        double invFLs2 = 0.5 / FLs2;

        double sum = 0.0;
        int kernelCount = 0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int Lcenter = cLambda + iL;
            double deltaL = vL - Lcenter * dL;
            double deltaL2 = deltaL * deltaL;

            // Pre-compute the lambda-width exponent.
            double L2exp = exp(-deltaL2 * invLs2);

            // Mirror condition for Lambda counts.
            int lcount = Lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                // The width of the first and last bins is dLambda_2,
                // so the mirror condition is to double their counts.
                mirrorFactor = 2.0;
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
                double deltaFL = vFL - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double rc = recursionKernel[lcount][FLcenter];
                if (rc > 0.0) {
                    kernelCount++;
                } else {
                    continue;
                }
                double weight = mirrorFactor * rc;
                if (weight > 0) {
                    double e = weight * biasMag * L2exp * exp(-deltaFL2 * invFLs2);
                    sum += e;
                }
            }
        }

        // Up-weight the 2D bias to account for bins outside the Wall who have no weight.
//        if (hardWallConstraint && kernelCount > 0) {
//            if (hardWallConstraint && propagateLambda) {
//                logger.severe(" The Hard Wall constraint is not supported for MD-OST.");
//            }
//            int maxCount = (2 * biasCutoff + 1);
//            maxCount *= maxCount;
//            double upweightFactor = (double) maxCount / (double) kernelCount;
//            sum *= upweightFactor;
//        }

        return sum;
    }

    /**
     * <p>evaluateHistogram.</p>
     *
     * @param lambda the lambda value.
     * @param dUdL   the dU/dL value.
     * @return The value of the Histogram.
     */
    double evaluateHistogram(double lambda, double dUdL) {
        int lambdaBin = binForLambda(lambda);
        int dUdLBin = binForFLambda(dUdL);
        try {
            return recursionKernel[lambdaBin][dUdLBin];
        } catch (Exception e) {
            // Catch an index out of bounds exception.
            return 0.0;
        }
    }

    /**
     * If necessary, allocate more space.
     */
    void checkRecursionKernelSize(double dEdLambda) {
        if (dEdLambda > maxFLambda) {
            logger.info(format(" Current F_lambda %8.2f > maximum histogram size %8.2f.",
                    dEdLambda, maxFLambda));

            double origDeltaG = updateFLambda(false, false);

            int newFLambdaBins = FLambdaBins;
            while (minFLambda + newFLambdaBins * dFL < dEdLambda) {
                newFLambdaBins += 100;
            }
            double[][] newRecursionKernel = new double[lambdaBins][newFLambdaBins];

            // We have added bins above the indeces of the current counts just copy them into the new array.
            for (int i = 0; i < lambdaBins; i++) {
                arraycopy(recursionKernel[i], 0, newRecursionKernel[i], 0, FLambdaBins);
            }
            recursionKernel = newRecursionKernel;
            FLambdaBins = newFLambdaBins;
            kernelValues = new double[FLambdaBins];
            maxFLambda = minFLambda + dFL * FLambdaBins;
            logger.info(format(" New histogram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            double newFreeEnergy = updateFLambda(false, false);
            assert (origDeltaG == newFreeEnergy);
        }
        if (dEdLambda < minFLambda) {
            logger.info(format(" Current F_lambda %8.2f < minimum histogram size %8.2f.",
                    dEdLambda, minFLambda));

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

            logger.info(format(" New histogram %8.2f to %8.2f with %d bins.\n",
                    minFLambda, maxFLambda, FLambdaBins));

            double newFreeEnergy = updateFLambda(false, false);
            assert (origDeltaG == newFreeEnergy);
        }
    }

    /**
     * Eqs. 7 and 8 from the 2012 Crystal Thermodynamics paper.
     */
    public double updateFLambda(boolean print, boolean save) {
        double freeEnergy = 0.0;
        double minFL = Double.MAX_VALUE;

        // Total histogram weight.
        double totalWeight = 0;
        double beta = 1.0 / (R * temperature);
        StringBuilder stringBuilder = new StringBuilder();

        for (int iL = 0; iL < lambdaBins; iL++) {
            int ulFL = -1;
            int llFL = -1;

            // Find the smallest FL bin.
            for (int jFL = 0; jFL < FLambdaBins; jFL++) {
                double count = recursionKernel[iL][jFL];
                if (count > 0) {
                    llFL = jFL;
                    break;
                }
            }

            // Find the largest FL bin.
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

                    double weight = exp(kernel * beta);
                    partitionFunction += weight;

                    double currentFLambda = minFLambda + jFL * dFL + dFL_2;
                    ensembleAverageFLambda += currentFLambda * weight;
                    lambdaCount += recursionKernel[iL][jFL];
                }
                if (minFL > maxBias) {
                    minFL = maxBias;
                }
                FLambda[iL] = ensembleAverageFLambda / partitionFunction;
                lla = minFLambda + llFL * dFL;
                ula = minFLambda + (ulFL + 1) * dFL;
            }

            // The first and last lambda bins are half size.
            double delta = dL;
            if (iL == 0 || iL == lambdaBins - 1) {
                delta = dL_2;
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

                double midLambda = (llL + ulL) / 2.0;
                double bias1D = current1DBiasEnergy(midLambda, false);
                double bias2D = computeBiasEnergy(midLambda, FLambda[iL]) - bias1D;

                stringBuilder.append(format(" %6.2e %7.5f %7.1f %7.1f %8.2f %8.2f %8.2f %8.2f %8.2f   %8.2f\n",
                        lambdaCount, midLambda, lla, ula, FLambda[iL], bias1D, bias2D, bias1D + bias2D,
                        freeEnergy, bias1D + bias2D + freeEnergy));
            }
        }

        if (tempering) {
            double temperEnergy = (minFL > temperOffset) ? temperOffset - minFL : 0;
            temperingWeight = exp(temperEnergy / deltaT);
        }

        freeEnergy = integrateNumeric(FLambda, integrationType);

        if (print && abs(freeEnergy - previousFreeEnergy) > 0.001) {
            logger.info("  Weight   Lambda      dU/dL Bins  <dU/dL>    g(L)  f(L,<dU/dL>) Bias    dG(L) Bias+dG(L)");
            logger.info(stringBuilder.toString());
            logger.info(format(" The free energy is %12.4f kcal/mol (Total Weight: %6.2e, Tempering: %6.4f, Counts: %12d).",
                    freeEnergy, totalWeight, temperingWeight, biasCount));
            logger.info(format(" Minimum Bias %8.3f", minFL));
            previousFreeEnergy = freeEnergy;
        } else if (!save && (print || biasCount % printFrequency == 0)) {
            logger.info(format(" The free energy is %12.4f kcal/mol (Total Weight: %6.2e, Tempering: %6.4f, Counts: %12d).",
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
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] gradient) {

        forceFieldEnergy = potential.energyAndGradient(x, gradient);

        // OSRW is propagated with the slowly varying terms.
        if (state == STATE.FAST) {
            return forceFieldEnergy;
        }

        gLdEdL = 0.0;
        dForceFieldEnergydL = lambdaInterface.getdEdL();
        dUdLambda = dForceFieldEnergydL;

        d2UdL2 = lambdaInterface.getd2EdL2();
        int lambdaBin = binForLambda(lambda);
        int FLambdaBin = binForFLambda(dUdLambda);

        // First derivative of the 2D bias with respect to Lambda.
        double dGdLambda = 0.0;
        // First derivative of the 2D bias with respect to  dU/dL.
        double dGdFLambda = 0.0;
        // Calculate recursion kernel G(L, F_L) and its derivatives with respect to L and F_L.
        double ls2 = (2.0 * dL) * (2.0 * dL);
        double FLs2 = (2.0 * dFL) * (2.0 * dFL);
        int kernelCount = 0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int lcenter = lambdaBin + iL;
            double deltaL = lambda - (lcenter * dL);
            double deltaL2 = deltaL * deltaL;
            // Mirror conditions for recursion kernel counts.
            int lcount = lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                mirrorFactor = 2.0;
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
                double deltaFL = dUdLambda - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double rc = recursionKernel[lcount][FLcenter];
                if (rc > 0.0) {
                    kernelCount++;
                } else {
                    continue;
                }
                double weight = mirrorFactor * rc;
                double bias = weight * biasMag
                        * exp(-deltaL2 / (2.0 * ls2))
                        * exp(-deltaFL2 / (2.0 * FLs2));
                gLdEdL += bias;
                dGdLambda -= deltaL / ls2 * bias;
                dGdFLambda -= deltaFL / FLs2 * bias;
            }
        }

        // Up-weight the 2D bias to account for bins outside the Wall who have no weight.
//        if (hardWallConstraint && kernelCount > 0) {
//            if (hardWallConstraint && propagateLambda) {
//                logger.severe(" The Hard Wall constraint is not supported for MD-OST.");
//            }
//            int maxCount = (2 * biasCutoff + 1);
//            maxCount *= maxCount;
//            double upweightFactor = (double) maxCount / (double) kernelCount;
//            dGdLambda *= upweightFactor;
//            dGdFLambda *= upweightFactor;
//        }

        // Lambda gradient due to recursion kernel G(L, F_L).
        dUdLambda += dGdLambda + dGdFLambda * d2UdL2;

        // Cartesian coordinate gradient due to recursion kernel G(L, F_L).
        fill(dUdXdL, 0.0);
        lambdaInterface.getdEdXdL(dUdXdL);
        for (int i = 0; i < nVariables; i++) {
            gradient[i] += dGdFLambda * dUdXdL[i];
        }

        // Compute the energy and gradient for the recursion slave at F(L) using interpolation.
        double bias1D = 0.0;
        if (include1DBias) {
            bias1D = current1DBiasEnergy(lambda, true);
        }
        biasEnergy = bias1D + gLdEdL;

        if (print) {
            logger.info(format(" %s %16.8f", "Bias Energy       ", biasEnergy));
            logger.info(format(" %s %16.8f  %s",
                    "OSRW Potential    ", forceFieldEnergy + biasEnergy, "(Kcal/mole)"));
        }

        if (propagateLambda) {
            energyCount++;

            // Log the current Lambda state.
            if (energyCount % printFrequency == 0) {
                double dBdL = dUdLambda - dForceFieldEnergydL;
                if (lambdaBins < 1000) {
                    logger.info(format(" L=%6.4f (%3d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                            lambda, lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, halfThetaVelocity));
                } else {
                    logger.info(format(" L=%6.4f (%4d) F_LU=%10.4f F_LB=%10.4f F_L=%10.4f V_L=%10.4f",
                            lambda, lambdaBin, dForceFieldEnergydL, dBdL, dUdLambda, halfThetaVelocity));
                }
            }

            // Metadynamics grid counts (every 'countInterval' steps).
            if (energyCount % countInterval == 0) {
                addBias(dForceFieldEnergydL, x, gradient);
            }

            langevin();
        }

        totalEnergy = forceFieldEnergy + biasEnergy;

        return totalEnergy;
    }

    /**
     * Add a Gaussian hill to the Histogram at (lambda, dEdU).
     */
    void addBias(double dEdU, double[] x, double[] gradient) {
        if (asynchronous) {
            asynchronousSend(lambda, dEdU);
        } else {
            synchronousSend(lambda, dEdU);
        }

        biasCount++;

        // Update F(L)
        fLambdaUpdates++;
        boolean printFLambda = fLambdaUpdates % fLambdaPrintInterval == 0;
        updateFLambda(printFLambda, false);

        if (osrwOptimization && lambda > osrwOptimizationLambdaCutoff) {
            if (gradient == null) {
                gradient = new double[x.length];
            }
            optimization(forceFieldEnergy, x, gradient);
        }

        // Write out restart files.
        if (energyCount > 0 && energyCount % saveFrequency == 0) {
            writeRestart();
        }
    }

    /**
     * Write histogram and lambda restart files.
     */
    void writeRestart() {
        if (algorithmListener != null) {
            algorithmListener.algorithmUpdate(molecularAssembly);
        }

        // Only the rank 0 process writes the histogram restart file.
        if (rank == 0) {
            StringBuilder stringBuilder = new StringBuilder(" Current Lambda Values:");
            for (double lambda : currentLambdaValues) {
                stringBuilder.append(format(" %6.4f", lambda));
            }
            logger.info(stringBuilder.toString());
            try {
                HistogramWriter ttOSRWHistogramRestart = new HistogramWriter(this,
                        new BufferedWriter(new FileWriter(histogramFile)));
                ttOSRWHistogramRestart.writeHistogramFile();
                ttOSRWHistogramRestart.flush();
                ttOSRWHistogramRestart.close();
                logger.info(format(" Wrote TTOSRW histogram restart file to %s.", histogramFile.getName()));
            } catch (IOException ex) {
                String message = " Exception writing TTOSRW histogram restart file.";
                logger.log(Level.INFO, Utilities.stackTraceToString(ex));
                logger.log(Level.SEVERE, message, ex);
            }
        }

        // All ranks write a lambda restart file.
        try {
            LambdaWriter ttOSRWLambdaRestart = new LambdaWriter(this, new BufferedWriter(new FileWriter(lambdaFile)));
            ttOSRWLambdaRestart.writeLambdaFile();
            ttOSRWLambdaRestart.flush();
            ttOSRWLambdaRestart.close();
            logger.info(format(" Wrote TTOSRW lambda restart file to %s.", lambdaFile.getName()));
        } catch (IOException ex) {
            String message = format(" Exception writing TTOSRW lambda restart file %s.", lambdaFile);
            logger.log(Level.INFO, Utilities.stackTraceToString(ex));
            logger.log(Level.SEVERE, message, ex);
        }
    }

    /**
     * Set a threshold to control writing of coordinate snapshots.
     *
     * @param lambdaWriteOut
     */
    void setLambdaWriteOut(double lambdaWriteOut) {
        this.lambdaWriteOut = lambdaWriteOut;
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Set lambda write out threshold to %6.3f lambda", lambdaWriteOut));
        }
    }

    /**
     * <p>checkRecursionKernelSize.</p>
     */
    public void checkRecursionKernelSize() {
        double[] x = new double[nVariables];
        x = potential.getCoordinates(x);
        double[] g = new double[nVariables];
        potential.energyAndGradient(x, g, false);
        double dudl = lambdaInterface.getdEdL();
        checkRecursionKernelSize(dudl);
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
     * Propagate Lambda using Langevin dynamics.
     */
    private void langevin() {
        // Compute the random force pre-factor (kcal/mol * psec^-2).
        double rt2 = 2.0 * Constants.R * temperature * thetaFriction / dt;

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
        lambda = sinTheta * sinTheta;
        lambdaInterface.setLambda(lambda);
    }

    /**
     * <p>binForLambda.</p>
     *
     * @param lambda a double.
     * @return a int.
     */
    int binForLambda(double lambda) {
        int lambdaBin = (int) floor((lambda - minLambda) / dL);
        if (lambdaBin < 0) {
            lambdaBin = 0;
        }
        if (lambdaBin >= lambdaBins) {
            lambdaBin = lambdaBins - 1;
        }
        return lambdaBin;
    }

    /**
     * <p>binForFLambda.</p>
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
     * <p>getForceFielddEdL.</p>
     *
     * @return a double.
     */
    double getForceFielddEdL() {
        return dForceFieldEnergydL;
    }

    /**
     * <p>getTotaldEdLambda.</p>
     *
     * @return a double.
     */
    public double getTotaldEdLambda() {
        return dUdLambda;
    }

    /**
     * <p>Getter for the field <code>forceFieldEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getForceFieldEnergy() {
        return forceFieldEnergy;
    }

    /**
     * <p>Getter for the field <code>biasEnergy</code>.</p>
     *
     * @return a double.
     */
    double getBiasEnergy() {
        return biasEnergy;
    }

    /**
     * <p>getPotentialEnergy.</p>
     *
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential getPotentialEnergy() {
        return potential;
    }

    /**
     * This calculates the 1D OSRW bias and its derivative with respect to Lambda.
     * <p>
     * See Equation 8 in http://doi.org/10.1021/ct300035u.
     *
     * @return a double.
     */
    private double current1DBiasEnergy(double currentLambda, boolean gradient) {
        if (!include1DBias) {
            return 0.0;
        }
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
     * <p>evaluate2DPMF.</p>
     *
     * @return A StringBuffer with 2D Bias PMF.
     */
    public StringBuffer evaluate2DPMF() {
        StringBuffer sb = new StringBuffer();
        for (int fLambdaBin = 0; fLambdaBin < FLambdaBins; fLambdaBin++) {
            double currentFL = minFLambda + fLambdaBin * dFL + dFL_2;
            sb.append(format(" %16.8f", currentFL));
            for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
                double bias = -evaluateKernel(lambdaBin, fLambdaBin);
                sb.append(format(" %16.8f", bias));
            }
            sb.append("\n");
        }
        return sb;
    }

    /**
     * <p>evaluateTotalPMF.</p>
     *
     * @return A StringBuffer the total 2D PMF.
     */
    public StringBuffer evaluateTotalPMF() {
        StringBuffer sb = new StringBuffer();
        for (int fLambdaBin = 0; fLambdaBin < FLambdaBins; fLambdaBin++) {

            double currentFL = minFLambda + fLambdaBin * dFL + dFL_2;
            sb.append(format(" %16.8f", currentFL));

            for (int lambdaBin = 0; lambdaBin < lambdaBins; lambdaBin++) {
                lambda = lambdaBin * dL + dL_2;
                double bias1D = -current1DBiasEnergy(lambda, false);
                double totalBias = bias1D - evaluateKernel(lambdaBin, fLambdaBin);
                sb.append(format(" %16.8f", totalBias));
            }
            sb.append("\n");
        }
        return sb;
    }

    /**
     * <p>Setter for the field <code>lambda</code>.</p>
     *
     * @param lambda a double.
     */
    public void setLambda(double lambda) {
        lambdaInterface.setLambda(lambda);
        this.lambda = lambda;
        theta = asin(sqrt(lambda));
    }

    /**
     * <p>Getter for the field <code>lambda</code>.</p>
     *
     * @return a double.
     */
    public double getLambda() {
        return lambda;
    }

    /**
     * <p>Getter for the field <code>lambdaInterface</code>.</p>
     *
     * @return a {@link ffx.potential.bonded.LambdaInterface} object.
     */
    public LambdaInterface getLambdaInterface() {
        return lambdaInterface;
    }

    /**
     * <p>Setter for the field <code>thetaMass</code>.</p>
     *
     * @param thetaMass a double.
     */
    public void setThetaMass(double thetaMass) {
        this.thetaMass = thetaMass;
    }

    /**
     * <p>setThetaFrication.</p>
     *
     * @param thetaFriction a double.
     */
    public void setThetaFrication(double thetaFriction) {
        this.thetaFriction = thetaFriction;
    }

    /**
     * Set the OSRW count interval. Every 'countInterval' steps the
     * recursionKernel will be incremented based on the current value of the
     * lambda state variable and the derivative of the energy with respect to
     * lambda (dU/dL).
     *
     * @param countInterval Molecular dynamics steps between counts.
     */
    public void setCountInterval(int countInterval) {
        if (countInterval > 0) {
            this.countInterval = countInterval;
        } else {
            logger.info(" OSRW count interval must be greater than 0.");
        }
    }

    /**
     * <p>getOSRWOptimumEnergy.</p>
     *
     * @return a double.
     */
    public double getOSRWOptimumEnergy() {
        if (osrwOptimum == Double.MAX_VALUE) {
            logger.info("Lambda optimization cutoff was not reached. Try increasing the number of timesteps.");
        }
        return osrwOptimum;
    }

    /**
     * <p>getOSRWOptimumCoordinates.</p>
     *
     * @return an array of {@link double} objects.
     */
    public double[] getOSRWOptimumCoordinates() {
        if (osrwOptimum < Double.MAX_VALUE) {
            return osrwOptimumCoords;
        } else {
            logger.info("Lambda optimization cutoff was not reached. Try increasing the number of timesteps.");
            return null;
        }
    }

    public void setMolecularAssembly(MolecularAssembly molecularAssembly) {
        this.molecularAssembly = molecularAssembly;
    }

    /**
     * <p>setOptimization.</p>
     *
     * @param osrwOptimization  a boolean.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     */
    public void setOptimization(boolean osrwOptimization, MolecularAssembly molecularAssembly) {
        this.osrwOptimization = osrwOptimization;
        this.molecularAssembly = molecularAssembly;
        File file = this.molecularAssembly.getFile();

        String fileName = FilenameUtils.removeExtension(file.getAbsolutePath());
        String ext = FilenameUtils.getExtension(file.getAbsolutePath());

        if (osrwOptimizationFilter == null) {
            if (ext.toUpperCase().contains("XYZ")) {
                osrwOptimizationFile = new File(fileName + "_opt.xyz");
                osrwOptimizationFilter = new XYZFilter(osrwOptimizationFile, this.molecularAssembly, null, null);
            } else {
                osrwOptimizationFile = new File(fileName + "_opt.pdb");
                osrwOptimizationFilter = new PDBFilter(osrwOptimizationFile, this.molecularAssembly, null, null);
            }
        }
    }

    /**
     * Returns the number of energy evaluations performed by this ttOSRW,
     * including those picked up in the lambda file.
     *
     * @return Number of energy steps taken by this walker.
     */
    public int getEnergyCount() {
        return energyCount;
    }

    /**
     * If this flag is true, (lambda, dU/dL) Monte Carlo samples that have no weight in the Histogram are rejected.
     */
    public void setHardWallConstraint(boolean hardWallConstraint) {
        this.hardWallConstraint = hardWallConstraint;
    }

    /**
     * If the dUdLHardWall flag is set to true, this method will
     * return false if the (lambda, dU/dL) sample is has not been seen.
     *
     * @param lambda The proposed lambda value.
     * @param dUdL   The proposed dU/dL value.
     * @return Returns false only if the dUdLHardWall flag is true, and the (lambda, dU/dL) sample has not been seen.
     */
    public boolean insideHardWallConstraint(double lambda, double dUdL) {
        if (hardWallConstraint) {
            double weight = evaluateHistogram(lambda, dUdL);
            return weight > 0.0;
        }
        return true;
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
    public double[] getScaling() {
        return potential.getScaling();
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
    public double[] getMass() {
        return potential.getMass();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Return a reference to each variables type.
     */
    @Override
    public Potential.VARIABLE_TYPE[] getVariableTypes() {
        return potential.getVariableTypes();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
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
    public void setVelocity(double[] velocity) {
        potential.setVelocity(velocity);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setAcceleration(double[] acceleration) {
        potential.setAcceleration(acceleration);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        potential.setPreviousAcceleration(previousAcceleration);
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
    public double[] getAcceleration(double[] acceleration) {
        return potential.getAcceleration(acceleration);
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
    public void setCrystal(Crystal crystal) {
        potential.setCrystal(crystal);
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
    public boolean destroy() {
        if (receiveThread != null && receiveThread.isAlive()) {
            double[] killMessage = new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};
            DoubleBuf killBuf = DoubleBuf.buffer(killMessage);
            try {
                logger.fine(" Sending the termination message.");
                world.send(rank, killBuf);
                logger.fine(" Termination message was sent successfully.");
                logger.fine(format(" Receive thread alive %b status %s", receiveThread.isAlive(), receiveThread.getState()));
            } catch (Exception ex) {
                String message = format(" Asynchronous Multiwalker OSRW termination signal " +
                        "failed to be sent for process %d.", rank);
                logger.log(Level.SEVERE, message, ex);
            }
        } else {
            logger.fine(" CountReceiveThread was either not initialized, or is not alive. This is the case for the Histogram script.");
        }
        return potential.destroy();
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
    public double getd2EdL2() {
        throw new UnsupportedOperationException(" Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        throw new UnsupportedOperationException(" Second derivatives of the bias are not implemented, as they require third derivatives of the potential.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean dEdLZeroAtEnds() {
        return false;
    }

}
