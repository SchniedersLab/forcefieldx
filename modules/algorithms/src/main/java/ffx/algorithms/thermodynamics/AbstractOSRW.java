/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms.thermodynamics;

import java.io.File;
import java.util.Random;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.Comm;
import edu.rit.pj.cluster.JobBackend;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.Barostat;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;

/**
 * An implementation of the Orthogonal Space Random Walk algorithm.
 *
 * @author Michael J. Schnieders, Wei Yang and Pengyu Ren
 */
public abstract class AbstractOSRW implements CrystalPotential {

    private static final Logger logger = Logger.getLogger(AbstractOSRW.class.getName());

    /**
     * The MolecularAssembly being simulated.
     */
    protected MolecularAssembly molecularAssembly;
    /**
     * A potential energy that implements the LambdaInterface.
     */
    final LambdaInterface lambdaInterface;
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
    final File lambdaFile;
    /**
     * Each walker reads the same histogram restart file. Only the walker of
     * rank 0 writes the histogram restart file.
     */
    final File histogramFile;
    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    protected double lambda;
    /**
     * Flag to indicate that the Lambda particle should be propagated.
     */
    boolean propagateLambda = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "propagateLambda" flag true.
     */
    int energyCount;
    int biasCount = 0;
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
     * Variable that dictates whether to use restart writing capabilites for Monte Carlo OSRW.
     */
    boolean mcRestart = false;
    /**
     * Parallel Java world communicator.
     */
    protected final Comm world;
    /**
     * Number of processes.
     */
    protected final int numProc;
    /**
     * Rank of this process.
     */
    protected final int rank;
    /**
     * A reference to the JobBackend to log messages for the Webserver.
     */
    final JobBackend jobBackend;
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
    double maxFLambda;
    /**
     * Force Field Potential Energy (i.e. with no bias terms added).
     */
    protected double forceFieldEnergy;
    /**
     * Partial derivative of the force field energy with respect to lambda.
     */
    double dForceFieldEnergydL;
    /**
     * OSRW Bias energy.
     */
    double biasEnergy;
    /**
     * Total partial derivative of the potential (U) being sampled
     * with respect to lambda.
     */
    double dUdLambda;
    /**
     * Second partial derivative of the potential being sampled with respect to lambda.
     */
    double d2UdL2;
    /**
     * Mixed second partial derivative with respect to coordinates and lambda.
     */
    double[] dUdXdL;
    /**
     * Magnitude of each hill (not including tempering).
     * <p>
     * The default biasMag = 0.05 (kcal/mol).
     */
    protected double biasMag;
    /**
     * 1D PMF with respect to lambda F(L).
     */
    double[] FLambda;
    /**
     * Magnitude of the 2D orthogonal space bias G(L,dE/dL).
     */
    double gLdEdL = 0.0;
    /**
     * First derivative of the 2D bias with respect to Lambda.
     */
    double dGdLambda;
    /**
     * First derivative of the 2D bias with respect to  dU/dL.
     */
    double dGdFLambda;
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
    int printFrequency;
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
    boolean osrwOptimization = false;
    /**
     * OSRW optimization only runs if Lambda is greater than the osrwOptimizationLambdaCutoff.
     * <p>
     * The default osrwOptimizationLambdaCutoff = 0.8.
     */
    double osrwOptimizationLambdaCutoff = 0.8;
    /**
     * Holds the lowest potential-energy parameters for loopBuilder runs from
     * all visits to lambda &gt; osrwOptimizationLambdaCutoff.
     */
    double[] osrwOptimumCoords;
    /**
     * The lowest energy found via optimizations.
     * <p>
     * The osrwOptimum is initially set to Double.MAX_VALUE.
     */
    double osrwOptimum = Double.MAX_VALUE;
    /**
     * The OSRW optimization frequency.
     */
    int osrwOptimizationFrequency = 10000;
    /**
     * The OSRW optimization convergence criteria.
     * <p>
     * The default osrwOptimizationEps = 0.1.
     */
    double osrwOptimizationEps = 0.1;
    /**
     * The OSRW optimization tolerance.
     */
    double osrwOptimizationTolerance = 1.0e-8;
    /**
     * The OSRW optimization energy window.
     */
    double osrwOptimizationEnergyWindow = 2.0;
    /**
     * File instance used for saving optimized structures.
     */
    File osrwOptimizationFile;
    /**
     * SystemFilter used to save optimized structures.
     */
    SystemFilter osrwOptimizationFilter;
    /**
     * Interval between how often the 1D histogram is printed to screen versus
     * silently updated in background.
     * <p>
     * The fLambdaPrintInterval is 25.
     */
    final int fLambdaPrintInterval = 25;
    /**
     * A count of FLambdaUpdates.
     */
    int fLambdaUpdates = 0;
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
    protected double totalEnergy;
    /**
     * Thermodynamic integration from Lambda=0 to Lambda=1.
     */
    double totalFreeEnergy;
    /**
     * Save the previous free energy, in order to limit logging to time points
     * where the free energy has changed.
     */
    double previousFreeEnergy = 0.0;
    /**
     * Equilibration counts.
     */
    int equilibrationCounts = 0;
    /**
     * Are FAST varying energy terms being computed, SLOW varying energy terms,
     * or BOTH. OSRW is not active when only FAST varying energy terms are being
     * propagated.
     */
    protected STATE state = STATE.BOTH;
    /**
     * Flag to indicate if OSRW should send and receive counts between processes
     * synchronously or asynchronously. The latter can be faster by ~40% because
     * simulation with Lambda &gt; 0.75 must compute two condensed phase
     * self-consistent fields to interpolate polarization.
     */
    protected final boolean asynchronous;
    /**
     * Write out structures only for lambda values greater than or equal to this threshold.
     */
    double lambdaWriteOut = 0.0;
    /**
     * Should the 1D OSRW bias be included in the target function.
     */
    boolean include1DBias;
    /**
     * Atom gradient array for use if "energy" is called.
     */
    private double[] grad = null;
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
     * OSRW Asynchronous MultiWalker Constructor.
     *
     * @param lambdaInterface    defines Lambda and dU/dL.
     * @param potential          defines the Potential energy.
     * @param lambdaFile         contains the current Lambda particle position and
     *                           velocity.
     * @param histogramFile      contains the Lambda and dU/dL histogram.
     * @param properties         defines System properties.
     * @param temperature        the simulation temperature.
     * @param dt                 the time step.
     * @param printInterval      number of steps between logging updates.
     * @param checkpointInterval number of steps between restart file updates.
     * @param asynchronous       set to true if walkers run asynchronously.
     * @param algorithmListener  the AlgorithmListener to be notified of
     *                           progress.
     */
    public AbstractOSRW(LambdaInterface lambdaInterface, CrystalPotential potential,
                        File lambdaFile, File histogramFile, CompositeConfiguration properties,
                        double temperature, double dt, double printInterval,
                        double checkpointInterval, boolean asynchronous,
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
        if (checkpointInterval >= this.dt) {
            saveFrequency = (int) (checkpointInterval / this.dt);
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
        FLambdaBins = 401;
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
        jobBackend = JobBackend.getJobBackend();

        include1DBias = properties.getBoolean("osrw-1D-bias", true);

        // Log OSRW parameters.
        logger.info("\n Orthogonal Space Random Walk Parameters");
        logger.info(format("  Gaussian Bias Magnitude:       %6.4f (kcal/mol)", biasMag));
        logger.info(format("  Gaussian Bias Cutoff:           %6d bins", biasCutoff));
        logger.info(format("  Print Interval:                 %6.3f psec", printInterval));
        logger.info(format("  Save Interval:                  %6.3f psec", checkpointInterval));
        logger.info(format("  Include 1D OSRW Bias:           %6b", include1DBias));
    }

    /**
     * <p>addBias.</p>
     *
     * @param dUdL     The current value of dU/dL
     * @param x        Atomic coordinate array.
     * @param gradient Gradient array.
     */
    public abstract void addBias(double dUdL, double[] x, double[] gradient);

    /**
     * Write out OSRW restart files.
     */
    public abstract void writeRestart();

    /**
     * Abstract method used to set the lambda threshold controlled write out of snapshots.
     *
     * @param lambdaWriteOut Lambda must be above this cut-off for writing out snapshots.
     */
    public abstract void setLambdaWriteOut(double lambdaWriteOut);

    /**
     * Compute the OSRW bias energy at (lambda, dU/dL).
     *
     * @param currentLambda The value of lambda.
     * @param currentdUdL   The value of dU/dL.
     * @return The value of the bias.
     */
    public abstract double computeBiasEnergy(double currentLambda, double currentdUdL);

    /**
     * <p>evaluateKernel.</p>
     *
     * @param cLambda   the current Lambda bin.
     * @param cF_Lambda the current dU/dL bin.
     * @return The value of the recursion kernel.
     */
    protected abstract double evaluateKernel(int cLambda, int cF_Lambda);

    /**
     * Update the 1D Bias based on the recursion kernel.
     *
     * @param print True requests verbose logging.
     * @param save  True requests saving the histogram.
     * @return The current free energy.
     */
    public abstract double updateFLambda(boolean print, boolean save);

    /**
     * {@inheritDoc}
     * <p>
     * Shuts down resources associated with this OSRW; primarily the receive
     * thread.
     */
    @Override
    public abstract boolean destroy();

    /**
     * <p>checkRecursionKernelSize.</p>
     *
     * @param dudl a double.
     */
    protected abstract void checkRecursionKernelSize(double dudl);

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
     * Specify if MC-OSRW restart files should be written.
     * TODO: move this method to MC-OSRW.
     *
     * @param mcRestart True to indicate writing MC OSRW restart files.
     */
    public void setMCRestartWriter(boolean mcRestart) {
        this.mcRestart = mcRestart;
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
    protected void langevin() {
        // Compute the random force pre-factor (kcal/mol * psec^-2).
        double rt2 = 2.0 * Thermostat.R * temperature * thetaFriction / dt;

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
    double current1DBiasEnergy(double currentLambda, boolean gradient) {
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
     * <p>lastFreeEnergy.</p>
     *
     * @return a double.
     */
    public double lastFreeEnergy() {
        return totalFreeEnergy;
    }

    /**
     * <p>evaluate2DPMF.</p>
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
     * <p>Setter for the field <code>resetStatistics</code>.</p>
     *
     * @param resetStatistics a boolean.
     */
    public void setResetStatistics(boolean resetStatistics) {
        this.resetStatistics = resetStatistics;
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
     * Set the OSRW Gaussian biasing potential magnitude (kcal/mol).
     *
     * @param biasMag Gaussian biasing potential magnitude (kcal/mol)
     */
    public void setBiasMagnitude(double biasMag) {
        this.biasMag = biasMag;
        logger.info(format(" Gaussian Bias Magnitude: %6.4f (kcal/mol)", biasMag));
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
    public VARIABLE_TYPE[] getVariableTypes() {
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
    public void setEnergyTermState(STATE state) {
        this.state = state;
        potential.setEnergyTermState(state);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {

        boolean origPropagateLambda = propagateLambda;

        propagateLambda = false;

        if (grad == null || grad.length != x.length) {
            grad = new double[x.length];
        }

        double energy = energyAndGradient(x, grad);

        propagateLambda = origPropagateLambda;

        return energy;
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

}
