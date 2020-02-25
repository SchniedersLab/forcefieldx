//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.cli;

import java.io.File;
import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.thermodynamics.MonteCarloOST;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.OptimizationParameters;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.WriteoutOptions;

import picocli.CommandLine;

/**
 * Represents command line options for scripts that utilize variants of the
 * Orthogonal Space Tempering (OST) algorithm. Metadynamics will be treated
 * as a special case of OST where there is no dU/dL axis.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class OSTOptions {

    private static final Logger logger = Logger.getLogger(OSTOptions.class.getName());

    /**
     * -c or --count Sets the number of time steps between OST counts.
     */
    @CommandLine.Option(names = {"-C", "--count"}, paramLabel = "10",
            description = "Time steps between MD Orthogonal Space counts.")
    private int countFreq = 10;

    /**
     * --bM or --biasMag sets the initial Gaussian bias magnitude in kcal/mol.
     */
    @CommandLine.Option(names = {"--bM", "--biasMag"}, paramLabel = "0.05", split = ",",
            description = "Orthogonal Space Gaussian bias magnitude (kcal/mol); repex OST uses a comma-separated list.")
    private double[] biasMag = new double[]{0.05};

    /**
     * --iW or --independentWalkers enforces that each walker maintains their own histogram.
     */
    @CommandLine.Option(names = {"--iW", "--independentWalkers"},
            description = "Enforces that each walker maintains their own histogram. ")
    private boolean independentWalkers = false;

    /**
     * --tp or --temperingParam sets the Dama et al tempering rate parameter, in
     * multiples of kBT.
     */
    @CommandLine.Option(names = {"--tp", "--temperingParam"}, paramLabel = "4.0", split=",",
            description = "Tempering rate parameter in multiples of kBT; repex OST uses a comma-separated list.")
    private double[] temperParam = new double[]{4.0};

    /**
     * --tth or --temperingThreshold sets the tempering threshold/offset in kcal/mol.
     */
    @CommandLine.Option(names = {"--tth", "--temperingThreshold"}, paramLabel = "20*bias", split=",",
            description = "Tempering threshold in kcal/mol; repex OST uses a comma-separated list.")
    private double[] temperThreshold = new double[]{-1};

    /**
     * --mc or --monteCarlo sets the Monte Carlo scheme for Orthogonal Space Tempering.
     */
    @CommandLine.Option(names = {"--mc", "--monteCarlo"},
            description = "Specify use of Monte Carlo OST")
    private boolean mc = false;

    /**
     * --mcHW or --monteCarloHardWall sets the Monte Carlo scheme to use a hard wall that rejects any sample
     * (Lambda, dU/dL) located in an empty histogram bin.
     */
    @CommandLine.Option(names = {"--mcHW", "--monteCarloHardWall"},
            description = "Monte Carlo OST hard wall constraint.")
    private boolean mcHW = false;

    /**
     * --mcmD or --mcTraj Sets the number of steps to take for each MD trajectory for MC-OST.
     */
    @CommandLine.Option(names = {"--mcMD", "--mcTraj"}, paramLabel = "100",
            description = "Number of dynamics steps to take for each MD trajectory for Monte Carlo OST")
    private int mcMD = 100;

    /**
     * --mcL or --mcLambdaStd Sets the standard deviation for lambda moves.
     */
    @CommandLine.Option(names = {"--mcL", "--mcLambdaStd"}, paramLabel = "0.01",
            description = "Standard deviation for lambda move.")
    private double mcL = 0.01;

    /**
     * --ts or --twoStep MC Orthogonal Space sampling using separate lambda and MD moves.
     */
    @CommandLine.Option(names = {"--ts", "--twoStep"},
            description = "MC Orthogonal Space sampling using separate lambda and MD moves.")
    private boolean ts = false;

    /**
     * --lw or --lambdaWritOut Only write out snapshots if lambda is greater than the value specified.
     */
    @CommandLine.Option(names = {"--lw", "--lambdaWritOut"}, paramLabel = "0.0",
            description = "Only write out snapshots if lambda is greater than the value specified.")
    private double lambdaWriteOut = 0.0;

    /**
     * <p>
     * Getter for the field <code>temperParam</code>.</p>
     *
     * @return a double.
     */
    public double getTemperParam() {
        return temperParam[0];
    }

    /**
     * Returns an array of all tempering parameters.
     *
     * @return All tempering parameters to be used.
     */
    public double[] getAllTemperingParameters() {
        return Arrays.copyOf(temperParam, temperParam.length);
    }

    /**
     * Returns an array of all tempering thresholds.
     * @return All tempering thresholds to be used.
     */
    public double[] getTemperThreshold() {
        return Arrays.copyOf(temperThreshold, temperThreshold.length);
    }

    /**
     * <p>
     * Getter for the field <code>lambdaWriteOut</code>.</p>
     *
     * @return a double.
     */
    public double getLambdaWriteOut() {
        return lambdaWriteOut;
    }

    /**
     * <p>
     * constructOST.</p>
     *
     * @param potential        a {@link ffx.crystal.CrystalPotential} object.
     * @param lambdaRestart    a {@link java.io.File} object.
     * @param histogramRestart a {@link java.io.File} object.
     * @param firstAssembly    a {@link ffx.potential.MolecularAssembly} object.
     * @param dynamics         a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param mdo              a {@link ffx.algorithms.cli.MultiDynamicsOptions} object.
     * @param thermo           a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param aListener        a {@link ffx.algorithms.AlgorithmListener} object.
     * @return a {@link OrthogonalSpaceTempering} object.
     */
    public OrthogonalSpaceTempering constructOST(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                 MolecularAssembly firstAssembly, DynamicsOptions dynamics,
                                                 MultiDynamicsOptions mdo, ThermodynamicsOptions thermo, AlgorithmListener aListener) {
        return constructOST(potential, lambdaRestart, histogramRestart,
                firstAssembly, null, dynamics, mdo, thermo, aListener);
    }

    /**
     * <p>
     * constructOST.</p>
     *
     * @param potential        a {@link ffx.crystal.CrystalPotential} object.
     * @param lambdaRestart    a {@link java.io.File} object.
     * @param histogramRestart a {@link java.io.File} object.
     * @param firstAssembly    a {@link ffx.potential.MolecularAssembly} object.
     * @param addedProperties  a {@link org.apache.commons.configuration2.Configuration} object.
     * @param dynamics         a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param mdo              a {@link ffx.algorithms.cli.MultiDynamicsOptions} object.
     * @param thermo           a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param aListener        a {@link ffx.algorithms.AlgorithmListener} object.
     * @return a {@link OrthogonalSpaceTempering} object.
     */
    public OrthogonalSpaceTempering constructOST(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                 MolecularAssembly firstAssembly, Configuration addedProperties,
                                                 DynamicsOptions dynamics, MultiDynamicsOptions mdo, ThermodynamicsOptions thermo,
                                                 AlgorithmListener aListener) {
        return constructOST(potential, lambdaRestart, histogramRestart, firstAssembly, addedProperties, dynamics, thermo, aListener, !mdo.isSynchronous(), 0);
    }

    /**
     * <p>
     * constructOST.</p>
     *
     * @param potential        a {@link ffx.crystal.CrystalPotential} object.
     * @param lambdaRestart    a {@link java.io.File} object.
     * @param histogramRestart a {@link java.io.File} object.
     * @param firstAssembly    a {@link ffx.potential.MolecularAssembly} object.
     * @param addedProperties  a {@link org.apache.commons.configuration2.Configuration} object.
     * @param dynamics         a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param thermo           a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param aListener        a {@link ffx.algorithms.AlgorithmListener} object.
     * @param async            If OST should use asynchronous communications.
     * @param replicateNum     Index of the replicate (0 for non-repex OST).
     * @return a {@link OrthogonalSpaceTempering} object.
     */
    public OrthogonalSpaceTempering constructOST(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                 MolecularAssembly firstAssembly, Configuration addedProperties,
                                                 DynamicsOptions dynamics, ThermodynamicsOptions thermo,
                                                 AlgorithmListener aListener, boolean async, int replicateNum) {

        LambdaInterface linter = (LambdaInterface) potential;
        CompositeConfiguration allProperties = new CompositeConfiguration(firstAssembly.getProperties());
        if (addedProperties != null) {
            allProperties.addConfiguration(addedProperties);
        }
        double temp = dynamics.getTemp();
        double dT = dynamics.getDt();
        double report = dynamics.getReport();
        double ckpt = dynamics.getCheckpoint();
        boolean resetNSteps = thermo.getResetNumSteps();
        OrthogonalSpaceTempering orthogonalSpaceTempering = new OrthogonalSpaceTempering(linter, potential, lambdaRestart,
                histogramRestart, allProperties, temp, dT, report, ckpt, async, resetNSteps, aListener);
        orthogonalSpaceTempering.setHardWallConstraint(mcHW);

        // Do NOT run applyOSTOptions here, because that can mutate the OST to a Barostat.
        return orthogonalSpaceTempering;
    }

    /**
     * Applies relevant options to a OST, and returns either
     * the OST object or something that wraps the OST (such as a
     * Barostat).
     *
     * @param orthogonalSpaceTempering Orthogonal Space Tempering.
     * @param firstAssembly            Primary assembly in OST.
     * @param dynamics                 MD options.
     * @param lpo                      Lambda particle options.
     * @param barostat                 NPT options.
     * @param histogramExists          If the histogram file exists already.
     * @return a {@link ffx.crystal.CrystalPotential} object.
     */
    public CrystalPotential applyAllOSTOptions(OrthogonalSpaceTempering orthogonalSpaceTempering,
                                               MolecularAssembly firstAssembly,
                                               DynamicsOptions dynamics, LambdaParticleOptions lpo,
                                               BarostatOptions barostat, boolean histogramExists) {

        applyHistogramOptions(orthogonalSpaceTempering, histogramExists);
        if (!histogramExists) {
            orthogonalSpaceTempering.setThetaFrication(lpo.getLambdaFriction());
            orthogonalSpaceTempering.setThetaMass(lpo.getLambdaMass());
        }
        if (dynamics.getOptimize()) {
            OptimizationParameters opt = orthogonalSpaceTempering.getOptimizationParameters();
            opt.setOptimization(true, firstAssembly);
        }
        return barostat.checkNPT(firstAssembly, orthogonalSpaceTempering);
    }

    /**
     * Assembles a MolecularDynamics wrapped around a Potential.
     *
     * @param topologies MolecularAssembly[]
     * @param potential  Potential to run on
     * @param dynamics   DynamicsOptions
     * @param aListener  AlgorithmListener
     * @return           MolecularDynamics
     */
    public MolecularDynamics assembleMolecularDynamics(MolecularAssembly[] topologies, CrystalPotential potential,
                                                       DynamicsOptions dynamics, AlgorithmListener aListener) {
        // Create the MolecularDynamics instance.
        MolecularAssembly firstTop = topologies[0];
        CompositeConfiguration props = firstTop.getProperties();

        dynamics.init();

        MolecularDynamics molDyn = MolecularDynamics.dynamicsFactory(firstTop, potential, props,
                aListener, dynamics.thermostat, dynamics.integrator, MolecularDynamics.DynamicsEngine.FFX);
        for (int i = 1; i < topologies.length; i++) {
            molDyn.addAssembly(topologies[i], topologies[i].getProperties());
        }
        molDyn.setRestartFrequency(dynamics.getCheckpoint());

        return molDyn;
    }

    /**
     * Begins MD-OST sampling from an assembled OST.
     *
     * @param orthogonalSpaceTempering The OST object.
     * @param topologies               All MolecularAssemblys.
     * @param potential                The top-layer CrystalPotential.
     * @param dynamics                 Dynamics options.
     * @param writeOut                 a {@link WriteoutOptions} object.
     * @param thermo                   Thermodynamics options.
     * @param dyn                      The .dyn dynamics restart file.
     * @param aListener                AlgorithmListener
     */
    public void beginMDOST(OrthogonalSpaceTempering orthogonalSpaceTempering,
                           MolecularAssembly[] topologies, CrystalPotential potential,
                           DynamicsOptions dynamics, WriteoutOptions writeOut, ThermodynamicsOptions thermo,
                           File dyn, AlgorithmListener aListener) {
        // Create the MolecularDynamics instance.
        MolecularAssembly firstTop = topologies[0];

        dynamics.init();

        MolecularDynamics molDyn = assembleMolecularDynamics(topologies, potential, dynamics, aListener);

        boolean initVelocities = true;
        long nSteps = dynamics.steps;
        // Start sampling.
        long nEquil = thermo.getEquilSteps();
        if (nEquil > 0) {
            logger.info("\n Beginning equilibration");
            orthogonalSpaceTempering.setPropagateLambda(false);
            runDynamics(molDyn, nEquil, dynamics, writeOut, true, dyn);
            logger.info(" Beginning OST sampling");
            orthogonalSpaceTempering.setPropagateLambda(true);
        } else {
            logger.info(" Beginning OST sampling without equilibration");
            if (!thermo.getResetNumSteps()) {
                long nEnergyCount = orthogonalSpaceTempering.getEnergyCount();
                if (nEnergyCount > 0) {
                    nSteps -= nEnergyCount;
                    logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
                    initVelocities = false;
                }
            }
        }
        if (nSteps > 0) {
            runDynamics(molDyn, nSteps, dynamics, writeOut, initVelocities, dyn);
        } else {
            logger.info(" No steps remaining for this process!");
        }
    }

    /**
     * <p>
     * setupMCOST.</p>
     *
     * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
     * @param topologies               an array of {@link ffx.potential.MolecularAssembly} objects.
     * @param dynamics                 a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param thermodynamics           a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param verbose                  Whether to print out additional information about MC-OST.
     * @param listener                 An AlgorithmListener
     * @return                         An assembled MonteCarloOST ready to run.
     */
    public MonteCarloOST setupMCOST(OrthogonalSpaceTempering orthogonalSpaceTempering, MolecularAssembly[] topologies,
                           DynamicsOptions dynamics, ThermodynamicsOptions thermodynamics, boolean verbose,
                           AlgorithmListener listener) {
        dynamics.init();

        MonteCarloOST monteCarloOST = new MonteCarloOST(orthogonalSpaceTempering.getPotentialEnergy(),
                orthogonalSpaceTempering, topologies[0], topologies[0].getProperties(), listener, dynamics, verbose, mcMD);

        long nEquil = thermodynamics.getEquilSteps();
        if (nEquil > 0) {
            monteCarloOST.setEquilibration(true);
        }
        return monteCarloOST;
    }

    /**G
     * Runs MC-OST.
     *
     * @param monteCarloOST MC-OST to run.
     */
    public void beginMCOST(MonteCarloOST monteCarloOST, DynamicsOptions dynamics, ThermodynamicsOptions thermo) {
        long nEquil = thermo.getEquilSteps();

        if (nEquil > 0) {
            logger.info("\n Beginning MC-OST equilibration.");
            monteCarloOST.setMDMoveParameters(nEquil);
            if (ts) {
                monteCarloOST.sampleTwoStep();
            } else {
                monteCarloOST.sampleOneStep();
            }
            monteCarloOST.setEquilibration(false);
            logger.info("\n Finished MC-OST equilibration.");
        }

        logger.info("\n Beginning MC-OST sampling.");
        monteCarloOST.setLambdaStdDev(mcL);
        monteCarloOST.setMDMoveParameters(dynamics.steps);
        if (lambdaWriteOut >= 0.0 && lambdaWriteOut <= 1.0) {
            monteCarloOST.setLambdaWriteOut(lambdaWriteOut);
        }

        if (ts) {
            monteCarloOST.sampleTwoStep();
        } else {
            monteCarloOST.sampleOneStep();
        }
    }

    /**
     * <p>
     * applyHistogramOptions.</p>
     *
     * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
     * @param histogramExists          a boolean.
     */
    public void applyHistogramOptions(OrthogonalSpaceTempering orthogonalSpaceTempering, boolean histogramExists) {
        applyHistogramOptions(orthogonalSpaceTempering, histogramExists, 0);
    }

    /**
     * Applies histogram-related options to a histogram.
     * @param orthogonalSpaceTempering OST containing the histogram.
     * @param histogramExists          Whether the .his restart file exists.
     * @param index                    Index of the histogram.
     */
    public void applyHistogramOptions(OrthogonalSpaceTempering orthogonalSpaceTempering, boolean histogramExists, int index) {
        assert index >= 0;
        Histogram histogram = orthogonalSpaceTempering.getHistogram(index);
        if (!histogramExists) {
            orthogonalSpaceTempering.setCountInterval(countFreq);
            histogram.setBiasMagnitude(getBiasMag(index));
        }
        histogram.setIndependentWalkers(independentWalkers);
        histogram.setTemperingThreshold(getTemperingThreshold(index));
        histogram.setTemperingParameter(getTemperingParameter(index));
        histogram.checkRecursionKernelSize();
    }

    private void runDynamics(MolecularDynamics molDyn, long numSteps, DynamicsOptions dynamics,
                             WriteoutOptions writeout, boolean initVelocities, File dyn) {
        molDyn.dynamic(numSteps, dynamics.dt, dynamics.report, dynamics.write, dynamics.temp,
                initVelocities, writeout.getFileType(), dynamics.getCheckpoint(), dyn);
    }

    /**
     * Checks if independent walkers has been specified.
     *
     * @return Walker independence.
     */
    public boolean getIndependentWalkers() {
        return independentWalkers;
    }

    /**
     * Checks if use of the Monte Carlo algorithm has been specified.
     *
     * @return Monte Carlo OST (as opposed to molecular dynamics OST).
     */
    public boolean isMc() {
        return mc;
    }

    /**
     * Returns true if the 2-step option is enabled (not guaranteed to
     * also mean that MC is enabled!).
     *
     * @return If --ts is enabled.
     */
    public boolean isTwoStep() {
        return ts;
    }

    /**
     * Returns the initial bias magnitude associated with a walker.
     *
     * @param i Index of a walker
     * @return  Its intended initial bias magnitude in kcal/mol.
     */
    public double getBiasMag(int i) {
        return biasMag.length > 1 ? biasMag[i] : biasMag[0];
    }

    /**
     * Returns the tempering threshold associated with a walker.
     *
     * @param i Index of a walker
     * @return  Its intended tempering threshold in kcal/mol.
     */
    public double getTemperingThreshold(int i) {
        return temperThreshold.length > 1 ? temperThreshold[i] : temperThreshold[0];
    }

    /**
     * Returns the tempering parameter associated with a walker.
     *
     * @param i Index of a walker
     * @return  Its intended tempering parameter in kBT.
     */
    public double getTemperingParameter(int i) {
        return temperParam.length > 1 ? temperParam[i] : temperParam[0];
    }
}
