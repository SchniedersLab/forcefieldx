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
package ffx.algorithms.cli;

import java.io.File;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.thermodynamics.MonteCarloOSRW;
import ffx.algorithms.thermodynamics.TransitionTemperedOSRW;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.AlchemicalOptions;

import picocli.CommandLine;

/**
 * Represents command line options for scripts that utilize variants of the
 * Orthogonal Space Random Walk (OSRW) algorithm. Metadynamics will be treated
 * as a special case of OSRW where there is no dU/dL axis.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class OSRWOptions {

    private static final Logger logger = Logger.getLogger(OSRWOptions.class.getName());

    /**
     * -c or --count sets the number of time steps between OSRW counts.
     */
    @CommandLine.Option(names = {"-C", "--count"}, paramLabel = "10", description = "Time steps between MD-OSRW counts.")
    private int countFreq = 10;

    /**
     * --bM or --biasMag sets the initial Gaussian bias magnitude in kcal/mol.
     */
    @CommandLine.Option(names = {"--bM", "--biasMag"}, paramLabel = "0.05", description = "OSRW Gaussian bias magnitude (kcal/mol).")
    private double biasMag = 0.05;

    /**
     * --tp or --temperingParam sets the Dama et al tempering rate parameter, in
     * multiples of kBT.
     */
    @CommandLine.Option(names = {"--tp", "--temperingParam"}, paramLabel = "8.0", description = "Dama et al tempering rate parameter in multiples of kBT")
    private double temperParam = 8.0;

    /**
     * --mc or --monteCarlo sets the Monte Carlo scheme for OSRW.
     */
    @CommandLine.Option(names = {"--mc", "monteCarlo"}, description = "Monte Carlo OSRW")
    private boolean mc = false;

    /**
     * --mcmD or --mcTraj sets the number of steps to take for each MD
     * trajectory for MCOSRW.
     */
    @CommandLine.Option(names = {"--mcMD", "--mcTraj"}, paramLabel = "100", description = "Number of dynamics steps to take for each MD trajectory for Monte Carlo OSRW")
    private int mcMD = 100;

    /**
     * --mcL or --mcLambdaStd sets the standard deviation for lambda.
     */
    @CommandLine.Option(names = {"--mcL", "--mcLambdaStd"}, paramLabel = "0.1", description = "Standard deviation for lambda move.")
    private double mcL = 0.1;

    /**
     * --ts or --twoStep Sample MC-OSRW using separate lambda and MD moves.
     */
    @CommandLine.Option(names = {"--ts", "--twoStep"}, description = "Sample MC-OSRW using separate lambda and MD moves.")
    private boolean ts = false;

    /**
     * --lw or --lambdaWritOut Only write out snapshots if lambda is greater than the value specified.
     */
    @CommandLine.Option(names = {"--lw", "--lambdaWritOut"}, paramLabel = "1.0",
            description = "Only write out snapshots if lambda is greater than the value specified.")
    private double lambdaWriteOut = 1.0;

    /**
     * <p>
     * Getter for the field <code>temperParam</code>.</p>
     *
     * @return a double.
     */
    public double getTemperParam() {
        return temperParam;
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
     * constructOSRW.</p>
     *
     * @param potential        a {@link ffx.crystal.CrystalPotential} object.
     * @param lambdaRestart    a {@link java.io.File} object.
     * @param histogramRestart a {@link java.io.File} object.
     * @param firstAssembly    a {@link ffx.potential.MolecularAssembly} object.
     * @param dynamics         a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param mdo              a {@link ffx.algorithms.cli.MultiDynamicsOptions} object.
     * @param thermo           a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param aListener        a {@link ffx.algorithms.AlgorithmListener} object.
     * @return a {@link TransitionTemperedOSRW} object.
     */
    public TransitionTemperedOSRW constructOSRW(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                MolecularAssembly firstAssembly, DynamicsOptions dynamics,
                                                MultiDynamicsOptions mdo, ThermodynamicsOptions thermo, AlgorithmListener aListener) {
        return constructOSRW(potential, lambdaRestart, histogramRestart, firstAssembly, null, dynamics, mdo, thermo, aListener);
    }

    /**
     * <p>
     * constructOSRW.</p>
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
     * @return a {@link TransitionTemperedOSRW} object.
     */
    public TransitionTemperedOSRW constructOSRW(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                MolecularAssembly firstAssembly, Configuration addedProperties,
                                                DynamicsOptions dynamics, MultiDynamicsOptions mdo, ThermodynamicsOptions thermo,
                                                AlgorithmListener aListener) {

        LambdaInterface linter = (LambdaInterface) potential;
        CompositeConfiguration allProperties = new CompositeConfiguration(firstAssembly.getProperties());
        if (addedProperties != null) {
            allProperties.addConfiguration(addedProperties);
        }
        double temp = dynamics.getTemp();
        double dT = dynamics.getDt();
        double report = dynamics.getReport();
        double ckpt = dynamics.getCheckpoint();
        boolean async = !mdo.isSynchronous();
        boolean resetNSteps = thermo.getResetNumSteps();
        TransitionTemperedOSRW ttOSRW = new TransitionTemperedOSRW(linter, potential, lambdaRestart,
                histogramRestart, allProperties, temp, dT, report, ckpt, async, resetNSteps, aListener);
        ttOSRW.checkRecursionKernelSize();

        // Do NOT run applyOSRWOptions here, because that can mutate the TT-OSRW to a Barostat.
        return ttOSRW;
    }

    /**
     * Applies relevant options to a TransitionTemperedOSRW, and returns either
     * the TTOSRW object or something that wraps the TTOSRW (such as a
     * Barostat).
     *
     * @param ttOSRW          Transition-Tempered Orthogonal Space Random Walk.
     * @param firstAssembly   Primary assembly in ttOSRW.
     * @param dynamics        MD options.
     * @param lpo             Lambda particle options.
     * @param alch            Alchemy options.
     * @param barostat        NPT options.
     * @param lamExists       If the lambda file exists for this walker.
     * @param histogramExists If the histogram file exists already.
     * @return a {@link ffx.crystal.CrystalPotential} object.
     */
    public CrystalPotential applyAllOSRWOptions(TransitionTemperedOSRW ttOSRW, MolecularAssembly firstAssembly,
                                                DynamicsOptions dynamics, LambdaParticleOptions lpo, AlchemicalOptions alch,
                                                BarostatOptions barostat, boolean lamExists, boolean histogramExists) {

        CrystalPotential cpot = ttOSRW;
        applyOSRWOptions(ttOSRW, histogramExists);
        if (histogramExists) {
            ttOSRW.setThetaFrication(lpo.getLambdaFriction());
            ttOSRW.setThetaMass(lpo.getLambdaMass());
        }

        if (dynamics.getOptimize()) {
            ttOSRW.setOptimization(true, firstAssembly);
            // TODO: Apply other minimization parameters.
        }

        if (!lamExists) {
            double lam = alch.getInitialLambda();
            logger.info(String.format(" Setting lambda to %5.3f", lam));
            ttOSRW.setLambda(lam);
        }
        cpot = barostat.checkNPT(firstAssembly, cpot);

        return cpot;
    }

    /**
     * <p>
     * applyOSRWOptions.</p>
     *
     * @param ttOSRW          a {@link TransitionTemperedOSRW} object.
     * @param histogramExists a boolean.
     */
    public void applyOSRWOptions(TransitionTemperedOSRW ttOSRW, boolean histogramExists) {
        ttOSRW.setTemperingParameter(temperParam);
        if (!histogramExists) {
            ttOSRW.setCountInterval(countFreq);
            ttOSRW.setBiasMagnitude(biasMag);
        }
    }

    /**
     * Begins MD-OSRW sampling from an assembled TT-OSRW.
     *
     * @param ttOSRW     The TT-OSRW object.
     * @param topologies All MolecularAssemblys.
     * @param potential  The top-layer CrystalPotential.
     * @param dynamics   Dynamics options.
     * @param writeout   a {@link ffx.algorithms.cli.WriteoutOptions} object.
     * @param thermo     Thermodynamics options.
     * @param dyn        The .dyn dynamics restart file.
     * @param aListener  AlgorithmListener
     */
    public void beginMDOSRW(TransitionTemperedOSRW ttOSRW, MolecularAssembly[] topologies, CrystalPotential potential,
                            DynamicsOptions dynamics, WriteoutOptions writeout, ThermodynamicsOptions thermo,
                            File dyn, AlgorithmListener aListener) {
        // Create the MolecularDynamics instance.
        MolecularAssembly firstTop = topologies[0];
        CompositeConfiguration props = firstTop.getProperties();

        dynamics.init();

        MolecularDynamics molDyn = MolecularDynamics.dynamicsFactory(firstTop, potential, props,
                aListener, dynamics.thermostat, dynamics.integrator, MolecularDynamics.DynamicsEngine.FFX);
        for (int i = 1; i < topologies.length; i++) {
            molDyn.addAssembly(topologies[i], topologies[i].getProperties());
        }

        boolean initVelocities = true;
        int nSteps = dynamics.steps;
        molDyn.setRestartFrequency(dynamics.getCheckpoint());
        // Start sampling.
        int nEquil = thermo.getEquilSteps();
        if (nEquil > 0) {
            logger.info(" Beginning equilibration");
            ttOSRW.setPropagateLambda(false);
            runDynamics(molDyn, nEquil, dynamics, writeout, true, dyn);
            logger.info(" Beginning Transition-Tempered OSRW sampling");
            ttOSRW.setPropagateLambda(true);
        } else {
            logger.info(" Beginning Transition-Tempered OSRW sampling without equilibration");
            if (!thermo.getResetNumSteps()) {
                int nEnergyCount = ttOSRW.getEnergyCount();
                if (nEnergyCount > 0) {
                    nSteps -= nEnergyCount;
                    logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
                    initVelocities = false;
                }
            }
        }
        if (nSteps > 0) {
            runDynamics(molDyn, nSteps, dynamics, writeout, initVelocities, dyn);
        } else {
            logger.info(" No steps remaining for this process!");
        }
    }

    /**
     * <p>
     * beginMCOSRW.</p>
     *
     * @param ttOSRW         a {@link TransitionTemperedOSRW} object.
     * @param topologies     an array of {@link ffx.potential.MolecularAssembly} objects.
     * @param potential      a {@link ffx.crystal.CrystalPotential} object.
     * @param dynamics       a {@link ffx.algorithms.cli.DynamicsOptions} object.
     * @param writeout       a {@link ffx.algorithms.cli.WriteoutOptions} object.
     * @param thermodynamics a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
     * @param dyn            a {@link java.io.File} object.
     * @param aListener      a {@link ffx.algorithms.AlgorithmListener} object.
     */
    public void beginMCOSRW(TransitionTemperedOSRW ttOSRW, MolecularAssembly[] topologies, CrystalPotential potential,
                            DynamicsOptions dynamics, WriteoutOptions writeout, ThermodynamicsOptions thermodynamics,
                            File dyn, AlgorithmListener aListener) {

        dynamics.init();

        MonteCarloOSRW mcOSRW = new MonteCarloOSRW(ttOSRW.getPotentialEnergy(), ttOSRW, topologies[0],
                topologies[0].getProperties(), null, ThermostatEnum.ADIABATIC, dynamics.integrator);

        int nEquil = thermodynamics.getEquilSteps();
        if (nEquil > 0) {
            logger.info("\n Beginning MC Transition-Tempered OSRW equilibration");
            mcOSRW.setEquilibration(true);
            mcOSRW.setMDMoveParameters(nEquil, mcMD, dynamics.dt);
            mcOSRW.sampleTwoStep();
            mcOSRW.setEquilibration(false);
            logger.info("\n Finished MC Transition-Tempered OSRW equilibration");
        }

        logger.info("\n Beginning MC Transition-Tempered OSRW sampling");
        mcOSRW.setLambdaStdDev(mcL);
        //lambdaWriteOut = getLambdaWriteOut();
        mcOSRW.setMDMoveParameters(dynamics.steps, mcMD, dynamics.dt);
        if (lambdaWriteOut >= 0.0 && lambdaWriteOut <= 1.0) {
            mcOSRW.setLambdaWriteOut(lambdaWriteOut);
        }
        if (ts) {
            mcOSRW.sampleTwoStep();
        } else {
            mcOSRW.sampleOneStep();
        }
    }

    private void runDynamics(MolecularDynamics molDyn, int numSteps, DynamicsOptions dynamics, WriteoutOptions writeout, boolean initVelocities, File dyn) {
        molDyn.dynamic(numSteps, dynamics.dt, dynamics.report, dynamics.write, dynamics.temp, initVelocities, writeout.getFileType(), dynamics.getCheckpoint(), dyn);
    }
}
