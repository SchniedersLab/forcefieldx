/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.algorithms.cli;

import java.io.File;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.TransitionTemperedOSRW;
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
     * --tp or --temperingParam sets the Dama et al tempering rate parameter,
     * in multiples of kBT.
     */
    @CommandLine.Option(names = {"--tp", "--temperingParam"}, paramLabel = "8.0", description = "Dama et al tempering rate parameter in multiples of kBT")
    private double temperParam = 8.0;

    public double getTemperParam() {
        return temperParam;
    }

    public TransitionTemperedOSRW constructOSRW(CrystalPotential potential, File lambdaRestart, File histogramRestart,
                                                MolecularAssembly firstAssembly, DynamicsOptions dynamics,
                                                MultiDynamicsOptions mdo, ThermodynamicsOptions thermo, AlgorithmListener aListener) {

        LambdaInterface linter = (LambdaInterface) potential;
        CompositeConfiguration firstProps = firstAssembly.getProperties();
        double temp = dynamics.getTemp();
        double dT = dynamics.getDt();
        double report = dynamics.getReport();
        double ckpt = dynamics.getCheckpoint();
        boolean async = !mdo.isSynchronous();
        boolean resetNSteps = thermo.getResetNumSteps();
        TransitionTemperedOSRW ttOSRW = new TransitionTemperedOSRW(linter, potential, lambdaRestart,
                histogramRestart, firstProps, temp, dT, report, ckpt, async, resetNSteps, aListener);

        // Do NOT run applyOSRWOptions here, because that can mutate the TT-OSRW to a Barostat.
        return ttOSRW;
    }

    /**
     * Applies relevant options to a TransitionTemperedOSRW, and returns either the TTOSRW
     * object or something that wraps the TTOSRW (such as a Barostat).
     *
     * @param ttOSRW Transition-Tempered Orthogonal Space Random Walk.
     * @param firstAssembly Primary assembly in ttOSRW.
     * @param dynamics MD options.
     * @param lpo Lambda particle options.
     * @param alch Alchemy options.
     * @param barostat NPT options.
     * @param lamExists If the lambda file exists for this walker.
     * @param histogramExists If the histogram file exists already.
     * @return TTOSRW, possibly wrapped.
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

    public void applyOSRWOptions(TransitionTemperedOSRW ttOSRW, boolean histogramExists) {
        ttOSRW.setTemperingParameter(temperParam);
        if (!histogramExists) {
            ttOSRW.setCountInterval(countFreq);
            ttOSRW.setBiasMagnitude(biasMag);
        }
    }

    /**
     * Begins MD-OSRW sampling from an assembled TT-OSRW.
     * @param ttOSRW The TT-OSRW object.
     * @param topologies All MolecularAssemblys.
     * @param potential The top-layer CrystalPotential.
     * @param dynamics Dynamics options.
     * @param thermo Thermodynamics options.
     * @param dyn The .dyn dynamics restart file.
     * @param aListener AlgorithmListener
     */
    public void beginMDOSRW(TransitionTemperedOSRW ttOSRW, MolecularAssembly[] topologies, CrystalPotential potential,
                            DynamicsOptions dynamics, WriteoutOptions writeout, ThermodynamicsOptions thermo,
                            File dyn, AlgorithmListener aListener) {
        // Create the MolecularDynamics instance.
        MolecularAssembly firstTop = topologies[0];
        CompositeConfiguration props = firstTop.getProperties();
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

    private void runDynamics(MolecularDynamics molDyn, int numSteps, DynamicsOptions dynamics, WriteoutOptions writeout, boolean initVelocities, File dyn) {
        molDyn.dynamic(numSteps, dynamics.dt, dynamics.report, dynamics.write, dynamics.temp, initVelocities, writeout.getFileType(), dynamics.getCheckpoint(), dyn);
    }
}
