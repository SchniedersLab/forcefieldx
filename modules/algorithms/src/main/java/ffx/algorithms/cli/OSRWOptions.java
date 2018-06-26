/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.RotamerOptimization;
import ffx.algorithms.TransitionTemperedOSRW;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.cli.AlchemicalOptions;
import org.apache.commons.configuration.CompositeConfiguration;
import picocli.CommandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
     * -b or --bias sets the initial Gaussian bias magnitude in kcal/mol.
     */
    @CommandLine.Option(names = {"-B", "--biasMag"}, paramLabel = "0.05", description = "OSRW Gaussian bias magnitude (kcal/mol).")
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
                                                LambdaParticleOptions lpo, MinimizeOptions min, MultiDynamicsOptions mdo,
                                                ThermodynamicsOptions thermo, AlchemicalOptions alch, AlgorithmListener aListener) {

        LambdaInterface linter = (LambdaInterface) potential;
        CompositeConfiguration firstProps = firstAssembly.getProperties();
        double temp = dynamics.getTemp();
        double dT = dynamics.getDt();
        double report = dynamics.getReport();
        double ckpt = dynamics.getCheckpoint();
        boolean async = ! mdo.isSynchronous();
        boolean resetNSteps = thermo.getResetNumSteps();
        TransitionTemperedOSRW ttOSRW = new TransitionTemperedOSRW(linter, potential, lambdaRestart,
                histogramRestart, firstProps, temp, dT, report, ckpt, async, resetNSteps, aListener);

        // Do NOT run applyOptionsToTTOSRW here, because that can mutate the TT-OSRW to a Barostat.
        return ttOSRW;
    }

    public CrystalPotential applyOptionsToTTOSRW(TransitionTemperedOSRW ttOSRW, MolecularAssembly firstAssembly, DynamicsOptions dynamics,
                                     LambdaParticleOptions lpo, MinimizeOptions min, MultiDynamicsOptions mdo,
                                     ThermodynamicsOptions thermo, AlchemicalOptions alch, boolean lamExists, boolean hisExists) {

        CrystalPotential cpot = ttOSRW;
        applyOptionsToTTOSRW(ttOSRW, hisExists);
        if (hisExists) {
            ttOSRW.setThetaFrication(lpo.getLambdaFriction());
            ttOSRW.setThetaMass(lpo.getLambdaMass());
        }

        if (dynamics.getOptimize()) {
            ttOSRW.setOptimization(true, firstAssembly);
            // TODO: Apply other minimization parameters.
        }

        if (lamExists) {
            ttOSRW.setLambda(alch.getInitialLambda());
        }

        // TODO: cpot = PressureOptions.applyPressureOptions(cpot);
        return cpot;
        /*
        if (options.mc) {
            MonteCarloOSRW mcOSRW = new MonteCarloOSRW(osrw.getPotentialEnergy(), osrw, topologies[0],
                topologies[0].getProperties(), null, ThermostatEnum.ADIABATIC, options.integrator);

            if (options.nEquil > 0) {
                logger.info("\n Beginning MC Transition-Tempered OSRW equilibration");
                mcOSRW.setEquilibration(true)
                mcOSRW.setMDMoveParameters(options.nEquil, options.mcMD, options.dt)
                mcOSRW.sample()
                mcOSRW.setEquilibration(false)
                logger.info("\n Finished MC Transition-Tempered OSRW equilibration");
            }

            logger.info("\n Beginning MC Transition-Tempered OSRW sampling");
            mcOSRW.setLambdaStdDev(options.mcL)
            mcOSRW.setMDMoveParameters(options.steps, options.mcMD, options.dt)
            mcOSRW.sample()
        } else {
            // Create the MolecularDynamics instance.
            // If we switch over to using the factory method, request the FFX Dynamics engine.
            MolecularDynamics molDyn = new MolecularDynamics(topologies[0], potential,
                topologies[0].getProperties(), null, options.tstat, options.integrator);
            for (int i = 1; i < topologies.size(); i++) {
                molDyn.addAssembly(topologies.get(i), properties.get(i));
            }

            boolean initVelocities = true;
            double restartInterval = options.write;
            String fileType = "XYZ";
            int nSteps = options.steps;
            molDyn.setRestartFrequency(options.checkpoint);
            // Start sampling.
            if (options.nEquil > 0) {
                logger.info(" Beginning equilibration");
                osrw.setPropagateLambda(false);
                molDyn.dynamic(options.nEquil, options.dt, options.report, options.write, options.temp, initVelocities, dyn);
                logger.info(" Beginning Transition-Tempered OSRW sampling");
                osrw.setPropagateLambda(true);
                molDyn.dynamic(nSteps, options.dt, options.report, options.write, options.temp, false,
                    fileType, restartInterval, dyn);
            } else {
                logger.info(" Beginning Transition-Tempered OSRW sampling without equilibration");
                boolean resetSteps = true;
                if (options.resetStepsString) {
                    if (options.resetStepsString.equalsIgnoreCase("false")) {
                        resetSteps = false;
                    }
                }
                if (!resetSteps) {
                    int nEnergyCount = osrw.getEnergyCount();
                    if (nEnergyCount > 0) {
                        nSteps -= nEnergyCount;
                        logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
                        initVelocities = false;
                    }
                }
                if (nSteps > 0) {
                    molDyn.dynamic(nSteps, options.dt, options.report, options.write, options.temp, initVelocities,
                        fileType, restartInterval, dyn);
                } else {
                    logger.info(" No steps remaining for this process!");
                }
            }
        }
         */

    }

    public void applyOptionsToTTOSRW(TransitionTemperedOSRW ttOSRW, boolean histogramExists) {
        ttOSRW.setTemperingParameter(temperParam);
        if (!histogramExists) {
            ttOSRW.setCountInterval(countFreq);
            ttOSRW.setBiasMagnitude(biasMag);
        }
    }

    /**
     * Distribute side-chain conformations of mola.
     *
     * @param mola To distribute
     * @param pot Potential to use
     */
    private void optStructure(MolecularAssembly mola, Potential pot, String[] distribRes, AlgorithmFunctions aFuncts, int rank, int worldSize) {
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        if (distribRes == null || distribRes.length == 0) {
            throw new IllegalArgumentException(" Programming error: Must have list of residues to split on!");
        }

        LambdaInterface linter = (pot instanceof LambdaInterface) ? (LambdaInterface) pot : null;
        double initLam = -1.0;
        if (linter != null) {
            initLam = linter.getLambda();
            linter.setLambda(0.5);
        }

        Pattern chainMatcher = Pattern.compile("^([a-zA-Z])?([0-9]+)$");

        List<Residue> residueList = new ArrayList<>(distribRes.length);

        for (String ts : distribRes) {
            Matcher m = chainMatcher.matcher(ts);
            Character chainID;
            int resNum;
            if (m.find()) {
                if (m.groupCount() == 2) {
                    chainID = m.group(1).charAt(0);
                    resNum = Integer.parseInt(m.group(2));
                } else {
                    chainID = ' ';
                    resNum = Integer.parseInt(m.group(1));
                }
            } else {
                logger.warning(String.format(" Could not parse %s as a valid residue!", ts));
                continue;
            }
            logger.info(String.format(" Looking for chain %c residue %d", chainID, resNum));

            for (Polymer p : mola.getChains()) {
                if (p.getChainID() == chainID) {
                    for (Residue r : p.getResidues()) {
                        if (r.getResidueNumber() == resNum && r.getRotamers(rLib) != null) {
                            residueList.add(r);
                        }
                    }
                }
            }
        }

        if (residueList.isEmpty()) {
            throw new IllegalArgumentException(" No valid entries for distWalkers!");
        }

        AlgorithmListener alist = aFuncts.getDefaultListener();
        RotamerOptimization ropt = new RotamerOptimization(mola, pot, alist);

        ropt.setThreeBodyEnergy(false);
        ropt.setVerboseEnergies(true);
        if (System.getProperty("ro-ensembleNumber") == null && System.getProperty("ro-ensembleEnergy") == null) {
            logger.info(String.format(" Setting ensemble to default of number of walkers %d", worldSize));
            ropt.setEnsemble(worldSize);
        }
        ropt.setPrintFiles(false);
        ropt.setResiduesIgnoreNull(residueList);

        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(false);
        RotamerLibrary.measureRotamers(residueList, false);

        String oldLazyMat = System.getProperty("ro-lazyMatrix");
        System.setProperty("ro-lazyMatrix", "true");

        ropt.optimize(RotamerOptimization.Algorithm.ALL);
        ropt.setCoordinatesToEnsemble(rank);

        // One final energy call to ensure the coordinates are properly set at the
        // end of rotamer optimization.
        double[] xyz = new double[pot.getNumberOfVariables()];
        pot.getCoordinates(xyz);
        logger.info(" Final Optimized Energy:");
        pot.energy(xyz, true);

        if (linter != null) {
            linter.setLambda(initLam);
        }

        if (oldLazyMat != null) {
            System.setProperty("ro-lazyMatrix", oldLazyMat);
        } else {
            System.clearProperty("ro-lazyMatrix");
        }
    }
}
