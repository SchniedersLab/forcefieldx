
package ffx.algorithms.groovy

import edu.rit.pj.ParallelTeam
import ffx.algorithms.AbstractOSRW
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.BarostatOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.cli.LambdaParticleOptions
import ffx.algorithms.cli.MultiDynamicsOptions
import ffx.algorithms.cli.OSRWOptions
import ffx.algorithms.cli.RandomSymopOptions
import ffx.algorithms.cli.ThermodynamicsOptions
import ffx.algorithms.cli.WriteoutOptions
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import org.apache.commons.configuration2.Configuration
import picocli.CommandLine
import org.apache.commons.io.FilenameUtils
import edu.rit.pj.Comm
import ffx.crystal.CrystalPotential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.parameters.ForceField

import java.util.stream.Collectors

/**
 * The Thermodynamics script uses the Transition-Tempered Orthogonal Space Random Walk
 * algorithm to estimate a free energy.
 * <br>
 * Usage:
 * <br>
 * ffxc NewThermodynamics [options] &lt;filename&gt [file2...];
 */
@CommandLine.Command(description = " Use the Transition-Tempered Orthogonal Space Random Walk algorithm to estimate a free energy.", name = "ffxc Thermodynamics")
class NewThermodynamics extends AlgorithmsScript {

    @CommandLine.Mixin
    private AlchemicalOptions alchemical

    @CommandLine.Mixin
    private TopologyOptions topology

    @CommandLine.Mixin
    private BarostatOptions barostat;

    @CommandLine.Mixin
    private DynamicsOptions dynamics;

    @CommandLine.Mixin
    private WriteoutOptions writeout;

    @CommandLine.Mixin
    private LambdaParticleOptions lambdaParticle;

    @CommandLine.Mixin
    private MultiDynamicsOptions multidynamics;

    @CommandLine.Mixin
    private OSRWOptions osrwOptions;

    @CommandLine.Mixin
    private RandomSymopOptions randomSymop;

    @CommandLine.Mixin
    private ThermodynamicsOptions thermodynamics;

    /**
     * The final argument(s) should be one or more filenames.
     */
    @CommandLine.Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null;

    private int threadsAvail = ParallelTeam.getDefaultThreadCount();
    private int threadsPer = threadsAvail;
    MolecularAssembly[] topologies;

    CrystalPotential potential;
    AbstractOSRW osrw;

    private Configuration additionalProperties;

    /**
     * Sets an optional Configuration with additional properties.
     * @param additionalProps
     */
    void setProperties(Configuration additionalProps) {
        this.additionalProperties = additionalProps;
    }

    @Override
    NewThermodynamics run() {

        // Begin boilerplate "make a topology" code.

        if (!init()) {
            return;
        }

        List<String> arguments = filenames;
        // Check nArgs should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1;
        nArgs = (nArgs < 1) ? 1 : nArgs;

        topologies = new MolecularAssembly[nArgs];

        int numParallel = topology.getNumParallel(threadsAvail, nArgs);
        threadsPer = threadsAvail / numParallel;

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
           The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nArgs == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        double initLambda = alchemical.getInitialLambda()

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        List<MolecularAssembly> topologyList = new ArrayList<>(nArgs);

        Comm world = Comm.world();
        int size = world.size();
        int rank = (size > 1) ? world.rank() : 0;

        // Segment of code for MultiDynamics and OSRW.
        List<File> structureFiles = arguments.stream().
                map{fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath())}.
                collect(Collectors.toList())

        File firstStructure = structureFiles.get(0);
        String baseFilename = FilenameUtils.removeExtension(firstStructure.getName());
        File histogramRestart = new File(baseFilename + ".his");

        // For a multi-process job, try to get the restart files from rank sub-directories.
        String withRankName = baseFilename;
        if (size > 1) {
            List<File> rankedFiles = new ArrayList<>(nArgs);
            for (File structureFile : structureFiles) {
                File rankDirectory = new File(structureFile.getParent() + File.separator
                        + Integer.toString(rank));
                if (!rankDirectory.exists()) {
                    rankDirectory.mkdir();
                }
                withRankName = rankDirectory.getPath() + File.separator + baseFilename;
                rankedFiles.add(new File(rankDirectory.getPath() + File.separator + structureFile.getName()));
            }
            structureFiles = rankedFiles;
        }

        File lambdaRestart = new File(withRankName + ".lam");
        File dyn = new File(withRankName + ".dyn");
        /**
         * Used for the obsolete traversals option.
         * File lambdaOneFile = new File(withRankName + ".lam1");
         * File lambdaZeroFile = new File(withRankName + ".lam0");
         */
        // End multidynamics-relevant code.

        /**
         * Read in files.
         */
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = algorithmFunctions.getActiveAssembly()
            if (mola == null) {
                return helpString()
            }
            arguments = new ArrayList<>()
            arguments.add(mola.getFile().getName())
            topologyList.add(alchemical.processFile(Optional.of(topology), mola, 0))
        } else {
            logger.info(String.format(" Initializing %d topologies...", nArgs))
            for (int i = 0; i < nArgs; i++) {
                topologyList.add(multidynamics.openFile(algorithmFunctions, Optional.of(topology), threadsPer, arguments.get(i), i, alchemical, structureFiles.get(i), rank));
            }
        }

        MolecularAssembly[] topologies = topologyList.toArray(new MolecularAssembly[topologyList.size()])

        /**
         * Configure the potential to test.
         */
        StringBuilder sb = new StringBuilder("\n Running Transition-Tempered Orthogonal Space Random Walk for ")
        potential = (CrystalPotential) topology.assemblePotential(topologies, threadsAvail, sb)

        logger.info(sb.toString());

        LambdaInterface linter = (LambdaInterface) potential;

        // End of boilerplate code.

        if (!dyn.exists()) {
            dyn = null;
        }
        boolean lamExists = lambdaRestart.exists();
        boolean hisExists = histogramRestart.exists();

        logger.info(" Starting energy (before .dyn restart loaded):");
        boolean updatesDisabled = topologies[0].getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);
        linter.setLambda(initLambda);

        potential.energy(x, true);

        if (nArgs == 1) {
            randomSymop.randomize(topologies[0], potential);
        }

        multidynamics.distribute(topologies, potential, algorithmFunctions, rank, size);

        osrw = osrwOptions.constructOSRW(potential, lambdaRestart, histogramRestart, topologies[0], additionalProperties, dynamics, multidynamics, thermodynamics, algorithmListener);

        // Can be either the TT-OSRW or a Barostat on top of it.
        // Cannot be the Potential underneath the TT-OSRW.
        CrystalPotential osrwPotential = osrwOptions.applyAllOSRWOptions(osrw, topologies[0], dynamics, lambdaParticle, alchemical, barostat, lamExists, hisExists);

        // Old code for obsolete options.

        /*osrw.setResetStatistics(options.reset);
        if (options.traversals) {
            if (nArgs == 1) {
                osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[0]);
            } else if (nArgs == 2) {
                osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[1]);
            }
        }*/

        if (false) { // TODO: MC-OSRW code.
            /*if (options.mc) {
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
            mcOSRW.sample()*/
        } else {
            osrwOptions.beginMDOSRW(osrw, topologies, osrwPotential, dynamics, writeout, thermodynamics, dyn, algorithmListener);
        }

        return this;
    }

    AbstractOSRW getOSRW() {
        return osrw;
    }

    public CrystalPotential getPotential() {
        return potential;
    }
}

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
