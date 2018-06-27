package ffx.realspace

import ffx.algorithms.cli.AlgorithmsScript
import ffx.realspace.cli.RealSpaceOptions
import ffx.realspace.parsers.RealSpaceFile
import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.algorithms.MolecularDynamics
import ffx.algorithms.TransitionTemperedOSRW
import ffx.algorithms.integrators.IntegratorEnum
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode

import ffx.algorithms.cli.DynamicsOptions
//import ffx.realspace.cli.RealSpaceOptions

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option

import java.util.logging.Logger

/**
 * The Alchemical Changes script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Alchemical [options] &lt;filename&gt;
 */

@Command(description = " Simulated annealing on a Real Space target.", name = "ffxc realspace.Alchemical")
class Alchemical extends AlgorithmsScript {

    @Mixin
    DynamicsOptions dynamicsOptions

//    @Mixin
//    RealSpaceOptions realSpaceOptions

    private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

    /**
     * -X or --wA The weight of the real space data (wA).
     */
    @Option(names = ["-X", "--wA"], paramLabel = "1.0",
            description = "The weight of the real space data (wA).")
    double wA = 1.0;

    /**
     * -y or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
     */
    @Option(names = ["-y", "--data"], split = ",",
            description = "Specify input data filename and its weight (wA) (e.g. -y filename,1.0).")
    String[] data = null;

    /**
     * -u or --mode sets the desired refinement mode
     * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
     */
    @Option(names = ["-u", "--mode"], paramLabel = "coordinates",
            description = "Refinement mode: coordinates, bfactors and/or occupancies.")
    String modeString = "coordinates";

    /**
     * The refinement mode to use.
     */
    RefinementMode refinementMode = RefinementMode.COORDINATES;
    /**
     * --ions or --doions Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.
     */
    @Option(names = ["--ions", "--doions"], paramLabel = '',
            description = '* --ion Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.')
    boolean opt_ions = false
    /**
     * --itype Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.
     */
    @Option(names = ["--itype", "--iontype"], paramLabel = '',
            description = 'Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.')
    String [] iontype = ''
    /**
     * --neut or --neutralize Adds more of the selected ion in order to neutralize the crystal's charge.
     */
    @Option(names = ["--neut", "--neutralize"], paramLabel = 'false',
            description = 'Add more of the selected ion to neutralize the crystal\'s charge')
    boolean neutralize = false
    /**
     * -W or --waters sets whether or not water positions are also optimized (default is false).
     */
    @Option(names = "--dowaters", paramLabel = 'true',
            description = 'Specify whether or not to optimize water positions.')
    boolean opt_waters = true

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames

    // Temperature in degrees Kelvin.
    double temperature = 298.15;

    // Time step in femtoseconds.
    double timeStep = 1.0;

    // Frequency to log thermodynamics information in picoseconds.
    double printInterval = 0.01;

    // Frequency to write out coordinates in picoseconds.
    double saveInterval = 100;

    // Value of Lambda.
    double lambda = 1.0;

    // ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
    ThermostatEnum thermostat = ThermostatEnum.ADIABATIC;

    // IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC ]
    IntegratorEnum integrator = IntegratorEnum.STOCHASTIC;

    // OSRW
    boolean runOSRW = true;

    // Number of molecular dynamics steps: default is 100 nanoseconds.
    int nSteps = 50000;

    // Reset velocities (ignored if a restart file is given)
    boolean initVelocities = true;

    // File type of coordinate snapshots to write out.
    String fileType = "PDB";

    // Frequency to write out restart information in picoseconds.
    double restartInterval = 1.0;

    def run() {

        if (!init()) {
            return
        }

        dynamicsOptions.init()

        String modelfilename
        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
            assemblies = { activeAssembly };
        }

        logger.info("\n Running Alchemical Changes on " + modelfilename)

        File structureFile = new File(FilenameUtils.normalize(modelfilename));
        structureFile = new File(structureFile.getAbsolutePath());
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
        File histogramRestart = new File(baseFilename + ".his");
        File lambdaRestart = new File(baseFilename + ".lam");
        File dyn = new File(baseFilename + ".dyn");

        Comm world = Comm.world();
        int size = world.size();
        int rank = 0;
        double[] energyArray = new double[world.size()];
        for (int i = 0; i < world.size(); i++) {
            energyArray[i] = Double.MAX_VALUE;
        }

        // For a multi-process job, try to get the restart files from rank sub-directories.
        if (size > 1) {
            rank = world.rank();
            File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank));
            if (!rankDirectory.exists()) {
                rankDirectory.mkdir();
            }
            lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
            dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
            structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
        }

        if (!dyn.exists()) {
            dyn = null;
        }

        // Set built atoms active/use flags to true (false for other atoms).
        Atom[] atoms = activeAssembly.getAtomArray();

        // Get a reference to the first system's ForceFieldEnergy.
        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy();
        forceFieldEnergy.setPrintOnFailure(false, false);

        // Configure all atoms to:
        // 1) be used in the potential
        // 2) be inactive (i.e. cannot move)
        // 3) not be controlled by the lambda state variable.
        for (int i = 0; i <= atoms.length; i++) {
            Atom ai = atoms[i - 1]
            ai.setUse(true)
            ai.setActive(false)
            ai.setApplyLambda(false)
        }

        double crystalCharge = activeAssembly.getCharge(true)
        logger.info(" Overall crystal charge: " + crystalCharge)
        ArrayList<MSNode> ions = assemblies[0].getIons()
        ArrayList<MSNode> waters = assemblies[0].getWaters()

        if (opt_ions) {
            if (ions == null || ions.size() == 0) {
                logger.info("\n Please add an ion to the PDB file to scan with.")
                return
            }

            for (MSNode msNode : ions) {
                for (Atom atom : msNode.getAtomList()) {
                    // Scan with the last ion in the file.
                    atom.setUse(true)
                    atom.setActive(true)
                    atom.setApplyLambda(true)
                    logger.info(" Alchemical atom: " + atom.toString())
                }
            }
        }

        // Lambdize waters for position optimization, if this option was set to true
        if (opt_waters) {
            for (MSNode msNode : waters) {
                for (Atom atom : msNode.getAtomList()) {
                    // Scan with the last ion in the file.
                    atom.setUse(true)
                    atom.setActive(true)
                    atom.setApplyLambda(true)
                    logger.info(" Water atom:      " + atom.toString())
                }
            }
        }

        // set up real space map data (can be multiple files)
        List mapfiles = new ArrayList();
        if (filenames.size() > 1) {
            RealSpaceFile realspacefile = new RealSpaceFile(filenames.get(1), 1.0);
            mapfiles.add(realspacefile);
        }
        if (data) {
            for (int i = 0; i < data.size(); i += 2) {
                double wA = Double.parseDouble(data[i + 1]);
                RealSpaceFile realspacefile = new RealSpaceFile(data[i], wA);
                mapfiles.add(realspacefile);
            }
        }

        if (mapfiles.size() == 0) {
            RealSpaceFile realspacefile = new RealSpaceFile(assemblies[0], 1.0);
            mapfiles.add(realspacefile);
        }

        RealSpaceData realSpaceData = new RealSpaceData(activeAssembly, activeAssembly.getProperties(),
                activeAssembly.getParallelTeam(), mapfiles.toArray(new RealSpaceFile[mapfiles.size()]));

        // Replace map file set-up above with commented lines below once RealSpaceOptions can be used (after conflicting flags are resolved)
//        List<RealSpaceFile> mapfiles = realSpaceOptions.processData(filenames, assemblies)
//
//        RealSpaceData realSpaceData = new RealSpaceData(activeAssembly, activeAssembly.getProperties(),
//                activeAssembly.getParallelTeam(), mapfiles.toArray(new RealSpaceFile[mapfiles.size()]))
//

        RefinementEnergy refinementEnergy = new RefinementEnergy(realSpaceData, RefinementMode.COORDINATES, null);
        refinementEnergy.setLambda(lambda);

        boolean asynchronous = true;
        Potential osrw;

        osrw = new TransitionTemperedOSRW(refinementEnergy, refinementEnergy, lambdaRestart, histogramRestart,
                assemblies[0].getProperties(), dynamicsOptions.temp, dynamicsOptions.dt, dynamicsOptions.report,
                dynamicsOptions.write, true, false, algorithmFunctions.getDefaultListener());

        osrw.setLambda(lambda);
        osrw.setThetaMass(5.0e-19);
        osrw.setOptimization(true, activeAssembly);
        // Create the MolecularDynamics instance.
        MolecularDynamics molDyn = new MolecularDynamics(assemblies[0], osrw, assemblies[0].getProperties(),
                null, dynamicsOptions.thermostat, dynamicsOptions.integrator);

        algorithmFunctions.energy(assemblies[0]);

        molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write, dynamicsOptions.temp, true,
                fileType, dynamicsOptions.write, dyn);

        logger.info(" Searching for low energy coordinates");
        double[] lowEnergyCoordinates = osrw.getLowEnergyLoop();
        double currentOSRWOptimum = osrw.getOSRWOptimum();
        if (lowEnergyCoordinates != null) {
            forceFieldEnergy.setCoordinates(lowEnergyCoordinates);
        } else {
            logger.info(" OSRW stage did not succeed in finding a minimum.");
        }
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
