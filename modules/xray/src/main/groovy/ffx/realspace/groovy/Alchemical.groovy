package ffx.realspace.groovy

import ffx.algorithms.AbstractOSRW
import ffx.algorithms.cli.AlgorithmsScript
import ffx.realspace.parsers.RealSpaceFile
import groovy.ui.SystemOutputInterceptor
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
import ffx.realspace.RealSpaceData

import ffx.algorithms.cli.DynamicsOptions
import ffx.realspace.cli.RealSpaceOptions

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

@Command(description = " Alchemical changes on a Real Space target.", name = "ffxc realspace.Alchemical")
class Alchemical extends AlgorithmsScript {

    @Mixin
    DynamicsOptions dynamicsOptions

    @Mixin
    RealSpaceOptions realSpaceOptions

    private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName())

    /**
     * -I or --onlyions sets whether or not ion positions are optimized (default is false; must set at least one of either '-W' or '-I') (only one type of ion is chosen).
     */
    @Option(names = ["-I", "--onlyions"], paramLabel = 'false',
            description = 'Set to only optimize ions (of a single type).')
    boolean onlyIons = false

    /**
     * --itype or --iontype Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.
     */
    @Option(names = ["--itype", "--iontype"], paramLabel = 'null',
            description = 'Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.')
    String [] iontype = null

    /**
     * --N or --neutralize Adds more of the selected ion in order to neutralize the crystal's charge.
     */
    @Option(names = ["-N", "--neutralize"], paramLabel = 'false',
            description = 'Neutralize the crystal\'s charge by adding more of the selected ion')
    boolean neutralize = false

    /**
     * -W or --onlywaters sets whether or not water positions are optimized (default is false; must set at least one of either '-W' or '-I').
     */
    @Option(names = ["-W", "--onlywaters"], paramLabel = 'false',
            description = 'Set to only optimize waters.')
    boolean onlyWaters = false

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames

    // Value of Lambda.
    double lambda = 0.0

    // ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
    ThermostatEnum thermostat = ThermostatEnum.ADIABATIC

    // IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC ]
    IntegratorEnum integrator = IntegratorEnum.STOCHASTIC

    // File type of coordinate snapshots to write out.
    String fileType = "PDB"

    // OSRW
    boolean runOSRW = true

    // Reset velocities (ignored if a restart file is given)
    boolean initVelocities = true

    private AbstractOSRW osrw;

    @Override
    Alchemical run() {

        if (!init()) {
            return
        }
        dynamicsOptions.init()
        System.setProperty("lambdaterm", "true")

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
            assemblies = { activeAssembly }
        }

        logger.info("\n Running Alchemical Changes on " + modelfilename)

        File structureFile = new File(FilenameUtils.normalize(modelfilename))
        structureFile = new File(structureFile.getAbsolutePath())
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName())
        File histogramRestart = new File(baseFilename + ".his")
        File lambdaRestart = new File(baseFilename + ".lam")
        File dyn = new File(baseFilename + ".dyn")

        Comm world = Comm.world()
        int size = world.size()
        int rank = 0
        double[] energyArray = new double[world.size()]
        for (int i = 0; i < world.size(); i++) {
            energyArray[i] = Double.MAX_VALUE
        }

        // For a multi-process job, try to get the restart files from rank sub-directories.
        if (size > 1) {
            rank = world.rank()
            File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank))
            if (!rankDirectory.exists()) {
                rankDirectory.mkdir()
            }
            lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam")
            dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn")
            structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName())
        }
        if (!dyn.exists()) {
            dyn = null
        }

        // Set built atoms active/use flags to true (false for other atoms).
        Atom[] atoms = activeAssembly.getAtomArray()

        // Get a reference to the first system's ForceFieldEnergy.
        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
        forceFieldEnergy.setPrintOnFailure(false, false)

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

//      Consider the option of creating a composite lambda gradient from vapor phase to crystal phase
        if (!onlyWaters) {
            logger.info("Doing ions.")
            if (ions == null || ions.size() == 0) {
                logger.info("\n Please add an ion to the PDB file to scan with.")
                return
            }
            for (MSNode msNode : ions) {
                try {
                    logger.info("Selecting ion.")
                    if (msNode.getAtomList().name == ionType) {
                        logger.info("Ion has been selected.")
                        for (Atom atom : msNode.getAtomList()) {
                            System.out.println("Activating ions")
                            atom.setUse(true)
                            atom.setActive(true)
                            atom.setApplyLambda(true)
                            logger.info(" Alchemical atom: " + atom.toString())
                        }
                    }
                } catch (MissingPropertyException e) {
                    logger.info("Ion has not been selected.")
                    if (neutralize) {
                        logger.info("Neutralizing crystal.")
                        double ionCharge = 0
                        for (Atom atom : msNode.getAtomList()) {
                            ionCharge += atom.multipoleType.getCharge()
                        }
                        logger.info("Ion charge is: " + ionCharge.toString())
                        int numIons = (int) -1*(Math.ceil(crystalCharge/ionCharge))
                        if (numIons > 0) {
                            logger.info(numIons + " " + msNode.getAtomList().name
                                    + " ions needed to neutralize the crystal.")
                            ionType = msNode.getAtomList().name
                            for (Atom atom : msNode.getAtomList()) {
                                atom.setUse(true)
                                atom.setActive(true)
                                atom.setApplyLambda(true)
                                logger.info(" Alchemical atom: " + atom.toString())
                            }
                        }
                    }
                    else {
                        ionType = msNode.getAtomList().name
                        for (Atom atom : msNode.getAtomList()) {
                            atom.setUse(true)
                            atom.setActive(true)
                            atom.setApplyLambda(true)
                            logger.info(" Alchemical atom: " + atom.toString())
                        }
                    }
                }
            }
        }

        // Lambdize waters for position optimization
        if (!onlyIons) {
            logger.info(waters.size + " water molecules in this PDB.")
            if (waters == null || waters.size() == 0) {
                logger.info("\n Please add water to the PDB file to scan with.")
                return
            }
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

        List<RealSpaceFile> mapfiles = realSpaceOptions.processData(filenames, assemblies)

        RealSpaceData realSpaceData = new RealSpaceData(activeAssembly, activeAssembly.getProperties(),
                activeAssembly.getParallelTeam(), mapfiles.toArray(new RealSpaceFile[mapfiles.size()]))

        RefinementEnergy refinementEnergy = new RefinementEnergy(realSpaceData, realSpaceOptions.refinementMode, null)
        refinementEnergy.setLambda(lambda)

        boolean asynchronous = true

        osrw = new TransitionTemperedOSRW(refinementEnergy, refinementEnergy, lambdaRestart, histogramRestart,
                assemblies[0].getProperties(), dynamicsOptions.temp, dynamicsOptions.dt, dynamicsOptions.report,
                dynamicsOptions.write, true, false, algorithmFunctions.getDefaultListener())

        osrw.setLambda(lambda);
        osrw.setThetaMass(5.0e-19);
        osrw.setOptimization(true, activeAssembly);
        // Create the MolecularDynamics instance.

        MolecularDynamics molDyn = new MolecularDynamics(assemblies[0], osrw, assemblies[0].getProperties(),
                null, thermostat, integrator)

        algorithmFunctions.energy(assemblies[0])
        molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write, dynamicsOptions.temp, true,
                fileType, dynamicsOptions.write, dyn)
        logger.info(" Searching for low energy coordinates")
        double[] lowEnergyCoordinates = osrw.getLowEnergyLoop()
        double currentOSRWOptimum = osrw.getOSRWOptimum()
        if (lowEnergyCoordinates != null) {
            forceFieldEnergy.setCoordinates(lowEnergyCoordinates)
            logger.info("\n Minimum coordinates found: " + lowEnergyCoordinates)
        } else {
            logger.info(" OSRW stage did not succeed in finding a minimum.")
        }
    }

    @Override
    public List<Potential> getPotentials() {
        return osrw == null ? new ArrayList<>() : Collections.singletonList(osrw);
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