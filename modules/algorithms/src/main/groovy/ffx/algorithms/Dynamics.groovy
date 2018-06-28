package ffx.algorithms

import ffx.algorithms.cli.WriteoutOptions
import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.BarostatOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.parameters.ForceField

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Dynamics script implements molecular and stochastic dynamics algorithms.
 * <br>
 * Usage:
 * <br>
 * ffxc Dynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run dynamics on a system.", name = "ffxc Dynamics")
class Dynamics extends AlgorithmsScript {

    @Mixin
    DynamicsOptions dynamics;

    @Mixin
    BarostatOptions barostatOpt;

    @Mixin
    WriteoutOptions writeout;

    /**
     * -r or --repEx to execute temperature replica exchange.
     */
    @Option(names = ['-x', '--repEx'], paramLabel = 'false',
            description = 'Execute temperature replica exchange')
    boolean repEx = false;

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    def run() {

        if (!init()) {
            return
        }

        dynamics.init()

        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
        }

        File structureFile = new File(FilenameUtils.normalize(modelfilename))
        structureFile = new File(structureFile.getAbsolutePath())
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName())

        Potential potential = activeAssembly.getPotentialEnergy()
        logger.info(" Starting energy (before .dyn restart loaded):")
        boolean updatesDisabled = activeAssembly.getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false)
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart")
        }
        double[] x = new double[potential.getNumberOfVariables()]
        potential.getCoordinates(x)

        potential.energy(x, true)

        if (barostatOpt.pressure > 0) {
            if (potential instanceof ffx.potential.ForceFieldEnergyOpenMM) {
                logger.warning(" NPT with OpenMM acceleration is still experimental and may not function correctly.")
            }
            logger.info(String.format(" Running NPT dynamics at pressure %7.4g", barostatOpt.pressure))
            CrystalPotential crystalPotential = (CrystalPotential) potential
            Barostat barostat = new Barostat(activeAssembly, crystalPotential)
            barostat.setPressure(barostatOpt.pressure)
            barostat.setMaxDensity(barostatOpt.maxDensity)
            barostat.setMinDensity(barostatOpt.minDensity)
            barostat.setMaxSideMove(barostatOpt.maxSideMove)
            barostat.setMaxAngleMove(barostatOpt.maxAngleMove)
            barostat.setMeanBarostatInterval(barostatOpt.meanInterval)
            potential = barostat
        }

        Comm world = Comm.world()
        int size = world.size()

        if (!repEx || size < 2) {
            logger.info("\n Running molecular dynamics on " + modelfilename)
            // Restart File
            File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn")
            if (!dyn.exists()) {
                dyn = null
            }

            MolecularDynamics molDyn = dynamics.getDynamics(writeout, potential, activeAssembly, algorithmListener)

            molDyn.dynamic(dynamics.steps, dynamics.dt,
                    dynamics.report, dynamics.write, dynamics.temp, true, dyn)

        } else {
            logger.info("\n Running replica exchange molecular dynamics on " + modelfilename)
            int rank = world.rank()
            File rankDirectory = new File(structureFile.getParent() + File.separator
                    + Integer.toString(rank))
            if (!rankDirectory.exists()) {
                rankDirectory.mkdir()
            }
            String withRankName = rankDirectory.getPath() + File.separator + baseFilename
            File dyn = new File(withRankName + ".dyn")
            if (!dyn.exists()) {
                dyn = null
            }

            MolecularDynamics molDyn = dynamics.getDynamics(writeout, potential, activeAssembly, algorithmListener)
            ReplicaExchange replicaExchange = new ReplicaExchange(molDyn, algorithmListener, dynamics.temp)

            int totalSteps = dynamics.steps
            int nSteps = 100
            int cycles = totalSteps / nSteps
            if (cycles <= 0) {
                cycles = 1
            }

            replicaExchange.sample(cycles, nSteps, dynamics.dt, dynamics.report, dynamics.write)
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
