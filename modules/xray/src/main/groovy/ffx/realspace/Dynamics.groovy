package ffx.realspace

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.xray.cli.RealSpaceOptions

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.MolecularDynamics
import ffx.potential.MolecularAssembly
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode

import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The Real Space Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Dynamics [options] &lt;filename&gt;
 */
@Command(description = " Molecular dynamics on a Real Space target.", name = "ffxc realspace.Dynamics")
class Dynamics extends AlgorithmsScript {

    /**
     * -l or --log sets the thermodynamics logging frequency in picoseconds (0.1 psec default).
     */
    @Option(names = "--log", paramLabel = "0.25",
            description = "Interval to report thermodynamics (psec).")
    double log = 0.25

    @Mixin
    DynamicsOptions dynamicsOptions

    @Mixin
    RealSpaceOptions realSpaceOptions

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames

    def run() {

        if (!init()) {
            return
        }

        dynamicsOptions.init()
        realSpaceOptions.init()

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

        logger.info("\n Running Real Space Dynamics on " + modelfilename)

        List<RealSpaceFile> mapfiles = realSpaceOptions.processData(filenames, assemblies)

        RealSpaceData realspacedata = new RealSpaceData(activeAssembly, activeAssembly.getProperties(),
                activeAssembly.getParallelTeam(),
                mapfiles.toArray(new RealSpaceFile[mapfiles.size()]))

        algorithmFunctions.energy(assemblies[0])

        RefinementEnergy refinementEnergy = new RefinementEnergy(realspacedata, RefinementMode.COORDINATES)

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn")
        if (!dyn.exists()) {
            dyn = null
        }
        MolecularDynamics molDyn = new MolecularDynamics(activeAssembly, refinementEnergy,
                activeAssembly.getProperties(), refinementEnergy, dynamicsOptions.thermostat, dynamicsOptions.integrator, null)
        refinementEnergy.setThermostat(molDyn.getThermostat())

        // Reset velocities (ignored if a restart file is given)
        boolean initVelocities = true
        molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt, log, dynamicsOptions.write, dynamicsOptions.temp, initVelocities, dyn)
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