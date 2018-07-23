package ffx.xray.groovy

import ffx.numerics.Potential
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import ffx.algorithms.MolecularDynamics
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.potential.MolecularAssembly
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The X-ray Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Dynamics [options] &lt;filename&gt;
 */
@Command(description = " Run Dynamics on an X-ray target.", name = "ffxc xray.Dynamics")
class Dynamics extends AlgorithmsScript {

    @Mixin
    XrayOptions xrayOptions

    @Mixin
    DynamicsOptions dynamicsOptions

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
    private List<String> filenames
    private RefinementEnergy refinementEnergy;

    @Override
    Dynamics run() {

        if (!init()) {
            return this
        }

        dynamicsOptions.init()
        xrayOptions.init()

        String modelfilename
        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath();
        }

        logger.info("\n Running xray.Dynamics on " + modelfilename)

        // Load parsed X-ray properties.
        CompositeConfiguration properties = assemblies[0].getProperties()
        xrayOptions.setProperties(parseResult, properties)

        // Set up diffraction data (can be multiple files)
        List<DiffractionData> diffractionFiles = xrayOptions.processData(filenames, assemblies)
        DiffractionData diffractionData = new DiffractionData(assemblies, properties,
                xrayOptions.solventModel, diffractionFiles.toArray(new DiffractionFile[diffractionFiles.size()]))

        diffractionData.scaleBulkFit()
        diffractionData.printStats()
        algorithmFunctions.energy(assemblies[0])

        RefinementEnergy refinementEnergy = new RefinementEnergy(diffractionData, xrayOptions.refinementMode)

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn")
        if (!dyn.exists()) {
            dyn = null
        }

        MolecularDynamics molDyn = dynamicsOptions.getDynamics(refinementEnergy, activeAssembly, algorithmListener)
        refinementEnergy.setThermostat(molDyn.getThermostat())
        boolean initVelocities = true
        molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt,
                dynamicsOptions.report, dynamicsOptions.write, dynamicsOptions.temp, initVelocities, dyn)

        return this
    }

    @Override
    public List<Potential> getPotentials() {
        return refinementEnergy == null ? new ArrayList<>() : Collections.singletonList(refinementEnergy);
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
