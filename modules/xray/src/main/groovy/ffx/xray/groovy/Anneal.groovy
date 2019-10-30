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
package ffx.xray.groovy

import ffx.algorithms.optimize.anneal.SimulatedAnnealing
import ffx.potential.cli.WriteoutOptions
import org.apache.commons.configuration2.CompositeConfiguration
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.AnnealOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.cli.XrayOptions
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The X-ray Annealing script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Anneal [options] &lt;filename&gt;
 */
@Command(description = " Simulated annealing on an X-ray target.", name = "ffxc xray.Anneal")
class Anneal extends AlgorithmsScript {

    @Mixin
    XrayOptions xrayOptions

    @Mixin
    DynamicsOptions dynamics

    @Mixin
    WriteoutOptions writeout;

    @Mixin
    AnnealOptions anneal

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
    private List<String> filenames

    private SimulatedAnnealing simulatedAnnealing = null;

    private Potential potential;
    private RefinementEnergy refinementEnergy;

    @Override
    Anneal run() {

        if (!init()) {
            return this
        }

        dynamics.init()
        // Added vs. regular Anneal script.
        xrayOptions.init();

        String modelFilename
        MolecularAssembly[] assemblies;
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        }

        modelFilename = activeAssembly.getFile().getAbsolutePath();

        logger.info("\n Running simulated annealing on X-ray target including " + modelFilename + "\n")

        // Restart File
        File dyn = new File(FilenameUtils.removeExtension(modelFilename) + ".dyn")
        if (!dyn.exists()) {
            dyn = null
        }

        // Differs between regular Anneal and x-ray Anneal.
        CompositeConfiguration properties = activeAssembly.getProperties();
        DiffractionData diffractionData = xrayOptions.getDiffractionData(filenames, assemblies, parseResult);
        potential = xrayOptions.toXrayEnergy(diffractionData, assemblies, algorithmFunctions);
        simulatedAnnealing = anneal.createAnnealer(dynamics, activeAssembly,
                potential, properties,
                algorithmListener, dyn);

        simulatedAnnealing.setPrintInterval(dynamics.report);
        simulatedAnnealing.setSaveFrequency(dynamics.write);
        simulatedAnnealing.setRestartFrequency(dynamics.checkpoint)
        simulatedAnnealing.setTrajectorySteps(dynamics.trajSteps);

        simulatedAnnealing.anneal();

        diffractionData.scaleBulkFit()
        diffractionData.printStats()

        double[] x = new double[potential.getNumberOfVariables()];
        x = potential.getCoordinates(x);
        potential.energy(x, true)

        diffractionData.writeData(FilenameUtils.removeExtension(modelFilename) + ".mtz")

        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        fileName = FilenameUtils.removeExtension(fileName)

        writeout.saveFile(String.format("%s%s", dirName, fileName), algorithmFunctions, activeAssembly);

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return refinementEnergy == null ? Collections.emptyList() : Collections.singletonList(refinementEnergy);
    }
}

