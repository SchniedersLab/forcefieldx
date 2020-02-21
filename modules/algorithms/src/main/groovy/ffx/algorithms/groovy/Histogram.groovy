//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Histogram script prints out a Orthogonal Space histogram from a *.his file.
 * <br>
 * Usage:
 * <br>
 * ffxc Histogram [options] &lt;filename&gt;
 */
@Command(description = " Evaluate the Orthogonal Space Histogram.", name = "ffxc Histogram")
class Histogram extends AlgorithmsScript {

    /**
     * -p or --pmf Save the histogram, PMF and 2D bias to files.
     */
    @Option(names = ['-p', '--pmf'], paramLabel = 'false',
            description = 'Save the bias histogram to a file.')
    boolean pmf = false

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "XYZ or PDB input files.")
    private List<String> filenames

    private OrthogonalSpaceTempering orthogonalSpaceTempering
    private File saveDir = null

    @Override
    Histogram run() {

        if (!init()) {
            return this
        }

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = [algorithmFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        println("\n Evaluating Histogram for " + modelFilename)

        File structureFile = new File(FilenameUtils.normalize(modelFilename))
        structureFile = new File(structureFile.getAbsolutePath())
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName())
        File histogramRestart = new File(baseFilename + ".his")
        File lambdaRestart = null

        // Get a reference to the active system's ForceFieldEnergy and atom array.
        ForceFieldEnergy energy = activeAssembly.getPotentialEnergy()

        // Print the current energy
        energy.energy(true, true)

        // These fields are needed for the OST constructor, but otherwise are not used.
        boolean asynchronous = false
        double timeStep = 1.0
        double printInterval = 1.0
        double saveInterval = 100.0
        double temperature = 298.15

        modelFilename = activeAssembly.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        orthogonalSpaceTempering = new OrthogonalSpaceTempering(energy, energy, lambdaRestart, histogramRestart,
                activeAssembly.getProperties(), temperature, timeStep, printInterval,
                saveInterval, asynchronous, algorithmListener)

        if (pmf) {
            orthogonalSpaceTempering.setMolecularAssembly(activeAssembly)
            OrthogonalSpaceTempering.Histogram histogram = orthogonalSpaceTempering.getHistogram()
            histogram.updateFLambda(false, true)
            StringBuffer sb = histogram.evaluateTotalPMF()

            String dirName = saveDir.toString() + File.separator
            String file = dirName + "pmf.txt"
            logger.info(" Writing " + file)
            FileWriter fileWriter = new FileWriter(file)
            fileWriter.write(sb.toString())
            fileWriter.close()

            sb = histogram.evaluate2DPMF()
            file = dirName + "pmf.2D.txt"
            logger.info(" Writing " + file)
            fileWriter = new FileWriter(file)
            fileWriter.write(sb.toString())
            fileWriter.close()
        }
        return this
    }

    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (orthogonalSpaceTempering == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList(orthogonalSpaceTempering)
        }
        return potentials
    }

}
