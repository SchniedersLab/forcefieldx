package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.thermodynamics.AbstractOSRW
import ffx.algorithms.thermodynamics.OSRW
import ffx.algorithms.thermodynamics.TransitionTemperedOSRW
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Histogram script prints out a (TT-)OSRW histogram from a *.his file.
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
     * -u or --untempered Histogram for untempered OSRW.
     */
    @Option(names = ['-u', '--untempered'], paramLabel = 'false',
            description = 'Histogram for untempered OSRW.')
    boolean untempered = false

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "XYZ or PDB input files.")
    private List<String> filenames

    private AbstractOSRW abstractOSRW
    private File saveDir = null

    @Override
    Histogram run() {

        if (!init()) {
            return this
        }

        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
        }

        println("\n Evaluating Histogram for " + modelfilename)

        File structureFile = new File(FilenameUtils.normalize(modelfilename))
        structureFile = new File(structureFile.getAbsolutePath())
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName())
        File histogramRestart = new File(baseFilename + ".his")
        File lambdaRestart = null

        // Get a reference to the active system's ForceFieldEnergy and atom array.
        ForceFieldEnergy energy = activeAssembly.getPotentialEnergy()

        // Print the current energy
        energy.energy(true, true)

        // These fields are needed for the OSRW constructor, but otherwise are not used.
        boolean asynchronous = false
        double timeStep = 1.0
        double printInterval = 1.0
        double saveInterval = 100.0
        double temperature = 298.15

        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator

        if (!untempered) {
            TransitionTemperedOSRW osrw = new TransitionTemperedOSRW(energy, energy, lambdaRestart, histogramRestart,
                    activeAssembly.getProperties(), temperature, timeStep, printInterval,
                    saveInterval, asynchronous, algorithmListener)
            if (pmf) {
                osrw.setMolecularAssembly(activeAssembly)
                osrw.updateFLambda(false, true)

                StringBuffer sb = osrw.evaluateTotalPMF()
                String file = dirName + "pmf.txt"
                logger.info(" Writing " + file)
                FileWriter fileWriter = new FileWriter(file)
                fileWriter.write(sb.toString())
                fileWriter.close()

                sb = osrw.evaluate2DPMF()
                file = dirName + "pmf.2D.txt"
                logger.info(" Writing " + file)
                fileWriter = new FileWriter(file)
                fileWriter.write(sb.toString())
                fileWriter.close()
            }
            abstractOSRW = osrw
        } else {
            OSRW osrw = new OSRW(energy, energy, lambdaRestart, histogramRestart,
                    activeAssembly.getProperties(), temperature, timeStep, printInterval,
                    saveInterval, asynchronous, algorithmListener)
            if (pmf) {
                osrw.setMolecularAssembly(activeAssembly)
                osrw.updateFLambda(false, true)

                StringBuffer sb = osrw.evaluateTotalPMF()
                String file = dirName + "pmf.txt"
                logger.info(" Writing " + file)
                FileWriter fileWriter = new FileWriter(file)
                fileWriter.write(sb.toString())
                fileWriter.close()

                sb = osrw.evaluate2DPMF()
                file = dirName + "pmf.2D.txt"
                logger.info(" Writing " + file)
                fileWriter = new FileWriter(file)
                fileWriter.write(sb.toString())
                fileWriter.close()
            }
            abstractOSRW = osrw
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return abstractOSRW == null ? Collections.emptyList() : Collections.singletonList(abstractOSRW)
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
