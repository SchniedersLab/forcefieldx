package ffx.xray.groovy

import org.apache.commons.configuration.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.potential.MolecularAssembly
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The X-ray Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Minimize [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Refine an X-ray/Neutron target.", name = "ffxc xray.Minimize")
class Minimize extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    @Mixin
    XrayOptions xrayOptions

    /**
     * -t or --threeStage Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true).
     */
    @Option(names = ['-t', '--threeStage'], paramLabel = 'false',
            description = 'Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true)')
    boolean threeStage = false
    /**
     * -E or --eps3 RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determine eps for each stage).
     */
    @Option(names = ['-E', '--eps3'], paramLabel = '-1.0', arity = '3',
            description = 'RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determines eps for each stage).')
    double[] eps3 = [-1.0, -1.0, -1.0]

    /**
     * -s or --suffix Specify the suffix to apply to output files. For example, for 1abc_refine.pdb, write out 1abc_refine_refine.[pdb|mtz] at the end.
     */
    @Option(names = ['--suffix'], paramLabel = '_refine',
            description = 'Suffix to apply to files written out by minimization.')
    String suffix = "_refine"

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
    private List<String> filenames

    @Override
    Minimize run() {

        if (!init()) {
            return this
        }

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
            modelfilename = activeAssembly.getFile().getAbsolutePath()
            assemblies = [activeAssembly]
        }

        logger.info("\n Running xray.Minimize on " + modelfilename)

        if (assemblies.length > 1) {
            logger.info(assemblies.toString())
        }

        // Load parsed X-ray properties.
        CompositeConfiguration properties = activeAssembly.getProperties()
        xrayOptions.setProperties(parseResult, properties)

        // Set up diffraction data (can be multiple files)
        List<ffx.xray.DiffractionData> diffractionFiles = xrayOptions.processData(filenames, assemblies)

        ffx.xray.DiffractionData diffractionData = new ffx.xray.DiffractionData(assemblies, properties,
                xrayOptions.solventModel, diffractionFiles.toArray(new DiffractionFile[diffractionFiles.size()]))

        diffractionData.scaleBulkFit()
        diffractionData.printStats()

        algorithmFunctions.energy(activeAssembly)

        // RMS gradient convergence criteria for three stage refinement
        double coordeps = eps3[0]
        double beps = eps3[1]
        double occeps = eps3[2]

        // Maximum number of refinement cycles.
        int maxiter = minimizeOptions.iterations

        if (threeStage) {
            ffx.xray.RefinementMinimize refinementMinimize = new ffx.xray.RefinementMinimize(diffractionData, RefinementMode.COORDINATES)
            if (coordeps < 0.0) {
                coordeps = refinementMinimize.getEps()
            }

            if (maxiter < Integer.MAX_VALUE) {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", coordeps, maxiter));
            } else {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", coordeps));
            }

            refinementMinimize.minimize(coordeps, maxiter)
            diffractionData.scaleBulkFit()
            diffractionData.printStats()
            algorithmFunctions.energy(activeAssembly)

            refinementMinimize = new ffx.xray.RefinementMinimize(diffractionData, RefinementMode.BFACTORS)
            if (beps < 0.0) {
                beps = refinementMinimize.getEps()
            }

            if (maxiter < Integer.MAX_VALUE) {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", beps, maxiter));
            } else {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", beps));
            }

            refinementMinimize.minimize(beps, maxiter)
            diffractionData.scaleBulkFit()
            diffractionData.printStats()

            if (diffractionData.getAltResidues().size() > 0
                    || diffractionData.getAltMolecules().size() > 0) {
                refinementMinimize = new ffx.xray.RefinementMinimize(diffractionData, RefinementMode.OCCUPANCIES)
                if (occeps < 0.0) {
                    occeps = refinementMinimize.getEps()
                }

                if (maxiter < Integer.MAX_VALUE) {
                    logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", occeps, maxiter));
                } else {
                    logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", occeps));
                }

                refinementMinimize.minimize(occeps, maxiter)
                diffractionData.scaleBulkFit()
                diffractionData.printStats()
            } else {
                logger.info("Occupancy refinement not necessary, skipping")
            }
        } else {
            // Type of refinement.
            RefinementMode refinementMode = xrayOptions.refinementMode
            ffx.xray.RefinementMinimize refinementMinimize = new ffx.xray.RefinementMinimize(diffractionData, refinementMode)
            double eps = minimizeOptions.eps
            if (eps < 0.0) {
                eps = refinementMinimize.getEps()
            }

            if (maxiter < Integer.MAX_VALUE) {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps, maxiter));
            } else {
                logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", eps));
            }

            refinementMinimize.minimize(eps, maxiter)
            diffractionData.scaleBulkFit()
            diffractionData.printStats()
        }

        algorithmFunctions.energy(activeAssembly)

        diffractionData.writeModel(FilenameUtils.removeExtension(modelfilename) + suffix + ".pdb")
        diffractionData.writeData(FilenameUtils.removeExtension(modelfilename) + suffix + ".mtz")

        return this
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