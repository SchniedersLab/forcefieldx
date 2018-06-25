package ffx.xray

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.cli.TimerOptions
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The X-ray Timer script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Timer [options] &lt;filename&gt;
 */
@Command(description = " Time calculation of the X-ray target.", name = "ffxc Timer")
class Timer extends AlgorithmsScript {

    @Mixin
    TimerOptions timerOptions

    @Mixin
    XrayOptions xrayOptions

    /**
     * -r or --mode sets the desired refinement mode
     * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
     */
    @Option(names = ['-r', '--mode'], paramLabel = "mode", description = 'Refinement mode: coordinates, bfactors and/or occupancies.')
    String modeString

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
    private List<String> filenames


    def run() {

        if (!init()) {
            return
        }

        xrayOptions.init()

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
            modelfilename = activeAssembly.getFile().getAbsolutePath();
        }

        logger.info("\n Running xray.Timer on " + modelfilename)

        // Set up diffraction data (can be multiple files)
        List<DiffractionData> diffractionfiles = xrayOptions.processData(filenames, assemblies);

        DiffractionData diffractiondata = new DiffractionData(assemblies, assemblies[0].getProperties(),
                xrayOptions.solventModel, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

        // Set the number of threads (needs to be done before opening the files).
        if (timerOptions.threads > 0) {
            System.setProperty("pj.nt", Integer.toString(timerOptions.threads))
        }

        diffractiondata.scaleBulkFit()
        diffractiondata.printStats()

        algorithmFunctions.energy(assemblies[0])

        // Type of refinement.
        RefinementMode refinementmode = RefinementMinimize.parseMode(modeString)

        RefinementEnergy refinementEnergy = new RefinementEnergy(diffractiondata, refinementmode)
        int n = refinementEnergy.getNumberOfVariables()
        double[] x = new double[n]
        double[] g = new double[n]
        refinementEnergy.getCoordinates(x)
        Potential energy = refinementEnergy.getDataEnergy()

        for (int i = 0; i < timerOptions.iterations; i++) {
            long time = -System.nanoTime()
            if (!timerOptions.gradient) {
                energy.energyAndGradient(x, g)
            } else {
                energy.energy(x)
            }
            time += System.nanoTime()
            logger.info(String.format(" Time %12.8f (sec)", time * 1.0e-9))
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