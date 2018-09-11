package ffx.xray.groovy

import java.util.stream.Collectors

import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.DiffractionData
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The X-ray ModelvsData script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.ModelvsData [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Compare the PDB model to the diffraction data.", name = "ffxc ModelvsData")
class ModelvsData extends AlgorithmsScript {

    @Mixin
    XrayOptions xrayOptions

    /**
     * -m or --maps Output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps.
     */
    @Option(names = ['-p', '--maps'], paramLabel = 'false', description = 'Output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps.')
    boolean maps = false
    /**
     * -t or --timings Perform FFT timings.
     */
    @Option(names = ['-t', '--timings'], paramLabel = 'false', description = 'Perform FFT timings.')
    boolean timings = false
    /**
     * -w or --mtz Write out MTZ containing structure factor coefficients.
     */
    @Option(names = ['-w', '--mtz'], paramLabel = 'false',
            description = 'write out MTZ containing structure factor coefficients.')
    boolean mtz = false

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
    private List<String> filenames
    private DiffractionData diffractiondata;
    private MolecularAssembly[] assemblies;

    @Override
    ModelvsData run() {

        if (!init()) {
            return this
        }

        xrayOptions.init()

        String modelfilename
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

        logger.info("\n Running xray.ModelvsData on " + modelfilename)

        // Load parsed X-ray properties.
        CompositeConfiguration properties = activeAssembly.getProperties()
        xrayOptions.setProperties(parseResult, properties)

        // Set up diffraction data (can be multiple files)
        List<DiffractionData> diffractionfiles = xrayOptions.processData(filenames, assemblies);

        diffractiondata = new DiffractionData(assemblies, properties,
                xrayOptions.solventModel, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

        diffractiondata.scaleBulkFit()
        diffractiondata.printStats()

        algorithmFunctions.energy(assemblies[0])

        if (mtz) {
            diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_ffx.mtz")
        }

        if (maps) {
            diffractiondata.writeMaps(FilenameUtils.removeExtension(modelfilename) + "_ffx")
        }

        if (timings) {
            diffractiondata.timings()
        }

        return this
    }

    @Override
    public List<Potential> getPotentials() {
        if (assemblies == null) {
            return new ArrayList<Potential>();
        } else {
            return Arrays.stream(assemblies).
                    filter { a -> a != null }.
                    map { a -> a.getPotentialEnergy() }.
                    filter { e -> e != null }.
                    collect(Collectors.toList());
        }
    }

    @Override
    public boolean destroyPotentials() {
        return diffractiondata == null ? true : diffractiondata.destroy();
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