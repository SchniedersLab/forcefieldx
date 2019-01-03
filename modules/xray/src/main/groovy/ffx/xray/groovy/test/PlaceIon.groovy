package ffx.xray.groovy.test

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.AlgorithmFunctions
import ffx.algorithms.AlgorithmUtils
import ffx.crystal.Crystal
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.xray.CrystalReciprocalSpace.SolventModel
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.parsers.DiffractionFile

/**
 * The X-ray PlaceIon script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.PlaceIon [options] &lt;filename [file2...]&gt;
 */
class PlaceIon extends Script {

    /**
     * Options for the PlaceIon Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.PlaceIon [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -S or --spacing Specify the spacing in fractional coordinates.
         */
        @Option(shortName = 'S', longName = 'spacing', defaultValue = '0.05', description = 'Search spacing in fractional coordinates.')
        double spacing;
        /**
         * -s or --suffix Specify the suffix to apply to output files. For example, for 1abc_refine.pdb, write out 1abc_refine_refine.[pdb|mtz] at the end.
         */
        @Option(shortName = 's', longName = 'suffix', defaultValue = '_refine', description = 'Suffix to apply to files written out by minimization.')
        String suffix;
        /**
         * -D or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
         */
        @Option(shortName = 'D', longName = 'data', defaultValue = '', numberOfArguments = 3, valueSeparator = ',',
                description = 'Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.')
        String[] data
        /**
         * The final arguments should be a PDB filename and data filename (CIF or MTZ).
         */
        @Unparsed(description = "PDB file and a CIF or MTZ file.")
        List<String> filenames
    }

    @Override
    PlaceIon run() {

        def cli = new CliBuilder()
        cli.name = "ffxc xray.PlaceIon"

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            cli.usage()
            return this
        }

        AlgorithmFunctions aFuncts
        try {
            // getAlgorithmUtils is a magic variable/closure passed in from ModelingShell
            aFuncts = getAlgorithmUtils()
        } catch (MissingMethodException ex) {
            // This is the fallback, which does everything necessary without magic names
            aFuncts = new AlgorithmUtils()
        }

        List<String> arguments = options.filenames

        String modelfilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelfilename = arguments.get(0)
        } else if (active == null) {
            cli.usage()
            return this
        } else {
            modelfilename = active.getFile()
        }

        logger.info("\n Running xray.PlaceIon on " + modelfilename)

        MolecularAssembly[] systems = aFuncts.openAll(modelfilename)

        // Set up diffraction data (can be multiple files)
        List diffractionfiles = new ArrayList()
        if (arguments.size() > 1) {
            DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        if (options.data) {
            for (int i = 0; i < options.data.size(); i += 3) {
                double wA = Double.parseDouble(options.data[i + 1])
                boolean neutron = Boolean.parseBoolean(options.data[i + 2])
                DiffractionFile diffractionfile = new DiffractionFile(options.data[i], wA, neutron)
                diffractionfiles.add(diffractionfile)
            }
        }

        if (diffractionfiles.size() == 0) {
            DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false)
            diffractionfiles.add(diffractionfile)
        }

        DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(),
                SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

        diffractiondata.scaleBulkFit()
        diffractiondata.printStats()

        aFuncts.energy(systems[0])

        MolecularAssembly molecularAssembly = systems[0]
        ArrayList<MSNode> ions = molecularAssembly.getIons()
        if (ions == null || ions.size() == 0) {
            logger.info(" Please add an ion to the PDB file to scan with.")
            return
        }

        // Scan with the last ion in the file.
        int nIon = ions.size()
        Atom ion = ions.get(nIon - 1).getAtomList().get(0)
        logger.info(" Preparing to scan with ion:\n " + ion.toString())

        RefinementMinimize refinementMinimize = new RefinementMinimize(diffractiondata, RefinementMode.COORDINATES)
        RefinementEnergy refinementEnergy = refinementMinimize.refinementEnergy
        int n = refinementEnergy.getNumberOfVariables()
        double[] x = new double[n]

        Crystal crystal = molecularAssembly.getCrystal()

        double spacing = options.spacing
        double[] fracCoords = new double[3]
        double[] realCoords = new double[3]
        double[] minCoords
        double eMin = Double.MAX_VALUE

        refinementEnergy.getCoordinates(x)
        double eStart = refinementEnergy.energy(x, true)
        double[] startingCoords = [ion.getX(), ion.getY(), ion.getZ()]

        // Scan the ion over the fractional coordinates of the unit cell.
        for (double ai = 0; ai <= 1.0; ai += spacing) {
            fracCoords[0] = ai
            for (double bi = 0; bi <= 1.0; bi += spacing) {
                fracCoords[1] = bi
                for (double ci = 0; ci <= 1.0; ci += spacing) {
                    fracCoords[2] = ci
                    // Get the real space coordinates.
                    crystal.toCartesianCoordinates(fracCoords, realCoords)
                    ion.setXYZ(realCoords)
                    if (crystal.isSpecialPosition(realCoords)) {
                        logger.info(" Ion special position skipped " + ion.toString())
                        continue
                    }
                    long time = -System.nanoTime()
                    refinementEnergy.getCoordinates(x)
                    double e = refinementEnergy.energy(x, true)
                    time += System.nanoTime()
                    if (e < eMin) {
                        eMin = e
                        minCoords = [realCoords[0], realCoords[1], realCoords[2]]
                        logger.info(String.format(" New minimum %16.8f at (%12.6f %12.6f %12.6f) in %6.3f (sec)",
                                eMin, realCoords[0], realCoords[1], realCoords[2], time * 1.0e-9))
                    } else {
                        logger.info(String.format(" Energy      %16.8f at (%12.6f %12.6f %12.6f) in %6.3f (sec)",
                                e, realCoords[0], realCoords[1], realCoords[2], time * 1.0e-9))
                    }
                }
            }
        }

        if (eMin < Double.MAX_VALUE) {
            logger.info(String.format(" Best minimum   %16.8f found at (%16.8f %16.8f %16.8f)", eMin, minCoords[0], minCoords[1], minCoords[2]))
        }

        if (eMin < eStart) {
            // Set the ion coordinates to the minimum found.
            ion.setXYZ(minCoords)

            // Report the force field energy.
            aFuncts.energy(systems[0])

            // Write out the results.
            diffractiondata.writeModel(FilenameUtils.removeExtension(modelfilename) + options.suffix + ".pdb")
            diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + options.suffix + ".mtz")
        } else {
            logger.info(String.format(" Initial position %16.8f kept at (%16.8f %16.8f %16.8f)",
                    eStart, startingCoords[0], startingCoords[1], startingCoords[2]))
        }

        return this
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