package ffx.xray

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.crystal.Crystal
import ffx.crystal.ReflectionList
import ffx.crystal.Resolution
import ffx.potential.MolecularAssembly
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils
import ffx.xray.parsers.CIFFilter
import ffx.xray.parsers.MTZWriter
import ffx.xray.parsers.MTZWriter.MTZType

/**
 * The CIF2MTZ script saves a CIF file to MTZ format.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.CIF2MTZ &lt;filename&gt;
 */
class CIF2MTZ extends Script {

    /**
     * Options for the CIF2MTZ Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.CIF2MTZ &lt;filename&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help;
        /**
         * The final argument should be an CIF filename.
         */
        @Unparsed(description = 'A PDB file and a CIF file.')
        List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {
        def cli = new CliBuilder();
        cli.name = "ffxc xray.CIF2MTZ";

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String datafilename
        String modelfilename
        if (arguments != null && arguments.size() > 1) {
            // Read in command line.
            modelfilename = arguments.get(0)
            datafilename = arguments.get(1)
        } else {
            return cli.usage()
        }

        logger.info("\n Running CIF2MTZ on " + modelfilename + " and " + datafilename + "\n")

        // This is an interface specifying the closure-like methods.
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsUtils()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }

        // Use PotentialsFunctions methods instead of Groovy method closures to do work.
        MolecularAssembly[] systems = functions.open(modelfilename)

        CIFFilter ciffilter = new CIFFilter();
        ReflectionList reflectionlist = ciffilter.getReflectionList(new File(datafilename), systems[0].getProperties());

        if (reflectionlist == null) {
            println(" Using crystal information from PDB to generate MTZ file.");

            Crystal crystal = systems[0].getCrystal().getUnitCell();
            double res = ciffilter.getResolution(new File(datafilename), crystal);
            if (res < 0.0) {
                println(" Resolution could not be determined from PDB and CIF file.");
                return;
            }

            Resolution resolution = new Resolution(res);
            reflectionlist = new ReflectionList(crystal, resolution, systems[0].getProperties());
        }

        DiffractionRefinementData refinementdata = new DiffractionRefinementData(systems[0].getProperties(), reflectionlist);
        ciffilter.readFile(new File(datafilename), reflectionlist, refinementdata, systems[0].getProperties());

        MTZWriter mtzwriter = new MTZWriter(reflectionlist, refinementdata, FilenameUtils.removeExtension(datafilename) + ".mtz", MTZType.DATAONLY);
        mtzwriter.write();
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
