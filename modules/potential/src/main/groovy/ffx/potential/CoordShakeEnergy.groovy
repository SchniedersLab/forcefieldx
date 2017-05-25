
package ffx.potential

import ffx.potential.bonded.Atom


// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;
import groovy.util.CliBuilder;

// FFX Imports
import ffx.numerics.Potential;

import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.potential.OpenMMForceFieldEnergy

import java.util.function.ToDoubleFunction;

/**
 * The CoordShakeEnergy script evaluates the energy of a system before and after moving coordinates by a fixed offset.
 * Intended primarily to test coordinate restraints. Could easily be extended to do randomized offsets, etc.
 * <br>
 * Usage:
 * <br>
 * ffxc CoordShakeEnergy &lt;filename&gt;
 */
class CoordShakeEnergy extends Script {

    /**
     * Options for the Energy Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Energy [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -d or --deltaX to set the distance by which the coordinates should be moved. Applies a constant offset in x,
         * reverts, in y, reverts, then z, and reverts again.
         */
        @Option(shortName='d', longName='deltaX', defaultValue='1.0', description='Distance to move coordinates (applies to all equally)') double dX;
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }


    /**
     * Execute the script.
     */
    def run() {

        def cli = new CliBuilder(usage:' ffxc Energy [options] <filename>', header:' Options:')

        // String args[] = getArgs();

        def options = new Options();
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
            //open(modelFilename)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Running Energy on " + modelFilename)

        // This is an interface specifying the closure-like methods.
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsUtils();
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }

        // Use PotentialsFunctions methods instead of Groovy method closures to do work.
        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        Potential thePotential = activeAssembly.getPotentialEnergy();

        def eFunct;

        if (thePotential instanceof OpenMMForceFieldEnergy) {
            OpenMMForceFieldEnergy ommE = (OpenMMForceFieldEnergy) thePotential;
            eFunct = { double[] coords -> return ommE.energyVsFFX(coords, true); };
        } else {
            eFunct = { double[] coords -> return thePotential.energy(coords, true) };
        }

        logger.info(" Starting energy: ");
        double[] x = thePotential.getCoordinates();
        eFunct(x);

        Atom[] atoms = activeAssembly.getAtomArray();

        def axes = ["X","Y","Z"];
        double[] atomXYZ = new double[3];
        for (int i = 0; i < 3; i++) {
            logger.info(String.format(" Moving atoms +1 unit in %s", axes[i]));
            for (Atom atom : atoms) {
                atom.getXYZ(atomXYZ);
                atomXYZ[i] += options.dX;
                atom.setXYZ(atomXYZ);
            }
            thePotential.getCoordinates(x);
            eFunct(x);

            logger.info(String.format(" Moving atoms -1 unit in %s", axes[i]));
            for (Atom atom : atoms) {
                atom.getXYZ(atomXYZ);
                atomXYZ[i] -= (2.0 * options.dX);
                atom.setXYZ(atomXYZ);
            }
            thePotential.getCoordinates(x);
            eFunct(x);

            for (Atom atom : atoms) {
                atom.getXYZ(atomXYZ);
                atomXYZ[i] += options.dX;
                atom.setXYZ(atomXYZ);
            }
        }

        logger.info(" After movement (may have small rounding errors): ");
        thePotential.getCoordinates(x);
        eFunct(x);
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
