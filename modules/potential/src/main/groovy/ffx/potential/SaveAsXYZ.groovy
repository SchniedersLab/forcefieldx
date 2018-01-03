package ffx.potential

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.crystal.Crystal
import ffx.crystal.SymOp
import ffx.potential.bonded.Atom
import ffx.potential.parameters.ForceField
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The SaveAsXYZ script saves a file as an XYZ file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsXYZ [options] &lt;filename&gt;
 */
class SaveAsXYZ extends Script {
    /**
     * Options for the SaveAsXYZ script
     * <br>
     * Usage:
     * <br>
     * ffxc SaveAsXYZ [options] &lt;filename&gt;
     */
    public class Options{
        /**
         * -h or --help to print a help message
         */
        @Option(longName='help', shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -p or --pos-offset to set the positive atom type offset
         */
        @Option(longName='pos-offset', shortName='p', defaultValue='0', description='Positive offset of atom types in the new file') int posOffset;
        /**
         * -n or --neg-offset to set the negative atom type offset
         */
        @Option(longName='neg-offset', shortName='n', defaultValue='0', description='Negative offset of atom types in the new file.') int negOffset
        /**
         * -r or --random to apply a random Cartesian symmetry operator with the specified translation range -X .. X (no default).
         */
        @Option(longName='random', shortName='r', defaultValue='-1.0', description='Apply a random Cartesian SymOp with translation range -X .. X.') double scalar
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }


    /**
     * Execute the script.
     */
    def run() {
        int offset = 0;

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc SaveAsXYZ [options] <filename>');
        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        // Positive offset atom types.
        if (options.posOffset > 0) {
            offset = options.posOffset;
        }

        // Negative offset atom types.
        if (options.negOffset > 0) {
            offset = options.negOffset;
            offset = -offset;
        }

        List<String> arguments = options.filenames;

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

        logger.info("\n Writing out XYZ for " + modelFilename);

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

        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        modelFilename = FilenameUtils.removeExtension(modelFilename) + ".xyz";

        //open(modelFilename);
        // Offset atom type numbers.
        if (offset != 0) {
            logger.info("\n Offset atom types by " + offset);
            ForceField forceField = activeAssembly.getForceField();
            forceField.renumberForceField(0,offset,0);
        }

        if (options.scalar > 0.0) {
            SymOp symOp = SymOp.randomSymOpFactory(options.scalar);

            logger.info(String.format("\n Applying random Cartesian SymOp\n: %s", symOp.toString()));

            Crystal crystal = activeAssembly.getCrystal();
            Atom[] atoms = activeAssembly.getAtomArray();
            double[] xyz = new double[3];
            for (int i=0; i<atoms.length; i++) {
                atoms[i].getXYZ(xyz);
                crystal.applyCartesianSymOp(xyz, xyz, symOp);
                atoms[i].setXYZ(xyz);
            }
        }

        functions.save(activeAssembly, new File(modelFilename));

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
