package ffx.potential.test

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Residue.ResidueType
import ffx.potential.parameters.ForceField
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The MultiResidue script evaluates the energy of a MultiResidue system.
 * <br>
 * Usage:
 * <br>
 * ffxc test.MultiResidue [options] &lt;filename&gt;
 */
class MultiResidue extends Script {

    /**
     * Options for the MultiResidue script.
     * <br>
     * Usage:
     * <br>
     * ffxc test.MultiResidue [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message.
         */
        @Option(shortName='h', description='Print this help message.') boolean help
        /**
         * -r or --resID to define the residue number (default is 1).
         */
        @Option(shortName='r', defaultValue='1', description='Residue number.') Integer resID
        /**
         * -c or --chain to set the single character chain name (default is \' \').'
         */
        @Option(shortName='c', defaultValue=' ', description='Single character chain name (default is \' \').') Character chain
        /**
         * -n or --name to set the name of residue to switch (default is ALA).
         */
        @Option(shortName='n', defaultValue='ALA', description='Name of residue to switch to.') String name
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc test.MultiResidue [options] <filename>', header:' Options:');
        def options = new Options()
        cli.parseFromInstance(options, args)
        if (options.help) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Running MultiResidue on " + modelFilename + "\n");

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
        ffx.potential.MolecularAssembly[] assemblies = functions.open(modelFilename)
        ForceField forceField = assemblies[0].getForceField();
        ffx.potential.ForceFieldEnergy forceFieldEnergy = assemblies[0].getPotentialEnergy();

        int resID = options.resID
        Character chain = options.chain
        String name = options.name

        ffx.potential.bonded.MultiResidue multiResidue;
        Residue residue;
        Polymer[] polymers = assemblies[0].getChains();
        for (int i = 0; i < polymers.length; i++) {
            Polymer polymer = polymers[i];
            if (chain.equals(polymer.getChainID())) {
                residue = polymer.getResidue(resID);
                if (residue != null) {
                    multiResidue = new ffx.potential.bonded.MultiResidue(residue, forceField, forceFieldEnergy);
                    polymer.addMultiResidue(multiResidue);
                }
            }
        }

        if (residue == null) {
            logger.info(" Chain " + chain + " residue " + resID + " was not found." );
            return;
        }

        ResidueType type = residue.getResidueType();
        int resNumber = residue.getResidueNumber();
        multiResidue.addResidue(new Residue(name, resNumber, type));



        int numResidues = multiResidue.getResidueCount();
        for (int i=0; i<numResidues; i++) {
            multiResidue.setActiveResidue(i);
            logger.info(" Active Residue: " + multiResidue.toString());
            forceFieldEnergy.energy(true,true);
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
