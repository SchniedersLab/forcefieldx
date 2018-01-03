
package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary

/**
 * The ApplyTorsions script is used to apply a custom set of torsions to a side 
 * chain.
 * <br>
 * Usage:
 * <br>
 * ffxc ApplyTorsions [options] &lt;filename&gt;
 * @author Jacob Litman
 * @author Michael Schnieders
 */
class ApplyTorsions extends Script {

    /**
     * Options for the ApplyTorsions Script.
     * <br>
     * Usage:
     * <br>
     * ffxc ApplyTorsions [options] &lt;filename&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', longName='help', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -c or --chain selects the chain name to use.
         */
        @Option(shortName='c', longName='chain', defaultValue=' ', description='Single character chain name (default is \' \').') String chain;
        /**
         * -r or --resid selects the residue to apply torsions to.
         */
        @Option(shortName='r', longName='resid', defaultValue='1', description='Residue number.') int resID;
        /**
         * -t or --torsionSets is a variable-length, colon-delimited set of comma-separated torsion sets; e.g. colons separate rotamers, commas individual values. Should not be left as final argument, as then it attempts to include the filename.
         */
        @Option(shortName='t', longName='torsionSets', defaultValue='0,0:180,0', numberOfArgumentsString='+', valueSeparator=':', description='Torsion sets to apply (torsions comma-separated, sets colon-separated). Do not leave as last argument.') String[] torSets;
        /**
         * -n or --nChi is the number of torsions available to this side chain; torsions unspecified in -t will be filled with 0.0.
         */
        @Option(shortName='n', longName='nChi', defaultValue='1', description='Number of torsions (will fill unspecified with 0.0') int nChi;
        /**
         * -vf or --videoFile is the name of the file to print torsion snapshots to; defaults to filename_rots.pdb
         */
        @Option(shortName='vf', longName='videoFile', description='File to print torsion snapshots to.') String vidFileName;
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames;
    }
    
    def run() {
        def cli = new CliBuilder(usage:' ffxc ApplyTorsions [options] <filename>', header:' Options:');
        
        AlgorithmFunctions afuncts;
        try {
            afuncts = getAlgorithmUtils();
        } catch (MissingMethodException ex) {
            afuncts = new AlgorithmUtils();
        }
        
        def options = new Options();
        cli.parseFromInstance(options, afuncts.getArguments());

        if (options.help == true) {
            return cli.usage();
        }
        List<String> arguments = options.filenames;
        
        String filename;
        MolecularAssembly activeAssembly;
        if (arguments && !arguments.isEmpty()) {
            filename = arguments.get(0);
            activeAssembly = afuncts.open(filename)[0];
        } else {
            // Last-ditch attempt to try to find an assembly loaded through the GUI, etc.
            activeAssembly = afuncts.getActiveAssembly();
            if (activeAssembly) {
                filename = activeAssembly.getFile().getName();
            } else {
                return cli.usage();
            }
        }
        
        String newName = FilenameUtils.getBaseName(filename);
        
        if (options.vidFileName) {
            videoFile = options.vidFileName;
        } else {
            videoFile = newName + "_rots.pdb";
        }
        
        File outFile = new File(newName + ".rotout.tmp");
        outFile.deleteOnExit();
        
        logger.info("\n Saving torsions for residue number " + options.resID + " of chain " + options.chain + ".");
        
        RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains());
        Polymer polymer = activeAssembly.getChain(options.chain);
        if (polymer == null) {
            logger.info(" Polymer + " + options.chain + " does not exist.");
            return;
        }
        Residue residue = polymer.getResidue(options.resID);
        if (residue == null) {
            logger.info(" Residue + " + options.resID + " does not exist.");
            return;
        }

        String ext = FilenameUtils.getExtension(filename);
        filename = FilenameUtils.removeExtension(filename);

        AlgorithmListener aList = afuncts.getDefaultListener();
        GenerateRotamers genr = new GenerateRotamers(activeAssembly, activeAssembly.getPotentialEnergy(), residue, outFile, options.nChi, aList);
        genr.setVideo(videoFile);
        genr.applyAndSaveTorsions(options.torSets);
        
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
