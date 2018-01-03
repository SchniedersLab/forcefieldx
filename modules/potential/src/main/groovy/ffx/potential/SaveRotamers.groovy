
package ffx.potentials

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The SaveRotamers script saves a file as a PDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveRotamers [options] &lt;filename&gt;
 */
class SaveRotamers extends Script {
    /**
     * Options for the SaveRotamers
     * <br>
     * Usage:
     * <br>
     * ffxc SaveRotamers [options] &lt;filename&gt;
     */
    public class Options{
        /**
         * -h or --help to print a help message
         */
        @Option(longName='help', shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -c or --chain to specify chain
         */
        @Option(longName='chain', shortName='c', description='Single character chain name (default is \' \').') char chain;
        /**
         * -l or --library to select rotamer library 
         */
        @Option(longName='library', shortName='l', description='Available rotamer libraries are Ponder and Richards (1) or Richardson (2).') int library;
        /**
         * -r or --resid to select residue number 
         */
        @Option(longName='resid', shortName='r', description='Residue number.') int resid;
        /**
         * -i or --independent to draw nucleic acid rotamers independently of chain context
         */ 
        @Option(longName='independent', shortName='i', defaultValue='false', description='Independent draws nucleic acid rotamers independently of chain context (true if flag is present).') boolean independent;
        /**
         * -s or --start to select first rotamer to draw. Indexed form rotamer 0
         */
        @Option(longName='start', shortName='s', description='First rotamer to draw. Indexed from rotamer 0.') int start;
        /**
         * -f or --finish to select last rotamer to draw. Indexed from rotamer 0
         */
        @Option(longName='finish', shortName='f', description='Last rotamer to draw. Indexed from rotamer 0.') int finish;
        /**
         * -x or --all to draw all rotamers beginning from the passed rotamer number (overrides other options). Indexed from rotamer 0.
         */
        @Option(longName='all', shortName='x', description='Draw all rotamers beginning from the passed rotamer number (overrides other options). Indexed from rotamer 0.') int all;
        /**
         * -u or --upstreamPucker to adjust the pucker of the 5\' residue to match the rotamer. Use flag to turn on
         */ 
        @Option(longName='upstreamPucker', shortName='u', defaultValue='false', description='Adjusts the pucker of the 5\' residue to match the rotamer (true if flag is present)') boolean upstreamPucker;
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames;
    }

    /**
     * Execute the script.
     */
    def run() {
        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc SaveRotamers [options] <PDB|XYZ>');
        def options = new Options()
        cli.parseFromInstance(options, args)
        if (options.help == true) {
            return cli.usage()
        }
    
        List<String> arguments = options.filenames;
  
        int resID = 1;
        String chain = " ";
        int library = 1;
        boolean independent = false;
        int start = -1;
        int finish = -1;
        int allStart = 0;
        // Will be checked for validity, and set false if invalid.
        boolean upstreamPucker = false;
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();

        // Residue number.
        if (options.resid) {
            resID = options.resid;
        }

        // Chain Name.
        if (options.chain) {
            chain = Character.toString(options.chain);
        }

        // Rotamer Library.
        if (options.library) {
            library = options.library;
        }

        // Rotamer independence
        if (options.independent) {
            independent = options.independent;
        }

        // Default behavior: draw all rotamers starting from 0 (as though set -x 0).
        boolean saveAllRotamers = true;

        // Start point
        if (options.start) {
            start = options.start;
            saveAllRotamers = false;
        }

        // Finish point
        if (options.finish) {
            finish = options.finish;
        }

        // Start from to all
        if (options.all) {
            allStart = options.all;
            saveAllRotamers = true;
            // Over-rides options.s
        }

        // Upstream pucker adjustment
        if (options.upstreamPucker) {
            upstreamPucker = options.upstreamPucker;
        }

        logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".");
        
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsFunctions()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }

        // Read in command line.
        //String modelFilename = arguments.get(0);

        // Open the file.
        //systems = open(modelFilename);
        // Read in command line.
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

        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains());
        Polymer polymer = activeAssembly.getChain(chain);
        if (polymer == null) {
            logger.info(" Polymer + " + chain + " does not exist.");
            return;
        }
        Residue residue = polymer.getResidue(resID);
        if (residue == null) {
            logger.info(" Residue + " + resID + " does not exist.");
            return;
        }

        if (library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        }

        Rotamer[] rotamers = residue.getRotamers(rLib);
        if (rotamers == null) {
            logger.severe(" There are no rotamers for residue + " + residue.toString());
        }
        int nrotamers = rotamers.length;

        boolean isDeoxy = false; // Applies to prevResidue.
        Residue prevResidue;
        // If upstreamPucker is specified true, ensure that it is valid to apply it
        // (residue is nucleic acid with an upstream partner), and set isDeoxy.
        // If it is invalid, set upstreamPucker false.
        if (upstreamPucker) {
            // Exception gets thrown if it's an amino acid, since "NA" is undefined.
            try {
                if (residue.getResidueType() == NA) {
                    prevResidue = (Residue) residue.getPreviousResidue();
                    // If no previous residue, set upstream pucker false.
                    // The method used will ensure prevResidue is a nucleic acid.
                    if (prevResidue == null) {
                        upstreamPucker = false;
                    } else {
                        Atom HOs = (Atom) prevResidue.getAtomNode("HO\'");
                        if (HOs == null) {
                            isDeoxy = true;
                        }
                    }
                }  else {
                    upstreamPucker = false;
                }
            } catch (Exception e) {
                upstreamPucker = false;
            }
        }

        String ext = FilenameUtils.getExtension(modelFilename);
        modelFilename = FilenameUtils.removeExtension(modelFilename);

        if (saveAllRotamers) {
            if (allStart >= nrotamers) {
                logger.info(" Specified start range is outside of rotamer range. No action taken.");
            } else {
                for (int i = allStart; i < nrotamers; i++) {
                    RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
                    if (upstreamPucker) {
                        double prevDelta = rotamers[i].chi1;
                        if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                            // North pucker
                            RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true);
                        } else {
                            RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true);
                        }
                    }
                    if (ext.toUpperCase().contains("XYZ")) {
                        logger.info(String.format("Saving rotamer %d", i))
                        functions.saveAsXYZ(activeAssembly,new File(modelFilename + ".xyz"));
                    } else {
                        //functions.saveAsPDB(systems, new File(modelFilename + ".pdb"));
                        logger.info(String.format("Saving rotamer %d", i))
                        functions.saveAsPDB(activeAssembly, new File(modelFilename + ".pdb"));
                    }
                }
            }
        } else {
            if (start >= nrotamers) {
                logger.info(" Specified start range is outside of rotamer range. No action taken.");
            } else {
                if (finish >= nrotamers) {
                    finish = nrotamers - 1;
                } else if (finish < start) {
                    logger.info(" Specified finish point is before the start point; drawing only rotamer " + start);
                    finish = start;
                }
                for (int i = start; i <= finish; i++) {
                    RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
                    if (upstreamPucker) {
                        double prevDelta = rotamers[i].chi1;
                        if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                            // North pucker
                            RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true);
                        } else {
                            RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true);
                        }
                    }
                    if (ext.toUpperCase().contains("XYZ")) {
                        logger.info(String.format("Saving rotamer %d", i))
                        functions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"));
                    } else {
                        //saveAsPDB(systems, new File(modelFilename + ".pdb"));
                        logger.info(String.format("Saving rotamer %d", i))
                        functions.saveAsPDB(activeAssembly, new File(modelFilename + ".pdb"));
                    }
                }
            }
        }
    }
}


/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
